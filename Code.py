import numpy as np
import matplotlib.pyplot as plt
from bioservices import KEGG
from scipy.integrate import solve_ivp
import networkx as nx
import plotly.graph_objects as go
import pandas as pd
import matplotlib.pyplot as plt



# Fetch reactions for a molecule
def fetch_reactions_for_molecule(molecule_id):
    """
    Fetch and neatly format reactions from KEGG involving the molecule of interest.
    """
    from tabulate import tabulate

    kegg = KEGG()
    results = kegg.find("reaction", molecule_id)
    reactions = results.split("\n") if results else []

    parsed_reactions = []
    skipped_lines = []

    for r in reactions:
        if "\t" in r:
            reaction_id, description = r.split("\t", 1)
            parsed_reactions.append((reaction_id, description))
        else:
            skipped_lines.append(r)

    print("=" * 50)
    print(f"Reactions involving {molecule_id}: {len(parsed_reactions)} found")
    print("=" * 50)

    if parsed_reactions:
        print(tabulate(parsed_reactions, headers=["Reaction ID", "Description"], tablefmt="grid"))
    else:
        print("No valid reactions found.")

    if skipped_lines:
        print("\nSkipped lines due to formatting issues:")
        for line in skipped_lines:
            print(f" - {line}")
    else:
        print("\nNo malformed lines encountered.")

    return parsed_reactions


# Build the reaction network
def build_reaction_network(reactions):
    """
    Create a reaction network using NetworkX for visualization.
    """
    G = nx.DiGraph()
    for reaction_id, description in reactions:
        try:
            substrates, products = description.split("<=>")
            substrates = substrates.split(" + ")
            products = products.split(" + ")
            for s in substrates:
                for p in products:
                    G.add_edge(s.strip(), p.strip(), reaction=reaction_id)
        except ValueError:
            print(f"Skipping malformed reaction description: {description}")
    return G


# Simulate reaction dynamics
def simulate_reaction_dynamics(reactions, initial_concentrations, time_span):
    """
    Simulate the system dynamics using ODEs, considering cascading reactions.
    """
    equilibrium_constants = {reaction[0]: np.random.uniform(0.1, 10) for reaction in reactions}

    # Parse species
    species = list(set(sum([
        [s.strip() for s in r[1].split("<=>")[0].split(" + ")] +
        [p.strip() for p in r[1].split("<=>")[1].split(" + ")]
        for r in reactions], [])))
    num_species = len(species)

    # Ensure initial_concentrations matches the number of species
    if len(initial_concentrations) < num_species:
        initial_concentrations += [0] * (num_species - len(initial_concentrations))

    def reaction_kinetics(t, y):
        dydt = np.zeros_like(y)
        for reaction_id, description in reactions:
            try:
                substrates, products = description.split("<=>")
                substrates = substrates.split(" + ")
                products = products.split(" + ")
                substrate_indices = [species.index(s.strip()) for s in substrates if s.strip() in species]
                product_indices = [species.index(p.strip()) for p in products if p.strip() in species]
                k_eq = equilibrium_constants[reaction_id]
                forward_flux = np.prod([y[idx] for idx in substrate_indices if idx < len(y)])
                reverse_flux = k_eq * np.prod([y[idx] for idx in product_indices if idx < len(y)])
                for idx in substrate_indices:
                    if idx < len(y):
                        dydt[idx] -= forward_flux - reverse_flux
                for idx in product_indices:
                    if idx < len(y):
                        dydt[idx] += forward_flux - reverse_flux
            except ValueError:
                print(f"Skipping malformed reaction: {description}")
        return dydt

    solution = solve_ivp(reaction_kinetics, time_span, initial_concentrations[:num_species],
                         t_eval=np.linspace(time_span[0], time_span[1], 100))
    return solution.t, solution.y, species



# Plot dynamics
def plot_dynamics(t, y, species, initial_concentration):
    """
    Plot the dynamics of the chemical reaction system.
    """
    plt.figure(figsize=(10, 6))
    for idx, sp in enumerate(species):
        plt.plot(t, y[idx], label=f"{sp} (Initial: {initial_concentration:.2f})")

    # Set legend to the right of the plot
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=10)
    plt.title(f"Reaction Dynamics (Initial Concentration: {initial_concentration:.2f})", fontsize=14)
    plt.xlabel("Time", fontsize=12)
    plt.ylabel("Concentration", fontsize=12)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout(pad=3)  # Add padding for layout adjustment
    plt.show()

def analyze_changes(y, species):
    """
    Analyze the changes in concentrations and output major and minor changes.
    """
    final_concentrations = y[:, -1]
    initial_concentrations = y[:, 0]

    changes = final_concentrations - initial_concentrations
    change_data = {
        "Species": species,
        "Initial Concentration": initial_concentrations,
        "Final Concentration": final_concentrations,
        "Change": changes,
    }

    df = pd.DataFrame(change_data)
    df["Magnitude"] = pd.cut(abs(df["Change"]), bins=[0, 0.1, 1, float("inf")], labels=["Minor", "Moderate", "Major"])
    print("\nSummary of Changes in Concentrations:")
    print(df.to_string(index=False, float_format="{:.2f}".format))

# Main script
if __name__ == "__main__":
    molecule_id = input("Enter the KEGG ID of the molecule (e.g., C00031 for Glucose): ")
    time_span = (0, 10)
    initial_concentration_range = [0.1, 0.5, 1.0, 2.0]

    reactions = fetch_reactions_for_molecule(molecule_id)
    if reactions:
        G = build_reaction_network(reactions)

        for initial_concentration in initial_concentration_range:
            initial_concentrations = [initial_concentration] + [0.1] * (len(reactions) - 1)
            t, y, species = simulate_reaction_dynamics(reactions, initial_concentrations, time_span)
            plot_dynamics(t, y, species, initial_concentration)
            analyze_changes(y, species)
