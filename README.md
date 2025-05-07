# KEGG Reaction Dynamics Simulator

This Python script fetches biochemical reactions from the **KEGG database** for a given molecule ID and simulates the **reaction dynamics** using ODEs (ordinary differential equations). The reactions are also visualized as a **network graph**, and each species' concentration change is analysed and plotted.

---

## 🔬 Features

- **Fetch Reactions**: Get all KEGG reactions involving a specified molecule.
- **Build Reaction Network**: Construct a directional graph of substrate-product pairs.
- **Simulate Reactions**: Solve a system of differential equations modeling the reaction kinetics.
- **Plot Concentration Dynamics**: Visualise how concentrations of each species evolve over time.
- **Analyse Changes**: Quantitatively analyze concentration changes (minor, moderate, major).

---

## 📦 Requirements

Install the required libraries using:

```bash
pip install numpy matplotlib bioservices scipy networkx plotly pandas tabulate
```

---
How It Works
- Input: User enters a KEGG molecule ID (e.g., C00031 for glucose).

- Fetch Reactions: All related reactions from KEGG are parsed and formatted.

- Build Network: A directed reaction graph is constructed.

- Simulate Dynamics:

- - Reactions are converted into differential equations.

- - The system is solved over a defined period using solve_ivp.

- Plot Results: Time-series plots show species concentrations.

- Analyse: Initial and final concentrations are compared and categorised.
