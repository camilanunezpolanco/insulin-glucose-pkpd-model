# Insulin-Glucose PK/PD Model Reconstruction

This project reconstructs and explores the insulin pharmacokinetic/pharmacodynamic framework described by Berger and Rodbard (1989). The goal is to examine whether the published equations and parameters are sufficient to reproduce the reported time courses for subcutaneous insulin absorption and plasma insulin.

## Project goals

- Reconstruct the published subcutaneous insulin absorption model.
- Simulate plasma insulin profiles for regular, NPH, lente, and ultralente insulin.
- Compare the timing and shape of the simulated profiles with the published figures.
- Explore parameter sensitivity and dynamical stability.
- Document ambiguities that affect reproducibility.

## Model overview

For dose `D`, the formulation-specific half-absorption time is

```text
T50(D) = aD + b.
```

The percentage of the injected dose remaining at the subcutaneous depot is

```text
A_depot(t) = 100 T50^s / (T50^s + t^s),
```

where `s` controls the absorption-curve shape. Plasma insulin is modeled as the balance between absorption from the depot and first-order elimination:

```text
dQ/dt = absorption_rate(t) - ke Q.
```

The clean implementation reports plasma concentration as `Q / Vi`. It also provides normalized profiles for shape-based comparison when the original publication does not provide enough scaling detail for an exact numerical reconstruction.

## Repository structure

```text
.
├── notebooks/
│   └── pkpd_model_exploration.ipynb
├── src/
│   └── pkpd_model.py
├── presentation/
│   └── PKPD_Model_Reproducibility.pdf
├── CITATION.cff
├── requirements.txt
├── .gitignore
└── README.md
```

## Running the clean simulation

Create a Python environment and install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the model:

```bash
python src/pkpd_model.py
```

The script creates:

- `outputs/insulin_profiles.png`
- `outputs/summary_metrics.csv`
- `outputs/local_sensitivity.csv`

## Notebook scope

The notebook contains the full exploratory workflow used during the project, including:

- Alternative scalings and replications of the published curves
- Local and Sobol sensitivity experiments
- Symbolic Jacobian and equilibrium calculations
- Exploratory insulin-glucose extensions

Some notebook sections are research experiments rather than a single finalized model. The script in `src/` is the recommended reproducible entry point for the core absorption and plasma-insulin simulation.

## Reproducibility findings

The reconstructed curves preserve the expected qualitative ordering:

- Regular insulin peaks earliest.
- NPH and lente produce intermediate profiles.
- Ultralente produces a flatter, longer-lasting profile.

Small numerical differences from the published figures can arise from rounded parameters, solver choice, normalization, and details not fully specified in the article. This project therefore distinguishes qualitative replication from exact numerical validation.

## Reference

Berger, M., & Rodbard, D. (1989). Computer simulation of plasma insulin and glucose dynamics after subcutaneous insulin injection. *Diabetes Care, 12*(10), 725-736. [https://doi.org/10.2337/diacare.12.10.725](https://doi.org/10.2337/diacare.12.10.725)

## Author

Camila Nunez Polanco, University of Massachusetts Amherst

