# Final Year Project: Adaptive MCMC Approach 

This repository contains the codebase, analytical pipelines, simulation models, and LaTeX thesis source for a Final Year Project (FYP) in **Data Science and Analytics (DSA)** at the **National University of Singapore (NUS)**.

The project was conducted at the **NUS Saw Swee Hock School of Public Health** under the supervision of **Prof. Alex R Cook**. The research focuses on the application of computational statistics to epidemiological modeling, specifically utilizing Bayesian inference via Adaptive MCMC and stochastic SIR simulations to evaluate public health interventions (e.g., school closures).

---

## - Repository Overview

```
.
├── Kiguchi_FYP.Rproj   # RStudio Project File
├── renv.lock                      # R package dependencies lockfile
├── renv/                          # renv environment configuration
├── data/
│   ├── raw/                       # Raw contact data
│   └── processed/                 # Stores MCMC chain results
├── scripts/
│   ├── src/                       # R helper functions & core MCMC code
│   ├── Data_*.qmd                 # EDA of POLYMOD dataset
│   ├── Network_01_*.qmd           # Network sanity checks & Bootstrap analyses
│   ├── Network_02-06_*.qmd/.R     # Bootstraping individuals from the POLYMOD participants
│   ├── Network_07-16_*.qmd/.R     # Adaptive MCMC estimation (Chains 1–3 for M1, M2A, M2B)
│   ├── Network_17-23_*.qmd        # Data prep, traceplots, diagnostics & model assessment
│   └── SIR_01-05_*.qmd            # Stochastic SIR simulations & visualization on the network
├── output/                        # Generated figures, tables, and simulation results
└── thesis/                        # LaTeX source code for the final thesis document
    ├── main.tex                   # Main LaTeX compiler entrypoint
    ├── references.bib             # Bibliography file
    └── sections/                  # Modular thesis chapter files
```

---

## - Analytical Pipeline & Workflow

The workflow is structured into three main phases:

### 1. Data Processing & Network Inference
- **Data Cleaning & EDA**: `Data_01_cleaning.qmd` and `Data_02_EDA.qmd` process raw contact matrices and demographic attributes.
- **Model Fitting (MCMC & Bootstrap)**:
  - Fits network models (**M1**, **M2A**, **M2B**) using custom Adaptive MCMC algorithms (`Network_07_M1_adaptive_mcmc.R`).
  - Evaluates convergence via multiple MCMC chains and produces traceplots.
  - Performs diagnostic checks and model comparison.

### 2. Stochastic Epidemic Simulation (SIR)
- **Synthetic Network Generation**: Constructs realistic 10,000-node networks incorporating demographic stratification.
- **Epidemic Dynamics**: Runs SIR simulation algorithms on inferred networks.
- **Targeted Interventions**: Simulates school-specific intervention strategies.
- **Visualization**: Generates comparative epidemic curves and infection trajectory visualizations.

### 3. Thesis Document Construction
- Source files located in `thesis/` can be compiled into the final PDF document using LaTeX.

---

## - Getting Started & Reproducibility

### Prerequisites
- **R** (>= 4.2.0)
- **RStudio** (Recommended as the primary IDE)
- **Quarto CLI** (for rendering `.qmd` notebooks)
- **TeX Live / MacTeX** (for LaTeX thesis compilation)

### Package Environment Setup
This project uses [`renv`](https://rstudio.github.io/renv/) to manage package dependencies and ensure full reproducibility.

1. Clone this repository:
   ```bash
   git clone https://github.com/yoheionigiri/final-year-project.git
   cd final-year-project
   ```

2. Open `Kiguchi_FYP_Submission.Rproj` in RStudio or launch R in the project root.

3. Restore the R environment dependencies:
   ```R
   renv::restore()
   ```

---

## - Author & Acknowledgments

- **Author**: Yohei Kiguchi (NUS College, Data Science and Analytics)
- **Supervisor**: Prof. Alex R Cook (NUS Saw Swee Hock School of Public Health)
- **Scope**: Epidemiological Network Modeling, Computational Statistics, and Global Health.