# Biolum

This repository contains all code and data used in our publication [A stochastic oscillator model simulates the entrainment of vertebrate cellular clocks by light](https://www.nature.com/articles/s41598-021-93913-2).

## Installation 
To install necessery dependencies, start [Julia](https://julialang.org/) in the project folder and then instantiate a package enviroment by running `] activate .` and `] instantiate`.

This project was written in [Julia for VSCode](https://www.julia-vscode.org/). Other popular editors include an Atom extension [Juno](https://junolab.org/) and [Jupyter notebooks](https://jupyter.org/).

## Luminescence recordings
Luminescence recordings are stored in [`data/biolumdata.xlsx`](data/biolumdata.xlsx). The file contains 8 sheets. Each sheet contains traces obtained from one cell culture plate. The first column on each sheet is called `Time` and contains time in hours from the beginning of the recording. The second column is called `Light` and decodes whether the light was turned on (1) or off (0) at the given time. The remaining columns contain luminescence recordings in counts per second (cps). The plate names correspond to the publication followingly
- `Plate U1`: Fitting data (untreated).
- `Plate U2`: Validation data (untreated).
- `Plate U3`: 15:15 LD cycle (untreated).
- `Plate U4`: 10:10 LD cycle (untreated).
- `Plate D1 A`: LD cycle (set A, pharmacological treatments FOR, DBC, U0126).
- `Plate D2 A`: Constant darkness (set A, pharmacological treatments FOR, DBC, U0126).
- `Plate D1 B`: LD cycle (set B, pharmacological treatments EGF, PMA, RO).
- `Plate D2 B`: Constant darkness (set B, pharmacological treatments EGF, PMA, RO).

## Figures
All files strating with `figure_` contain scripts to generate figures from the paper. The generated figures are automatically stored in folder `figures` as SVG files.
- `figure_population_desynchronization.jl`: Figure 1C.
- `figure_C_validation.jl`: Figure 2A, B, C, D.
- `figure_prc.jl`: Figure 2E.
- `figure_drug_fit.jl`: Figure 3A, B and Supplementary Figure S4 and Supplementary Table S2.
- `figure_pca.jl`: Figure 3C and Supplementary Figure S3.
- `figure_parameters.jl`: Figure 3D.
- `figure_fig4_addition.jl`: Figure 4A, B.
- `figure_fig4.jl`: Figure 4C.
- `figure_entrainment_phase.jl`: Figure 4D and Supplementary Figure S5.
- `figure_pacing_protocols.jl`: Supplementary Figure S1.
- `figure_noise.jl`: Supplementary Figure S2.
- `figure_data_norm.jl`: Supplementary Figure S6.
- `figure_drug_norm.jl`: Supplementary Figures S7 and S8.
- `figure_K.jl`: Supplementary Figure S9.
- `figure_trajectories.jl`: Supplementary Figure S10.

## Scripts
Scripts to generate results from the paper. Output files are stored as CSV files in folder `outputs`.
- `script_A_optim.jl`: The first step of the optimization routine. Fitting properties of a deterministic model.
- `script_B_optim.jl`: The second step of the optimization routine. Fitting the noise intensity based on the damping ratio in constant dark.
- `script_C_optim.jl`: The third step of the optimization routine. Fitting data from untreated cells.
- `script_D_initial_conditions.jl`: Generate initial conditions for the drug data.
- `script_D_optim_A.jl`: Fit traces from drug plate set A.
- `script_D_optim_B.jl`: Fit traces from drug plate set B.
- `script_damping_ratio.jl`: Fit damped sine to the data.
- `script_drug_norm.jl`: Normalize trajectories of pharmacological treatments.
- `script_entrainment_phase.jl`: Generate phase angle curves.
- `script_prc.jl`: Generate a phase response curve.

## Miscellaneous
- `BCModel.jl`: Routines for easy manipulation with the models.
- `functions_A_optim.jl`: Support functions for `script_A_optim.jl`.
- `functions_B_optim.jl`: Support functions for `script_B_optim.jl`.
- `functions_C_optim.jl`: Support functions for `script_C_optim.jl`.
- `functions_D_optim.jl`: Support functions for `script_D_optim.jl`.
- `functions_model.jl`: Equations for the stochastic and deterministic models.
- `functions.jl`: Miscellaneous functions.
- `Peaks.jl`: Routines for peak detection.
