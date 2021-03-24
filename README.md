# Biolum

This repository contains all code and data used in our publication.

## Installation 
Julia can be downloaded at https://julialang.org/. To install necessery dependencies, start Julia, navigate yourself into the project folder by running `cd("my/path/to/biolum")`, and then either instantiate a package enviroment by running `] activate .` and `] instantiate` or just run `include("install.jl")` which will run `install.jl` script that will install all needed packages in the global enviroment.

This project was written in Julia for VSCode (https://www.julia-vscode.org/). Other popular editors include an Atom extension Juno (https://junolab.org/) and Jupyter notebooks (https://jupyter.org/).

## Luminescence recordings
Luminescence recordings are stored in `data/biolumdata.xlsx`. The file contains 8 sheets. Each sheet contains traces obtained from one cell culture plate. The first column on each sheet is called `Time` and contains time in hours from the beginning of the recording. The second column is called `Light` and decodes whether the light was turned on (1) or off (0) at the given time. The remaining columns contain luminescence recordings in counts per second (CPS). The plate names correspond to the publication followingly
- `Plate U1`: fitting data.
- `Plate U2`: validation data.
- `Plate U3`: 15:15 LD cycle.
- `Plate U4`: 10:10 LD cycle.
- `Plate D1 A`: LD cycle (set A).
- `Plate D2 A`: Constant darkness (set A).
- `Plate D1 B`: LD cycle (set B).
- `Plate D2 B`: Constant darkness (set B).

## Figures
All files strating with `figure_` contain scripts to generate figures from the paper. The generated figures are automatically stored in folder `figures` as SVG files.
- `figure_C_validation.jl`: Fig 2.
- `figure_data_norm.jl`: Fig S4.
- `figure_drug_fit.jl`: Fig3A, Fig3B, Fig S3. This script also generates the values for Table S2.
- `figure_drug_norm.jl`: Fig S5 and Fig S6.
- `figure_fig4.jl`: Fig 4.
- `figure_pacing_protocols.jl`: Fig S1.
- `figure_pca.jl`: Fig 3C and Fig S2.
- `figure_population_desynchronization.jl`: Fig 1C.

## Scripts
Scripts to generate results presented in the paper. Output files are stored as CSV files in folder `outputs`.
- `script_A_optim.jl`: The first step of the optimization routine. Fitting properties of a deterministic model.
- `script_B_optim.jl`: The second step of the optimization routine. Fitting the noise intensity based on the damping ratio in constant dark.
- `script_C_optim.jl`: The third step of the optimization routine. Fitting data from untreated cells.
- `script_D_initial_conditions.jl`: Generate initial conditions for the drug data.
- `script_D_optim_A.jl`: Fit traces from drug plate set A.
- `script_D_optim_B.jl`: Fit traces from drug plate set B.
- `script_damping_ratio.jl`: Fit damped sine to the data.
- `script_drug_norm.jl`: Normalize drug data.
