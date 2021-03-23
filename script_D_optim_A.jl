include("functions.jl")

MODEL_NAME = "model05"
RUN_NAME = "run 1"

## Control ====================================================================
params = Dict(
    :model_name => MODEL_NAME,
    :showplots => false,
    :plate_set => "A",
    :drug => "Control",
    :ode_alg => AutoTsit5(Rosenbrock23()),  # ODE solver
    :sde_alg => DRI1(),  # SDE solver
    :trajectories => 1000,
    :bbo_SearchRange => (0.0, 100.0),
    :bbo_MaxSteps => 25_000,
    :n_reval => 100,  # for fitness reevaluation
    :osc_con => true,  # constraint on sustained period
    
    # input files
    :in_initial_population_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "C_parameters.csv"),
    :in_initial_conditions_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_initial_conditions.csv"),
    :in_normalization_parameters_filename => "",  # empty string for Control
    
    # outout files
    :folder_progress => joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_control_a"),
    :out_results_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv"),
    :out_normalization_parameters_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_norm_parameters_a.csv")  # empty string for no saving
)
D_run_all!(deepcopy(params))

## Drugs ======================================================================
params[:osc_con] = false  # sustained period not forced for the drug treatments
params[:in_normalization_parameters_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_norm_parameters_a.csv")
params[:out_normalization_parameters_filename] = ""

# DBC -------------------------------------------------------------------------
params[:drug] = "DBC 0.5"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_1.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_dbc_1")
D_run_all!(deepcopy(params))

params[:drug] = "DBC 1"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_1.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_2.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_dbc_2")
D_run_all!(deepcopy(params))

params[:drug] = "DBC 3"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_2.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_3.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_dbc_3")
D_run_all!(deepcopy(params))

# U0126 -----------------------------------------------------------------------
params[:drug] = "U0126 10"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_1.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_U0126_1")
D_run_all!(deepcopy(params))

params[:drug] = "U0126 20"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_1.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_2.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_U0126_2")
D_run_all!(deepcopy(params))

params[:drug] = "U0126 40"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_2.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_3.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_U0126_3")
D_run_all!(deepcopy(params))

# Forskolin -------------------------------------------------------------------
params[:drug] = "Forskolin 5"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_1.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_forskolin_1")
D_run_all!(deepcopy(params))

params[:drug] = "Forskolin 10"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_1.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_2.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_forskolin_2")
D_run_all!(deepcopy(params))

params[:drug] = "Forskolin 15"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_2.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_3.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_forskolin_3")
D_run_all!(deepcopy(params))