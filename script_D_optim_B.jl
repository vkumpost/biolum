include("functions.jl")

MODEL_NAME = "model05"
RUN_NAME = "run 1"

## Control ====================================================================
params = Dict(
    :model_name => MODEL_NAME,
    :showplots => false,
    :plate_set => "B",
    :drug => "Control",
    :ode_alg => AutoTsit5(Rosenbrock23()),  # ODE solver
    :sde_alg => DRI1(),  # SDE solver
    :trajectories => 1000,
    :bbo_SearchRange => (0.0, 100.0),
    :bbo_MaxSteps => 25_000,
    :n_reval => 100,  # for firtness reevaluation
    :osc_con => true,  # constraint on sustained period
    
    # input files
    :in_initial_population_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "C_parameters.csv"),
    :in_initial_conditions_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_initial_conditions.csv"),
    :in_normalization_parameters_filename => "",  # empty string for Control
    
    # outout files
    :folder_progress => joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_control_b"),
    :out_results_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv"),
    :out_normalization_parameters_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_norm_parameters_b.csv")  # empty string for no saving
)
D_run_all!(deepcopy(params))

## Drugs ======================================================================
params[:osc_con] = false  # sustained period not forced for the drug treatments
params[:in_normalization_parameters_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_norm_parameters_b.csv")
params[:out_normalization_parameters_filename] = ""

# EGF -------------------------------------------------------------------------
params[:drug] = "EGF 30"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_1.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_egf_1")
D_run_all!(deepcopy(params))

params[:drug] = "EGF 50"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_1.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_2.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_egf_2")
D_run_all!(deepcopy(params))

params[:drug] = "EGF 80"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_2.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_3.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_egf_3")
D_run_all!(deepcopy(params))

# PMA -------------------------------------------------------------------------
params[:drug] = "PMA 0.5"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_1.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_pma_1")
D_run_all!(deepcopy(params))

params[:drug] = "PMA 1"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_1.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_2.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_pma_2")
D_run_all!(deepcopy(params))

params[:drug] = "PMA 3"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_2.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_3.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_pma_3")
D_run_all!(deepcopy(params))

# Ro-318220 -------------------------------------------------------------------
params[:drug] = "Ro-318220 2"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_1.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_ro_1")
D_run_all!(deepcopy(params))

params[:drug] = "Ro-318220 5"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_1.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_2.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_ro_2")
D_run_all!(deepcopy(params))

params[:drug] = "Ro-318220 8"
params[:in_initial_population_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_2.csv")
params[:out_results_filename] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_3.csv")
params[:folder_progress] = joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_ro_3")
D_run_all!(deepcopy(params))