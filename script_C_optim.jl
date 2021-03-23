include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 1"

params = Dict(
    :model_name => MODEL_NAME,
    :showplots => false,  # for testing

    :ode_alg => AutoTsit5(Rosenbrock23()),
    :sde_alg => DRI1(),
    :trajectories => 1000,

    :SearchRange => (0.0, 100.0),  # parameter range
    :MaxSteps => 100_000,  # number of optimization steps

    :n_reval => 100,  # for firtness reevaluation

    :in_initial_population_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)",
        "B_parameters.csv"),
    :out_final_population_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)",
        "C_parameters.csv"),
    :folder_progress => joinpath("$(MODEL_NAME) $(RUN_NAME)", "progress_c_optim")
)

## ==
C_load_initial_population!(params)
C_load_fitting_data!(params)
C_add_ode_model!(params)
C_add_sde_model!(params)
C_generate_cost_function!(params)
C_optimize!(params)
C_reval_fitness!(params)
C_save_results(params)