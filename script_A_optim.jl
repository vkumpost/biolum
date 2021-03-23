include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 1"

params = Dict(
    :model_name => MODEL_NAME,  # model name
    :showplots => false,  # for testing
    :alg => AutoTsit5(Rosenbrock23()),  # ODE solver
    :entrainment_days => 3,  # number of days allowed for entrainment
    :MaxFuncEvals => 300_000,  # max number of f evals during optimization
    :TraceMode => :compact,  # compact, silent
    :num_iter => 50,  # how many parameter sets to estimate
    
    # search range for parameters ---------------------------------------------
    :SearchRange => [(0.0, 100.0), (0.0, 1.0), (0.0, 1.0)],  # model04

    # output file
    :out_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "A_parameters.csv")
)

## ==
A_add_pacing_regimes!(params)
A_create_models!(params)
A_create_cost_function!(params)
# A_optimize!(params)  # for testing
A_find_sets!(params)
A_save_results(params)