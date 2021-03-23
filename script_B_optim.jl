include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 1"

params = Dict{Symbol, Any}(
    :model_name => MODEL_NAME,
    :showplots => false,  # for testing
    :ndays => 5,  # how many days to use to calculate damping ratio
    :alg => DRI1(),  # SDE solver
    :trajectories => 1000,  # number of SDE model evaluations
    :noise_range => (0.0, 0.1),  # range for noise values
    :PopulationSize => 10,  # init pop. size for the evol. alg.
    :MaxFuncEvals => 1000,  # max. number of f. eval. during optimization
    :i => 1,  # iterator counter for B_optimize!
    :in_A => joinpath("$(MODEL_NAME) $(RUN_NAME)", "A_parameters.csv"),
    :out_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "B_parameters.csv")
)

## ==
B_load_damping_ratio!(params)
B_load_A_results!(params)
# B_create_model!(params)  # for testing
# B_create_cost_function!(params)  # for testing
# B_optimize!(params)  # for testing
B_optimize_all!(params)
B_save_results(params)