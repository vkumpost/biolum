include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 2"

alg = DRI1()
input_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "C_parameters.csv")
output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_initial_conditions.csv")
trajectories = 30000

## load parameter values and model function ===================================
df = loadcsv(input_filename)
df_bestfitness = df[1, Not("fitness")]
pnames = [Symbol(x) for x in names(df_bestfitness)]
pvalues = Vector(df_bestfitness)

transientevents = createevents([(:DD, 24), (:LD, 10, 12, 12)])
model = loadmodel("$(MODEL_NAME)sde")
model = setparameters(model, pnames, pvalues)
model = addcallback(model, transientevents, :I)
model = remakemodel(model, tspan=(0.0, transientevents[end]))
printparameters(model)

# test simulation -------------------------------------------------------------
t, x, _ = simulatepopulation(model, alg, 10, saveat=0.1)
figure()
subplot(211)
plot(t, x)
plotevents(transientevents)
subplot(212)
for k = 1:5
    t, x, _ = solvemodel(model, alg, saveat=0.1)
    plot(t, x)
end
plotevents(transientevents)

## generate intial conditions =================================================
icarr = Array{Float64}(undef, trajectories, 3)  # store all initial conditions here
for i = 1:trajectories
    println(i)
    _, _, sol = solvemodel(model, alg; save_everystep=false)
    icarr[i, :] = sol.u[2]
end

## convert the initial conditions into DataFrame and save as csv ==============
df_icarr = DataFrame(icarr, [model.species...])
savecsv(df_icarr, output_filename)
