include("functions.jl")

model_name = "model04"
run_name = "run 1"
n_cells = 5

in_filename = joinpath("$(model_name) $(run_name)", "B_parameters.csv")
out_filename = joinpath("population_desynchronization.svg")

df = loadcsv(in_filename)
p = DataFrame(df[13, :])

model = loadmodel(model_name * "sde")
model = setparameters(model, p)
printparameters(model)

events = createevents([(:DD, 24), (:LD, 10, 12, 12)])
model = addcallback(model, events, :I)
model = remakemodel(model, tspan = (0.0, events[end] + 4 * 24.0),
    saveat = events[end-2, 1]:0.1:events[end] + 4 * 24.0)

alg = DRI1()
@time (tout, xout, X, sim) = simulatepopulation(model, alg, 1000)

## Figure =====================================================================
fig, (ax2, ax1) = subplots(1, 2; figsize = (7, 2))

# Mean of traces --------------------------------------------------------------
ax1.plot(tout, xout; color = :black)
max_x = maximum(xout)
min_x = minimum(xout)
plotevents(ax1, events; ylims=[0.0, max_x + min_x])

ax1.set_xlim(tout[1], tout[end])
ax1.set_xticks((tout[1]:48.0:tout[end]))
ax1.set_xticklabels(convert.(Int, (tout[1]:48.0:tout[end]) .- tout[1]))
ax1.set_xlabel("Time (hours)", labelpad=0)

ax1.set_yticks(0.0:0.1:max_x + min_x)
ax1.set_ylim(0.0, max_x + min_x)
ax1.set_ylabel("Lumin. (au)", labelpad=0)
ax1.set_title("Mean of traces", pad=0, loc="left")

# Individual traces -----------------------------------------------------------
for i = 1:n_cells
    ax2.plot(tout, X[:, i]; color = :black)
end
max_x = maximum(X[:, 1:n_cells])

plotevents(ax2, events; ylims=[0.0, max_x])

ax2.set_xticks(tout[1]:48.0:tout[end])
ax2.set_xticklabels(convert.(Int, (tout[1]:48.0:tout[end]) .- tout[1]))
ax2.set_xlabel("Time (hours)", labelpad=0)

ax2.set_yticks(0.0:0.1:max_x)
ax2.set_ylim(0.0, max_x)
ax2.set_ylabel("Lumin. (au)", labelpad=0)

ax2.set_title("Example of 5 single-cell traces", pad=0, loc="left")
ax2.set_xlim(tout[1], tout[end])

fig.tight_layout()

## Save figure ================================================================
savefigure(fig, out_filename)
