include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 2"
TRAJECTORIES = 30000

df = loadcsv("$(MODEL_NAME) $(RUN_NAME)", "C_parameters.csv")
p = DataFrame(df[1, 1:end-1])  # the best set of parameters without fitness

output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "C_validation.svg")
alg = DRI1()  # sde solver

## Prepare figure layout ======================================================
fig = figure(figsize=(7, 5), constrained_layout=true)  # 7.5, 8.75
gs = fig.add_gridspec(3, 2)
ax1 = fig.add_subplot(get(gs, (0, pyslice(0, 2))))
ax2 = fig.add_subplot(get(gs, (1, pyslice(0, 2))))
ax3 = fig.add_subplot(get(gs, (2, 0)))
ax4 = fig.add_subplot(get(gs, (2, 1)))

## Create model ===============================================================
model = loadmodel(MODEL_NAME * "sde")
model = setparameters(model, p)

## Precalculate initial conditions ============================================
function initial_conditions(model, alg; pacing=12)
    println("Calculating initial conditions wiht pacing $pacing...")
    events = createevents([(:DD, 24), (:LD, 10, pacing, pacing)])
    model = addcallback(model, events, :I)
    model = remakemodel(model, tspan=(0.0, events[end]), saveat=[])
    initial_conditions = Array{Float64}(undef, TRAJECTORIES, 3)
    for i = 1:TRAJECTORIES
        _, _, sol = solvemodel(model, alg; save_everystep=false)
        initial_conditions[i, :] = sol.u[2]
    end
    return initial_conditions
end
ic = initial_conditions(model, alg)

## Plate U1 ===================================================================
data = loaddata("Plate U1")
data = selectwells(data, r"WT [A|B]")
events = detectevents(data)
bounds = [events[4, 2], events[4, 2] + 300.0]
idx = bounds[1] .< data.Time .< bounds[2]
data = data[idx, :]
data.Time .-= bounds[1]
data = transformwells(data, zscore)
traces = Matrix(data[:, r"WT"])
data = aggregatewells(data, mean)
events = detectevents(data)
model = addcallback(model, events, :I)
model = remakemodel(model, tspan=(0.0, data.Time[end]), saveat=data.Time)
tout, xout = simulatepopulation2(model, alg, TRAJECTORIES; u0arr=ic)

Rs = Float64[]
Ls = Float64[]
println("Fitting data")
for i = 1:16
    traces[:, i] = zscore(traces[:, i])
    R = cor(traces[:, i], zscore(xout))^2
    R = rsquared(traces[:, i], zscore(xout))
    println("R = $R")
    push!(Rs, R)
end
rsquared_mean = round(mean(Rs), digits=2)
rsquared_std = round(std(Rs), digits=2)
println("R² = $rsquared_mean ± $rsquared_std")

h_data = ax1.plot(data.Time, traces, color=:gray)
h_model = ax1.plot(tout, zscore(xout), color=:black)
plotevents(ax1, events)
ax1.autoscale(enable=true, tight=true)
ax1.set_xlabel("Time (hours)", labelpad=0)
ax1.set_ylabel("Lumin. (au)", labelpad=0)
ax1.set_title("Fitting data", pad=0, loc="left")
ax1.set_xticks(0:24:maximum(tout))
ax1.legend([h_model[1], h_data[1]], ["Model", "Data"];
    fancybox=false, edgecolor=:black, framealpha=1.0)

## Plate U2 ===================================================================
data = loaddata("Plate U2")
data = selectwells(data, r"DAP")
events = detectevents(data)
bounds = [events[2, 2], events[2, 2] + 300.0]
idx = bounds[1] .< data.Time .< bounds[2]
data = data[idx, :]
data.Time .-= bounds[1]
data = transformwells(data, zscore)
traces = Matrix(data[:, r"DAP"])
events = detectevents(data)

model = addcallback(model, events, :I)
model = remakemodel(model, tspan=(0.0, data.Time[end]), saveat=data.Time)
tout, xout, _ = simulatepopulation(model, alg, TRAJECTORIES; u0arr=ic)

Rs = Float64[]
println("Plate U2")
for i = 1:4
    traces[:, i] = zscore(traces[:, i])
    R = cor(traces[:, i], zscore(xout))^2
    R = rsquared(traces[:, i], zscore(xout))
    println("R = $R")
    push!(Rs, R)
end
rsquared_mean = round(mean(Rs), digits=2)
rsquared_std = round(std(Rs), digits=2)
println("R² = $rsquared_mean ± $rsquared_std")

ax2.plot(data.Time, traces, color=:gray)
ax2.plot(tout, zscore(xout), color=:black)
plotevents(ax2, events)
ax2.autoscale(enable=true, tight=true)
ax2.set_xlabel("Time (hours)", labelpad=0)
ax2.set_ylabel("Lumin. (au)", labelpad=0)
ax2.set_title("Validation data", pad=0, loc="left")
ax2.set_xticks(0:24:maximum(tout))

## Plate U3 ===================================================================
data = loaddata("Plate U3")
data = selectwells(data, r"Pool")
events = detectevents(data)

bounds = [events[1, 2], events[1, 2] + 150.0]
idx = bounds[1] .< data.Time .< bounds[2]
data = data[idx, :]
data.Time .-= bounds[1]
data = transformwells(data, zscore)
traces = Matrix(data[:, r"Pool"])
events = detectevents(data)

model = addcallback(model, events, :I)
model = remakemodel(model, tspan=(0.0, data.Time[end]), saveat=data.Time)
ic = initial_conditions(model, alg; pacing=15)
tout, xout, _ = simulatepopulation(model, alg, TRAJECTORIES; u0arr=ic)

Rs = Float64[]
println("Plate U3")
for i = 1:size(traces, 2)
    traces[:, i] = zscore(traces[:, i])
    R = cor(traces[:, i], zscore(xout))^2
    R = rsquared(traces[:, i], zscore(xout))
    println("R = $R")
    push!(Rs, R)
end
rsquared_mean = round(mean(Rs), digits=2)
rsquared_std = round(std(Rs), digits=2)
println("R² = $rsquared_mean ± $rsquared_std")

ax3.plot(data.Time, traces, color=:gray)
ax3.plot(tout, zscore(xout), color=:black)
plotevents(ax3, events)
ax3.autoscale(enable=true, tight=true)
ax3.set_xlabel("Time (hours)", labelpad=0)
ax3.set_ylabel("Lumin. (au)", labelpad=0)
ax3.set_title("15:15 LD cycle", pad=0, loc="left")
ax3.set_xticks(0:24:maximum(tout))

## Plate U4 ===================================================================
data = loaddata("Plate U4")
data = selectwells(data, r"Pool")
events = detectevents(data)

bounds = [events[1, 2], events[1, 2] + 150.0]
idx = bounds[1] .< data.Time .< bounds[2]
data = data[idx, :]
data.Time .-= bounds[1]
data = transformwells(data, zscore)
traces = Matrix(data[:, r"Pool"])
events = detectevents(data)

model = addcallback(model, events, :I)
model = remakemodel(model, tspan=(0.0, data.Time[end]), saveat=data.Time)
ic = initial_conditions(model, alg; pacing=10)
tout, xout, _ = simulatepopulation(model, alg, TRAJECTORIES; u0arr=ic)

Rs = Float64[]
println("Plate U4")
for i = 1:size(traces, 2)
    traces[:, i] = zscore(traces[:, i])
    R = cor(traces[:, i], zscore(xout))^2
    R = rsquared(traces[:, i], zscore(xout))
    println("R = $R")
    push!(Rs, R)
end
rsquared_mean = round(mean(Rs), digits=2)
rsquared_std = round(std(Rs), digits=2)
println("R² = $rsquared_mean ± $rsquared_std")

ax4.plot(data.Time, traces, color=:gray)
ax4.plot(tout, zscore(xout), color=:black)
plotevents(ax4, events)
ax4.autoscale(enable=true, tight=true)
ax4.set_xlabel("Time (hours)", labelpad=0)
ax4.set_ylabel("Lumin. (au)", labelpad=0)
ax4.set_title("10:10 LD cycle", pad=0, loc="left")
ax4.set_xticks(0:24:maximum(tout))

## Save figure ================================================================
savefigure(fig, output_filename)
