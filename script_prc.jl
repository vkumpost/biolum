include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 2"

params = Dict(
    :model_name => MODEL_NAME,
    :alg => DRI1(),  # solver
    :trajectories => 30000,
    :pulse_length => 12,
    
    # input files
    :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "C_parameters.csv"),
    
    # outout files
    :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "PRC_C_parameters.csv")
)

# Get parameters
filename = params[:filename]
output_filename = params[:output_filename]
model_name = params[:model_name]
alg = params[:alg]
trajectories = params[:trajectories]
pulse_length = params[:pulse_length]

# Load parameters
df = loadcsv(filename)
p = DataFrame(df[1, Not(:fitness)])

# Load model
events = createevents([(:LD, 10, 12, 12)])
tspan = (0.0, events[end] + 9*24)
model = loadmodel(model_name * "sde")
model = addcallback(model, events, :I)
model = remakemodel(model, tspan=tspan, saveat=(tspan[1]+120.0):0.1:tspan[end])
model = setparameters(model, p)

@time tout, xout = simulatepopulation2(model, alg, trajectories)

Aref = ptp(xout)
peaks = findpeaks(xout, tout; minp = Aref * 0.1)

## ==
idx = peaks.locs .< events[end]
C0 = median(peaks.locs[idx] .- events[end-4:end, 1])

idx = peaks.locs .> events[end]
P1 = peaks.locs[idx][3]  # reference peak 1 (3rd peak after LD cycle)
P2 = peaks.locs[idx][4]  # reference peak 2 (4th peak after LD cycle)
FRP = P2 - P1  # free-running period
C1 = P1 - C0  # the beginning of the subjective daytime
C2 = P2 - C0  # the end of the subjective daytime
pulse_times = range(C1, C2; length = 10)  # 25
# window_start = C2+24  # window to estimate PCR  1000
window_start = C2+48  # window to estimate PCR

figure()
subplot(211)
plot(tout, xout; color = :black)
plot(peaks.locs, peaks.pks, "x"; color = :red)

subplot(212)
plot(tout, xout; color = :black)
plotevents(events)
vlines.(pulse_times, minimum(xout), maximum(xout))
vlines(C1, minimum(xout), maximum(xout), color = :red, zorder = 1000)
vlines(C2, minimum(xout), maximum(xout), color = :orange, zorder = 1000)
vlines(window_start, minimum(xout), maximum(xout), color = :green, zorder = 1000)

## ==
n = length(pulse_times)
fig, axarr = subplots(n+1)
for ax in axarr
    ax.axis("off")
end

window_signals = Vector{Vector{Float64}}(undef, n+1)
axarr[1].plot(tout, xout; color = :black)
plotevents(axarr[1], events)
axarr[1].vlines(window_start, minimum(xout), maximum(xout), color = :red)
# axarr[1].set_xlim(240.0, tspan[end])  # 1000
axarr[1].set_xlim(120.0, tspan[end])

window_time = tout[tout .>= window_start] .- tout[tout .>= window_start][1]
window_signals[1] = xout[tout .>= window_start]

## ==
for i = 1:n
    global prob, model

    pulse_time = pulse_times[i]
    pulse_events = vcat(events, [pulse_time pulse_time + pulse_length])
    model = addcallback(model, pulse_events, :I)
    @time tout, xout = simulatepopulation2(model, alg, trajectories)
    peaks = findpeaks(xout, tout; minp = 0.01)
    window_signals[i + 1] = xout[tout .>= window_start]

    axarr[i+1].plot(tout, xout; color = :black)
    plotevents(axarr[i+1], pulse_events)
    axarr[i+1].vlines(window_start, minimum(xout), maximum(xout); color = :red)
    # axarr[i+1].set_xlim(240.0, tspan[end])  # 1000
    axarr[i+1].set_xlim(120.0, tspan[end])  # 1000
end

## ==
close(:all)
peak_PCR = []
for i = 2:n+1
    t = window_time
    x1 = zscore(window_signals[1])
    x2 = zscore(window_signals[i])
    peaks1 = findpeaks(x1; minp = 0.1 * (maximum(x1) - minimum(x1)), showplot = false)
    ref_loc = peaks1.locs[2]  # 1
    peaks2 = findpeaks(x2; minp = 0.1 * (maximum(x2) - minimum(x2)), showplot = false)
    diff_locs = abs.(peaks1.locs[2] .- peaks2.locs)  # 1
    idx = argmin(diff_locs)
    shift_loc = peaks2.locs[idx]
    Δ = t[ref_loc] - t[shift_loc]
    push!(peak_PCR, Δ)
    
    figure()
    plot(t, x1; color = :black)
    plot(t[ref_loc], x1[ref_loc], "o"; color = :red)
    plot(t, x2)
    plot(t[shift_loc], x2[shift_loc], "o"; color = :orange)
    title("Δt = $(round(Δ, digits = 2)) hours")
end

circadian_time = round.(abs.(pulse_times .- C1) ./ FRP * 24.0, digits = 2)
circadian_shift = peak_PCR ./ FRP * 24.0

PRC = DataFrame(times = circadian_time, shift = circadian_shift)
savecsv(PRC, output_filename)