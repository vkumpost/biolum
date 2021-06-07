include("functions.jl")

MODEL_NAME = "model06"
RUN_NAME = "run 2 x"

close(:all)

params = Dict(

    :model_name => MODEL_NAME,
    :alg => DRI1(),  # solver
    :trajectories => 30000,
    :T_cycles => [8, 12, 16, 20, 24, 28, 32, 36, 40],
    
    # input/output files =====================================================

    # Control
    :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv"),
    :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_control_a.csv")

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_control_b.csv")

    # DBC
    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_1.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_forskolin_1.csv")

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_2.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_forskolin_2.csv")

    # DBC
    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_1.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_dbc_1.csv")

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_2.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_dbc_2.csv")

    # U0126
    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_1.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_u0126_1.csv"),

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_2.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_u0126_2.csv")

    # EGF
    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_1.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_egf_1.csv")

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_2.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_egf_2.csv")

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_3.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_egf_3.csv")

    # PMA
    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_2.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_pma_2.csv")

    # :filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_3.csv"),
    # :output_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "T_pma_3.csv")

)

# Get parameters
filename = params[:filename]
output_filename = params[:output_filename]
model_name = params[:model_name]
alg = params[:alg]
trajectories = params[:trajectories]
T_cycles = params[:T_cycles]

# Load parameters
df = loadcsv(filename)
p = DataFrame(df[1, Not(:fitness)])

# Load model
events = createevents([(:LD, 20, 12, 12)])
tspan = (0.0, events[end] + 24.0)
model = loadmodel(model_name * "sde")
model = addcallback(model, events, :I)
model = remakemodel(model, tspan=tspan, saveat=events[end-9, 1]:0.1:(events[end] + 12.0))
model = setparameters(model, p)

@time tout, xout = simulatepopulation2(model, alg, trajectories)

Aref = ptp(xout)
peaks = findpeaks(xout, tout; minp = Aref * 0.01)

n_T_cycles = length(T_cycles)

figure()
subplot(211)
plot(tout, xout; color = :black)
plot(peaks.locs, peaks.pks, "x"; color = :red)

subplot(212)
plot(tout, xout; color = :black)
plotevents(events)

## ==
fig, axarr = subplots(n_T_cycles+1)
for ax in axarr
    ax.axis("off")
end

Aref = maximum(xout) - minimum(xout)
peaks = findpeaks(xout, tout; minp = Aref * 0.01)

if length(peaks) == 10
    phase_angle = mean(events[end-9:end, 1] .- peaks.locs)
    println(phase_angle)
end

axarr[1].plot(tout, xout; color = :black)
plotevents(axarr[1], events)
axarr[1].set_xlim(240.0, tspan[end])
axarr[1].plot(peaks.locs, peaks.pks, "x"; color = :red)
phase_angles = []
## ==
i_T_cycle = 1
repeats = 0
while i_T_cycle <= n_T_cycles
    global prob, model, i_T_cycle, T_cycles, repeats

    T_cycle = T_cycles[i_T_cycle]
    println(T_cycle)
    n_transition_cycles = round(240.0 / T_cycle)
    println(n_transition_cycles)
    T_events = createevents([(:LD, n_transition_cycles + 10, T_cycle/2, T_cycle/2)])

    T_model = remakemodel(model,
            tspan = (0.0, T_events[end]),
            saveat = T_events[end-9, 1]-(T_cycle/2):0.1:(T_events[end])
        )
    T_model = addcallback(T_model, T_events, :I)
    @time tout, xout = simulatepopulation2(T_model, alg, trajectories)

    axarr[i_T_cycle+1].plot(tout, xout; color = :black)
    plotevents(axarr[i_T_cycle+1], T_events)
    axarr[i_T_cycle+1].set_xlim(tout[1], tout[end])
    axarr[i_T_cycle+1].set_ylabel(T_cycle)

    Aref = maximum(xout) - minimum(xout)
    peaks = findpeaks(xout, tout; minp = Aref * 0.1)
    axarr[i_T_cycle+1].plot(peaks.locs, peaks.pks, "x"; color = :red)

    if length(peaks) == 10
        phase_angle = median(T_events[end-9:end, 1] .- peaks.locs) / T_cycle * 24.0
        push!(phase_angles, phase_angle)
        println(phase_angle)
        i_T_cycle += 1
        repeats = 0
    elseif repeats < 5
        repeats += 1
        println("Repeat $repeats")
    else
        println("No success after $repeats repeats")
        push!(phase_angles, NaN)
        i_T_cycle += 1
        repeats = 0
    end

end

## ==

PRC = DataFrame(T_cycles = T_cycles, phase_angles = phase_angles)
savecsv(PRC, output_filename)