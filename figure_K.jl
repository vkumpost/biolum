include("functions.jl")

close(:all)

small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

MODEL_NAME = "model04b"
RUN_NAME = "run 2"

# Load model parameters
filename = joinpath("model04 $(RUN_NAME)", "C_parameters.csv")
df = loadcsv(filename)
p = DataFrame(df[1, 1:end-1])

fig, axs = subplots(3, figsize=(8, 4))
for (i, K) = enumerate([10^-4, 10^-8, 0.0])

    # Load model
    model = loadmodel(MODEL_NAME * "sde")
    model = setparameters(model, p)

    # Set up model parameters
    events = createevents([(:DD, 24), (:LD, 10, 12, 12)])
    tspan = (0.0, events[end] + 12.0 + 24.0*5.0)
    model = setparameters(model, [:K], [K])
    model = addcallback(model, events, :I)
    model = remakemodel(model, tspan=tspan, saveat=0.1)

    # Solve model
    @time tout, xout = simulatepopulation2(model, DRI1(), 30000)

    events_plot = events[end-5:end, :]
    offset = events_plot[1]
    events_plot .-= offset
    idx = tout .> offset
    tout_plot = tout[idx] .- offset
    xout_plot = xout[idx]

    axs[i].plot(tout_plot, xout_plot, color = :black)
    plotevents(axs[i], events_plot)
    if K == 0
        axs[i].set_title("K = 0")
    else
        axs[i].set_title("K = 1e$(round(Int, log10(K)))")
    end
    axs[i].set_ylabel("Lumin. (au)")
    axs[i].set_xlabel("Time (hours)")
    axs[i].set_xticks(0:24:tout_plot[end])
    axs[i].set_xlim(tout_plot[1], tout_plot[end])
    axs[i].set_ylim(minimum(xout_plot), maximum(xout_plot))

end

fig.tight_layout()

savefigure(fig, "fig_K.svg")