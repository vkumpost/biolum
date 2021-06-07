include("functions.jl")

close(:all)

small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

MODEL_NAME = "model04"
RUN_NAME = "run 2"

# Load model parameters
filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "C_parameters.csv")
df = loadcsv(filename)
p = DataFrame(df[1, 1:end-1])

fig, axs = subplots(4, figsize=(8, 5))
for (i, σ) = enumerate([0.002, 0.02, 0.04, 0.08])

    # Load model
    model = loadmodel(MODEL_NAME * "sde")
    model = setparameters(model, p)

    # Set up model parameters
    events = createevents([(:DD, 24), (:LD, 10, 12, 12)])
    tspan = (0.0, events[end] + 12.0 + 24.0*5.0)
    model = setparameters(model, [:noise], [σ])
    model = addcallback(model, events, :I)
    model = remakemodel(model, tspan=tspan, saveat=0.1)

    # Solve model
    @time tout, xout = simulatepopulation2(model, DRI1(), 30000)
    @time _, _, X, _ = simulatepopulation(model, DRI1(), 100)

    events_plot = events[end-5:end, :]
    offset = events_plot[1]
    events_plot .-= offset
    idx = tout .> offset
    tout_plot = tout[idx] .- offset
    xout_plot = xout[idx]
    X_plot = X[idx, :]

    h2 = nothing
    for j = 1:10
        h2 = axs[i].plot(tout_plot, X_plot[:, j], color = :gray, alpha = 0.5)
    end
    h1 = axs[i].plot(tout_plot, xout_plot, color = :black)
    plotevents(axs[i], events_plot; ylims = [0.0, 0.5])
    axs[i].text(0, 0.52, "σ = $σ", fontsize = 10)
    axs[i].set_ylabel("Lumin. (au)")
    axs[i].set_xlabel("Time (hours)")
    axs[i].set_xticks(0:24:tout_plot[end])
    axs[i].set_xlim(tout_plot[1], tout_plot[end])
    axs[i].set_ylim(0.0, 0.5)
    if i == 1
       axs[i].legend([h1[1], h2[1]], ["Population-level mean", "Individual cells"], ncol=2, fancybox=true, edgecolor=:black, framealpha=1.0)
    end

end

fig.tight_layout()

savefigure(fig, "fig_noise.svg")