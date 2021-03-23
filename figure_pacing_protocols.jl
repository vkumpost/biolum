include("functions.jl")

## Plate names ================================================================
platenames = ["Plate U1", "Plate U2", "Plate U3", "Plate U4", "Plate D1 A",
    "Plate D2 A", "Plate D1 B", "Plate D2 B"]

## Create figure layout =======================================================
small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

fig = figure(figsize=(6, 5), constrained_layout=true)
gs = fig.add_gridspec(5, 2)
ax1 = fig.add_subplot(get(gs, (0, pyslice(0, 2))))
ax2 = fig.add_subplot(get(gs, (1, pyslice(0, 2))))
ax3 = fig.add_subplot(get(gs, (2, 0)))
ax4 = fig.add_subplot(get(gs, (2, 1)))
ax5 = fig.add_subplot(get(gs, (3, 0)))
ax6 = fig.add_subplot(get(gs, (3, 1)))
ax7 = fig.add_subplot(get(gs, (4, 0)))
ax8 = fig.add_subplot(get(gs, (4, 1)))
axarr = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

## Plot plates ================================================================
for (platename, ax) in zip(platenames, axarr)

    plate = loaddata(platename)
    events_all = detectevents(plate)

    if platename == "Plate U1"
        offset = events_all[4, 2]
        events = events_all[5:8, :] .- offset
        t = [0.0; sort(vec(events))]
        title_str = "Fitting data"
    elseif platename == "Plate U2"
        offset = events_all[2, 2]
        events = events_all[3:9, :] .- offset
        t = [0.0; sort(vec(events))]
        title_str = "Validation data"
    elseif platename in ["Plate U3", "Plate U4"]
        offset = events_all[1, 2]
        events = events_all[2:3, :] .- offset
        tend = events_all[1, 2] + 130.0
        t = [0.0; sort(vec(events)); tend]
        if platename == "Plate U3"
            title_str = "15:15 LD cycle"
        elseif platename == "Plate U4"
            title_str = "10:10 LD cycle"
        end
    elseif platename in ["Plate D1 A", "Plate D1 B"]
        offset = events_all[1, 2]
        events = events_all[2:4, :] .- offset
        tend = events_all[4, 2] - offset + 5
        t = [0.0; sort(vec(events)); tend]
        if platename == "Plate D1 A"
            title_str = "LD cycle (set A)"
        elseif platename == "Plate D1 B"
            title_str = "LD cycle (set B)"
        end
    elseif platename in ["Plate D2 A", "Plate D2 B"]
        offset = events_all[1, 2]
        tend = events_all[2, 1] - offset
        t = [0.0, tend]
        if platename == "Plate D2 A"
            title_str = "Constant darkness (set A)"
        elseif platename == "Plate D2 B"
            title_str = "Constant darkness (set B)"
        end
    end

    x = fill(0.0, length(t))
    for i = 2:2:length(t) - 1
        x[i] = 0
        x[i + 1] = 1
    end

    ax.step(t, x; color = :black)
    ax.set_yticks([0.0, 1.0])
    ax.set_yticklabels(["Dark", "Light"])
    ax.set_xticks(0:24:t[end])
    ax.set_xlabel("Time (hours)")
    ax.set_title(title_str)
    ax.autoscale(tight = true, axis = "x")
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    
end

## Save figure ================================================================
savefigure(fig, "pacing_protocols.svg")