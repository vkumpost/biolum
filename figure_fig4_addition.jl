include("functions.jl")

close(:all)

small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

params = Dict(
    :model_name => "model06",
    :dir => "model06 run 2 x",
    :plate_sets => [
        "A",
        "A",
        "A"
    ],
    :drugs => [
        "Control",
        "DBC 0.5",
        "U0126 20",
    ],
    :drug_filenames => [
        "D_control_a.csv",
        "D_dbc_1.csv",
        "D_u0126_2.csv"
    ],
    :titles => [
        "Control A", 
        "$(translatelabel("DBC 0.5"))",
        "$(translatelabel("U0126 20"))",
    ],
    :colors => [:black, CB_COLORS[2], CB_COLORS[3]],
    :savefigs => true,
    :alg => DRI1(),  # solver
    :trajectories => 30000
)

## load data ==================================================================
function load_data!(params)
    dir = params[:dir]
    drugs = params[:drugs]
    drug_filenames = params[:drug_filenames]
    ndrug = length(drugs)

    dfparr = [
        loadcsv(dir, drug_filenames[x]) for x in 1:ndrug
    ]

    params[:dfparr] = dfparr
end

load_data!(params)

## Setup a model ==============================================================
function create_models!(params)

    model_name = params[:model_name]
    dfparr = params[:dfparr]

    models = BCModel[]
    events = nothing
    for i = 1:length(dfparr)
        
        p = DataFrame(dfparr[i][1, Not("fitness")])
        events = createevents([(:DD, 24), (:LD, 10, 12, 12)])
        tspan = (0.0, events[end] + 12 + 24*5)
        model = loadmodel(model_name*"sde")
        model = setparameters(model, p)  
        model = addcallback(model, events, :I)
        model = remakemodel(model, tspan = tspan, saveat = 1.0)
        push!(models, model)
        
    end

    params[:models] = models
    params[:events] = events

end

create_models!(params)


## Run all simulations =======================================================
function run_all_simulations!(params)

    models = params[:models]
    trajectories = params[:trajectories]
    alg = params[:alg]
      
    toutarr = []
    xoutarr = []
    Xarr = []
    simarr = []
    
    n_models = length(models)

    for i = 1:n_models
        
        println("Model $(i)")

        @time tout, xout, X, sim = simulatepopulation(models[i], alg, trajectories)

        push!(toutarr, tout)
        push!(xoutarr, xout)
        push!(Xarr, X)
        push!(simarr, sim)

    end

    params[:toutarr] = toutarr
    params[:xoutarr] = xoutarr
    params[:Xarr] = Xarr
    params[:simarr] = simarr

end

run_all_simulations!(params)

## Show results ==============================================================

colors = params[:colors]
events = params[:events][end-4:end, :]
labels = params[:titles]

fig, axarr = subplots(4, figsize=(8, 4.5))

# first subfigure
max_value = -Inf
min_value = Inf
for i = 1:3
    global max_value, min_value
    
    xout = params[:xoutarr][i]
    tout = params[:toutarr][i]
    X = params[:Xarr][i]

    idx = events[1] .<= tout .<= events[end] + 12 + 24*5
    offset = tout[idx][1]
    tout = tout[idx] .- offset
    xout = xout[idx]

    if maximum(xout) > max_value
        max_value = maximum(xout)
    end

    if minimum(xout) < min_value
        min_value = minimum(xout)
    end

    axarr[1].plot(tout, xout, color = colors[i])

    if i == 3
        plotevents(axarr[1], events .- offset; ylims = [min_value-1, max_value+1])
        # axarr[1].set_title("Simulated time series")
        axarr[1].set_ylabel("Lumin. (au)")
        axarr[1].set_xlabel("Time (hours)")
        axarr[1].set_xlim(tout[1], tout[end])
        axarr[1].set_xticks(0:24:tout[end])
        axarr[1].set_ylim(min_value-0.001, max_value+0.001)
        axarr[1].legend(labels, ncol = 3, edgecolor=:black, framealpha=1.0, loc="upper right")
        axarr[1].set_title("Population-level simulated luminescence", loc = "left", y = 0.0, x = 0.0)
    end

end

# other subfigures (individual cells)
for i = 1:3
    
    xout = params[:xoutarr][i]
    tout = params[:toutarr][i]
    X = params[:Xarr][i]

    idx = events[1] .<= tout .<= events[end] + 12 + 24*5
    offset = tout[idx][1]
    tout = tout[idx] .- offset
    xout = xout[idx]

    max_value = 0.75
    min_value = 0.0
    # if maximum(xout) > max_value
    #     max_value = maximum(xout)
    # end

    # if minimum(xout) < min_value
    #     min_value = minimum(xout)
    # end

    for j = 1:10
        axarr[i+1].plot(tout, X[idx, j], color = colors[i], alpha = 1.0)
    end
    # axarr[1].plot(tout, xout, color = colors[i])

  
    plotevents(axarr[i+1], events .- offset; ylims = [min_value-1, max_value+1])
    # axarr[i+1].set_title("Simulated time series")
    axarr[i+1].set_ylabel("Lumin. (au)")
    axarr[i+1].set_xlabel("Time (hours)")
    axarr[i+1].set_xlim(tout[1], tout[end])
    axarr[i+1].set_xticks(0:24:tout[end])
    axarr[i+1].set_ylim(min_value-0.001, max_value+0.001)
    axarr[i+1].set_title("Individual simulated cells: $(labels[i])", loc = "left", y = 0.0, x = 0.0)



end

fig.tight_layout()
savefigure(fig, "fig4_simulations.svg")

## ==
#
locs_LD_arr = []
locs_DD_arr = []
for i = 1:3
    xout = params[:xoutarr][i]
    tout = params[:toutarr][i]
    X = params[:Xarr][i]
    events = params[:events]
    trajectories = params[:trajectories]

    idx_LD = events[end-1, 1] .<= tout .<= events[end, 1]
    idx_DD = events[end] + 12 + 2 .<= tout .<= events[end] + 12 + 24 + 2

    t_LD = tout[idx_LD]
    t_DD = tout[idx_DD]

    locs_LD = []
    locs_DD = []

    for j = 1:trajectories
        peaks = findpeaks(X[:, j], tout, minp = 0.01, sortstr = "descend")
        
        locs = peaks.locs
        idx = t_LD[1] .<= peaks.locs .< t_LD[end]
        if sum(idx) > 0
            imax = argmax(peaks.pks[idx])
            push!(locs_LD, peaks.locs[idx][imax])
        end
        
        locs = peaks.locs
        idx = t_DD[1] .<= peaks.locs .< t_DD[end]
        if sum(idx) > 0
            imax = argmax(peaks.pks[idx])
            push!(locs_DD, peaks.locs[idx][imax])
        end
        
    end

    push!(locs_LD_arr, locs_LD)
    push!(locs_DD_arr, locs_DD)

end

fig, ax = subplots(1, 2, figsize=(8, 2))
for i = 1:3
    color = colors[i]
    locs_LD = locs_LD_arr[i] .- minimum(locs_LD_arr[i])
    locs_DD = locs_DD_arr[i] .- minimum(locs_DD_arr[i])
    if i == 1
        alpha = 1.0
    else
        alpha = 1.0
    end
    h1 = ax[1].hist(locs_LD, bins = minimum(locs_LD):(maximum(locs_LD)+1), color = color, alpha = alpha, density=true, histtype=:step, linewidth=2)
    ax[2].hist(locs_DD, bins = minimum(locs_DD):(maximum(locs_DD)+1), color = color, alpha = alpha, density=true, histtype=:step, linewidth=2)

    n = length(h1[1])
    ax[1].hlines(1/n, 0, 24, color = :gray, linewidth = 2, zorder = 1000)
    ax[2].hlines(1/n, 0, 24, color = :gray, linewidth = 2, zorder = 1000)

    ax[1].set_xlim(0, 24)
    ax[2].set_xlim(0, 24)
    ax[1].set_title("LD cycle")
    ax[2].set_title("Constant darkness")
    ax[1].set_xlabel("Time (hours)")
    ax[2].set_xlabel("Time (hours)")
    ax[1].set_ylabel("Peak pdf (-)")
    ax[2].set_ylabel("Peak pdf (-)")
    ax[1].set_xticks(0:4:24)
    ax[2].set_xticks(0:4:24)
    ax[1].set_ylim(0, 0.08)
    ax[2].set_ylim(0, 0.08)
end
ax[2].legend([labels...; "Uniform distribution"], edgecolor=:black, framealpha=1.0, loc = "lower left", ncol=2)
fig.tight_layout()
savefigure(fig, "fig4_pdf.svg")