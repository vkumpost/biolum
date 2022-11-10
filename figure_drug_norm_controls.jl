include("functions.jl")

colors = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

# Load data
data_dict = Dict{String, Dict}()

plate_names = ["norm_plate_d1_a", "norm_plate_d1_b", "norm_plate_d2_a",
    "norm_plate_d2_b", "raw_plate_d1_a", "raw_plate_d1_b", "raw_plate_d2_a",
    "raw_plate_d2_b"]

for plate_name in plate_names
    df = loadcsv(plate_name * ".csv")
    data_dict[plate_name] = Dict(
        "time" => df[:, 1] .- df[1, 1],
        "events" => detectevents(df) .- df[1, 1],
        "mean" => mean(Matrix(df[:, end-3:end]), dims=2)
    )
end

# Create figure raw
for norm_name in ["raw", "norm"]

    fig = figure(figsize=(6.5, 3.2), constrained_layout=true)
    gs = fig.add_gridspec(2, 2; width_ratios=[80, 1])
    ax1 = fig.add_subplot(get(gs, (0, pyslice(0, 2))))
    ax2 = fig.add_subplot(get(gs, (1, 0)))
    axarr = [ax1, ax2]

    min_x = Inf
    max_x = -Inf

    for (i, control_name) = enumerate(["a", "b"])    

        for (k, plate_number) = enumerate(["d1", "d2"])
            t = data_dict["$(norm_name)_plate_$(plate_number)_$(control_name)"]["time"]
            x = data_dict["$(norm_name)_plate_$(plate_number)_$(control_name)"]["mean"]
            events = data_dict["$(norm_name)_plate_$(plate_number)_$(control_name)"]["events"]
            
            min_x = min(min_x, minimum(x))
            max_x = max(max_x, maximum(x))
            
            if i == 1
                axarr[k].plot(t, x, color="black", label="Control A", linewidth=2)
            elseif i == 2
                axarr[k].plot(t, x, color="gray", label="Control B", linewidth=2)
                plotevents(axarr[k], events, ylims=[min_x-0.1, max_x+0.1])
            end
            if norm_name == "raw"
                axarr[k].set_ylabel("Lumin. (cps)", labelpad=0)
            else
                axarr[k].set_ylabel("Lumin. (au)", labelpad=0)
            end
            axarr[k].set_xlabel("Time (hours)", labelpad=0)
            if i == k
                axarr[k].set_title("Control - constant darkness - $norm_name", loc="left", pad=0)
            else
                axarr[k].set_title("Control - jet lag - $norm_name", loc="left", pad=0)
            end


            axarr[k].set_xticks(0:24:maximum(t))
            if k == 2
                axarr[k].legend(
                         loc="upper left", bbox_to_anchor=(1, 1.0), fancybox=false,
                         edgecolor=:black)
            end
        
            axarr[k].autoscale(tight=true)
            axarr[k].set_ylim([min_x-0.1, max_x+0.1])

        end        

    end

    # Safe figure
    filename = lowercase("controls_$(norm_name).svg")
    savefigure(fig, "drug_plots", filename)
    
end
