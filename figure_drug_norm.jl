include("functions.jl")

druglists = [
    ["Control", "Forskolin 5", "Forskolin 10", "Forskolin 15"],
    ["Control", "DBC 0.5", "DBC 1", "DBC 3"],
    ["Control", "U0126 10", "U0126 20", "U0126 40"],
    ["Control", "EGF 30", "EGF 50", "EGF 80"],
    ["Control", "PMA 0.5", "PMA 1", "PMA 3"],
    ["Control", "Ro-318220 2", "Ro-318220 5", "Ro-318220 8"],
]

platelists = [
    ["Plate D1 A", "Plate D2 A"],
    ["Plate D1 A", "Plate D2 A"],
    ["Plate D1 A", "Plate D2 A"],
    ["Plate D1 B", "Plate D2 B"],
    ["Plate D1 B", "Plate D2 B"],
    ["Plate D1 B", "Plate D2 B"],
]

colors = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

for (druglist, platelist) = zip(druglists, platelists)

    nplate = length(platelist)

    # Load data
    datarawarr = [DataFrame(), DataFrame()]
    datanormarr = [DataFrame(), DataFrame()]
    for i = 1:nplate
        reg = Regex(join(druglist, "|"))
        
        filename = "raw_$(lowercase(replace(platelist[i], " " => "_"))).csv"
        dataraw = loadcsv(filename)
        datarawarr[i] = selectwells(dataraw, reg)

        filename = "norm_$(lowercase(replace(platelist[i], " " => "_"))).csv"
        datanorm = loadcsv(filename)
        datanormarr[i] = selectwells(datanorm, reg)
    end

    drugname = split(druglist[2], " ")[1]

    # Create figure raw
    for (normname, dataarr) in zip(["original", "normalized", "naive"], [datarawarr, datanormarr, datarawarr])

        fig = figure(figsize=(6.5, 3.2), constrained_layout=true)
        gs = fig.add_gridspec(2, 2; width_ratios=[80, 1])
        ax1 = fig.add_subplot(get(gs, (0, pyslice(0, 2))))
        ax2 = fig.add_subplot(get(gs, (1, 0)))
        axarr = [ax1, ax2]

        ndrug = length(druglist)
        for i = 1:nplate
            platename = platelist[i]

            handles = []
            for j = 1:ndrug
                X = wells2matrix(selectwells(dataarr[i], Regex(druglist[j])))
                if normname == "naive"
                    for i_column = 1:size(X, 2)
                        X[:, i_column] = zscore(X[:, i_column])
                    end
                end
                h = axarr[i].plot(dataarr[i].Time, X, color=colors[j])
                push!(handles, h[1])
            end
            plotevents(axarr[i], detectevents(dataarr[i]))
            if normname == "original"
                axarr[i].set_ylabel("Lumin. (cps)", labelpad=0)
            else
                axarr[i].set_ylabel("Lumin. (au)", labelpad=0)
            end
            axarr[i].set_xlabel("Time (hours)", labelpad=0)
            if i == 1
                axarr[i].set_title("$drugname - jet lag - $normname", loc="left", pad=0)
            else
                axarr[i].set_title("$drugname - constant darkness - $normname", loc="left", pad=0)
            end
            axarr[i].set_xticks(0:24:maximum(dataarr[i].Time))
            if i == 2
                axarr[i].legend(handles, translatelabel.(druglist),
                    loc="upper left", bbox_to_anchor=(1, 1.0), fancybox=false,
                    edgecolor=:black)
            end
            axarr[i].autoscale(tight=true)
        end

        # Safe figure
        filename = lowercase("$(drugname)_$(normname).svg")
        savefigure(fig, "drug_plots", filename)
        
    end
end