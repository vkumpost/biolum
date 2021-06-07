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
    :n_population => 50,
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
    :savefigs => true
)

## load data ==================================================================
function load_data!(params)
    model_name = params[:model_name]
    dir = params[:dir]
    drugs = params[:drugs]
    drug_filenames = params[:drug_filenames]
    ndrug = length(drugs)

    dfparr = [
        loadcsv(dir, drug_filenames[x]) for x in 1:ndrug
    ]

    u0arr = loadcsv(dir, "D_initial_conditions.csv")

    params[:u0arr] = Matrix(u0arr)
    params[:dfparr] = dfparr
end

load_data!(params)


## Setup a model ==============================================================
function create_models!(params)

    model_name = params[:model_name]
    dfparr = params[:dfparr]

    models = Vector{Vector{BCModel}}(undef, length(dfparr))
    for i = 1:length(dfparr)
        models[i] = BCModel[]
        for j = 1:51
            p = DataFrame(dfparr[i][j, Not("fitness")])
            events = createevents([(:DD, 1)])
            tspan = (0.0, 72.0)
            model = loadmodel(model_name*"sde")
            model = setparameters(model, p)  
            model = addcallback(model, events, :I)
            model = remakemodel(model, tspan=tspan,
                saveat=1.0)
            push!(models[i], model)
        end
    end

    params[:models] = models

end

create_models!(params)

## Run all simulations =======================================================
function run_all_simulations!(params)

    models = params[:models]
    dir = params[:dir]
    n_population = params[:n_population]
    drugs = params[:drugs]
    u0arr = params[:u0arr]

    solutions = []
    times = []
    Xarr = []
    simarr = []
    
    n_drugs = length(drugs)

    noisearr = Vector{Vector{Float64}}(undef, n_drugs)
    Iarr = Vector{Vector{Float64}}(undef, n_drugs)
    Aarr = Vector{Vector{Float64}}(undef, n_drugs)
    meanRarr = Vector{Vector{Float64}}(undef, n_drugs)
    Larr = Vector{Vector{Float64}}(undef, n_drugs)
    Lmaxarr = Vector{Vector{Float64}}(undef, n_drugs)
    
    for i = 1:n_drugs
        println(drugs[i])

        noisearr[i] = []
        Iarr[i] = []
        Aarr[i] = []
        meanRarr[i] = []
        Larr[i] = []
        Lmaxarr[i] = []

        for j = 1:n_population+1
            println("Step $j")

            tout, xout, X, sim = simulatepopulation(models[i][j], DRI1(), 1000;
                u0arr = u0arr)
            if j == 1
                push!(solutions, xout)
                push!(times, tout)
                push!(Xarr, X)
                push!(simarr, sim)
                continue
            end
            Rarr = []
            R = fill(0.0, length(xout))
            for k = 1:length(sim)
                R .+= sim[k][3, :]
            end
            R ./= length(sim)

            A, noise, I = getparameters(models[i][j], [:A, :noise, :I])
            

            push!(noisearr[i], noise)
            push!(Iarr[i], I)
            push!(Aarr[i], A)
            push!(meanRarr[i], mean(R))
            push!(Larr[i], mean(xout))
            push!(Lmaxarr[i], mean(xout))
        end

    end

    params[:noisearr] = noisearr
    params[:Iarr] = Iarr
    params[:Aarr] = Aarr
    params[:meanRarr] = meanRarr
    params[:Larr] = Larr

    
    params[:solutions] = solutions
    params[:times] = times
    params[:Xarr] = Xarr
    params[:simarr] = simarr

end

run_all_simulations!(params)

## ==
function plot_parameters!(params)
    dir = params[:dir]
    labels = params[:titles]
    meanRarr = params[:meanRarr]
    noisearr = params[:noisearr]
    Iarr = params[:Iarr] 
    Aarr = params[:Aarr]
    Larr = params[:Larr]

    n = length(Iarr[1])

    x = 1:length(labels)  # the label locations

    fig, ax = subplots(1, 2; figsize = (4, 2))


    for (i, x) in enumerate(Larr)
        ax[1].hlines(mean(x), i - 0.5 / 1.5, i + 0.5 / 1.5; color = :red,
        zorder = 1000)
    end
    ax[1].plot(
        [[0.5 .* rand(n).+ x .- 0.5 ./ 2 for x in 1:length(Larr)]...;],
        [Larr...;]; color = :black, marker = ".", linewidth = 0.0,
        markeredgewidth = 0.0, markersize = 4.0
    )

    println("L = $(mean.(Larr)) ± $(std.(Larr))")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax[1].set_ylabel("Transcription activation")
    ax[1].set_title("Transcription activation")
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(labels; rotation=20, ha="right")
    
    ax[1].legend(["Parameter sets", "Mean"]; fancybox=false, edgecolor=:black)

    for (i, x) in enumerate(noisearr)
        ax[2].hlines(mean(x), i - 0.5 / 1.5, i + 0.5 / 1.5; color = :red,
        zorder = 1000)
    end
    ax[2].plot(
        [[0.5 .* rand(n).+ x .- 0.5 ./ 2 for x in 1:length(noisearr)]...;],
        [noisearr...;]; color = :black, marker = ".", linewidth = 0.0,
        markeredgewidth = 0.0, markersize = 4.0 #, alpha = 0.2
        )

    println("σ = $(mean.(noisearr)) ± $(std.(noisearr))")

    ax[2].set_ylabel("Noise intensity")
    ax[2].set_title("Noise intensity")
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(labels; rotation=20, ha="right")

    fig.tight_layout()

    println("DRAWN")
    savefigure(fig, "fig4_parameters.svg")

end

plot_parameters!(params)


## plot_simulations ===========================================================
function plot_simulations(params)

    colorarr = params[:colors]
    dir = params[:dir]
    solutions = params[:solutions]
    times = params[:times]
    drugs = params[:drugs]
    Xarr = params[:Xarr]
    ndrug = length(drugs)
    titles = params[:titles]

    fig, axarr = subplots(3; figsize=(4, 4))

    for i = 1:ndrug
        y = solutions[i]
        ty = times[i]
        y_cells = Xarr[i]
        axarr[i].plot(ty, y; color = colorarr[i], linewidth = 2)
        axarr[i].autoscale(tight = true)
        axarr[i].text(0, 0.8, titles[i], fontsize = big_font)
        axarr[i].autoscale(tight = true)
        axarr[i].set_xticks(minimum(ty):24:maximum(ty))
        axarr[i].set_xticklabels(round.(Int, (minimum(ty):24:maximum(ty)) .- minimum(ty)))
        axarr[i].plot(ty, y_cells[:, 1:5]; color = colorarr[i], linewidth = 1)

        axarr[i].set_xlabel("Time (hours)")
        axarr[i].set_ylabel("Lumin. (au)")
        axarr[i].set_ylim(0.0, 0.8)
        axarr[i].legend(["Mean", "Individual traces"]; ncol = 2, fancybox=false, edgecolor=:black)
    end

    fig.tight_layout()

    savefigure(fig, dir, "fig4a.svg")
end

plot_simulations(params)