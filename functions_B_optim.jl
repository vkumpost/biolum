## Load reference damping ration ==============================================
function B_load_damping_ratio!(params)
    df_dref = loadcsv("damped_sine_parameters.csv")
    dref = df_dref[1, :d]
    params[:dref] = dref
end

## Load the parameters form A =================================================
function B_load_A_results!(params)
    filename = params[:in_A]
    showplots = params[:showplots]

    df = loadcsv(filename)
    pnames = [Symbol(x) for x in names(df)]
    population = Matrix(df)

    if showplots
        fig, ax = plt.subplots()
        ax.imshow(transpose(population))
        ax.set_yticks(0:length(pnames)-1)
        ax.set_yticklabels(pnames)
        ax.set_ylabel("Parameter")
        ax.set_xlabel("Parameter set")
        ax.set_title("Parameters from A")
    end

    params[:pnames] = pnames
    params[:population] = population
end

## Set up the model =========================================================
function B_create_model!(params)
    ndays = params[:ndays]
    alg = params[:alg]
    trajectories = params[:trajectories]
    model_name = params[:model_name]
    pnames = params[:pnames]
    population = params[:population]
    i = params[:i]
    showplots = params[:showplots]

    events = createevents([(:DD, 1*24), (:LD, 10, 12, 12)])
    A = events[end, 2] + 12.0  # first point from which to estimate daping ratio
    B = events[end, 2] + 12.0 + ndays*24.0  # last point

    model = loadmodel(model_name * "sde")
    model = setparameters(model, pnames, population[i, :])
    model = addcallback(model, events, :I)
    model = remakemodel(model, tspan=(0.0, B), saveat=0.1)

    if showplots
        printparameters(model)
        tout, xout, X, sim = simulatepopulation(model, alg, trajectories)
        sol = dropdims(mean(sim, dims=3), dims=3)
        figure()
        subplot(211)
        for j = 1:length(sol[:, 1])
            plot(sim[1].t, sol[j, :])
        end
        plotevents(events)
        subplot(212)
        plot(tout, xout)
        plotevents(events)
    end

    model = remakemodel(model, saveat=A:0.1:B)

    if showplots
        tout, xout, X, sim = simulatepopulation(model, alg, trajectories)
        sol = dropdims(mean(sim, dims=3), dims=3)
        figure()
        subplot(211)
        for j = 1:length(sol[:, 1])
            plot(sim[1].t, sol[j, :])
        end
        plotevents(events)
        subplot(212)
        plot(tout, xout)
        plotevents(events)
    end

    params[:model] = model

end


## Set up the cost function =================================================
function B_create_cost_function!(params)
    showplots = params[:showplots]
    model = params[:model]
    alg = params[:alg]
    trajectories = params[:trajectories]
    dref = params[:dref]

    function costfunction(x; showplots=false)

        if x isa Array
            x = x[1]
        end

        optmodel = setparameters(model, [:noise], [x])

        t, x, _ = simulatepopulation(optmodel, alg, trajectories)
        xx = zscore(x)
        tt = t .- t[1]

        p = fitdampedsine(tt, xx, [2.5, 0.01, 25.0, 0.0])
        xxest = dampedsine(tt, p)
        R = cor(xx, xxest) ^ 2
        if showplots
            println()
            println("Damped sine fit:")
            println("- Amplitude: $(p[1])")
            println("- Damping ration: $(p[2])")
            println("- Period: $(p[3])")
            println("- Phase: $(p[4])")
            println()
            println("Quality of fit:")
            println("- R squared: $R")
        end
        if showplots
            plt.plot(tt, xx, color=:gray)
            plt.plot(tt, xxest, color=:black)
        end
        if R > 0.8
            return (p[2] - dref)^2
        else
            return 1000.0
        end
    end

    if showplots
        costfunction(0.03; showplots=showplots)
    end

    params[:costfunction] = costfunction

end


## Set up optimizer ===========================================================
function B_optimize!(params)
    noise_range = params[:noise_range]
    i = params[:i]
    population = params[:population]
    showplots = params[:showplots]
    MaxFuncEvals = params[:MaxFuncEvals]
    PopulationSize = params[:PopulationSize]
    dref = params[:dref]

    if !(:res in keys(params))
        params[:res] = fill(NaN, size(population, 1))
    end

    B_create_model!(params)
    B_create_cost_function!(params)

    cf = params[:costfunction]

    optimizer_kwargs = Dict(
        :method => :adaptive_de_rand_1_bin_radiuslimited,
        :SearchRange => (noise_range[1], noise_range[2]),
        :NumDimensions => 1,
        :PopulationSize => PopulationSize,
        :MaxFuncEvals => MaxFuncEvals,
        :TraceMode => :compact,
        :TargetFitness => 0.0,
        :FitnessTolerance => (0.1*dref)^2
    )
    res = BBO.bboptimize(cf; optimizer_kwargs...)

    if showplots
        cf(BBO.best_candidate(res), showplots=showplots)
    end

    params[:res][i] = BBO.best_candidate(res)[1]
end


## Optimize all ===========================================================
function B_optimize_all!(params)
    population = params[:population]
    for i = 1:size(population, 1)
        println("Iteration $i")
        params[:i] = i
        B_optimize!(params)
    end
end


## Save results =============================================================
function B_save_results(params)
    filename = params[:out_filename]

    population = params[:population]
    pnames = params[:pnames]
    res = params[:res]

    labels = [pnames...; :noise]
    data = hcat(population, res)

    df = DataFrame(data, labels)
    savecsv(df, filename)

end