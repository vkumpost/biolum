function D_load_initial_population!(params)
    initial_population_filename = params[:in_initial_population_filename]

    df = loadcsv(initial_population_filename)
    population = Matrix(df[2:end, Not(:fitness)])
    parameters = df[1, Not(:fitness)]
    pnames = [Symbol(x) for x in names(parameters)]
    p = Vector(parameters)

    params[:p] = p
    params[:pnames] = pnames
    params[:population] = population
end

function D_load_normalization_parameters!(params)
    filename = params[:in_normalization_parameters_filename]
    if isempty(filename)
        params[:normpars] = [NaN, NaN, NaN, NaN]
    else
        df = loadcsv(filename)
        params[:normpars] = Vector(df[1, :])
    end
end


function D_load_training_data!(params)
    plate_set = params[:plate_set]
    drug = params[:drug]
    showplot = params[:showplots]

    df = loadcsv(lowercase("norm_plate_d1_$(plate_set).csv"))
    df_sel = selectwells(df, Regex(drug))
    tdata1 = df_sel[:, "Time"]
    X = wells2matrix(df_sel)
    n = size(X, 1)
    xdata1 = vec(mean(X, dims=2))
    events1 = detectevents(df_sel)

    df = loadcsv(lowercase("norm_plate_d2_$(plate_set).csv"))
    df_sel = selectwells(df, Regex(drug))
    tdata2 = df_sel[:, "Time"]
    X = wells2matrix(df_sel)
    n = size(X, 1)
    xdata2 = vec(mean(X, dims=2))
    events2 = detectevents(df_sel)

    tdataarr = [tdata1, tdata2]
    xdataarr = [xdata1, xdata2]
    eventsarr = [events1, events2]

    if showplot
        figure()
        for i = 1:2
            subplot(2, 1, i)
            plot(tdataarr[i], xdataarr[i], color=:black)
            plotevents(eventsarr[i])
        end
    end

    params[:tdataarr] = tdataarr
    params[:xdataarr] = xdataarr
    params[:eventsarr] = eventsarr
end

function D_add_ode_model!(params)
    
    model_name = params[:model_name]
    showplots = params[:showplots]
    alg = params[:ode_alg]
    p = params[:p]

    model = loadmodel(model_name)
    events = createevents([(:DD, 10*24)])
    model = addcallback(model, events, :I)
    model = remakemodel(model, saveat = 240.0:0.1:480.0, tspan = (0.0, 480.0))
    model = setparameters(model, model.parameters, p[1:end-1])

    if showplots
        figure()
        _, _, sol = solvemodel(model, alg)
        for j = 1:length(sol[1])
            plot(sol.t, sol[j, :])
        end
    end

    params[:ode_model] = model

end

function D_add_sde_models!(params)
    model_name = params[:model_name]
    p = params[:p]
    pnames = params[:pnames]
    eventsarr = params[:eventsarr]
    tdataarr = params[:tdataarr]

    models = [loadmodel(model_name * "sde"), loadmodel(model_name * "sde")]
    for i = 1:2
        models[i] = setparameters(models[i], pnames, p)
        models[i] = addcallback(models[i], eventsarr[i], :I)
        tspan = (0.0, maximum(tdataarr[i]))
        models[i] = remakemodel(models[i], tspan=tspan,
            saveat=tdataarr[i])
    end

    params[:models] = models
end

function D_load_initial_conditions!(params)
    filename = params[:in_initial_conditions_filename]
    models = params[:models]

    df = loadcsv(filename)
    state_names = [Symbol(x) for x in names(df)]
    u0arr = Matrix(df)
    @assert all(state_names .== models[1].species)
    params[:u0arr] = u0arr
end

function D_generate_cost_function!(params)

    osc_con = params[:osc_con]
    ode_model = params[:ode_model]
    ode_alg = params[:ode_alg]
    models = params[:models]
    pnames = params[:pnames]
    sde_alg = params[:sde_alg]
    trajectories = params[:trajectories]
    u0arr = params[:u0arr]
    normpars = params[:normpars]
    xdataarr = params[:xdataarr]
    tdataarr = params[:tdataarr]
    folder_progress = params[:folder_progress]

    function costfunction(x; showplot=false, figfilename=nothing)

        # ODE PART ------------------------------------------------------------
        optmodel = setparameters(ode_model, ode_model.parameters, x[1:end-1])
        tdd, xdd, sol = solvemodel(optmodel, ode_alg)
        if any((s.retcode != :Success for s in sol))
            return Inf
        end
        τ = sustainedperiod(xdd, tdd)
        if showplot
            fig, ax = subplots(3)
            ax[1].plot(tdd, xdd; color = :black)
            ax[1].set_title("τ = $τ")
        end
        if osc_con && isnan(τ)
            println("No-no")
            return Inf
        end   

        # SDE PART ------------------------------------------------------------
        optmodel01 = setparameters(models[1], pnames, x)
        tsim01, xsim01 = simulatepopulation2(optmodel01, sde_alg,
            trajectories; u0arr=u0arr)

        if any(isnan.(tsim01))
            return Inf
        end

        optmodel02 = setparameters(models[2], pnames, x)
        tsim02, xsim02 = simulatepopulation2(optmodel02, sde_alg,
            trajectories; u0arr=u0arr)

        if any(isnan.(tsim01))
            return Inf
        end

        if isnan(normpars[1])
            xsim01 = zscore(xsim01)
            xsim02 = zscore(xsim02) 
        else
            xsim01 = zscore(xsim01, normpars[1], normpars[2])
            xsim02 = zscore(xsim02, normpars[3], normpars[4]) 
        end

        L01 = squarederror(xdataarr[1], xsim01, fun=mean)
        L02 = squarederror(xdataarr[2], xsim02, fun=mean)

        totloss = L01 + L02

        # show plots ----------------------------------------------------------
        if showplot
            ax[2].plot(tdataarr[1], xdataarr[1], color=:black)
            ax[2].plot(tsim01, xsim01)
            ax[2].set_title("L01 = $(round(L01, digits=2))")
            ax[3].plot(tdataarr[2], xdataarr[2], color=:black)
            ax[3].plot(tsim02, xsim02)
            ax[3].set_title("L02 = $(round(L02, digits=2))")
            fig.tight_layout()
            if !isnothing(figfilename)
                savefigure(fig, folder_progress, figfilename)
            end
        end

        if isnan(totloss)
            return Inf
        end

        return totloss

    end

    if params[:showplots]
        costfunction(params[:p]; showplot=params[:showplots])
    end

    params[:costfunction] = costfunction

end

function D_optimize!(params)
    pnames = params[:pnames]
    population = params[:population]
    costfunction = params[:costfunction]
    MaxSteps = params[:bbo_MaxSteps]
    SearchRange = params[:bbo_SearchRange]

    function createoptcallback()
        best_fitness = Inf
        counter = 0
        function optcallback(oc)
            current_individuals = BBO.best_candidate(oc)
            current_fitness = BBO.best_fitness(oc)
            if current_fitness < best_fitness
                close(:all)
                counter += 1
                best_fitness = current_fitness
                costfunction(current_individuals; showplot=true,
                    figfilename="progress_$(counter + 10000000).png")
                println("NEW FITNESS: $current_fitness")
            end
        end
        return optcallback
    end

    optimizer_kwargs = Dict(
        :method => :adaptive_de_rand_1_bin_radiuslimited,
        :SearchRange => SearchRange,
        :NumDimensions => length(pnames),
        :MaxSteps => MaxSteps,
        :TraceMode => :compact,  # compact, silent
        :PopulationSize => size(population, 1),
        :Population => transpose(population),
        :CallbackFunction => createoptcallback(),
        :CallbackInterval => 0.0
    )

    res = BBO.bboptimize(costfunction; optimizer_kwargs...)

    if params[:showplots]
        costfunction(BBO.best_candidate(res); showplot=params[:showplots])
    end

    params[:optimization_results] = res

end


function D_reval_fitness!(params)

    n_reval = params[:n_reval]
    res = params[:optimization_results]
    pnames = params[:pnames]
    costfunction = params[:costfunction]


    population = transpose(BBO.population(res).individuals)
    fitness = BBO.population(res).fitness

    fitness_reval = []
    n = size(population, 1)
    for i = 1:n
        earr = []
        for j = 1:n_reval
            e = costfunction(Vector(population[i, :]))
            push!(earr, e)
        end
        push!(fitness_reval, median(earr))
        println("$i: $(fitness[i]) | $(median(earr))")
    end

    i = argmin(fitness_reval)
    best_candidate = population[i, :]
    best_fitness = fitness_reval[i]

    final_population = [reshape(best_candidate, 1, length(best_candidate)); population]
    final_fitness = [best_fitness; fitness_reval]
    df_results = DataFrame([final_population final_fitness], [pnames...; :fitness])
    params[:df_results] = df_results
end

function D_save_results(params)
    df = params[:df_results]
    filename = params[:out_results_filename]
    savecsv(df, filename)
end

function D_save_normalization_parameters(params)

    showplot = params[:showplots]
    df_results = params[:df_results]
    models = deepcopy(params[:models])
    alg = params[:sde_alg]
    trajectories = params[:trajectories]
    u0arr = params[:u0arr]
    filename = params[:out_normalization_parameters_filename]

    if isempty(filename)
        return
    end

    p = DataFrame(df_results[1, 1:end-1])
    for i = 1:2
        models[i] = setparameters(models[i], p)
    end

    if showplot
        figure()
    end

    normpars = [0.0, 0.0, 0.0, 0.0]
    for i = 1:1000
        tsim01, xsim01 = simulatepopulation2(models[1], alg, trajectories;
            u0arr=u0arr)
        normpars[1] += mean(xsim01)
        normpars[2] += std(xsim01)

        tsim02, xsim02 = simulatepopulation2(models[2], alg, trajectories;
            u0arr=u0arr)
        normpars[3] += mean(xsim02)
        normpars[4] += std(xsim02)

        if showplot
            subplot(2, 1, 1)
            plot(tsim01, xsim01)
            subplot(2, 1, 2)
            plot(tsim02, xsim02)
        end
    end

    normpars ./= 1000
    df = DataFrame([[x] for x in normpars], ["μ1", "σ1", "μ2", "σ2"])
    savecsv(df, filename)

end

function D_run_all!(params)
    D_load_initial_population!(params)
    D_load_normalization_parameters!(params)
    D_load_training_data!(params)
    D_add_ode_model!(params)
    D_add_sde_models!(params)
    D_load_initial_conditions!(params)
    D_generate_cost_function!(params)
    D_optimize!(params)
    D_reval_fitness!(params)
    D_save_results(params)
    D_save_normalization_parameters(params)
end