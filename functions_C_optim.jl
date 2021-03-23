## load_initial_population ====================================================
function C_load_initial_population!(params)
    model_name = params[:model_name]
    filename = params[:in_initial_population_filename]
    df = loadcsv(filename)
    params[:pnames] = Symbol.(names(df))
    params[:population] = Matrix(df)
end


## Load data ==================================================================
function C_load_fitting_data!(params)

    showplots = params[:showplots]

    data = loaddata("Plate U1")
    data = selectwells(data, r"WT [A|B]")
    data_events = detectevents(data)
    it = 10  # number of transient days
    events, offset = addtransientevents(data_events, 12, it; i=4)
    data.Time .+= offset
    idx = events[it, 2] .< data.Time .< events[it + 4, 2]
    data_sel = data[idx, :]
    
    if showplots
        figure()
        subplot(211)
        plot(data[:, "Time"], data[:, "WT A (5)"], color=:black)
        plot(data_sel[:, "Time"], data_sel[:, "WT A (5)"], color=:C00)
        plotevents(events)
    end

    idx = events[it, 2] .< data_sel.Time .< events[it+5, 1]
    tdata = Vector(data_sel[:, "Time"])
    df_rawdata = transformwells(data_sel, zscore)
    xdata = aggregatewells(df_rawdata, mean).mean
    X = wells2matrix(df_rawdata)

    if showplots
        subplot(212)
        plot(tdata, X, color=:black)
        plot(tdata, xdata, color=:red)
        plotevents(events)
    end

    params[:tdata] = tdata 
    params[:xdata] = xdata
    params[:events] = events

end


## create_models_ode! ========================================================
function C_add_ode_model!(params)
    
    model_name = params[:model_name]
    showplots = params[:showplots]
    alg = params[:ode_alg]
    population = params[:population]
    p = population[1, :]

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


## Prepare model ==============================================================
function C_add_sde_model!(params)
    model_name = params[:model_name]
    events = params[:events]
    tdata = params[:tdata]
    alg = params[:sde_alg]
    trajectories = params[:trajectories]
    showplots = params[:showplots]
    population = params[:population]
    p = population[1, :]
    
    model = loadmodel(model_name * "sde")
    model = addcallback(model, events, :I)
    model = remakemodel(model, tspan=(0.0, tdata[end]), saveat=0.1)
    model = setparameters(model, model.parameters, p)

    if showplots
        @time tout, xout = simulatepopulation2(model, alg, trajectories)
        @time tout2, xout2 = simulatepopulation2(model, alg, trajectories)

        figure()
        subplot(211)
        plot(tout, xout)
        plot(tout2, xout2)
        plotevents(events)
    end

    model = remakemodel(model, tspan=(0.0, tdata[end]), saveat=tdata)

    if showplots
        @time tout, xout = simulatepopulation2(model, alg, trajectories);
        @time tout2, xout2 = simulatepopulation2(model, alg, trajectories);
        subplot(212)
        plot(tout, xout)
        plot(tout2, xout2)
        plotevents(events)
    end

    params[:sde_model] = model
end


## Cost function ==============================================================
function C_generate_cost_function!(params)

    sde_model = params[:sde_model]
    tdata = params[:tdata]
    xdata = params[:xdata]
    sde_alg = params[:sde_alg]

    ode_model = params[:ode_model]
    ode_alg = params[:ode_alg]

    trajectories = params[:trajectories]
    showplots = params[:showplots]

    folder_progress = params[:folder_progress]

    # Cost function
    function costfunction(x; showplot=false, figfilename=nothing)
        
        # ODE PART ============================================================
        # == DD ==
        optmodel = setparameters(ode_model, ode_model.parameters, x[1:end-1])
        tdd, xdd, sol = solvemodel(optmodel, ode_alg)
        if any((s.retcode != :Success for s in sol))
            return Inf
        end
        τ = sustainedperiod(xdd, tdd)
        if showplot
            fig, ax = subplots(2)
            ax[1].plot(tdd, xdd; color = :black)
            ax[1].set_title("τ = $τ")
        end
        if isnan(τ)
            println("No-no")
            return Inf
        end

        # SDE PART ============================================================
        optmodel = setparameters(sde_model, sde_model.parameters, x)

        tout, xout = simulatepopulation2(optmodel, sde_alg, trajectories)

        if any(isnan.(tout))
            totloss = Inf
        else
            xout_norm = zscore(xout)
        end

        totloss = squarederror(xdata, xout_norm, fun=mean)
        
        if showplot
            ax[2].plot(tout, xout_norm, color=:red)
            ax[2].plot(tdata, xdata, color=:black)
            ax[2].set_title("L = $(round(totloss, digits=2))")
            fig.tight_layout()
            if !isnothing(figfilename)
                savefigure(fig, folder_progress, figfilename)
            end
        end

        return totloss

    end

    if showplots
        costfunction(sde_model.prob.p; showplot=showplots)
    end

    params[:costfunction] = costfunction
end

## optimalization =============================================================
function C_optimize!(params)
    SearchRange = params[:SearchRange]
    MaxSteps = params[:MaxSteps]
    model = params[:sde_model]
    costfunction = params[:costfunction]
    population = params[:population]

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
        :Method => :adaptive_de_rand_1_bin_radiuslimited,
        :SearchRange => SearchRange,
        :NumDimensions => length(model.parameters),
        :MaxSteps => MaxSteps,
        :TraceMode => :compact,  # compact, silent
        :PopulationSize => size(population, 1),
        :Population => transpose(population),
        :CallbackFunction => createoptcallback(),
        :CallbackInterval => 0.0
    )

    res = BBO.bboptimize(costfunction; optimizer_kwargs...)

    best_candidate = BBO.best_candidate(res)
    best_fitness = BBO.best_fitness(res)
    
    final_population = transpose(BBO.population(res).individuals)
    fitness = BBO.population(res).fitness

    population_save = [reshape(best_candidate, 1, length(best_candidate));
        final_population]
    fitness_save = [best_fitness; fitness]

    colnames = [[x for x in model.parameters]...; :fitness]
    final_population = DataFrame(cat(population_save, fitness_save, dims=2),
        colnames)

    params[:final_population] = final_population

end

## C_reeval_fitness ==========================================================
function C_reval_fitness!(params)
    costfunction = params[:costfunction]
    n_reval = params[:n_reval]
    final_population = deepcopy(params[:final_population])

    population = final_population[2:end, 1:end-1]
    fitness = final_population[2:end, end]

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

    # edit the original population based on the reval fitness
    final_population[2:end, end] = fitness_reval
    i = argmin(fitness_reval)
    final_population[1, :] = final_population[1+i, :]

    params[:final_population] = final_population
end

## Save sesutls ===============================================================
function C_save_results(params)
    model_name = params[:model_name]
    final_population = params[:final_population]
    showplots = params[:showplots]
    filename_population = params[:out_final_population_filename]

    if showplots
        costfunction = params[:costfunction]
        costfunction(Vector(final_population[1, 1:end-1]); showplot=showplots)
    end

    savecsv(final_population, filename_population)

end