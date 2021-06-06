## Create light-pacing regimes ================================================
function A_add_pacing_regimes!(params)
    eventsarr = Vector{Array{Float64,2}}(undef, 3)
    
    eventsarr[1] = createevents([(:DD, 10*24)])
    eventsarr[2] = createevents([(:DD, 10*24), (:LD, 150, 12, 12)])
    eventsarr[3] = createevents([(:DD, 10*24), (:DD, 12), (:LD, 150, 12, 12)])
    
    if params[:showplots]
        N = length(eventsarr)
        fig, ax = subplots(N, 1)
        for i = 1:N
            plotevents(ax[i], eventsarr[i])
            ax[i].set_xlim(0, 600.0)
        end
    end

    params[:eventsarr] = eventsarr
    params[:n] = 3

end


## Create models ==============================================================
function A_create_models!(params)
    
    model_name = params[:model_name]
    eventsarr = params[:eventsarr]
    entrainment_days = params[:entrainment_days]

    modelarr = Vector{BCModel}(undef, 3)

    model = loadmodel(model_name)
    model = addcallback(model, eventsarr[1], :I)
    model = remakemodel(model, saveat = 240.0:0.1:480.0, tspan = (0.0, 480.0))
    modelarr[1] = model

    model = loadmodel(model_name)
    model = addcallback(model, eventsarr[2], :I)
    model = remakemodel(model, saveat =
        eventsarr[2][entrainment_days + 1, 1]:0.1:600, tspan = (0.0, 600.0))
    modelarr[2] = model

    model = loadmodel(model_name)
    model = addcallback(model, eventsarr[3], :I)
    model = remakemodel(model, saveat =
        eventsarr[3][entrainment_days + 1, 1]:0.1:612, tspan = (0.0, 612.0))
    modelarr[3] = model

    if params[:showplots]
        for i = 1:3
            subplot(3, 1, i)
            model = setparameters(modelarr[i], [:I], [10.0])
            _, _, sol = solvemodel(model, params[:alg])
            for j = 1:length(sol[1])
                plot(sol.t, sol[j, :])
            end
            plotevents(eventsarr[i])
        end
    end

    params[:modelarr] = modelarr

end

## Cost function ==============================================================
function A_create_cost_function!(params)
   
    modelarr = params[:modelarr]
    eventsarr = params[:eventsarr]
    alg = params[:alg]
    xnames = params[:modelarr][1].parameters

    function costfunction(x; showplots = false)

        toterr = 0.0
    
        # == DD ==
        optmodel = setparameters(modelarr[1], xnames, x)
        tdd, xdd, sol = solvemodel(optmodel, alg)
        # check if the integration was success
        if any((s.retcode != :Success for s in sol))
            return Inf
        end
        τ = sustainedperiod(xdd, tdd)
        if showplots
            println("DD period is $(τ) hours.")
        end
        if isnan(τ)
            return 1000.0
        else
            toterr = toterr + (τ - 24.0)^2
        end
    
        # == LD ==
        optmodel = setparameters(modelarr[2], xnames, x)
        tld, xld, sol = solvemodel(optmodel, alg)
        # check if the integration was success
        if any((s.retcode != :Success for s in sol))
            return Inf
        end
        τ = sustainedperiod(xld, tld)
        if showplots
            println("LD period is $(τ) hours.")
        end
        if isnan(τ)
            return 1000.0
        else
            toterr = toterr + (τ - 24.0)^2
        end
    
        # == LD 2 ==
        optmodel = setparameters(modelarr[3], xnames, x)
        tld2, xld2, sol = solvemodel(optmodel, alg)
        # check if the integration was success
        if any((s.retcode != :Success for s in sol))
            return Inf
        end
    
        peaks01 = findpeaks(xld, tld)
        locs01 = peaks01.locs
        peaks02= findpeaks(xld2, tld2)
        locs02 = peaks02.locs
        if isempty(locs01) || isempty(locs02)
            return 1000.0
        end
    
        Δ = locs02[1] - locs01[1] - 12.0
        toterr = toterr + Δ^2
        if showplots
            println("LD phase difference is $(Δ) hours.")
        end
    
        if showplots
            figure()
            plt.subplot(311)
            plt.plot(tdd, xdd)
            plt.subplot(312)
            plt.plot(tld, xld)
            plotevents(eventsarr[2])
            plt.subplot(313)
            plt.plot(tld2, xld2)
            plotevents(eventsarr[3])
        end
    
        return toterr
    
    end

    if params[:showplots]
        costfunction(params[:modelarr][1].prob.p; showplots=params[:showplots])
    end

    params[:costfunction] = costfunction

end


## Run optimizer ==============================================================
function A_optimize!(params)
    
    model_parameters = params[:modelarr][1].parameters
    NumDimensions = length(model_parameters)
    costfunction = params[:costfunction]
    
    optimizer_kwargs = Dict(
        :Method => :adaptive_de_rand_1_bin_radiuslimited,
        :SearchRange => params[:SearchRange],
        :NumDimensions => NumDimensions,
        :MaxFuncEvals => params[:MaxFuncEvals],
        :TraceMode => params[:TraceMode],
        :TargetFitness => 0.0,
        :FitnessTolerance => 0.001
    )

    res = BBO.bboptimize(costfunction; optimizer_kwargs...) 

    if params[:showplots]
        estimated_parameters = BBO.best_candidate(res)
        costfunction(estimated_parameters; showplots=params[:showplots])
    end

    if BBO.best_fitness(res) < 0.001
        if :results in keys(params)
            push!(params[:results], BBO.best_candidate(res))
        else
            params[:results] = DataFrame([[x] for x in BBO.best_candidate(res)],
                [x for x in model_parameters])
        end
    end
end


## Find parameter sets ========================================================
function A_find_sets!(params)
    i = 0
    while i < params[:num_iter]
        println("Iteration $i")
        A_optimize!(params)
        if :results in keys(params)
            i = size(params[:results], 1)
        end
    end
end


## Save results ===============================================================
function A_save_results(params)
    filename = params[:out_filename]
    df = params[:results]
    savecsv(df, filename)
end