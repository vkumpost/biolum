# This file contains functions for the manipulation and simulation of the model.

"""
`BCModel`

Composite type representing bioluminescence circadian model.

**Fields**
- `species`: Tuple of symbols holding the variable names.
- `parameters`: Tuple of symbols holding the parameter names.
- `outfun`: Output function of the form `x = outfun(sol)`.
- `prob`: `ODEProblem` or `SDEProblem`.
"""
struct BCModel
    species
    parameters
    outfun
    prob
end


"""
`solvemodel(model::BCModel, alg=nothing; kwargs)`

Numerically solve differential equation model. Returns `(tout, xout, sol)`
where `tout` is a time vector, `xout` is a vector of output values and `sol` is
a `DifferentialEquations` solution.
"""
function solvemodel(model::BCModel, alg=nothing; kwargs...)
    sol = solve(model.prob, alg; kwargs...)
    tout = sol.t
    xout = model.outfun(sol)
    return tout, xout, sol
end


"""
`simulatepopulation(model::BCModel, alg, trajectories; u0arr=[], kwargs...)`

Run simulation of `model` multiple times. Returns `(tout, xout, X, sim)` where
`tout` is time time vector, `xout` is mean of indiviudal traces that are also
row-wise stored in `X`. `sim` is `EnsembleSolution`.

Important: You must specify `saveat` or other keyword that ensures that all the
traces will be saved at the same time point!
"""
function simulatepopulation(model::BCModel, alg, trajectories;
    u0arr=[], kwargs...)

    if isempty(u0arr)
        prob_func = (prob, i, repeat) -> (prob)
    else
        prob_func = (prob, i, repeat) -> remake(prob,
            u0=u0arr[i, :])
    end

    ensembleprob = EnsembleProblem(model.prob, prob_func=prob_func)
    sim = solve(ensembleprob, alg, EnsembleThreads();
        trajectories=trajectories, kwargs...)

    nspecies, nsamples, ntrajectories = size(sim)

    if any((s.retcode != :Success for s in sim))
        X = fill(NaN, nsamples, ntrajectories)
        tout = fill(NaN, nsamples)
        xout = fill(NaN, nsamples)
    else
        X = fill(NaN, nsamples, ntrajectories)
        tout = fill(0.0, nsamples)
        for i = 1:length(sim)
            X[:, i] = model.outfun(sim[i])
            tout .+= sim[i].t
        end
        tout ./= ntrajectories
        xout = vec(mean(X, dims=2))
    end

    return (tout, xout, X, sim)

end


"""
`simulatepopulation2(model::BCModel, alg, trajectories; u0arr=[], kwargs...)`

Works exactly like `simulatepopulation2` but returns only `(tout, xout)` where
`tout` is time time vector, `xout` is mean of indiviudal traces.

Important: You must specify `saveat` or other keyword that ensures that all the
traces will be saved at the same time point!
"""
function simulatepopulation2(model::BCModel, alg, trajectories;
    u0arr=[], kwargs...)

    if !isempty(u0arr)
        model = remakemodel(model, u0=u0arr[1, :])
    end
    (tout, xout, sol) = solvemodel(model, alg; kwargs...)

    Threads.@threads for i = 2:trajectories

        if !isempty(u0arr)
            model = remakemodel(model, u0=u0arr[i, :])
        end
        (_, xout_tmp, sol) = solvemodel(model, alg; kwargs...)

        if any((s.retcode != :Success for s in sol))
            tout = fill(NaN, length(xout))
            xout = fill(NaN, length(xout))
            return (tout, xout)
        else
            xout .+= xout_tmp
        end

    end

    xout ./= trajectories

    return (tout, xout)

end


"""
`argparameters(model::BCModel, pnames)`

Find parameter indices specified by `pnames` array.
"""
function argparameters(model::BCModel, pnames)
    indices = Int[]
    for pname in pnames
        idx = findfirst(model.parameters .== pname)
        push!(indices, idx)
    end
    return indices
end


"""
`setparameters(model::BCModel, pnames, pvalues)`

Set parameter values `pvalues` to parameters `pnames`.
"""
function setparameters(model::BCModel, pnames, pvalues)

    newmodel = deepcopy(model)
    idx = argparameters(model, pnames)
    newmodel.prob.p[idx] = pvalues

    return newmodel
end


"""
`setparameters(model::BCModel, df::DataFrame)`

Set parameter values from a `DataFrame`.
"""
function setparameters(model::BCModel, df::DataFrame)
    pnames = [Symbol(x) for x in names(df)]
    pvalues = Matrix(df)
    @assert size(pvalues, 1) == 1 "Dataframe has more than one row!"
    return setparameters(model, pnames, pvalues)
end


"""
`getparameters(model::BCModel, pnames)`

Get a parameter value.
"""
function getparameters(model::BCModel, pnames)
    idx = argparameters(model, pnames)
    return model.prob.p[idx]
end


"""
`printparameters(model::BCModel)`

Print parameter values to the standard output.
"""
function printparameters(model::BCModel)
    println("========================\nModel parameters")
    for (name, value) in zip(model.parameters, model.prob.p)
        println("------------------------")
        println("$(name) \t| $(value)")
    end
    println("========================")
end


"""
`remakemodel(model::BCModel; kwargs...)`

Apply `remake` to the model's problem.
"""
function remakemodel(model::BCModel; kwargs...)

    prob = remake(model.prob; kwargs...)
    newmodel = BCModel(model.species, model.parameters, model.outfun, prob)

    return newmodel

end


"""
`addcallback(model::BCModel, events::Matrix, pname)`

Add a callback to a `model` that applies `events` on the parameter `pname`.
"""
function addcallback(model::BCModel, events::Matrix, pname)

    ip = findfirst(model.parameters .== pname)
    cb = createcallback(events, ip)
    return remakemodel(model, callback=cb)

end