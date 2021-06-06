import CSV
import XLSX
import LsqFit
import Distributions
import MultivariateStats
import BlackBoxOptim

using StatsBase
using Statistics
using DataFrames
using PyPlot
using DifferentialEquations

const MS = MultivariateStats
const BBO = BlackBoxOptim
const PyCall = PyPlot.PyCall

include("BCModel.jl")
include("Peaks.jl")
include("functions_model.jl")
include("functions_A_optim.jl")
include("functions_B_optim.jl")
include("functions_C_optim.jl")
include("functions_D_optim.jl")

# Colorblind-friendly colors
const CB_COLORS = ["#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628",
    "#984ea3", "#999999", "#e41a1c", "#dede00"]

"""
`detectevents(t, x)`

Extract boundaries for the light-on events.
"""
function detectevents(t, x)
    if all(x .== 0)
        return Matrix{Float64}(undef, 0, 2)
    end
    event_starts = Array{Float64, 1}()
    event_ends = Array{Float64, 1}()
    for i = 1:length(x)-1
        if x[i] < 0.5 && x[i+1] > 0.5
            push!(event_starts, (t[i] + t[i+1])/2)
        elseif x[i] > 0.5 && x[i+1] < 0.5
            push!(event_ends, (t[i] + t[i+1])/2)
        end
    end
    if isempty(event_starts)
        event_starts = [t[1]; event_starts]
    end
    if isempty(event_ends)
        event_ends = [event_ends; t[end]]
    end
    if event_ends[1] < event_starts[1]
        event_starts = [t[1]; event_starts]
    end
    if event_ends[end] < event_starts[end]
        event_ends = [event_ends; t[end]]
    end
    events = [event_starts event_ends]
    return events
end


"""
`detectevents(df::DataFrame)`

Detect events in DataFrame `df` using columns `Time` and `Light`.
"""
function detectevents(df::DataFrame)
    return detectevents(df[:, "Time"], df[:, "Light"])
end


"""
`loaddata(plate; eventspan = [0, -1])`

Load plate as a `DataFrame`. `offset` specifies how many days to skip and then
aligns the data to the end of last skipped day.
"""
function loaddata(plate; eventspan = [0, -1])

    filename = joinpath(@__DIR__, "data", "biolumdata.xlsx")
    table, header = XLSX.readtable(filename, plate)
    floatarr = convert(Array{Array{Float64, 1}, 1}, table)
    df = DataFrame(floatarr, header)

    if eventspan[1] > 0

        first_event = eventspan[1]
        last_event = eventspan[2]

        events = detectevents(df)
        
        t_start = events[eventspan[1], 2]
        if last_event < first_event
            t_end = Inf
        else
            t_end = events[eventspan[2], 2]
        end
        
        idx = t_start .< df[:, "Time"] .< t_end
        df = df[idx, :]
        df[:, "Time"] .-= t_start
    end

    return df

end


"""
`plotdata(df)`

Plot `Biolum` data stored in a `DataFrame`.
"""
function plotdata(df)
    fig, ax = subplots()
    df2 = df[:, Not(["Time", "Light"])]
    ax.plot(df.Time, Matrix(df2), color=:black)
    plotevents(ax, detectevents(df))
end


"""
`plotevents(ax, events; ylims=[], color=:red, alpha=.25)`

Plot events.
"""
function plotevents(ax, events; ylims=[], color=:red, alpha=.25)
    if isempty(ylims)
        lines = ax.get_lines()
        if isempty(lines)
            B = 0.0
            T = 1.0
        else
            B = Inf
            T = -Inf
            for line in lines
                y = line.get_ydata()
                B = min(B, minimum(y))
                T = max(T, maximum(y))
            end
        end
    else
        B = ylims[1]
        T = ylims[2]
    end

    n = size(events, 1)
    for i = 1:n
        L = events[i, 1]
        R = events[i, 2]
        ax.fill_between([L, R], [B, B], [T, T], linewidth=0.0, color=color,
            alpha=alpha, zorder=0)
    end

end


function plotevents(events; kwargs...)
    plotevents(plt.gca(), events; kwargs...)
end


"""
`savefigure(fig, args...; kwargs...)`

Save a figure `fig` to the `figures` folder.
"""
function savefigure(fig, args...; kwargs...)
    filename = joinpath(@__DIR__, "figures", args...)
    @assert !isfile(filename) "File already exists!"
    dir, _ = splitdir(filename)
    if !isdir(dir)
        mkpath(dir)
    end
    return fig.savefig(filename; kwargs...)
end


"""
`savecsv(df::DataFrame, args...)`

Save a `DataFrame` as a csv file to the `outputs` folder.
"""
function savecsv(df::DataFrame, args...)
    filename = joinpath(@__DIR__, "outputs", args...)
    @assert !isfile(filename) "File already exists!"
    dir, _ = splitdir(filename)
    if !isdir(dir)
        mkpath(dir)
    end
    return CSV.write(filename, df)
end

"""
`selectwells(df::DataFrame, name)`

Select wells specified by `name`.
"""
function selectwells(df::DataFrame, name)
    df1 = df[:, r"Time|Light"]
    df2 = df[:, name]
    return hcat(df1, df2)
end


"""
`wells2matrix(df::DataFrame)`

Transform wells into from a `DataFrame` into a `Matrix`.
"""
function wells2matrix(df::DataFrame)
    return Matrix(df[:, Not(["Time", "Light"])])
end


"""
`transformwells(df::DataFrame, fun::Function)`

Apply function `fun` to each column of DataFrame `df` except for columns `Time`
and `Light`.
"""
function transformwells(df::DataFrame, fun::Function)
    df = deepcopy(df)
    for name in names(df)
        if name ∉ ("Time", "Light")
            df[:, name] = fun(df[:, name])
        end
    end
    return df
end


"""
`loadcsv(args...)`

Load a `DataFrame` from a csv file located in the `outputs` folder.
"""
function loadcsv(args...)
    filename = joinpath(@__DIR__, "outputs", args...)
    df = DataFrame(CSV.File(filename))
    return df
end


"""
`translatelabel(label)`

Translate column name from the excel file to drug abbreviation and concentration
with units.
"""
function translatelabel(label)
    labels = Dict(
        "Forskolin 5" => "[FOR] = 5 μM",
        "Forskolin 10" => "[FOR] = 10 μM",
        "Forskolin 15" => "[FOR] = 15 μM",
        "DBC 0.5" => "[DBC] = 0.5 mM",
        "DBC 1" => "[DBC] = 1 mM",
        "DBC 3" => "[DBC] = 3 mM",
        "U0126 10" => "[U0126] = 10 μM",
        "U0126 20" => "[U0126] = 20 μM",
        "U0126 40" => "[U0126] = 40 μM",
        "PMA 0.5" => "[PMA] = 0.5 μM",
        "PMA 1" => "[PMA] = 1 μM",
        "PMA 3" => "[PMA] = 3 μM",
        "EGF 30" => "[EGF] = 30 ng/ml",
        "EGF 50" => "[EGF] = 50 ng/ml",
        "EGF 80" => "[EGF] = 80 ng/ml",
        "Ro-318220 2" => "[RO] = 2 μM",
        "Ro-318220 5" => "[RO] = 5 μM",
        "Ro-318220 8" => "[RO] = 8 μM"
    )
    if label in keys(labels)
        new_label = labels[label]
    else
        new_label = label
    end
    return new_label
end

# helper function for creating figures with complex layouts
pyslice(i1, i2) = PyPlot.pycall(PyPlot.pybuiltin("slice"), PyPlot.PyObject, i1, i2)


"""
`dampedsine(t, p)`

Evaluated damped sine at time points `t` for parameters `p = [A, d, T, θ]`
where where `A` is amplitude, `d` damping ratio, `T` period and `θ` phase of
the damped sine.
"""
function dampedsine(t, p)
    A, d, T, θ = p
    x = A .* exp.(-d.*t) .* sin.( (2*pi*t) ./ T .+ θ)
    return x
end


"""
`fitdampedsine(tdata, xdata, p₀)`

Fit damped sine parameters `p = [A, d, T, θ]` to data.
"""
function fitdampedsine(tdata, xdata, p₀)
    lsqfit = LsqFit.curve_fit(dampedsine, tdata, xdata, p₀)
    return lsqfit.param
end


"""
`createevents(X)`

Create event matrix from specifications `X` defined as `X = [(:LL, 2*24),
(:LD, 5, 12, 12), (:DD, 5*24), (:LL, 5*24)]`
"""
function createevents(X)
    current_time = 0.0
    events = Matrix{Float64}(undef, 0, 2)
    for x in X
        code = Symbol(x[1])
        if code == :LD
            for iday = 1:x[2]
                if !isempty(events) && events[end, 2] == current_time
                    events[end, 2] = current_time + x[3]
                else
                    events = [events; [current_time current_time + x[3]]]
                end
                current_time = current_time + x[3] + x[4]
            end
        elseif code == :DD
            duration = x[2]
            current_time += duration
        elseif code == :LL
            if !isempty(events) && events[end, 2] == current_time
                events[end, 2] = current_time + x[2]
            else
                events = [events; [current_time current_time + x[2]]]
            end
            current_time += x[2]
        else
            throw("Unknown code!")
        end
    end
    return events
end

"""
`createcallback(events, i)`

Generate callback for events modyfing `i`-th parameter.
"""
function createcallback(events, i)

    p = NaN  # holds the original value of the parameter

    initialize = function (c, u, t, integrator)
        integrator.p = copy(integrator.p)  # protect from overwriting
        p = integrator.p[i]
        integrator.p[i] = 0.0
    end

    tstops_on = events[:, 1]
    affect_on! = (integrator) -> integrator.p[i] = p
    cb_on = PresetTimeCallback(tstops_on, affect_on!,
        save_positions=(false, false), initialize=initialize)

    tstops_off = events[:, 2]
    affect_off! = (integrator) -> integrator.p[i] = 0.0
    cb_off = PresetTimeCallback(tstops_off, affect_off!,
        save_positions=(false, false))

    return CallbackSet(cb_on, cb_off)

end

"""
`argpeaks(x)`

Return indicies of local maxima of the input vector.
"""
function argpeaks(x)

    # find peak indices
    idx = Array{Int64, 1}()
    for i = 2:length(x)-1
        if x[i-1] < x[i] >= x[i+1]
            push!(idx, i)
        end
    end

    return idx

end


function sustainedperiod(x, t; minamp = 0.01, minlocs = 3, minlocdiffr = 0.02, minpksr = 0.02)

    # find peaks and calculate statistics
    peaks = findpeaks(x, t)
    locs = peaks.locs
    pks = peaks.pks
    amp = maximum(x) - minimum(x)
    m, s = mean_and_std(diff(locs))

    if amp < minamp || length(locs) < minlocs
        return NaN
    end

    pks .-= minimum(x)
    ref_peak = pks[1]
    locs_diff = diff(locs)
    ref_locs_diff =  locs_diff[1]
    loc_condition = !all((ref_locs_diff * (1 - minlocdiffr)) .< locs_diff .< (ref_locs_diff * (1 + minlocdiffr)))
    pks_condition = !all((ref_peak * (1 - minpksr)) .< pks .< (ref_peak * (1 + minpksr)))
    if pks_condition || loc_condition
        return NaN
    end

    return m

end


"""
`addtransientevents(events, pacing, ndays; i=1)`

Add transient events to already existing events. `pacing` is in hours. Returns
`(newevents, offset)` where `offset` is the relative position of the new events
to the old ones.
"""
function addtransientevents(events, pacing, ndays; i=1)
    transientevents = createevents([(:DD, pacing),
        (:LD, ndays, pacing, pacing)])
    if i == 0
        offset_events = 0.0
    else
        offset_events = events[i, 2]
    end
    data_events_select = events[i+1:end, :] .- offset_events .+
    transientevents[end]
    offset = transientevents[end] - offset_events
    newevents = [transientevents; data_events_select]
    return (newevents, offset)
end


"""
`aggregatewells(df::DataFrame, fun::Function, args...)`

Aggregate columns of DataFrame `df` into a single column by applying function
`fun` to each row of the DataFrame. Columns `Time` and `Light` are no affected.
By passing `args` the user can specify which columns to aggregate. E.g.
`r"A|B|C" => "D"` will aggregate columns A, B and C into a new column D.
"""
function aggregatewells(df::DataFrame, fun::Function, args...)

    function mat2vec(df::DataFrame, fun::Function)
        mat = Matrix(df)
        n = size(df, 1)
        vec = fill(NaN, n)
        for i = 1:n
            vec[i] = fun(mat[i, :])
        end
        return vec
    end

    df2 = df[:, ["Time", "Light"]]
    df = df[:, Not(["Time", "Light"])]
    if isempty(args)
        df2[:, Symbol(fun)] = mat2vec(df, fun)
    else
        for arg in args
            df3 = df[:, arg[1]]
            @assert !isempty(df3) "No columns matching the input regex!"
            df2[:, arg[2]] = mat2vec(df3, fun)
        end
    end
    return df2
end


"""
`squarederror(x, y; fun=sum)`

Squared error.
"""
function squarederror(x, y; fun=sum)
    return fun((x .- y).^2)
end


"""
`rsquared(x, y)`

Calculate the coefficient of determination for `x` and `y` where `x` is the real
data value and `y` is the predicted (fitted, modeled) value.
"""
function rsquared(x, y)
    μ = mean(x)
    SSres = sum((x .- y).^2)
    SStot = sum((x .- μ).^2)
    R = 1 - SSres/SStot
    return R
end


"""
ptp(x)

Compute peak to peak amplitude of array `x`.
"""
function ptp(x)
    return maximum(x) - minimum(x)
end


"""
cmean(x)

Compute circular mean of array `x` that represents angles in radians.
"""
function cmean(x)
    S = mean(sin.(x))
    C = mean(cos.(x))
    m = atan(S, C)
    return m
end


"""
cstd(x)

Compute circular standard deviation of array `x` that represents angles in radians.
"""
function cstd(x)
    S = mean(sin.(x))
    C = mean(cos.(x))
    R = sqrt(S^2 + C^2)
    m = sqrt(-2*log(R))
    return m
end