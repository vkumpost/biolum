using PyPlot
import Base.isempty, Base.length, Base.getindex

struct Peaks
    pks::Vector  # Local maxima
    locs::Vector  # Peak locations
    w::Vector  # Peak widths
    p::Vector  # Peak prominences
end

function isempty(peaks::Peaks)
    return isempty(peaks.locs)
end

function length(peaks::Peaks)
    return length(peaks.locs)
end

function getindex(peaks::Peaks, idx)
    pks = peaks.pks[idx]
    locs = peaks.locs[idx]
    w = peaks.w[idx]
    p = peaks.p[idx]
    return Peaks(pks, locs, w, p)
end

"""
`findpeaks`

**Arguments**
- `data`: input signal vector.

**Optional Arguments**
- `x`: location vector, `locs` and `w` are expressed in terms of `x`.

**Keyword Arguments**
- `sortstr`: Peak sorting. The possible values are
    - `"none"`: no sorting (default).
    - `"ascend"`: increasing order, from the smallest to the largest value.
    - `"descend"`: descending order, from the largest to the smallest value.
- `minp`: minimum peak prominence. Default is `0.0`.
- `showplot`: show found peaks and their properties in a plot. Default is `false`.

**Returns**
- Struct `Peaks` holding the array of all found local maxima and their 
    properties. The struct contains fields
    - `pks`: loacal maxima.
    - `locs`: peak locations.
    - `w`: Peak widths.
    - `p`: peak prominence.
"""
function findpeaks(data; minp = 0.0, sortstr = "none", showplot = false)
    @assert all(.!isnan.(data)) "Data contains NaN values."
    @assert all(.!isinf.(data)) "Data contains Inf values."  # TO-DO

    locs = findlocations(data)

    if isempty(locs)
        return Peaks([], [], [], [])  # to-do correct type
    end

    pks = data[locs]
    p = findprominences(data, locs)
    w = findwidths(data, locs, p)

    idx = p .> minp
    pks = pks[idx]
    locs = locs[idx]
    p = p[idx]
    w = w[idx, :]

    if sortstr != "none"
        if sortstr == "ascend"
            idx = sortperm(pks)
        elseif sortstr == "descend"
            idx = sortperm(pks; rev = true)
        end
        locs = locs[idx]
        pks = pks[idx]
        p = p[idx]
        w = w[idx, :]
    end

    if showplot
        fig, ax = subplots()
        x = 1:length(data)
        ax.plot(x, data; color = :black, label = "Data")
        ax.plot(locs, pks, "v"; color = :red, label = "Peaks")
        for i = 1:length(locs)
            line_p = ax.vlines(locs[i], pks[i] - p[i], pks[i]; color = :orange)
            line_w = ax.hlines(pks[i] - p[i]/2, w[i, 1], w[i, 2]; color = :gray)
            if i == 1
                line_p.set_label("Prominence")
                line_w.set_label("Width")
            end
        end
        ax.legend()
    end

    return Peaks(pks, locs, diff(w; dims = 2)[:], p)
end


function findpeaks(data, x; args...)
    peaks = findpeaks(data; args...)
    return Peaks(peaks.pks, x[peaks.locs], peaks.w, peaks.p)
end


function findlocations(data)
    n = length(data)  # length of data
    locs = Int[]
    loc_candidate = 0
    for i = 2:n-1
        if data[i-1] < data[i]
            loc_candidate = i
        end
        if loc_candidate > 0 && data[i] > data[i+1]
            push!(locs, loc_candidate)
            loc_candidate = 0
        end
    end
    return locs
end


# function findwidths(data, locs)
#     nlocs = length(locs)
#     w = Vector{Int}(undef, 5)
#     for iloc = 1:nlocs

# end


function findprominences(data, locs)
    
    nlocs = length(locs)
    p = typeof(data)(undef, nlocs)
    for iloc = 1:nlocs

        loc = locs[iloc]
        marker = data[loc]

        # run left
        left_index = 1
        for i = iloc-1:-1:1
            if marker < data[locs[i]]
                left_index = locs[i]
                break                
            end
        end
        left_minimum = minimum(data[left_index:locs[iloc]])

        # run right
        right_index = length(data)
        for i = iloc+1:1:nlocs
            if marker < data[locs[i]]
                right_index = locs[i]
                break                
            end
        end
        right_minimum = minimum(data[locs[iloc]:right_index])

        reference_level = max(left_minimum, right_minimum)
        p[iloc] = marker - reference_level

    end    

    return p
    
end


function findwidths(data, locs, p)

    nlocs = length(locs)
    w = Matrix{Float64}(undef, nlocs, 2)
    for iloc = 1:nlocs

        loc = locs[iloc]
        pk = data[loc]  # peak
        wr = pk - p[iloc] * 0.5  # reference height for width measurements
        # println(wr)

        # run left
        left_index = 1
        for i = loc:-1:2
            x1 = data[i]
            x2 = data[i-1]
            if x2 < wr < x1
                # left_index = i
                left_index = i - (i - i + 1) / (x2 - x1) * (wr - x1)
                break                
            end
        end

        # run right
        right_index = length(data)
        for i = loc:1:length(data)-1
            x1 = data[i]
            x2 = data[i+1]
            if x2 < wr < x1
                # right_index = i
                right_index = (i - i + 1) / (x2 - x1) * (wr - x1) + i
                break
            end
        end

        w[iloc, :] = [left_index, right_index]

    end

    return w

end


# x = 0:0.1:8π
# data = sin.(x)
# plot(x, data)
# locs = findlocations(data)
# plot(x[locs], data[locs], "o")

# p = findprominences(data, locs)
# data = rand(100)
# data = [1, 2, 3, 4, 4, 4, 4, 2, 3, 2, 1, 1, 0, 0]
# peaks = findpeaks(data; sortstr = "ascend", showplot = true)
# peaks = findpeaks(data; sortstr = "descend", showplot = true)