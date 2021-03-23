include("functions.jl")

platenames = ["Plate D1 A", "Plate D2 A", "Plate D1 B", "Plate D2 B"]

for platename in platenames

    # Select wells corresponding to the desired drugs
    if platename[end] == 'A'
        drugs = ["Control", "Forskolin 5", "Forskolin 10", "Forskolin 15",
            "DBC 0.5", "DBC 1", "DBC 3", "U0126 10", "U0126 20", "U0126 40"]
    elseif platename[end] == 'B'
        drugs = ["Control", "PMA 0.5", "PMA 1", "PMA 3", "Ro-318220 2",
            "Ro-318220 5", "Ro-318220 8", "EGF 30", "EGF 50", "EGF 80"]
    end

    ndrug = length(drugs)
    reg = Regex(join(drugs, "|"))
    
    dataraw = loaddata(platename)
    data = selectwells(dataraw, reg)

    # Select times corresponding to the desired time interval
    if occursin("D1", platename)
        tspan = [13.0, 106.0]
    elseif occursin("D2", platename)
        tspan = [13.0, 79.0]
    end

    offset = detectevents(data)[1, 2]
    idx = tspan[1] .< data.Time .< tspan[2]
    datasel = data[idx, :]
    datasel.Time .-= offset

    # Normalize data
    params = Dict()
    datanorm = deepcopy(datasel)
    for i = 1:ndrug
        datatmp = selectwells(datasel, Regex(drugs[i]))
        X = wells2matrix(datatmp)
        params[drugs[i]] = [mean(mean(X, dims=1)), mean(std(X, dims=1))]

        mean_control = params["Control"][1]
        std_control = params["Control"][2]
        mean_drug = params[drugs[i]][1]
        std_drug = params[drugs[i]][2]

        mean_ratio = (mean_drug - mean_control) / std_control
        std_ratio = std_drug / std_control
        fun = x -> zscore(x) .* std_ratio .+ mean_ratio

        datanorm[:, Regex(drugs[i])] = wells2matrix(transformwells(datatmp, fun))
    end

    # Save data
    outputplatename = lowercase(replace(platename, " " => "_"))
    savecsv(datasel, "raw_$(outputplatename).csv")
    savecsv(datanorm, "norm_$(outputplatename).csv")
end