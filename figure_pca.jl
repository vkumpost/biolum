include("functions.jl")

MODEL_NAME = "model06"
RUN_NAME = "run 2 x"

params = Dict(
    :out_folder => "$MODEL_NAME $RUN_NAME",
    :drugs => [
        Dict(
            :name => "Control A",
            :dfarr => [
                loadcsv("$MODEL_NAME $RUN_NAME", "D_control_a.csv")[2:end, 1:end-1]
            ],
            :labels => [
                "Control A"
            ],
            :color => :black,
            :control => "Control A"
        ),
        Dict(
            :name => "Forskolin",
            :dfarr => [
                loadcsv("$MODEL_NAME $RUN_NAME", "D_forskolin_1.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_forskolin_2.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_forskolin_3.csv")[2:end, 1:end-1]
            ],
            :labels => [
                "[FOR] = 5 μM",
                "[FOR] = 10 μM",
                "[FOR] = 15 μM",
                ],
            :color => CB_COLORS[1],
            :control => "Control A"
        ),
        Dict(
            :name => "DBC",
            :dfarr => [
                loadcsv("$MODEL_NAME $RUN_NAME", "D_dbc_1.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_dbc_2.csv")[2:end, 1:end-1],
                # loadcsv("$MODEL_NAME $RUN_NAME", "D_dbc_3.csv")[2:end, 1:end-1]
            ],
            :labels => [
                "[DBC] = 0.5 mM",
                "[DBC] = 1 mM",
                # "[DBC] = 3 mM",
                ],
            :color => CB_COLORS[2],
            :control => "Control A"
        ),
        Dict(
            :name => "U0126",
            :dfarr => [
                loadcsv("$MODEL_NAME $RUN_NAME", "D_u0126_1.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_u0126_2.csv")[2:end, 1:end-1],
                # loadcsv("$MODEL_NAME $RUN_NAME", "D_u0126_3.csv")[2:end, 1:end-1]
            ],
            :labels => [
                "[U0126] = 10 μM",
                "[U0126] = 20 μM",
                # "[U0126] = 40 μM",
                ],
            :color => CB_COLORS[3],
            :control => "Control A"
        ),
        Dict(
            :name => "Control B",
            :dfarr => [
                loadcsv("$MODEL_NAME $RUN_NAME", "D_control_b.csv")[2:end, 1:end-1]
            ],
            :labels => [
                "Control B"
                ],
            :color => :gray,
            :control => "Control B"
        ),
        Dict(
            :name => "EGF",
            :dfarr => [
                loadcsv("$MODEL_NAME $RUN_NAME", "D_egf_1.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_egf_2.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_egf_3.csv")[2:end, 1:end-1]
            ],
            :labels => [
                "[EGF] = 30 ng/ml",
                "[EGF] = 50 ng/ml",
                "[EGF] = 80 ng/ml",
                ],
            :color => CB_COLORS[4],
            :control => "Control B"
        ),
        Dict(
            :name => "PMA",
            :dfarr => [
                # loadcsv("$MODEL_NAME $RUN_NAME", "D_pma_1.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_pma_2.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_pma_3.csv")[2:end, 1:end-1]
            ],
            :labels => [
                # "[PMA] = 0.5 μM",
                "[PMA] = 1 μM",
                "[PMA] = 3 μM",
                ],
            :color => CB_COLORS[5],
            :control => "Control B"
        ),
        Dict(
            :name => "RO",
            :dfarr => [
                # loadcsv("$MODEL_NAME $RUN_NAME", "D_ro_1.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_ro_2.csv")[2:end, 1:end-1],
                loadcsv("$MODEL_NAME $RUN_NAME", "D_ro_3.csv")[2:end, 1:end-1]
            ],
            :labels => [
                # "[RO] = 2 μM",
                "[RO] = 5 μM",
                "[RO] = 8 μM",
                ],
            :color => CB_COLORS[6],
            :control => "Control B"
        ),
    ]
)

## pca_fit! ===================================================================
function pca_fit!(params)

    drugs = params[:drugs]

    master_df = DataFrame()
    
    for drug in drugs
        dfarr = drug[:dfarr]
        labels = drug[:labels]
        for (df, label) in zip(dfarr, labels)
            master_df = vcat(master_df, df)
        end
    end
    M = Matrix(master_df)
    M_zscore = mapslices(zscore, M; dims=1)
    pca = MS.fit(MS.PCA, transpose(M_zscore); pratio = 1.0)

    params[:master_df] = master_df
    params[:pca] = pca
    params[:principal_components] = MS.transform(pca, transpose(M_zscore))
    params[:projection] = MS.projection(pca)

    for i = 1:length(pca.prinvars)
        exp_var = round(Int64, pca.prinvars[i]/pca.tvar*100)
        println("PC$i: $exp_var % explained variance")
    end

end
pca_fit!(params)

## pca_plot ===================================================================
function pca_plot(params)

    drugs = params[:drugs]
    principal_components = params[:principal_components]
    pca = params[:pca]
    out_folder = params[:out_folder]

    fig, ax = subplots(figsize=(7, 7))
    cluster_counter = 0
    i1 = 1  # principal component on the x-axis
    i2 = 2  # principal component on the y-axis

    for drug in drugs
        dfarr = drug[:dfarr]
        labels = drug[:labels]
        color = drug[:color]
        markers = ["o", "x", "+"]

        for i = 1:length(dfarr)        

            df = dfarr[i]
            label = labels[i]
            marker = markers[i]
            markersize = 3

            start_idx = cluster_counter + 1
            end_idx = cluster_counter + size(df, 1)
            cluster_counter += size(df, 1)
            
            ax.plot(
                principal_components[i1, start_idx:end_idx],
                principal_components[i2, start_idx:end_idx];
                linestyle = "",
                marker = marker,
                markersize = markersize,
                color = color,
                label = label,
                fillstyle = "none")
        end

    end    

    ax.legend(ncol=2, fancybox=false, edgecolor=:black)

    aaa = round(Int64, pca.prinvars[i1]/pca.tvar*100)
    bbb = round(Int64, pca.prinvars[i2]/pca.tvar*100)
    ax.set_xlabel("PC$i1 ($aaa % explained variance)")
    ax.set_ylabel("PC$i2 ($bbb % explained variance)")
    ax.set_title("PCA plot")
    fig.show()
    savefigure(fig, joinpath(out_folder, "pca_clusters.svg"))
    
end

pca_plot(params)


## pca_paths ==================================================================
function pca_paths(params)

    drugs = params[:drugs]
    principal_components = params[:principal_components]
    pca = params[:pca]
    out_folder = params[:out_folder]

    fig, ax = subplots(figsize=(5, 5))
    cluster_counter = 0
    i1 = 1  # principal component on the x-axis
    i2 = 2  # principal component on the y-axis

    control_memory = Dict()
    for drug in drugs
        name = drug[:name]
        dfarr = drug[:dfarr]
        color = drug[:color]
        labels = drug[:labels]
        control = drug[:control]
        path = []

        for i = 1:length(dfarr)

            df = dfarr[i]
            marker = "o"
            label = labels[i]

            start_idx = cluster_counter + 1
            end_idx = cluster_counter + size(df, 1)
            cluster_counter += size(df, 1)

            x_m = median(principal_components[i1, start_idx:end_idx])
            y_m = median(principal_components[i2, start_idx:end_idx])
            x_mad = mad(principal_components[i1, start_idx:end_idx])
            y_mad = mad(principal_components[i2, start_idx:end_idx])

            if control == label
                control_memory[control] = [x_m, y_m]
                zorder = 100
            else
                push!(path, [x_m, y_m])
                zorder = 1
            end

            ax.scatter(x_m, y_m, color = color, zorder = zorder, marker = marker)
            ax.hlines(y_m, x_m-x_mad, x_m+x_mad, color = color, alpha=0.5, zorder = zorder)
            ax.vlines(x_m, y_m-y_mad, y_m+y_mad, color = color, alpha=0.5, zorder = zorder)
        end

        x = [control_memory[control][1], [x[1] for x in path]...]
        y = [control_memory[control][2], [x[2] for x in path]...]
        ax.plot(x, y, color = color, zorder = -1, label = name)


    end    
    ax.legend(ncol=1, fancybox=false, edgecolor=:black)

    aaa = round(Int64, pca.prinvars[i1]/pca.tvar*100)
    bbb = round(Int64, pca.prinvars[i2]/pca.tvar*100)
    ax.set_xlabel("PC$i1 ($aaa % explained variance)", labelpad=0)
    ax.set_ylabel("PC$i2 ($bbb % explained variance)", labelpad=0)
    ax.set_title("PCA plot", pad=0, loc="left")
    fig.tight_layout()

    fig.show()
    savefigure(fig, joinpath(out_folder, "pca_paths.svg"))

end

pca_paths(params)


## pca_projections ============================================================
function pca_projections(params)

    pca = params[:pca]
    master_df = params[:master_df]
    projection = params[:projection]
    out_folder = params[:out_folder]

    fig, ax = subplots(figsize=(3, 2))
    i1 = 1
    i2 = 2

    parameter_names = names(master_df)
    for i = 1:size(projection, 1)
        x = projection[i, i1]
        y = projection[i, i2]
        ax.arrow(0.0, 0.0, x, y, head_width=0.05, color=:black, zorder = -100)
        ax.text(x, y, parameter_names[i], color=:black, zorder = -100)
    end

    e = 0.1
    ax.set_xlim(minimum(projection[:, i1]) - e, maximum(projection[:, i1]) + e)
    ax.set_ylim(minimum(projection[:, i2]) - e, maximum(projection[:, i2]) + e)

    aaa = round(Int64, pca.prinvars[i1]/pca.tvar*100)
    bbb = round(Int64, pca.prinvars[i2]/pca.tvar*100)
    ax.set_title("Model parameter projections", pad=0, loc="left")

    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.85, 0.85)

    fig.tight_layout()

    savefigure(fig, joinpath(out_folder, "pca_projections.svg"))

end

pca_projections(params)
