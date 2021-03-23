include("functions.jl")

small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

MODEL_NAME = "model06"
RUN_NAME = "run 2 x"

## Fig 3A =====================================================================
drug_fullname = ""
drug_names = ["Control", "DBC 0.5",  "U0126 20"]
drug_filenames = [
    joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv"),
    joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_1.csv"),
    joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_2.csv"),
]
plate_set = "A"
output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "fit_example_a.svg")
colorarr = [:black, CB_COLORS[2], CB_COLORS[3]]

## Fig 3B =====================================================================
# drug_fullname = ""
# drug_names = ["Control", "EGF 30",  "Ro-318220 5"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_2.csv"),
# ]
# plate_set = "B"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "fit_example_b.svg")
# colorarr = [:gray, CB_COLORS[4], CB_COLORS[6]]

## S3 Fig A ===================================================================
# drug_fullname = "Forskolin"
# drug_names = ["Control", "Forskolin 5",  "Forskolin 10", "Forskolin 15"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_2.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_forskolin_3.csv"),
# ]
# plate_set = "A"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "drug_fit_for.svg")
# colorarr = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

## S3 Fig B ===================================================================
# drug_fullname = "DBC"
# drug_names = ["Control", "DBC 0.5", "DBC 1", "DBC 3"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_2.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_dbc_3.csv"),
# ]
# plate_set = "A"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "drug_fit_dbc.svg")
# colorarr = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

## S3 Fig C ===================================================================
# drug_fullname = "U0126"
# drug_names = ["Control", "U0126 10",  "U0126 20", "U0126 40"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_a.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_2.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_u0126_3.csv"),
# ]
# plate_set = "A"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "drug_fit_u0126.svg")
# colorarr = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

## S3 Fig D ===================================================================
# drug_fullname = "EGF"
# drug_names = ["Control", "EGF 30",  "EGF 50", "EGF 80"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_2.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_egf_3.csv"),
# ]
# plate_set = "B"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "drug_fit_egf.svg")
# colorarr = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

## S3 Fig E ===================================================================
# drug_fullname = "PMA"
# drug_names = ["Control", "PMA 0.5",  "PMA 1", "PMA 3"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_2.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_pma_3.csv"),
# ]
# plate_set = "B"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "drug_fit_pma.svg")
# colorarr = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

## S3 Fig F ===================================================================
# drug_fullname = "Ro-318220"
# drug_names = ["Control", "Ro-318220 2",  "Ro-318220 5", "Ro-318220 8"]
# drug_filenames = [
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_control_b.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_1.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_2.csv"),
#     joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_ro_3.csv"),
# ]
# plate_set = "B"
# output_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "drug_fit_ro.svg")
# colorarr = [:black, CB_COLORS[1], CB_COLORS[2], CB_COLORS[3]]

params = Dict(
    :model_name => MODEL_NAME,
    :showplots => false,
    :plate_set => plate_set,
    :sde_alg => DRI1(),  # SDE solver
    :trajectories => 1000,
    :colorarr => colorarr,
    
    # input files
    :in_initial_conditions_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_initial_conditions.csv"),
    :in_normalization_parameters_filename => joinpath("$(MODEL_NAME) $(RUN_NAME)", "D_norm_parameters_$(lowercase(plate_set)).csv"),  # empty string for Control
    
)

paramsarr = []
for (name, filename) in zip(drug_names, drug_filenames)
    global paramsarr, params
    tmp = deepcopy(params)
    tmp[:drug] = name
    tmp[:in_initial_population_filename] = filename
    push!(paramsarr, tmp)
end

# Create the figure
fig = figure(figsize=(4, 3), constrained_layout=true)
gs = fig.add_gridspec(2, 2; width_ratios=[100, 1])
ax1 = fig.add_subplot(get(gs, (0, pyslice(0, 2))))
ax2 = fig.add_subplot(get(gs, (1, 0)))
axarr = [ax1, ax2]

handles_data = []
handles_model = []
println("Drug & R2D1 & R2D2")
for i = 1:length(paramsarr)
    global axarr, handles_data,handles_model
    local params
    params = paramsarr[i]
    D_load_initial_population!(params)
    D_load_normalization_parameters!(params)
    D_load_training_data!(params)
    D_add_sde_models!(params)
    D_load_initial_conditions!(params)

    drug = params[:drug]
    colorarr = params[:colorarr]
    models = params[:models]
    pnames = params[:pnames]
    p = params[:p]
    trajectories = params[:trajectories]
    u0arr = params[:u0arr]
    sde_alg = params[:sde_alg]
    normpars = params[:normpars]
    xdataarr = params[:xdataarr]
    eventsarr = params[:eventsarr]
    optmodel01 = setparameters(models[1], pnames, p)
    tsim01, xsim01 = simulatepopulation2(optmodel01, sde_alg, trajectories; u0arr=u0arr)
    optmodel02 = setparameters(models[2], pnames, p)
    tsim02, xsim02 = simulatepopulation2(optmodel02, sde_alg, trajectories; u0arr=u0arr)

    xsim01 = zscore(xsim01, normpars[1], normpars[2])
    xsim02 = zscore(xsim02, normpars[3], normpars[4])

    axarr[1].plot(tsim01, xsim01; color = colorarr[i])
    axarr[1].plot(tsim01, xdataarr[1]; color = colorarr[i], linestyle = "--")
    h_model = axarr[2].plot(tsim02, xsim02; color = colorarr[i])
    h_data = axarr[2].plot(tsim02, xdataarr[2]; color = colorarr[i], linestyle = "--")

    axarr[1].set_xticks(0:24:tsim01[end])
    axarr[2].set_xticks(0:24:tsim01[end])

    if i == length(paramsarr)
        plotevents(axarr[1], eventsarr[1])
    end
    
    push!(handles_data, h_data[1])
    push!(handles_model, h_model[1])

    # R2D1 = round(cor(xsim01, xdataarr[1])^2; digits = 2)
    # L2D1 = round(squarederror(xsim01, xdataarr[1]; fun=mean); digits = 2)
    # R2D2 = round(cor(xsim02, xdataarr[2])^2; digits = 2)
    # L2D2 = round(squarederror(xsim02, xdataarr[2]; fun=mean); digits = 2)
    R2D1 = round(rsquared(xsim01, xdataarr[1]); digits = 2)
    R2D2 = round(rsquared(xsim02, xdataarr[2]); digits = 2)

    println("$drug & $R2D1 & $R2D2")
end

axarr[1].autoscale(tight=true)
axarr[2].autoscale(tight=true)
if isempty(drug_fullname)
    axarr[1].set_title("LD cycle")
    axarr[2].set_title("Constant darkness")
else
    axarr[1].set_title(drug_fullname * " - LD cycle")
    axarr[2].set_title(drug_fullname * " - constant darkness")
end
axarr[1].set_xlabel("Time (hours)")
axarr[2].set_xlabel("Time (hours)")
axarr[1].set_ylabel("Lumin. (au)")
axarr[2].set_ylabel("Lumin. (au)")


handles = [handles_model[1], handles_data[1], handles_model...]
labels = ["Model", "Data", translatelabel.(drug_names)...]
axarr[2].legend(handles, labels,
                    loc="upper left", bbox_to_anchor=(1, 1.0), fancybox=false,
                    edgecolor=:black)
 
savefigure(fig, output_filename)