include("functions.jl")

dir = "model06 run 2 x"

filenames = [
    "D_control_a.csv", "D_control_b.csv",
    "D_forskolin_1.csv", "D_forskolin_2.csv", "D_forskolin_3.csv",
    "D_dbc_1.csv", "D_dbc_2.csv", # "D_dbc_3.csv",
    "D_u0126_1.csv", "D_u0126_2.csv", # "D_u0126_3.csv",
    "D_egf_1.csv", "D_egf_2.csv", "D_egf_3.csv",
    "D_pma_2.csv", "D_pma_3.csv",  # "D_pma_1.csv", 
    "D_ro_2.csv", "D_ro_3.csv"  # "D_ro_1.csv", 
]

labels = [
    "Control A", "Control B",
    "Forskolin 5", "Forskolin 10", "Forskolin 15",
    "DBC 0.5", "DBC 1", # "DBC 3",
    "U0126 10", "U0126 20", # "U0126 40",
    "EGF 30", "EGF 50", "EGF 80",
    "PMA 1", "PMA 3", # "PMA 0.5", 
    "Ro-318220 5", "Ro-318220 8" # "Ro-318220 2", 
]

colors = [
    "black", "gray",
    CB_COLORS[1], CB_COLORS[1], CB_COLORS[1],
    CB_COLORS[2], CB_COLORS[2],
    CB_COLORS[3], CB_COLORS[3],
    CB_COLORS[4], CB_COLORS[4], CB_COLORS[4],
    CB_COLORS[5], CB_COLORS[5],
    CB_COLORS[6], CB_COLORS[6]
]

n_filenames = length(filenames)
df_arr = Vector{DataFrame}(undef, n_filenames)

for i_filename = 1:n_filenames
    global filenames, df_arr, dir

    filename = filenames[i_filename]
    df_arr[i_filename] = loadcsv(dir, filename)

end

close(:all)
small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

fig, axs = subplots(2, 3; figsize=(8, 3))

axs = [axs...;]

parameter_names = names(df_arr[1])[1:6]

box_labels = translatelabel.(labels)
for i = 1:6
    x = [x[2:end, parameter_names[i]] for x in df_arr]
    x = [log.(x) for x in x]
    medians = [median(x) for x in x]
    x_max = [maximum(x) for x in x]
    x_min = [minimum(x) for x in x]
    for ind = 1:length(medians)
        axs[i].scatter(ind, medians[ind], marker="o", color=colors[ind], s=20)
        axs[i].vlines(ind, x_min[ind], x_max[ind], color=colors[ind], linestyle="-", lw=2)
    end

    yticks_arr = [round(maximum([x...;]), digits=1), round(minimum([x...;]), digits=1)]
    vlines_arr = [2.5, 5.5, 7.5, 9.5, 12.5, 14.5]
    for j = 1:6
        axs[i].vlines(vlines_arr[j], yticks_arr[1], yticks_arr[2], color=:gray, linestyle="--", lw = 1)
    end
    
    axs[i].set_title(parameter_names[i])
    axs[i].set_xticks([1.5, 4, 6.5, 8.5, 11, 13.5, 15.5])
    axs[i].set_xticklabels(["Control", "FOR", "DBC", "U0126", "EGF", "PMA", "RO"], rotation=25, ha="right")
    axs[i].set_yticks(yticks_arr)

end

fig.tight_layout()
savefigure(fig, dir, "parameter_values.svg")
