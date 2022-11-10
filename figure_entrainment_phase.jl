include("functions.jl")

filenames = [
    ("model06 run 2 x", ["T_control_a.csv"], 
    "Control", "black"),
    ("model06 run 2 x", ["T_dbc_1.csv"], 
    "[DBC] = 0.5 mM", CB_COLORS[2]),
    ("model06 run 2 x", ["T_U0126_2.csv"], 
    "[U0126] = 20 μM", CB_COLORS[3])
]
output_filename = "fig4_T.svg"

# filenames = [
#     ("model06 run 2 x", ["T_control_a.csv"], 
#     "Control A", "black"),
#     ("model06 run 2 x", ["T_dbc_1.csv"], 
#     "[DBC] = 0.5 mM", CB_COLORS[1]),
#     ("model06 run 2 x", ["T_dbc_2.csv"], 
#     "[DBC] = 1 mM", CB_COLORS[2]),
# ]
# output_filename = joinpath("model06 run 2 x", "T_dbc.svg")

# filenames = [
#     ("model06 run 2 x", ["T_control_a.csv"], 
#     "Control A", "black"),
#     ("model06 run 2 x", ["T_U0126_1.csv"], 
#     "[U0126] = 10 μM", CB_COLORS[1]),
#     ("model06 run 2 x", ["T_U0126_2.csv"], 
#     "[U0126] = 20 μM", CB_COLORS[2]),
# ]
# output_filename = joinpath("model06 run 2 x", "T_u0126.svg")

# filenames = [
#     ("model06 run 2 x", ["T_control_a.csv"], 
#     "Control A", "black"),
#     ("model06 run 2 x", ["T_forskolin_1.csv"], 
#     "[FOR] = 5 μM", CB_COLORS[1]),
#     ("model06 run 2 x", ["T_forskolin_2.csv"], 
#     "[FOR] = 10 μM", CB_COLORS[2]),
#     ("model06 run 2 x", ["T_forskolin_3.csv"], 
#     "[FOR] = 15 μM", CB_COLORS[3]),
# ]
# output_filename = joinpath("model06 run 2 x", "T_forskolin.svg")

# filenames = [
#     ("model06 run 2 x", ["T_control_b.csv"], 
#     "Control B", "black"),
#     ("model06 run 2 x", ["T_egf_1.csv"], 
#     "[EGF] = 30 ng/ml", CB_COLORS[1]),
#     ("model06 run 2 x", ["T_egf_2.csv"], 
#     "[EGF] = 50 ng/ml", CB_COLORS[2]),
#     ("model06 run 2 x", ["T_egf_3.csv"], 
#     "[EGF] = 80 ng/ml", CB_COLORS[3]),
# ]
# output_filename = joinpath("model06 run 2 x", "T_egf.svg")

# filenames = [
#     ("model06 run 2 x", ["T_control_b.csv"], 
#     "Control B", "black"),
#     ("model06 run 2 x", ["T_pma_2.csv"], 
#     "[PMA] = 1 μM", CB_COLORS[2]),
#     ("model06 run 2 x", ["T_pma_3.csv"], 
#     "[PMA] = 3 μM", CB_COLORS[3]),
# ]
# output_filename = joinpath("model06 run 2 x", "T_pma.svg")

## ==

N = length(filenames)
PRC_arr = Vector{DataFrame}(undef, N)
label_arr = Vector{String}(undef, N)
color_arr = Vector{String}(undef, N)
for i = 1:length(filenames)
    folder = filenames[i][1]
    filename_arr = filenames[i][2]
    PRC_tmp_arr = [loadcsv(joinpath(folder, x)) for x in filename_arr]
    PRC = DataFrame(times = Float64[], shift_mean = Float64[], shift_std = Float64[])
    for j = 1:size(PRC_tmp_arr[1], 1)
        global PRC_mean, PRC_std
        t = PRC_tmp_arr[1][j, 1]
        m = mean([x[j, 2] for x in PRC_tmp_arr])
        v = std([x[j, 2] for x in PRC_tmp_arr])
        push!(PRC, [t, m, v])
    end
    PRC_arr[i] = PRC
    label_arr[i] = filenames[i][3]
    color_arr[i] = filenames[i][4]
end


## ==
close(:all)
fig, ax = subplots(figsize = (4, 3))
markerarr = ["o", "v", "^", "s"]
for i = 1:N
    PRC = PRC_arr[i]
    label = label_arr[i]
    X = PRC[:, :times]
    Y = PRC[:, :shift_mean]
    ax.plot(X, Y; marker = markerarr[i], color = color_arr[i], label = label)
end

ax.set_xlabel("LD cycle period (hours)", labelpad=0)
ax.set_ylabel("Phase angle (hours)", labelpad=0)
ax.set_title("Phase of entrainment", pad=0, loc="left")
ax.legend(fancybox=true, edgecolor=:black, framealpha=1.0)
ax.set_xticks(PRC_arr[1][:, :times])
fig.tight_layout(pad=0)

savefigure(fig, output_filename)