#=
First, you must run `script_prc.jl` to estimate the PRC.
=#

include("functions.jl")

small_font = 8
medium_font = 9
big_font = 10
rc("font", family="arial", size=small_font)
rc("axes", titlesize=big_font, labelsize=medium_font)
rc("xtick", labelsize=small_font) 
rc("ytick", labelsize=small_font)

MODEL_NAME = "model04"
RUN_NAME = "run 2"

filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "PRC_C_parameters.csv")
fig_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "PRC_C_parameters.svg")    

df = loadcsv(filename)

fig, ax = subplots(figsize = (3.5, 2.5))
ax.plot(PRC[:, :times], PRC[:, :shift], marker = ".", color = :black)
ax.set_xlabel("Phase (hours)")
ax.set_ylabel("Phase shift (hours)")
ax.set_title("Phase response curve")
ax.hlines(0, -1, 25; color = :gray, label = "no light", zorder = -1000)
ax.set_xticks(0:4.0:24)
ax.set_yticks(-12:4.0:12)
ax.set_xlim(-1.0, 25.0)
ax.set_ylim(-13, 13)
fig.tight_layout()

savefigure(fig, fig_filename)