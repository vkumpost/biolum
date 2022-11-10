include("functions.jl")

MODEL_NAME = "model04"
RUN_NAME = "run 2"

filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "PRC_C_parameters.csv")
fig_filename = joinpath("$(MODEL_NAME) $(RUN_NAME)", "PRC_C_parameters.svg")    

PRC = loadcsv(filename)

fig, ax = subplots(figsize = (4, 3))
ax.plot(PRC[:, :times], PRC[:, :shift], marker = "o", color = :black)
ax.set_xlabel("Phase (hours)", labelpad=0)
ax.set_ylabel("Phase shift (hours)", labelpad=0)
ax.set_title("Phase response curve", pad=0, loc="left")
ax.hlines(0, -1, 25; color = :gray, label = "no light", zorder = -1000)
ax.set_xticks(0:4.0:24)
ax.set_yticks(-12:4.0:12)
ax.set_xlim(-1.0, 25.0)
ax.set_ylim(-13, 13)
fig.tight_layout()
fig.show()
savefigure(fig, fig_filename)
