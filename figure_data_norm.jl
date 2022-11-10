include("functions.jl")

data_raw = loaddata("Plate U1")
data_sel = selectwells(data_raw, r"WT [A|B]")
events = detectevents(data_sel)

X_all = wells2matrix(data_sel)
t_all = data_sel[:, "Time"]

first_event = 4
last_event = 8
idx = events[first_event, 2] .< t_all .< events[last_event, 2]
X = X_all[idx, :]
t = vec(t_all[idx, :] .- t_all[idx, :][1])
events = events[first_event+1:last_event, :] .- t_all[idx, :][1]

fig, (ax1, ax2, ax3) = subplots(3; figsize = (7, 4.5))
ax1.plot(t, X; color = :black)
plotevents(ax1, events)
ax1.autoscale(tight = true)
ax1.set_ylabel("Lumin. (cps)", labelpad=0)
ax1.set_xlabel("Time (hours)", labelpad=0)
ax1.set_title("Raw bioluminescence data", pad=0, loc="left")
ax1.set_xticks(0:48:t[end])

X_norm = mapslices(zscore, X; dims = 1)
ax2.plot(t, X_norm; color = :black)
plotevents(ax2, events)
ax2.autoscale(tight = true)
ax2.set_ylabel("Lumin. (au)", labelpad=0)
ax2.set_xlabel("Time (hours)", labelpad=0)
ax2.set_title("Normalized well traces", pad=0, loc="left")
ax2.set_xticks(0:48:t[end])

n = size(X_norm, 1)
distributions = Vector{Distributions.UnivariateDistribution}(undef, n)
for j = 1:n
    distributions[j] = Distributions.fit(Distributions.Normal, X_norm[j, :])
end

ax3.plot(t, [d.μ for d in distributions], color=:black)
ax3.fill_between(t, [d.μ + 3*d.σ for d in distributions],
    y2=[d.μ - 3*d.σ for d in distributions], color = :gray)
plotevents(ax3, events;  ylims=[-1000, 1000])
ax3.autoscale(tight = true)
ax3.set_ylim(minimum([d.μ - 3*d.σ for d in distributions]),
    maximum([d.μ + 3*d.σ for d in distributions]))
ax3.set_ylabel("Lumin. (au)", labelpad=0)
ax3.set_xlabel("Time (hours)", labelpad=0)
ax3.set_title("Mean of wells", pad=0, loc="left")
ax3.legend(["Mean", "3*SD"], edgecolor=:black, framealpha=1, ncol=2)
ax3.set_xticks(0:48:t[end])

fig.tight_layout(pad=0.3)

fig.show()
savefigure(fig, "data_norm.svg")
