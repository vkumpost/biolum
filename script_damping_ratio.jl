include("functions.jl")

data = loaddata("Plate U1")

t = data[:, "Time"]
df = select(data, r"WT [AB]")
events = detectevents(data)

idx = events[4, 2] + 12.0 .< t .< events[5, 1] - 24.0

df2 = DataFrame()
for name in names(df)
    col = df[idx, name]
    df2[:, name] = zscore(Float64.(col))
end
t2 = t[idx] .- t[idx][1]
x2 = dropdims(mean(Matrix(df2), dims=2), dims=2)

p = fitdampedsine(t2, x2, [2.5, 0.01, 25.0, 0.0])
xest = dampedsine(t2, p)

R = cor(x2, xest) ^ 2

plot(t2, Matrix(df2), color=:gray)
plot(t2, x2, color=:black, linewidth=2)
plot(t2, xest, color=:red, linewidth=2)

df = DataFrame([[p[1]], [p[2]], [p[3]], [p[4]], [R]],
    ["A", "d", "Period", "Phase", "R2"])
savecsv(df, "damped_sine_parameters.csv")

println()
println("Damped sine fit:")
println("- Amplitude: $(p[1])")
println("- Damping ration: $(p[2])")
println("- Period: $(p[3])")
println("- Phase: $(p[4])")
println()
println("Quality of fit:")
println("- Coefficient of determination: $R")
