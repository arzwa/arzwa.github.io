# This file was generated, do not modify it. # hide
using Distributions, Plots, StatsPlots # hide
α = 1.0
K = 3
d = Gamma(α, 1/α)  # the Gamma distribution object
q = quantile(d, 1/K:1/K:1)  # the discretization points (defining the classes)
m = quantile(d, 1/2K:1/K:1-1/2K)  # the medians in each class
@show round.(m, digits=2)
plot(d, grid=false, legend=false, size=(500,200), xlim=(0,5), xlabel="relative rate", ylabel="probability density", guidefont=8)
vline!(q, color=:black)
vline!(m, linestyle=:dot)

savefig("_assets/phylocourse/mliqtree/gamma.svg") # hide