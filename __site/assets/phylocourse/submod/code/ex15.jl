# This file was generated, do not modify it. # hide
using Plots
plot(first.(l), last.(l),
    xlabel="p", ylabel="P(data|p)",
    grid=false, legend=false, color=:black)
savefig("_assets/phylocourse/submod/lhood1.svg") # hide