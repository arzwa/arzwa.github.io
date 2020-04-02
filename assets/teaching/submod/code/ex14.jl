# This file was generated, do not modify it. # hide
using Plots
plot(probs[1,:], probs[2,:],
    xlabel="p", ylabel="P(data|p)",
    grid=false, legend=false, color=:black)
savefig("_assets/teaching/submod/lhood1.svg") # hide