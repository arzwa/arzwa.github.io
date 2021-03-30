# This file was generated, do not modify it. # hide
plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
savefig("_assets/phylocourse/submod/lhood3.svg") # hide