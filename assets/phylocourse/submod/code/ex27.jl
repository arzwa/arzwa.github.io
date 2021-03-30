# This file was generated, do not modify it. # hide
plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
dist = distance_JC(p)
vline!([dist])
hline!([ctmc_probability(seqa, seqb, dist, 1/3)])
savefig("_assets/phylocourse/submod/lhood4.svg") # hide