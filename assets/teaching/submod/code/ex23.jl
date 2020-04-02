# This file was generated, do not modify it. # hide
seqa = "TTTATCGACCTATTC"
seqb = "TAAAACGAACTATAC"
p = heatmap(0.1:0.02:2, 0.1:0.02:2, size=(400,350),
    (λ,t)->ctmc_probability(seqa, seqb, t, λ),
    fill=true, xlabel="\\lambda", ylabel="t", title="log-likelihood")
for a=0.1:0.25:5
    plot!(p, x->a/x, ylim=(0.1,2), xlim=(0.1,2), legend=false, color=:black)
end
savefig("_assets/teaching/submod/lhood2.svg") # hide