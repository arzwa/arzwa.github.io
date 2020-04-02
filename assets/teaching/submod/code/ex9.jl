# This file was generated, do not modify it. # hide
using Distributions
ntoi = Dict(n=>i for (i,n) in enumerate("ATCG"))
iton = Dict(i=>n for (i,n) in enumerate("ATCG"))
translate(seq::String) = [ntoi[n] for n in seq]
translate(seq::Vector{Int}) = join([iton[i] for i in seq])

function simulate(seq, P, n)
    x = translate(seq)
    Pn = P^n
    x = map(i->rand(Categorical(Pn[:,i])), x)
    return translate(x)
end

original = "ATCGGGCGGGATTATTACGGAT"
evolved  = simulate(original, Pmatrix(0.05), 10)
diffs = join([original[i] == evolved[i] ? "|" : " " for i=1:length(original)])
println(original, "\n", diffs, "\n", evolved)