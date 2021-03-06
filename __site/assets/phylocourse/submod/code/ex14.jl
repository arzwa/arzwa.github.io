# This file was generated, do not modify it. # hide
function test_different_ps(x, y, n)
    l = []
    for p=0.0:0.001:1.0
        Pn = Pmatrix(p)^n
        site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
        loglikelihood = sum(log.(site_probabilities))
        push!(l, (p, loglikelihood))
    end
    return l
end

l = test_different_ps(x, y, 10);