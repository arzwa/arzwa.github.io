# This file was generated, do not modify it. # hide
function test_different_ps(x, y, n)
    probs = []
    for p=0.:0.001:1.
        Pn = Pmatrix(p)^n
        site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
        push!(probs, [p,sum(log.(site_probabilities))])
    end
    return hcat(probs...)
end

probs = test_different_ps(x, y, 10)