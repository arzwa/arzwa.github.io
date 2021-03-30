# This file was generated, do not modify it. # hide
function ctmc_probability(seqa, seqb, t, λ)
    x = translate(seqa)
    y = translate(seqb)
    Pt = exp(Q(λ)*t)
    sum(log.([Pt[j,i] for (i,j) in zip(x,y)]))
end

ctmc_probability("TTAT", "TTGG", 0.1, 1.)