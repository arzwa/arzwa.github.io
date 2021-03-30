# This file was generated, do not modify it. # hide
function test_different_λs(seqa, seqb, t)
    l = []
    for λ=0.:0.01:10.
        push!(l, (λ, ctmc_probability(seqa, seqb, t, λ)))
    end
    return l
end

l = test_different_λs("TTAT", "TTGG", 1.4)

themax, index = findmax(last.(l))
println("Maximum likelihood value: P(data|̂λ) = $themax")
println("ML estimate: ̂λ = $(l[index])")