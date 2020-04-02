# This file was generated, do not modify it. # hide
function test_different_λs(seqa, seqb, t)
    probs = []
    for λ=0.:0.01:10.
        push!(probs, [λ, ctmc_probability(seqa, seqb, t, λ)])
    end
    return hcat(probs...)
end

probs = test_different_λs("TTAT", "TTGG", 1.4)

themax, index = findmax(probs[2,:])
println("Maximum likelihood value: P(data|̂λ) = $themax")
println("ML estimate: ̂λ = $(probs[1,index])")