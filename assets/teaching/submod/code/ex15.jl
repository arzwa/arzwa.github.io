# This file was generated, do not modify it. # hide
themax, index = findmax(probs[2,:])
println("Maximum likelihood value: P(data|̂p) = $themax")
println("ML estimate: ̂p = $(probs[1,index])")