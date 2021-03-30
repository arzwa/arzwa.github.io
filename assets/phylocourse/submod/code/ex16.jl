# This file was generated, do not modify it. # hide
themax, index = findmax(last.(l))
println("Maximum likelihood value: P(data|p=p̂) = $themax")
println("ML estimate: ̂p = $(l[index])")