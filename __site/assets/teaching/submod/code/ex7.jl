# This file was generated, do not modify it. # hide
for i in [1, 2, 5, 10, 20, 50, 100]
    @show round.(P^i*p₀, digits=3)
end