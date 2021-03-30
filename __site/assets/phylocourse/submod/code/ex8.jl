# This file was generated, do not modify it. # hide
for n in [1, 2, 5, 10, 20, 50, 100, 200]
    fn = P^n*f0
    println(round.(fn, digits=3))
end