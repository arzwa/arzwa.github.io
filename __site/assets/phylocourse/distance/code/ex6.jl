# This file was generated, do not modify it. # hide
using Distributions
p = plot(title="The Gamma distribution with mean 1")
for α in [0.1, 0.25, 0.5, 1.0, 5.0, 10., 100.]
    plot!(p, Gamma(α, 1/α), label="\\alpha = $α", xlim=(0,5), ylim=(0,5))
end
savefig(p, "_assets/phylocourse/distance/gamma.svg") # hide