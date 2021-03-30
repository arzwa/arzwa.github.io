# This file was generated, do not modify it. # hide
using Printf
randexp(λ) = -log(rand())/λ  # generate a random number from the exponential distribution with rate λ

function simulate(site, λ, t)
    @printf "X(%.2f) = %s\n" 0. site
    curr_t = randexp(λ)  # time of the first substitution event
    while curr_t < t
        site = rand(setdiff("ATCG", site))  # change state
        @printf "X(%.2f) = %s\n" curr_t site
        curr_t += randexp(λ)  # get time of next substitution event
    end
    @printf "X(%.2f) = %s\n" t site
end

simulate('A', 10., 1.)