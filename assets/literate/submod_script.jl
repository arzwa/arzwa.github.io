# This file was generated, do not modify it.

function simulate(site, n, p)
    print(site, " ")
    for i=1:n
        site = rand() < p ? rand(setdiff("ATCG", site)) : site
        print("⟶  $site ")
    end
end

simulate('A', 10, 0.4)

function simulate(site, n, p)
    print("(X(0) = $site")
    for i=1:n
        site = rand() < p ? rand(setdiff("ATCG", site)) : site
        print(", X($i) = $site")
    end
    print(")")
end

simulate('A', 5, 0.4)

Pmatrix(p) = [1-p p/3 p/3 p/3 ;
              p/3 1-p p/3 p/3 ;
              p/3 p/3 1-p p/3 ;
              p/3 p/3 p/3 1-p ]
P = Pmatrix(0.2)

p₀ = [0. 1. 0. 0.]'

p₁ = P*p₀

P^10*p₀

for i in [1, 2, 5, 10, 20, 50, 100]
    @show round.(P^i*p₀, digits=3)
end

P*[0.25 0.25 0.25 0.25]'

using Distributions
ntoi = Dict(n=>i for (i,n) in enumerate("ATCG"))
iton = Dict(i=>n for (i,n) in enumerate("ATCG"))
translate(seq::String) = [ntoi[n] for n in seq]
translate(seq::Vector{Int}) = join([iton[i] for i in seq])

function simulate(seq, P, n)
    x = translate(seq)
    Pn = P^n
    x = map(i->rand(Categorical(Pn[:,i])), x)
    return translate(x)
end

original = "ATCGGGCGGGATTATTACGGAT"
evolved  = simulate(original, Pmatrix(0.05), 10)
diffs = join([original[i] == evolved[i] ? "|" : " " for i=1:length(original)])
println(original, "\n", diffs, "\n", evolved)

Pn = Pmatrix(0.05)^10  # get the 10-step transition probabilities
x = translate(original)
y = translate(evolved)
site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
sequence_probability = prod(site_probabilities)

sum(log.(site_probabilities))

println(join(x), "\n", join(y))

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

using Plots
plot(probs[1,:], probs[2,:],
    xlabel="p", ylabel="P(data|p)",
    grid=false, legend=false, color=:black)
savefig("_assets/teaching/submod/lhood1.svg") # hide

themax, index = findmax(probs[2,:])
println("Maximum likelihood value: P(data|̂p) = $themax")
println("ML estimate: ̂p = $(probs[1,index])")

using Printf
randexp(λ) = -log(rand())/λ  # generate a random number from the exponentil ditribution with rate λ

function simulate(site, λ, t)
    @printf "  X(%.2f) = %s\n" 0. site
    curr_t = randexp(λ)  # time of the first substitution event
    while curr_t < t
        site = rand(setdiff("ATCG", site))  # change state
        @printf "↪ X(%.2f) = %s\n" curr_t site
        curr_t += randexp(λ)  # get time of next substitution event
    end
    @printf "↪ X(%.2f) = %s\n" t site
end

simulate('A', 10., 1.)

Q(λ) = [-3λ   λ   λ   λ ;
          λ -3λ   λ   λ ;
          λ   λ -3λ   λ ;
          λ   λ   λ -3λ ]

theQ = Q(0.2)

theP = exp(theQ*1.2)

exp(theQ*0.8)[1,4]

exp(theQ*10.)

function ctmc_probability(seqa, seqb, t, λ)
    x = translate(seqa)
    y = translate(seqb)
    Pt = exp(Q(λ)*t)
    sum(log.([Pt[j,i] for (i,j) in zip(x,y)]))
end

ctmc_probability("TTAT", "TTGG", 0.1, 1.)

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

seqa = "TTTATCGACCTATTC"
seqb = "TAAAACGAACTATAC"
p = heatmap(0.1:0.02:2, 0.1:0.02:2, size=(400,350),
    (λ,t)->ctmc_probability(seqa, seqb, t, λ),
    fill=true, xlabel="\\lambda", ylabel="t", title="log-likelihood")
for a=0.1:0.25:5
    plot!(p, x->a/x, ylim=(0.1,2), xlim=(0.1,2), legend=false, color=:black)
end
savefig("_assets/teaching/submod/lhood2.svg") # hide

plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
savefig("_assets/teaching/submod/lhood3.svg") # hide

distance_JC(p) = -0.75 * log(1. - 4p/3)
p = mapreduce(!=, +, seqa, seqb)/length(seqa)
distance_JC(p)

plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
dist = distance_JC(p)
vline!([dist])
hline!([ctmc_probability(seqa, seqb, dist, 1/3)])
savefig("_assets/teaching/submod/lhood4.svg") # hide

