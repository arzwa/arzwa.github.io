# This file was generated, do not modify it.

function simulate(site, n, p)
    print(site)
    for i=1:n
        if rand() < p
            site = rand(setdiff("ATCG", site))
        end
        print(" --> ", site)
    end
end

@show setdiff("ATCG", 'A')
@show rand(setdiff("ATCG", 'A'))

simulate('A', 10, 0.4)

Pmatrix(p) = [1-p p/3 p/3 p/3 ;
              p/3 1-p p/3 p/3 ;
              p/3 p/3 1-p p/3 ;
              p/3 p/3 p/3 1-p ]
P = Pmatrix(0.2)

f0 = [0, 1, 0, 0.];

f1 = P*f0

P^10*f0

for n in [1, 2, 5, 10, 20, 50, 100, 200]
    fn = P^n*f0
    println(round.(fn, digits=3))
end

P*[0.25, 0.25, 0.25, 0.25]

using Distributions
ntoi = Dict(n=>i for (i,n) in enumerate("ATCG"))
iton = Dict(i=>n for (i,n) in enumerate("ATCG"))
translate(seq::String) = [ntoi[n] for n in seq]
translate(seq::Vector{Int}) = join([iton[i] for i in seq])

function simulate(seq, P, n)
    x = translate(seq)   # translate nucleotides to integers
    Pn = P^n
    x = map(i->rand(Categorical(Pn[:,i])), x)
    return translate(x)  # translate integers back to nucleotides
end

original = "ATCGGGCGGGATTATTACGG"
evolved  = simulate(original, Pmatrix(0.05), 10)
diffs = join([original[i] == evolved[i] ? "|" : " " for i=1:length(original)])
println(original, "\n", diffs, "\n", evolved)

x = translate(original)
y = translate(evolved);
Pn = Pmatrix(0.05)^10  # get the 10-step transition probabilities
site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
sequence_probability = prod(site_probabilities)

sum(log.(site_probabilities))

println(join(x), "\n", join(y))

function test_different_ps(x, y, n)
    l = []
    for p=0.0:0.001:1.0
        Pn = Pmatrix(p)^n
        site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
        push!(l, (p,sum(log.(site_probabilities))))
    end
    return l
end

l = test_different_ps(x, y, 10);

using Plots
plot(first.(l), last.(l),
    xlabel="p", ylabel="P(data|p)",
    grid=false, legend=false, color=:black)
savefig("_assets/phylocourse/submod/lhood1.svg") # hide

themax, index = findmax(last.(l))
println("Maximum likelihood value: P(data|p=p̂) = $themax")
println("ML estimate: ̂p = $(l[index])")

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

seqa = "TTTATCGACCTATTC"
seqb = "TAAAACGAACTATAC"
p = heatmap(0.1:0.02:2, 0.1:0.02:2, size=(400,350),
    (λ,t)->ctmc_probability(seqa, seqb, t, λ),
    fill=true, xlabel="\\lambda", ylabel="t", title="log-likelihood")
for a=0.1:0.25:5
    plot!(p, x->a/x, ylim=(0.1,2), xlim=(0.1,2), legend=false, color=:black)
end
savefig("_assets/phylocourse/submod/lhood2.svg") # hide

plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
savefig("_assets/phylocourse/submod/lhood3.svg") # hide

distance_JC(p) = -0.75 * log(1. - 4p/3)
p = mapreduce(!=, +, seqa, seqb)/length(seqa)
distance_JC(p)

plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
dist = distance_JC(p)
vline!([dist])
hline!([ctmc_probability(seqa, seqb, dist, 1/3)])
savefig("_assets/phylocourse/submod/lhood4.svg") # hide

