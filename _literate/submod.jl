# [back](/phylocourse/)
# 
# \newcommand{\albet}{\mathcal{A}}
# \newcommand{\Dt}{\Delta t}
#
# \toc
# >**Exercises** will be marked in blocks like this, try to formulate a brief
# >answer to them.

# ## Substitution models for molecular phylogenetics
#
# By now you should have some idea of the major approaches towards the
# inference of phylogenies from molecular sequence data. As discussed in the
# course notes by prof. Van de Peer a rough partitioning of the different
# methods is the following

# |                     | Model-based                            | Not model-based          |
# |---------------------|----------------------------------------|--------------------------|
# | Character-based     | Maximum likelihood, Bayesian inference | Maximum parsimony        |
# | Not character-based | Distance methods                       | (no established methods) |

# Nowadays, most methods in phylogenetics -- be it for tree construction (which
# is the focus of the course), inference of adaptive protein evolution or
# macroevolutionary analyses of trait evolution and diversification -- are
# model-based approaches that employ a **probabilistic model of evolution** and
# seek to perform **statistical inference of the parameters of this model**. In
# these notes I will attempt to make clear what these models of evolution
# amount to and I will try to convey some intuition on how modern phylogenetic
# methods are largely *statistical* in nature. 

# Concerning the problem of phylogeny reconstruction from sequence data, the
# probabilistic model of sequence evolution is usually a **substitution model**
# which describes how a single site in an alignment evolves 'along' a tree.  The
# goal in statistical phylogenetics is then to somehow find the phylogenetic
# tree and parameters of the substitution model that 'best explain' the
# observed sequence data. Note that this is deliberately formulated in a vague
# way, as the specifics of this procedure depend strongly on your philosophy of
# statistical inference. In particular two major approaches are used in
# phylogenetics, based on either the principle of [Maximum
# likelihood](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) or
# [Bayesian inference](https://en.wikipedia.org/wiki/Bayesian_inference).  In
# addition, methods based on pairwise distances go somewhat halfway a full
# probabilistic modeling approach, using substitution models (and hence,
# statistical methods) to compute *pairwise* distances among sequences but not
# explicitly modeling the evolution of a group of sequences along a phylogeny. 

# Using substitution models to model the evolutionary process
# probabilistically, we can answer questions that become progressively more
# interesting.

# > Q1: What is the probability that a given nucleotide substitutes into another
# > over some time $t$ given the substitution model?

# > Q2: What is the probability that the sequence `ATCGAATCG` evolves into
# >`ATCCCAATC` over some time $t$ given the substitution model?

# > Q3: Assuming some substitution model, what is the expected number of
# > substitutions that happened during the course of evolution when we observe
# > $x$ differences between two sequences?

# > Q4: What is the probability of observing a set of $n$ sequences given that they
# > evolved along a phylogenetic tree $\mathcal{T}$ under a given substitution
# > model?

# > Q5: Which phylogenetic tree connecting a set of $n$ sequences is the most
# > likely to have resulted in the observed sequences given the substitution
# > model?

# I hope these questions already given you an idea of why such models of
# sequence evolution are interesting! The last question -- which asks for the
# tree that 'explains' the data best under the model -- is essentially the
# problem of maximum likelihood inference of phylogenies. The goal of these
# notes is to give you a feel of how these questions can be answered and how
# they relate to phylogenetic inference. 
#
# > **Exercise.** Consider the last question (Q5) in the above list. How does this 
# > compare to the approach taken in maximum parsimony tree inference? Under
# > what condition does the maximum parsimony tree answer the question asked?

# ## Discrete-time Markov models of nucleotide substitution

# ### Simulating the simplest model 

# The simplest non-trivial probabilistic model of DNA substitution could be
# specified in words as follows:

# > In a given time step of length $Δt$, a given nucleotide changes into random
# > other nucleotide with probability $p$, it stays the same nucleotide with
# > probability $1-p$.

# This is a simple probabilistic model of sequence evolution with a single
# parameter $p$. We could for instance use this model to **simulate** the
# evolution of a given site in a sequence over a series of $n$ time steps:

function simulate(site, n, p)
    print(site)
    for i=1:n
        if rand() < p 
            site = rand(setdiff("ATCG", site)) 
        end
        print(" --> ", site)
    end
end

# Note that the string `"ATCG"` is equivalent to a list of characters `['A',
# 'T', 'C', 'G']`, that `setdiff` takes the set difference and that `rand`
# picks a random element from a collection. To be clear:
@show setdiff("ATCG", 'A')
@show rand(setdiff("ATCG", 'A'))

# Now we do a *simulation* of 10 time steps assuming an initial state `A` and 
# with $p$ (the single parameter of the model) set to $0.4$:

simulate('A', 10, 0.4)

# The model of sequence evolution we specified is a simple kind of [stochastic
# process](https://en.wikipedia.org/wiki/Stochastic_process) belonging to the
# class of **discrete-time Markov chains (DTMC)**. A stochastic process can be
# thought of intuitively as something that evolves randomly in time. Here this
# 'something' is the state of the nucleotide of interest, and we denote the state
# of the nucleotide at time point $n$ by $X(n)$, taking values in the state
# space $\{A, T, C, G\}$.

# Note that the 'discrete-time' label refers to the fact that the evolutionary
# process is assumed to proceed in discrete time steps which can be indexed by
# integers, i.e. $n \in \mathbb{N}$ (continuous-time models, where we consider
# $X(t)$ with $t \in \mathbb{R}^+$, will be discussed below). In the simulation
# function above, we effectively iterate through the time steps, at each time
# point flipping a biased coin which lands heads with probability $p$ and tails
# with probability $1-p$. If the virtual coin lands heads, the nucleotide
# substitutes, if it lands tails it stays put.

# A sequence $(X(1), X(2), X(3), \dots)$ is called a *realization* of the
# stochastic process. The simulation function above prints such a realization
# given the initial state `A`.

# > **Exercise.** Use the simulation function above to estimate the probability
# > that a nucleotide `A` substitutes to a nucleotide `T` after 10 time steps
# > when $p = 0.1$

# ### The Markov property

# We called the above model a discrete-time **Markov** chain. The key feature
# of all stochastic processes with the **Markov property** is that at any given
# time, the future evolution of the process only depends on the current state
# and not the history of the process. In probabilistic notation[^pr] :

# $$ P\big[X(n+1)|X(n),X(n-1),X(n-2),\dots,X(0)\big] = P\big[X(n+1)|X(n)\big] $$

# The probability $P[X(n+1)|X(n)]$ is known as the **transition probability**.

# This property is what makes these models particularly mathematically tractable[^sp].
# Using the transition probabilities and the Markov property, we can compute the
# probability of a particular realization under the model. Consider for instance
# the simple model above with the parameter $p$ and an observed realization
# `A ⟶  G ⟶  G ⟶  T`. The probability of this realization under the model
# can be obtained as:

# \begin{align}
# P(&A \rightarrow G \rightarrow G \rightarrow T) = P(X_0=A, X_1=G, X_2=G, X_3=T) \\
# &= P(X_3=T|X_0=A,X_1=G,X_2=G)P(X_0=A,X_1=G,X_2=G) \\
# &= P(X_3=T|X_0=A,X_1=G,X_2=G) P(X_2=G|X_0=A,X_1=G) P(X_0=A,X_1=G) \\
# &= P(X_3=T|X_0=A,X_1=G,X_2=G) P(X_2=G|X_0=A,X_1=G) P(X_1=G|X_0=A) P(X_0=A)\\
# &= P(X_3=T|X_2=G)P(X_2=G|X_1=G)P(X_1=G|X_0=A)P(X_0=A) \\
# &= \frac{p}{3} \times (1-p) \times \frac{p}{3} \times 1 \\
# &= \Big(\frac{p}{3}\Big)^2(1-p)
# \end{align}

# Where I have written $X_n$ instead of $X(n)$ to avoid all those parentheses.
# While the equations are a bit convoluted, this is just a simple exercise in
# conditional probabilities (see also note [^pr]).

# But *why is being able to compute such a probability interesting?* Consider a
# question of the sort:

# >For a given site, what is the probability to observe a nucleotide `A` at the
# >present when it was an `T` in an ancestor $n$ generations ago, assuming the
# >probability of a substitution in a single generation is $p$ (with all kinds
# >of substitutions equally likely)?

# Clearly this starts to sound like something that can be of scientific
# interest (if you don't see why, please take a moment to ponder why being able
# to answer such a question is interesting!). We cannot quite apply the
# calculations we did above directly though, since we have not specified the
# entire sequence of transitions like we did above. We are looking for the
# probability

# $$P(T \rightarrow \underbrace{? \rightarrow ? \rightarrow \dots \rightarrow ?}_{n-1 \text{ intermediate states}} \rightarrow A)$$

# or the **$n$-step transition probability** $P(X_n=A |X_0=T)$ asociated with
# the Markov chain. If we write the transition probability $P[X(n+1)=A|X(n)=T]$
# (i.e. the probability of a `T --> A` transition) more concisely as $p_{AT}$.
# If we write the **alphabet** of the sequence data as $\albet = \{A,T,C,G\}$,
# the probability we are querying in the above question is of the form

# \begin{align}
# P(X_n=A |X_0=T) &= \sum_{x \in \albet} p_{Ax} P(X_{n-1}=x|X_0=T) \\
#   &= \sum_{x \in \albet} p_{Ax} \Big(
#           \sum_{y \in \albet} p_{xy} P(X_{n-2}=y|X(0)=T) \Big) \\
#   &= \sum_{x \in \albet} p_{Ax} \Big(
#           \sum_{y \in \albet} p_{xy} \Big( \dots \sum_{z \in \albet} p_{zT} \Big) \Big)
# \end{align}

# Nothing stops us from implementing this recursive formulation[^recursion] of
# the $n$-step transition probability in a computer program, however we can
# formulate this elegantly and concisely using a little bit of matrix algebra.
#
# > **Exercise.** Explain in words the reasoning behind the first equality in
# > the above equation. Can you think of a different way to start the derivation?

# ### The transition probability matrix

# *Any* DTMC model of nucleotide substitution is completely determined by its
# so-called **transition probability matrix** $P$:

# $$P = \begin{bmatrix}
#     p_{AA} & p_{AT} & p_{AC} & p_{AG} \\
#     p_{TA} & p_{TT} & p_{TC} & p_{TG} \\
#     p_{CA} & p_{CT} & p_{CC} & p_{CG} \\
#     p_{GA} & p_{GT} & p_{GC} & p_{GG}
# \end{bmatrix}$$

# You can fill this in for the simple model we defined above to get
# $$P_{simple} = \begin{bmatrix}
#     (1-p) & p/3 & p/3 & p/3 \\
#     p/3 & (1-p) & p/3 & p/3 \\
#     p/3 & p/3 & (1-p) & p/3 \\
#     p/3 & p/3 & p/3 & (1-p)
# \end{bmatrix}$$

# At this point, it's fairly easy to see that we could easily specify more
# complicated Markov models by introducing more parameters. For instance we could
# have different probabilities of a purine (`A` and `G`) to purine substitution
# versus a purine to pyrimidine substitution. But note that *however one
# specifies the model, each column of the transition probability matrix has to
# sum to one*. This simply means that the chain has to go *somewhere*, i.e. the
# total probability to be in *any* state in the next time step has to be one.
# To make this completely clear, consider the current state of the process is 
# $C$, then the third column of $P$ gives the probability distribution over
# the state at the next time step, i.e. with probability $p_{AC} = p/3$ will the
# chain move to state $A$, with probability $p_{TC} =p/3$ will the chain move to
# state $T$, with probability $p_{CC} = (1-p)$ there will be no change and again
# with probability $p/3$ the chain will move to state $G$.

# For a given initial distribution over states $f(0) = [f_A, f_T, f_C, f_G]$,
# it is not hard to verify that the following matrix multiplication

# $$f(1) = P f(0) = \begin{bmatrix} 
#   p_{AA} f_A + p_{AT} f_T +p_{AC} f_C + p_{AG} f_G \\
#   p_{TA} f_A + p_{TT} f_T +p_{TC} f_C + p_{TG} f_G \\
#   p_{CA} f_A + p_{CT} f_T +p_{CC} f_C + p_{CG} f_G \\
#   p_{GA} f_A + p_{GT} f_T +p_{GC} f_C + p_{GG} f_G \end{bmatrix}
#   $$

# gives the probability distribution over states at the first time point. 
 
# >**Exercise.** Explain in words what the above equation means.

# For instance, consider the transition probability matrix
## Note: in julia one can define short functions using the notation 
## `f(x) = <expression in x>`, without using the `function` keyword
Pmatrix(p) = [1-p p/3 p/3 p/3 ;
              p/3 1-p p/3 p/3 ;
              p/3 p/3 1-p p/3 ;
              p/3 p/3 p/3 1-p ]
P = Pmatrix(0.2)

# and if we assume an initial state `T`, the initial distribution is (if we
# keep the `ATCG` order):

f0 = [0, 1, 0, 0.];

# i.e. the probability the nucleotide is in state `A`, `T`, `C` or `G` at time
# point 0 is 0, 1, 0 and 0 respectively (this is just another way of
# expressing that the initial state is `T`). We obtain the probability
# distribution over states at the first time point as
f1 = P*f0

# This is simply the vector of transition probabilities $[p_{AT}, p_{TT},
# p_{CT}, p_{GT}]$. Now, of course we can repeat this

# $$f(2) = Pf(1) = P(Pf(0)) = P^2f(0)$$

# and generally we get

# $$f(n) = P^nf(0)$$

# which provides us with an easy way to compute the $n$-step transition
# probabilities. For instance, the vector

P^10*f0

# gives the 10-step transition probabilities for the initial state `T`.

# >**Exercise.** Why is it the case that $f(2) = P f(1) = P(Pf(0))$? What
# >assumption are we making here?

# Clearly, we can now answer our question from before:

# >For a given site, what is the probability to observe a nucleotide `A` at the
# >present when it was a `T` in an ancestor $n$ generations ago, assuming the
# >probability of a substitution in a single time step is $p$ (with all kinds of
# >substitutions equally likely)? 
#
# For $n = 10$ and $p = 0.2$, this is $0.239$ (see above).

# > **Exercise**. Use the code above (or write your own) to compute the
# > probability that a `T` substitutes to a `G` when $p = 0.1$ and $n=5$.

# Note that as the number of time steps $n$ grows, the probability distribution
# over states converges to an equilibrium
for n in [1, 2, 5, 10, 20, 50, 100, 200]
    fn = P^n*f0
    println(round.(fn, digits=3))
end

# This distribution is called the **stationary distribution** of the Markov
# chain, and in the simple model we defined above it is simply $[1/4, 1/4, 1/4,
# 1/4]$. This of course could be expected, since nothing in the model seems to
# systematically prefer one nucleotide over the other. Clearly, this entails
# that, given enough time, it will be impossible to determine what the initial
# state of any given nucleotide position was based on the observed state any
# better then a random guess would. 
#
# If we write the stationary distribution as $\pi$, it has the property that

# $$\pi = P\pi$$

P*[0.25, 0.25, 0.25, 0.25]

# Which suggests why 'stationary' is really the right word to use -- once the
# Markov chain is in the stationary state (in this case this means that once
# the probability a given nucleotide is in any state is 1/4), it will stay
# there. 
#
# The stationary distribution of a Markov chain is not only an interesting
# thing in Markov chain theory, but is also relevant for molecular evolution.
# For instance, if some genome has a strong GC bias, a model with stationary
# distribution `[0.25, 0.25, 0.25, 0.25]` may seem unrealistic. Substition 
# models with stationary distributions that are more realistic for the data
# at hand have therefore been developed -- as we will see later.

# Lastly, we note that for some Markov models the $n$-step transition
# probabilities have an analytical solution. For the simple one parameter DTMC,
# the matrix $P^n$ with the $n$-step transition probabilities is of the form

# $$P_{simple}^n = \begin{bmatrix}
#     p_0(n) & p_1(n) & p_1(n) & p_1(n) \\
#     p_1(n) & p_0(n) & p_1(n) & p_1(n) \\
#     p_1(n) & p_1(n) & p_0(n) & p_1(n) \\
#     p_1(n) & p_1(n) & p_1(n) & p_0(n) \end{bmatrix} $$

# with

# $$ p_0(n) = \frac{1}{4} + \frac{3}{4}\Big(1 - \frac{4}{3}p\Big)^n\ \text{  and  }\ p_1(n) = \frac{1}{4} - \frac{1}{4}\Big(1 - \frac{4}{3}p\Big)^n $$

# In other words, for a given site, the probability of observing a different
# nucleotide after $n$ time steps is $\frac{1}{4} - \frac{1}{4}\Big(1 -
# \frac{4}{3}p\Big)^n$. However for most Markov models such a simple formula
# does not exist, so it is generally more important to understand the approach
# using the transition probability matrix directly.

# ### From individual sites to sequences

# If we make the assumption that sites in a sequence evolve independently
# (clearly a terrible assumption from a biological point of view, but a very
# convenient one statistically), the DTMC for a single site quickly generalizes
# to a model for entire sequences. Consider the following simulation of a
# sequence evolving according to the simple DTMC:

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

# >**Exercise**. Try to understand the code above and explain in words how the
# >simulation algorithm works. 

# Now to obtain the probability that one sequence evolves into another we can
# simply multiply the probabilities of the corresponding site-wise
# probabilities because of our independence assumption. 
#

x = translate(original)
y = translate(evolved);
Pn = Pmatrix(0.05)^10  # get the 10-step transition probabilities
site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
sequence_probability = prod(site_probabilities)

# If $\ell_i$ denotes the probability that the nucleotide at position $i$ in
# sequence $x$ evolved into the nucleotide at position $i$ in sequence $y$ over
# 10 time steps, then the probability `sequence_probability` in the example
# above could be written formally as 
#
# $$P\big(\mathrm{seq}(10)=y|\mathrm{seq}(0)=x, p=0.05\big) = \ell_1 \ell_2 \dots \ell_n
#     = \prod_{i=1}^n \ell_i$$
#
# Make sure you understand what this means!
#  
# > **Exercise.** Note that these probabilities are *very* small values. Is this
# > surprising? Does this mean evolution is unlikely? Can you explain why (not)?
# > (Hint: how many possible DNA sequences are there of length 20?)

# When working with small values in numerical applications, its almost
# always better to work on a log scale, i.e. the log probability of this
# particular evolutionary scenario under the substitution model is

sum(log.(site_probabilities))

# Of course you recall that $\log ab = \log a + \log b$, so that multiplying a
# bunch of probabilities corresponds to summing them on a logarithmic scale.

# >**Exercise.** What is the probability that the sequence in the variable `x` 
# >in the code example above evolves into a sequence "AAAAAAAAAAAAAAAAAAAA"
# >assuming $p = 0.05$ and $n=10$? Compare this probability to the probability 
# >computed above for the evolution of `x` into `y`, explain why this probability
# >is different.

# ### A note on time-reversibility

# The simple DTMC model of sequence evolution we discussed above is a
# **time-reversible model**. Mathematically, what this means is that

# $$\pi_i p_{ji}(n) = \pi_j p_{ij}(n) $$

# For instance $\pi_A p_{TA}(n) = \pi_T p_{AT}(n)$. Here $\pi_i$ is the
# probability for state $i$ under the stationary distribution, which in turn we
# can think of as the probability that a typical nucleotide site evolved under
# the assumed model is in state $i$. What the above equation means is that the
# 'flow' from `A` to `T` is the same as the flow from `T` to `A` for instance.
# Another way of saying that the amount of evolutionary change from `T` to `A`
# is the same as the amount of change from `A` to `T` (this does *not* mean
# that the *rates* of change are the same). An important consequence of this
# time-reversibility is that the following two scenario's will have the same
# probability under the substitution model:

# 1. A sequence `ATTTCATG` evolving over $2n$ time steps into the sequence
#    `TTGGGATG`
# 2. An unknown ancestral sequence evolving independently into the sequences
#    `ATTTCATG` and `TTGGGATG` over a time $t$

# ![](/assets/phylocourse/submod/rev.png)

# In other words, from the point of view of the model, both scenario's are
# indistinguishable. Almost all commonly used substitution models in
# phylogenetics are reversible.

# ### An example of maximum likelihood estimation

# If you're unsure what all this stuff has to do with statistics, consider the
# following scenario: you observed the two sequences above, i.e. (in integer
# representation `A=1, T=2, C=3, G=4`):
println(join(x), "\n", join(y))

# This is the **data**. Now let's say we believe the simple DTMC model is a
# reasonable model for these sequences, and assume we know they are derived
# from a common ancestral sequence $5$ time steps in the past.  You have no clue
# however about the $p$ parameter of the model, and you would like to *infer
# this parameter from the observed data*.  One thing you could do is look for
# the value of $p$ that makes the data *most likely*, i.e. results in the
# highest probability under the model. In other words, we could look for the
# value of $p$ that makes $P(\text{data}|p)$ attain a maximum. A naive approach
# do this is the following:
function test_different_ps(x, y, n)
    l = []
    for p=0.0:0.001:1.0
        Pn = Pmatrix(p)^n
        site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
        loglikelihood = sum(log.(site_probabilities))
        push!(l, (p, loglikelihood))
    end
    return l 
end

l = test_different_ps(x, y, 10);

# We can plot the results
using Plots
plot(first.(l), last.(l),
    xlabel="p", ylabel="P(data|p)",
    grid=false, legend=false, color=:black)
savefig("_assets/phylocourse/submod/lhood1.svg") # hide

# ![](/assets/phylocourse/submod/lhood1.svg)

# Incidentally, the value $P(\text{data}|p)$ viewed as a function of $p$ for a
# fixed data set is called the **likelihood function** of $p$. The maximum
# likelihood estimate (MLE) using this naive brute force approach is

themax, index = findmax(last.(l))
println("Maximum likelihood value: P(data|p=p̂) = $themax")
println("ML estimate: ̂p = $(l[index])")

# If you understand this, you understand the essence of maximum likelihood
# estimation in statistics.

# >**Exercise.** Make sure you understand the above code. Is the MLE as
# >expected? Where did we use the time-reversibility property? Why does the
# >likelihood curve has this long plateau?

# ## The real deal: continuous-time Markov models of sequence evolution

# If you managed to get through the above sections, now comes a serious
# disappointment. DTMCs are generally *not* used in molecular phylogenetics.
# However, their closely related continuous-time counterparts are. In this
# section we will take a look at the **continuous-time Markov chain (CTMC)**
# analog of the simple model considered above. As is often the case when going
# from discrete to continuous things in mathematics, the former can be seen as
# a limit of the latter, and therefore our knowledge of DTMCs will be very
# helpful here.

# ### CTMC's and the Poisson process

# We will develop the continuous-time analog of the DTMC described above.  In
# the continuous time models, we no longer model the evolution of a nucleotide
# site in discrete time steps, but over the positive real line $\mathbb{R}$. We
# will denote the state at a time point $t$ as $X(t) \in \{A,T,C,G\}$, and our
# goal will be to develop a stochastic model for the evolution of $X$ over
# time.

# Whereas in the DTMC case we specified a probability of substitution per time
# step (the $p$ parameter), now we consider a **substitution rate** $\lambda$.
# The meaning of the $\lambda$ parameter is best understood as follows: the
# probability that a nucleotide substiutes in a small time interval $\Delta t$
# will be approximately $\lambda \Delta t$, with the approximation becoming
# more correct as $\Delta t \rightarrow 0$. The expected number of
# substitutions over a time length $t$ will be 

# $$\mathbb{E}[\text{number of substitutions}] = \lambda t$$

# Note that as in basic physics, multiplying a 'rate' times a 'time' variable
# (e.g. as in $\lambda t$) results in a **distance**. This single-parameter
# CTMC is closely related to a **[Poisson
# process](https://en.wikipedia.org/wiki/Poisson_point_process)**. The number
# of substitution events under the model happening over a time $t$ is a Poisson
# distributed random variable with mean $\lambda t$. It is a property of the
# Poisson process that the time between two events, the so-called **waiting
# time** $t_w$, is distributed according to an exponential distribution
# $\lambda e^{-\lambda t_w}$. The Poisson process and related CTMCs are very
# general models, used to model random phenomena as varied as the number of car
# accidents, telephone calls, radioactive decay events or, indeed, nucleotide
# substitutions in a certain time period.

# We can now do the sequence simulation again like we did in the very first
# section. We assume substitutions happen according to the Poisson process
# model, and that when a substitution happens, all kinds of substitutions are
# equally likely:

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

# > **Exercise.** Take some time to understand how the simulation works. Which
# > property of the Poisson process are we using to drive the simulation? How 
# > many substitutions did you expect to occur during the simulation above? 

# ### The rate matrix $Q$ and transition probabilities
#
# While the DTMC was directly defined in terms of a transition probability
# matrix, this is not the case for the CTMC.
# The continuous-time model we specified and simulated above is again a Markov
# model: both the time until the next substitution event, and the type of this
# substitution event are only dependent on the present state, and not on the
# past. The Markov chain is fully specified by its **rate matrix**
# $Q$[^generator]. The rate matrix for the simple model can be written as

# $$Q = \begin{bmatrix}
#     -3\lambda & \lambda & \lambda & \lambda \\
#     \lambda & -3\lambda & \lambda & \lambda \\
#     \lambda & \lambda & -3\lambda & \lambda \\
#     \lambda & \lambda & \lambda & -3\lambda \end{bmatrix} $$

# Where I have assumed that the total substitution rate (what we called
# $\lambda$ above) is given by $3\lambda$ (to avoid writing all those
# fractions). This continuous-time Markov chain (CTMC) model of nucleotide
# substitution is known as the **Jukes & Cantor model** (JC69, after Jukes &
# Cantor 1969). The entries of this matrix could be interpreted as rates of
# 'flow' between states in the model, with the diagonal elements the negative
# rates of flow *away* from some state.

# Now we wish to be able to compute transition probabilities of the form

# $$p_{AG}(\Delta t) = P\big(X(\Delta t)=A|X(0)=G\big)$$

# Note that $p_{ij}(0) = 0$ if $i \ne j$ and $p_{ii}(0) = 1$. In other words
# the transition probability matrix at $t=0$ is

# $$P(0) = I$$

# For a very small time interval $\Delta t$, we have approximately that
# $p_{ij}(\Delta t) = \lambda \Delta t$ if $i \ne j$ (this is the definition of
# the rate parameter $\lambda$) and $p_{ii}(\Delta t) = 1 - 3\lambda \Delta t$.
# In other words if we would take very small time steps, we could approximate
# the CTMC model by the DTMC model with $p=\lambda \Delta t$. In matrix notation:

# $$ P(\Dt) \approx I + Q \Dt = \begin{bmatrix}
#     1-3\lambda \Dt & \lambda \Dt & \lambda \Dt & \lambda \Dt \\
#     \lambda \Dt & 1-3\lambda \Dt & \lambda \Dt & \lambda \Dt \\
#     \lambda \Dt & \lambda \Dt & 1-3\lambda \Dt & \lambda \Dt \\
#     \lambda \Dt & \lambda \Dt & \lambda \Dt & 1-3\lambda \Dt \end{bmatrix} $$

# Divding by $\Delta t$ and taking the limit as $\Delta t \rightarrow 0$ (i.e.
# approximating the CTMC by a DTMC with smaller and smaller time steps), we find
# that

# $$\frac{dP(t)}{dt} \Bigg|_{t=0} = \lim_{\Dt \rightarrow 0} \frac{P(\Dt) - P(0)}{\Delta t} = 
# \lim_{\Dt \rightarrow 0}\frac{I + Q\Dt - I}{\Dt} = Q$$

# The rate matrix is therefore the derivative of the transition probability
# matrix at $t = 0$. The Markov property further implies that, in general
#
# $$\frac{dP(t)}{dt} = \frac{P(t + \Dt) - P(t)}{\Dt} = QP(t)$$  
#
# Which has the solution[^solution]

# $$ P(t) = e^{Qt} = I + Qt + Q^2\frac{t^2}{2!} + Q^3\frac{t}{3!} + \dots $$

# Which shows how the approximation relates to the exact solution (in the
# approximation we simply ignored all higher powers of $Q$). I agree having a
# matrix in the exponent looks scary, but don't worry, modern programming
# languages aren't scared of matrix exponentials:

Q(λ) = [-3λ   λ   λ   λ ;
          λ -3λ   λ   λ ;
          λ   λ -3λ   λ ;
          λ   λ   λ -3λ ]

theQ = Q(0.2)

# and now compute the transition probability matrix $P(t)$ at time $t = 1.2$
# using the matrix exponential formula:

theP = exp(theQ*1.2)

# The entry at $P[i,j]$ gives the transition probability $p_{ij}(t)$. So for
# instance $p_{AG}(0.8)$ is

exp(theQ*0.8)[1,4]

# By the way, note again how the matrix converges for large $t$

exp(theQ*10.)

# So now that we can compute $P(t)$, we can do everything that we did with the
# DTMC above, but in continuous time, which is much more convenient. For
# instance the log-probability that a sequence `TTAT` evolves into a sequence
# `TTGG` over a time $t=0.1$ million years given a substitution rate of
# $\lambda = 1$ substitution per million years can be computed as

function ctmc_probability(seqa, seqb, t, λ)
    x = translate(seqa)
    y = translate(seqb)
    Pt = exp(Q(λ)*t)
    sum(log.([Pt[j,i] for (i,j) in zip(x,y)]))
end

ctmc_probability("TTAT", "TTGG", 0.1, 1.)

# Or we can use the naive approach used above to obtain a maximum likelihood
# estimate for $\lambda$ (assuming we know $t$)

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

# To finish this section, again, we note that the simple Jukes-Cantor CTMC
# model has an analytical solution for the transition probabilities.
# Specifically the solution of $P(t) = e^{Qt}$ for the Jukes-Cantor rate matrix
# is

# $$P_{JC}(t) = \begin{bmatrix}
#     p_0(t) & p_1(t) & p_1(t) & p_1(t) \\
#     p_1(t) & p_0(t) & p_1(t) & p_1(t) \\
#     p_1(t) & p_1(t) & p_0(t) & p_1(t) \\
#     p_1(t) & p_1(t) & p_1(t) & p_0(t) \end{bmatrix} $$

# with

# $$ p_0(t) = \frac{1}{4} + \frac{3}{4} e^{-4 \lambda t}\ \text{  and  }\ p_1(t) = \frac{1}{4} - \frac{1}{4}e^{-4 \lambda t} $$

# which unsuprisingly shares a lot of structure with the analogous DTMC model.
# More general, parameter-rich substitution models do not admit such an
# expression in terms of simple formula's.

# ## So what are molecular distances anyway?

# The transition probbilities for the simplest CTMC model -- i.e. the
# Jukes-Cantor model -- are functions of two generally unknown parameters,
# being the substitution rate $\lambda$ and the time $t$ (in the DTMC model we
# had $p$ and $n$ serving similar roles). The transition probabilities involve
# only the product of these two parameters $\lambda t$, which is as we already
# noted above a **distance** (rate $\times$ time = distance). As you can easily
# see, the  transition probabilities will be the same when (1) the rate is
# $\lambda$ and the time is $t$ or when (2) the rate is $\lambda/2$ and the
# time is $2t$ or more generally (3) when the rate is $a\lambda$ and the time
# is $t/a$.

# So what happens if we would want to estimate both $\lambda$ and $t$ from the
# data using a maximum likelihood approach like the one used above?

seqa = "TTTATCGACCTATTC"
seqb = "TAAAACGAACTATAC"
p = heatmap(0.1:0.02:2, 0.1:0.02:2, size=(400,350),
    (λ,t)->ctmc_probability(seqa, seqb, t, λ),
    fill=true, xlabel="\\lambda", ylabel="t", title="log-likelihood")
for a=0.1:0.25:5
    plot!(p, x->a/x, ylim=(0.1,2), xlim=(0.1,2), legend=false, color=:black)
end
savefig("_assets/phylocourse/submod/lhood2.svg") # hide

# ![](/assets/phylocourse/submod/lhood2.svg)

# One can clearly see a boomerang shaped area in the likelihood surface
# corresponding to parameter region with high likelihood. 
#
# >**Exercise.** Interpret the likelihood surface plotted above, what can 
# >you say about the most likely values of the parameters of the model of
# >evolution?
#
# We are unable to obtain a unique maximum likelihood estimate (MLE) for both
# $t$ and $\lambda$.  However, we could try to obtain an estimate of the
# product of $\lambda$ and $t$ (i.e. the distance).

plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
savefig("_assets/phylocourse/submod/lhood3.svg") # hide

# ![](/assets/phylocourse/submod/lhood3.svg)

# What we've done here is fix the *total* substitution rate for each state at
# $1$ (i.e. $\lambda = 1/3$, giving a total substitution rate of $3\lambda =
# 1$) and then estimate the time $t$. In other words, we are estimating time,
# but on a *different time scale*, namely a scale of expected number of
# substitutions per site (i.e. one unit of this rescaled time corresponds to
# one expected substitution per site). This amounts to the same as estimating
# the product $\lambda t$ or in other words the **distance**.

# This example show how we can estimate distances for *any substitution model*
# (so also more complicated one than the JC model introduced above) given a rate
# matrix $Q$ that is scaled such that the average total substitution rate for
# each state is 1. We simply obtain the maximum likelihood estimate for the
# time parameter of the CTMC model keeping the rate matrix fixed using the
# transition probabilities of the CTMC and the observed data. However, as you
# might expect (or recall from the course notes), there is a simple formula for
# the distance under the JC model

# $$\hat{d} = -\frac{3}{4} \log \Big(1 - \frac{4}{3}p\Big)$$

# Where $p$ is the proportion of different sites in the two aligned sequences.
# We can use this formula to verify our graphical maximum likelihood approach
# above

distance_JC(p) = -0.75 * log(1. - 4p/3)
p = mapreduce(!=, +, seqa, seqb)/length(seqa)
distance_JC(p)

# It seems to work, let's verify more closely:

plot(0:0.01:1, d->ctmc_probability(seqa, seqb, d, 1/3),
    color=:black, legend=false, xlabel="distance", ylabel="log-likelihood")
dist = distance_JC(p)
vline!([dist])
hline!([ctmc_probability(seqa, seqb, dist, 1/3)])
savefig("_assets/phylocourse/submod/lhood4.svg") # hide

# ![](/assets/phylocourse/submod/lhood4.svg)

# Of course it worked. More on distances and distance-based phylogenetic
# inference in the [next section](../distance)!

# ## Exercises

# 1. Consider the figure below showing evolutionary histories for a single
#    nucleotide site for a set of three taxa. Compute numerically or write down
#    symbolically the probability of the observed evolutionary history in (1) and
#    (2).
# ![](/assets/phylocourse/submod/ex1.png)

# 2. I guess most of you have knowledge of some programming language, but
#    likely not julia. A good exercise -- for those who feel like it -- would be
#    to implement some of the bits of simulation and inference code above in your
#    programming language of choice (Python, R, Perl, ...).


# ------------------------------------------------------------------------------


# [^ternary]: While most julia code should be very easy to read even if you're unfamiliar, the [*ternary operator*](https://en.wikipedia.org/wiki/%3F:) used in this snippet might need some explanation (as Python does not have it for instance). It is simply a very concise `if/else` statement you can read as follows: `(if condition) ? (do this) : (else do that)`. It's a programming structure inherited from the `C` programming language.

# [^pr]: Some basic familiarity with probability theory is assumed here. Recall that $P(A|B)$ denotes the conditional probability of the event $A$ given the event $B$ (e.g. $P(\text{I have COVID}|\text{I cough all the time})$ can be read as 'the probability that I have COVID if it is the case that I cough all the time'). $P(A,B)$ is the joint probability of event $A$ and $B$, i.e. the probability that both $A$ and $B$ 'occur'. The product rule of probability theory relates these as $P(A,B) = P(A|B)P(B) = P(B|A)P(A)$, which also gives Bayes' rule $P(A|B) = P(B|A)P(A)/P(B)$. Note that this also works for more than two events, e.g. $P(A,B,C) = P(A,B|C)P(C) = P(A|B,C)P(B|C)P(C)$.

# [^sp]: Open any textbook on stochastic processes, and you will find that the mathematics of stochastic processes quickly becomes very involved.

# [^recursion]: Note that this recursive formulation is general, it will work for any discrete Markov model of DNA substitution, not ony the simple single parameter one. For the simple model we could probably find an easier formula for the $n$-step transition probability if we would do some maths here.

# [^allen]: A book on some basic stochastic processes theory which I think is quite nice is [Allen (2003)](https://www.worldcat.org/title/introduction-to-stochastic-processes-with-applications-to-biology/oclc/52092121).

# [^generator]: In the Markov chain literature this is often referred to as the *infinitesimal generator*, but in phylogenetics it's commonly reffered to simply as the rate matrix.

# [^solution]: Note that if you forget for a while that $P(t)$ and $Q$ represent matrices, and just look at this as a simple ordinary differential equation (ODE) of the form $$\frac{df(t)}{dt} = af(t)$$ the solution is very straightforward (you might recognize this form of ODE from population growth or radioactive decay models). If we bring $f(t)$ from the right hand side to the left, and integrate on both sides from $t_0$ to $t$, we get the very familiar solution $f(t) = f(t_0)e^{at}$. If we would naively apply this to our case we would have $P(t) = P(0)e^{Qt}$ and since we have $P(0) = I$ this happens to be the actual form of the solution.

# using Literate #src
# Literate.markdown(@__FILE__, "phylocourse", documenter=false) #src
