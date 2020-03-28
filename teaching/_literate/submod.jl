# \toc
# [back](/teaching/)
# # Substitution models for molecular phylogenetics
# If you've had a look at prof Van de Peer's course notes, you should have an
# idea already of the main approaches people have used to infer phylogenies
# using molecular data. The course notes provide a detailed discussion of (1)
# parsimony methods, which are character-based but do not involve a model of
# sequence evolution; and (2) distance-based methods, which are not character-
# based, but do involve a model of sequence evolution. Both approaches have
# been superseded by an increasingly statistical paradigm in molecular
# phylogenetics, where we adopt an approach that is both character- and
# model-based.

# In statistical phylogenetics, the approach is to assume a probabilistic model
# of sequence evolution, and to then somehow find the phylogenetic tree and
# parameters (of the model of evolution) that explain the sequence data. Note that
# this is deliberately formulated in a vague way, as the specifics of this procedure
# depend strongly on your philosophy of statistical inference. In particular
# two major approaches are used in phylogenetics, based on either the principle of [Maximum
# likelihood](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) or
# [Bayesian inference](https://en.wikipedia.org/wiki/Bayesian_inference).
# The cornerstone of both modern statistical phylogenetics and distance-based
# methods is the notion of a **substitution model**, and a good understanding of what
# such a model means is vital to develop a feel for modern phylogenetic methods.

# Using substitution models to model the evolutionary process probabilistically, we can answer questions that become progressively more interesting.

# > What is the probability that a given nucleotide substitutes into another over some time $t$ given the substitution model?

# > What is the probability that the sequence `ATCGAATCG` evolves into `ATCCCAATC` over some time $t$ given the substitution model?

# > What is the expected number of substitutions that happened during the course of evolution when we observe $x$ differences between two sequences given the substitution model?

# > What is the probability of observing a set of $n$ sequences given that they evolved along a phylogenetic tree $\mathcal{T}$ under a given substitution model?

# > Which phylogenetic tree connecting a set of $n$ sequences is the most likely to have resulted in the observed sequences given the substitution model?

# The goal of these notes is to give you a feel of how these questions can be answered and how they relate to phylogenetic inference.

# ## Discrete-time Markov models of nucleotide substitution

# ### The simplest model
# The simplest possible probabilistic model of DNA substitution could be specified in words as follows:
# > In a given time step of length $Δt$, a given nucleotide changes into random other nucleotide with probability $p$, it stays the same nucleotide with probability $1-p$.

# This is a simple probabilistic model of sequence evolution. We could for instance use this model to simulate the evolution of a given site in a sequence over a series of $n$ time steps: [^julia] [^ternary]

function simulate(site, n, p)
    print(site, " ")
    for i=1:n
        site = rand() < p ? rand(setdiff("ATCG", site)) : site
        print("⟶  $site ")
    end
end

simulate('A', 10, 0.4)

# The model of sequence evolution we specified is a simple kind of [stochastic process](https://en.wikipedia.org/wiki/Stochastic_process) belonging to the class of **discrete-time Markov chains (DTMC)**. A stochastic process can be thougth of intuitively as something that evolves randomly in time. Here this something is the state of the nucleotide of interest, and we denote the state of the nucleotide at time point $n$ by $X(n)$, taking values in the state space $\{A, T, C, G\}$.

# Note that the 'discrete-time' label refers to the fact that the evolutionary process is assumed to proceed in discrete time steps which can be indexed by integers, i.e. $n \in \mathbb{N}$ (continuous-time models, where we consider $X(t)$ with $t \in \mathbb{R}^+$, will be discussed below). In the simulation function above, we effectively iterate through the time steps, at each time point flipping a biased coin which lands heads with probability $p$ and tails with probability $1-p$. If the virtual coin lands heads, the nucleotide substitutes, if it lands tails it stays put.

# A sequence $(X(1), X(2), X(3), \dots)$ is called a *realization* of the stochastic process. To be clear:

function simulate(site, n, p)
    print("(X(0) = $site")
    for i=1:n
        site = rand() < p ? rand(setdiff("ATCG", site)) : site
        print(", X($i) = $site")
    end
    print(")")
end

simulate('A', 5, 0.4)

# shows a realization over a finite number of time steps of the Markov chain
# in our notation (given a particular initial state).

# ### Markov chains

# We called the above model a discrete-time **Markov** chain. The key feature of all stochastic processes with the **Markov property** is that at any given time, the future evolution of the process only depends on the current state and not the history of the process. In probabilistic notation[^pr] :
# $$ P\big[X(n+1)|X(n),X(n-1),X(n-2),\dots,X(0)\big] = P\big[X(n+1)|X(n)\big] $$

# The probability $P[X(n+1)|X(n)]$ is known as the **transition probability**.

# This property is what makes these models particularly mathematically tractable[^sp].
# Using the transition probabilities and the Markov property, we can compute the
# probability of a particular realization under the model. Consider for instance
# the simple model above with the parameter $p$ and an observed realization
# `A ⟶  G ⟶  G ⟶  T`. The probability of this realization under the model
# can be obtained as:

# \begin{align}
# P(&X_0=A, X_1=G, X_2=G, X_3=T) \\
# &= P(X_3=T|X_0=A,X_1=G,X_2=G)P(X_0=A,X_1=G,X_2=G) \\
# &= P(X_3=T|X_0=A,X_1=G,X_2=G) P(X_2=G|X_0=A,X_1=G) P(X_0=A,X_1=G) \\
# &= P(X_3=T|X_0=A,X_1=G,X_2=G) P(X_2=G|X_0=A,X_1=G) P(X_1=G|X_0=A) P(X_0=A)\\
# &= P(X_3=T|X_2=G)P(X_2=G|X_1=G)P(X_1=G|X_0=A)P(X_0=A) \\
# &= p \times (1-p) \times p \times 1 \\
# &= p^2(1-p)
# \end{align}

# Where I have written $X_n$ instead of $X(n)$ to avoid all those parentheses. While the equations are a bit convoluted, this is just a simple exercise in conditional probabilities (see also note [^pr]).

# But *why is being able to compute such a probability interesting?* Consider a question of the sort:

# >For a given homologous site, what is the probability to observe a nucleotide `A` in chimpanzee and nucleotide `T` in human, given that chimps and humans are separated by $n$ generations and the probability of a substitution in a single generation is $p$ (with all kinds of substitutions equally likely)?

# Clearly this starts to sound like something that can be of scientific interest. We cannot quite apply the calculations we did above directly though, since we have not specified the entire sequence of transitions like we did above. What we are effectively asking is the **$n$-step transition probability** $P(X_n=A |X_0=T)$. If we write the transition probability $P[X(n+1)=A|X(n)=T]$ (i.e. the probability of a `T ⟶ A` transition) more concisely as $p_{AT}$. The probability we are querying in the above question is of the form

# \begin{align}
# P(X_n=A &|X_0=T) = \sum_{x \in \{A,T,C,G\}} p_{Ax} P(X_{n-1}=x|X_0=T) \\
#   &= \sum_{x \in \{A,T,C,G\}} p_{Ax}
#           \sum_{y \in \{A,T,C,G\}} p_{xy} P(X_{n-2}=y|X(0)=T) \\
#   &= \sum_{x \in \{A,T,C,G\}} p_{Ax}
#           \sum_{y \in \{A,T,C,G\}} p_{xy} \sum_{z \in \{A,T,C,G\}} \dots
# \end{align}

# Nothing stops us from implementing this recursive formulation of the $n$-step transition probability in a computer program, however we can formulate this elegantly and concisely using a little bit of matrix algebra.

# ### The transition probability matrix

# Any DTMC model of nucleotide substitution is completely determined by its so-called **transition probability matrix** $P$:

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

# At this point, it's fairly easy to see that we could easily specify more complicated Markov models by introducing more parameters. For instance we could have different probabilities of a purine (`A` and `G`) to purine substitution versus a purine to pyrimidine substitution. But note that *however one specifies the model, each column of the transition probability matrix has to sum to one*. This simply means that the chain has to go *somewhere*, i.e. the total probability to be in *any* state in the next time step has to be one.

# For a given initial distribution over states $p(0)$, it is not hard to verify that the following matrix multiplication
# $$p(1) = Pp(0)$$
# gives probability the distribution over states at the first time point. For instance, consider the transition probability matrix
Pmatrix(p) = [1-p p/3 p/3 p/3 ;
              p/3 1-p p/3 p/3 ;
              p/3 p/3 1-p p/3 ;
              p/3 p/3 p/3 1-p ]
P = Pmatrix(0.2)

# and if we assume an initial state `T`, the initial distribution is (if we keep the `ATCG` order):

p₀ = [0. 1. 0. 0.]'

# i.e. the probability the nucleotide is in state `T` at time point 0 is one (this is just another way of expressing that the initial state is `T`). We obtain the probability distribution over states at the first time point as
p₁ = P*p₀

# This is simply the vector of transition probabilities $[p_{AT}, p_{TT}, p_{CT}, p_{GT}]^T$. Now, of course we can repeat this
# $$p(2) = Pp(1) = P(Pp(0)) = P^2p(0)$$

# and generally we get

# $$p(n) = P^np(0)$$

# which provides us with an easy way to compute the $n$-step transition probabilities. For instance, the vector

P^10*p₀

# gives the 10-step transition probabilities for the initial state `T`.

# Note that as time grows, the probability distribution over states converges to an equilibrium
for i in [1, 2, 5, 10, 20, 50, 100]
    @show round.(P^i*p₀, digits=3)
end

# This distribution is called the stationary distribution of the Markov chain, and in the simple model we defined above it is simply $[1/4, 1/4, 1/4, 1/4]$.

# ### A note on reversibility

# ### From individual sites to sequences
# If we make the assumption that sites in a sequence evolve independently
# (clearly a terrible assumption from a biological point of view, but a very
# convenient one statistically), the DTMC for a single site quickly generalizes
# to a model for entire sequences.

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

simulate("ATCGGGCGGGATTATTACGGAT", P(0.05), 10)

# Now to obtain the probability that one sequence evolves into another we can simply multiply the probabilities of the correspoding site-wise probabilities because of our independence assumption.

# ## The real deal: continuous-time Markov models of sequence evolution

# ### The rate matrix $Q$

Q(λ) = [-3λ   λ   λ   λ ;
          λ -3λ   λ   λ ;
          λ   λ -3λ   λ ;
          λ   λ   λ -3λ ]

# ### Transition probabilities

# ## So what are distances anyway?

# [^julia]: Code examples are in the julia programming language, which is a great language for scientific computing (which I prefer over Python or R, which you are probably familiar with). If you want to follow along, download the latest `julia` version from [https://julialang.org/](https://julialang.org/).

# [^ternary]: While most julia code should be very easy to read even if you're unfamiliar, the [*ternary operator*](https://en.wikipedia.org/wiki/%3F:) used in this snippet might need some explanation (as Python does not have it for instance). It is simply a very concise `if/else` statement you can read as follows: `(if condition) ? (do this) : (else do that)`. It's a programming structure inherited from the `C` programming language.

# [^pr]: Some basic familiarity with probability theory is assumed here. Recall that $P(A|B)$ denotes the conditional probability of the event $A$ given the event $B$ (e.g. $P(\text{I have COVID}|\text{I cough all the time})$ can be read as 'the probability that I have COVID if it is the case that I cough all the time'). $P(A,B)$ is the joint probability of event $A$ and $B$, i.e. the probability that both $A$ and $B$ 'occur'. The product rule of probability theory relates these as $P(A,B) = P(A|B)P(B) = P(B|A)P(A)$, which also gives Bayes' rule $P(A|B) = P(B|A)P(A)/P(B)$. Note that this also works for more than two events, e.g. $P(A,B,C) = P(A,B|C)P(C) = P(A|B,C)P(B|C)P(C)$.

# [^sp]: Open any textbook on stochastic processes, and you will find that the mathematics of stochastic processes quickly becomes very involved.
