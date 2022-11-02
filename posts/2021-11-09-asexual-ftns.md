---
title: Asexual FTNS
---
\newcommand{\ex}{\mathbb{E}}
\newcommand{\var}{\mathbb{V}}

We derive a form of the fundamental theorem of natural selection (FTNS) for
asexual organisms. Fisher's FTNS, which is not *that* fundamental, says that,
for a population consisting of $n$ types with constant relative fitnesses and
no mutation, (1) mean fitness is non-decreasing and (2) increases at a rate
proportional to the genetic variance in fitness (is this in fact what
*Fisher*'s version says?). We will here further assume generations are
non-overlapping and that there is asexual reproduction.

We refer to the expected number of offspring individuals of an individual of
type $i$ as the *fitness* of type $i$.  Let the fitnesses for the $n$ types be
$w_1, w_2, \dots, w_n$ and let $n_1(t), \dots, n_i(t)$ be the number of
individuals of each type in generation $t$, with $N(t) = \sum_i n_i(t)$ the
total population size.  The type frequencies are denoted by $p_i = n_i / N$.
Note that types could be alleles in haploids, genotypes in diploids and
polyploids, or something else, the only assumption at present being that there
is no recombination or mutation, so that types are faithfully inherited.  Our
goal is to find a solution for the change in mean fitness over time when
fitnesses are constant. We will have three steps:

(1) We first find a recursion for the (geno)type frequencies. We have
$$\Delta p_i = p_i' - p_i 
    = \frac{w_i n_i}{\sum_j n_j w_j} - p_i 
    = \frac{p_i (w_i - \bar{w})}{\bar{w}}$$ 
Where $\bar{w} = \ex[w] = \sum_j w_j p_j$ is the mean fitness. This
also holds when fitnesses are time-dependent (in which case we would
write $w_i(t)$).

(2) Now we find a recursion for mean fitness, we have
$$\Delta \bar{w} = \bar{w^\prime} - \bar{w} = \sum_i (w_i'p_i' - w_i p_i)$$
The goal now is to partition the change in mean fitness in a part due
to evolutionary change, and a part due to changes in fitness over
time. We can accomplish this by noting that
$$
\begin{align*}
    \Delta \bar{w} &= \sum_i ((w_i + \Delta w_i)p_i' - w_ip_i) \\
    &= \sum_i \Delta w_i p_i' + \sum_i \Delta p_i w_i \\
    &= \ex[\Delta w] + \sum_i \Delta p_i w_i
\end{align*}
$$
This shows that the change in mean fitness is equal to the mean of the
fitness changes and a component from evolutionary change (i.e. change
in genotype frequencies).

(3) We now plug our result from (1) in (2) to find
$$
\begin{align*}
\Delta \bar{w} &= 
    \ex[\Delta w] + \frac{1}{\bar{w}}\sum_i p_i w_i (w_i - \bar{w}) \\
    &= \ex[\Delta w] + \frac{1}{\bar{w}}(\sum_i w_i^2 p_i -
    \bar{w}\sum_i w_i p_i) \\
    &= \ex[\Delta w] + \frac{\var(w)}{\bar{w}}
\end{align*}
$$
Which is our desired result. Since the variance $\var(w)$ is
non-negative, we have that, when relative fitnesses are constant so
that $\ex[\Delta w] = 0$, mean fitness is non-decreasing. Furthermore,
mean fitness increases at a rate proportional to the genetic variance
in fitness.

**Some additional remarks.** 
We note that in the above, genotype frequencies $p_i$ and changes
therein only depend on ratio's of fitnesses of the form $w_i/\bar{w}$,
and not on absolute fitnesses, so we can always scale fitnesses when
dealing with genotype frequencies (and not absolute numbers of each
type in the population).
We can decompose (absolute) fitness in a viability and fertility
component. Viability $v_i$ will be the probability that an individual
of genotype $i$ reaches reproductive age. Fertility $f_i$ will be the
expected number of offspring from a type $i$ adult. We then have $w_i
= v_i f_i$, i.e. formally fitness (the expected contribution of an
individual of type $i$ to the next generation) is given by the law of
total expectation as
$$\ex[\text{#offspring}_i]=\ex[\ex[\text{#offspring}_i|\text{survival}_i]]
= f_i v_i$$ 
Assume that $w_1 > w_i$ for $i > 1$. We then always have $\bar{w} \le
w_1$ and it can be seen that, when we assume constant fitnesses
$$
\begin{align*}
p_1 &\rightarrow 1 \\
\bar{w} &\rightarrow w_1 \\
\var(w) &\rightarrow 0 
\end{align*}
$$
During evolution, the quantity $(w_1 - \bar{w})/w_1 = 1 - \bar{w}/w_1$
is minimized. This quantity is the *relative reproductive excess* of
genotype 1. With pure viability selection (i.e. all $f_i = f$),
$\bar{w}/w_1$ can be interpreted as the relative survival probability
of a random individual compared to the most fit genotype, and
the reproductive excess $1-\bar{w}/w_1$ can be interpreted as the
relative amount of selective mortality in the population compared to
the most fit genotype. This quantity is known as the *genetic load*.

---

Acknowledgement: the treatment here closely follows Nagylaki (1992)
