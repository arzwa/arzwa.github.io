---
title: Neutral evolution
author: Arthur Zwaenepoel
---

So I bought Lynch & Walsh ginormous book 'Evolution & Selection of Quantitative
Traits', and as I would like to study it *thoroughly*, I'll keep something of a
study diary on it.  This book, with it's bible-like impression, a book of
revelation so it seems, a work of titans deserves a profound study, despite
it's overwhelming magnitude. Indeed, the question is probably rather whether I
deserve the book or have the brains for it. Anyway, I'm not discussing the
first chapter, although it is very interesting, and here I will be delving in
the second chapter, which gives an impression of being quite basic, but
nevertheless touches on some advanced matters (for me at least). The chapter is
on neutral evolution in one and two-locus systems, and is evidently of great
importance for evolutionary genetics in general. It is concerned with genetic
drift, mutation and migration.

## The Wright-Fisher model

The well-known Wright-Fisher model is obviously a good place to start. It
assumes

- a constant population size $N$
- monoecious individuals
- random mating
- discrete generations

and obviously for now, no selection, mutation or migration. Consider one
allele/gene in this population, present in $i$ copies, the probability of
having $j$ copies in the next generation is given by the binomial distribution
$Binom(n, p)$ with $n = 2N$ (that is the number of genes to sample to get the
next generation) and $p = i/2N$ (the frequency of the gene we are following).
Specifically, we have

$$ P_{ij} = \binom{2N}{j} \Big(\frac{i}{2N}\Big)^{j}\Big(1 -
\frac{i}{2N}\Big)^{2N-j} $$

This $[P_{ij}]$ values expressing the probability of transition from copy
number $i$ to $j$ form a transition probability matrix $P$ with dimensions
$(2N+1) \times (2N+1)$. Note that all $P_{0.}$ values except $P_{00}$ are 0 and
likewise for $P_{1.}$. Also note that all rows sum to one. The probability
distribution over the copy number of the allele of interest at time $t$ is
given by

$$ x(t) = x(0)P^t $$

with $x$ a probability vector for the number of copies of the allele under
consideration. As indicated, there are two absorbing states in the model,
since there is no mutation or migration there is no source of novel variation
and an allele will either get lost or fixed eventually. In terms of the probability
vector over copy number values $x$, this translates to saying that the sum
$x_{2N} + x_{0}$ converges to one as $t \rightarrow \infty$.

The variance in alele frequencies attained conditional on the allel frequency
in the parental population can serve as a measure for the dispersion in allele
frequency induced by random sampling in the Wright-Fisher model. After one
generation of random mating the expected number of allel copies in the
population is given by the Hardy-Weinberg proportions, *i.e.*

$$ \mathbb{E}(N) = p^2 \times 2 + 2p(1-p) \times 1 + (1-p)^2\times 0 = 2p$$

and likewise

$$ \mathbb{E}(N^2) = p^2 \times 4 + 2p(1-p)\times 1 = 2p(p+1)$$

giving

$$ Var(n) = 2p(p+1) - 4p^2 = -2p^2 +2p = 2p(1-p) $$

where $n$ is for copy number. To be clear, this is the variance in the number
of copies of the allele of interest in an individual. To get the variance for
the allele frequency in the population, we need the variance in the total
number of copies in the population divde by the number of alleles in the
population, that is

$$ NVar\Big(\frac{n}{2N}\Big) = \frac{2Np(1-p)}{(2N)^2} = \frac{p(1-p)}{2N} $$

Which is the sampling variance of a binomial random variable $Binom(2N, p)$.

So why is this interesting? It shows that the expected fraction of
heterozygotes is a measure forthe variance in allele frequency, and
that we cn track heterozygosity to get  grip on the effects of drift
on the variation in the population.

## Loss of heterozygosity in a Wright-Fisher population

Now that we know that heterozygosity can serve as a measure for the effect of
drift in a Wright-Fisher population, it might be interesting to study this in
more detail. Considering an initially fully unrelated population, the probability
that two genes are identical by descent in the next generation is given by

$$ f(1) = \frac{1}{2N} + \Big(1 - \frac{1}{2N}\Big) $$

or in an arbitrary generation

$$ f(t) = \frac{1}{2N} + \Big(1 - \frac{1}{2N}\Big)f(t-1)$$

This is of course the expected inbreeding coefficient in generation $t$ under
random-mating in a Wright-Fisher population of size N. This is closely related to the heterozygosity in the population. Indeed, an individual can only be a heterozygous for a locus if the alleles are not IBD, which occurs with probability $1-f$. In a random mating population, the fraction of heterozygotes $2p(1-p)$ is reduced proportional to the amount of inbreeding in the population. This fractional reduction $1-f$ we term the heterozygosity, for which we can find a nice recursion:

$$\begin{eqnarray}
H(t) &=& 1 - \frac{1}{2N} - \Big(1-\frac{1}{2N}\Big)f(t-1) \\
     &=& \Big(1 - \frac{1}{2N}\Big)H(t-1) \\  
     &=& \Big(1 - \frac{1}{2N}\Big)^tH(0)
\end{eqnarray}$$

giving an expression for the *expected heterozygosity* at time $t$ relative to some base population with heterozygosity $H(0)$.

$$ H(t) \approx e^{\frac{-t}{2N}}H(0) $$

This expression can be used to calculate expected times for particular reductions etc.
