---
layout: post
title:  Recap - finite populations
categories: [qgen popgen evolution recap]
---

## Finite idealized population

A bit of a recap is never a bad thing to do. So considering a locus with two alleles $p$ and $q$ in a finite diploid population of size $N$, initially at frequencies $q_0$ and $p_0$. Assuming an idealized population, following all Hardy-Weinberg conditions except the infinite population size, the expected frequency of $q_1$ is then $q_0$. The variance of this quantity can be calculated as a binomial sampling variance with sample size of $2N$ and is given by:

$$ \sigma^2_{q_1} = \frac{q_0 (1-q_0)}{2N} $$

This can be interpreted as the variance in $q_1$ among different samples of size $N$ drawn from a base population with $q = q_0$. This is also the variance of the *change* in allele frequency $\sigma^2_{\Delta q}$.

We can also consider **genotype frequencies**, where the frequency of homozygotes in generation $t$ is of course given by $q_{t-1}^2$ and $p_{t-1}^2$. Consider the mean frequency of homozygotes among lines sampled from a base population with minor allele frequency $q_0$, which can be denoted $\overline{q^2}$. Recall that the expression for the variance of a random variable $X$ is given by $E(X^2) - E(X)^2$. Thus we have:

$$ \begin{eqnarray}
\sigma^2_{q} &=& \overline{q^2} - \overline{q}^2 \\
\overline{q^2} &=& \overline{q}^2 + \sigma^2_{q}
\end{eqnarray} $$

The consequences of this expression are clear. Indeed, $\overline{q}$ is just given by $q_0$, *i.e.* the allele frequency in the base population. So the expected number of homozygotes $\overline{q^2}$ is always in excess over the original fraction of homozygotes $q_0^2$. Likewise, the fraction of heterozygotes is reduced by $2\sigma^2_{q}$. Note that we have not considered the process in function of time. We have just observed a snapshot situation in a finite population and have related it to the base population from which it was originally drawn. By doing this we observed that the fraction of homozygotes is increased with the variance of the allele frequency $\sigma_q^2$. How $\sigma_q^2$ evolves with time is something we will consider next.

## Inbreeding

We will now consider the famous coefficient of inbreeding, which is a measure of the dispersive process happening in finite populations. The coefficient of inbreeding $F$ is the probability that two alleles at a randomly sampled locus are *identical by descent* (IBD). The average $F$ gives a measure of genetic drift in a population. Note that $F$ is only meaningful when compared to a base population!

Considering a base population with $2N$ genes which we consider all to be independent. We can easily see that $ F_1 = \frac{1}{2N} $. Similarly, we can reason that

$$ F_2 = \frac{1}{2N} + (1 - \frac{1}{2N})F_1 $$

This gives a recursion for $F$. However this recursion is as such not very easy to work with. We can easily see that there is always an addition of $\frac{1}{2N}$ to the inbreeding coefficient in subsequent generations, originating from the probability to sample identical alleles in the base population. Call this $\Delta F = \frac{1}{2N}$. Now we can write the following expression

$$ F_{t} = \Delta F + (1 -\Delta F)F_{t-1} $$

which can be rearranged to give

$$ \Delta F = \frac{F_{t} - F_{t-1}}{1-F_{t-1}} $$

which clearly shows the nature of $\Delta F$ as a measure for the rate of inbreeding. Indeed, it expresses the increase in the inbreeding coefficient normalized by the proportion of the population that was not yet inbred.

We now define the panmictic index $P_t = 1-F_t$ and we rearrange the last expression to give

$$ \frac{P_t}{P_{t-1}} = 1-\Delta F $$

which gives rise to a lovely recursion rendering

$$ \begin{eqnarray}
\frac{P_t}{P_{t-2}} = \frac{P_t}{P_{t-1}} \frac{P_{t-1}}{P_{t-2}} &=& (1-\Delta F)^2 \\
&\vdots& \\
P_t &=& (1 - \Delta F)^t
\end{eqnarray} $$

Which allows us to get a general expression for $F$ in function of $t$

$$ F_t = 1 - (1-\Delta F)^t $$

which holds of course only when the F_0 is considered as 0.

We can now revisit the expressions in the previous sections which were attained by considering properties of binomial samples. We have shown that $\sigma_{\Delta q}^2 = \frac{q_0(1-q_0)}{2N}$ which we can now express as $q_0(1-q_0)\Delta F = q_0p_0\Delta F$. Similarly we can consider the variance in gene frequencies among different hypothetical lines of a same base population in any generation as $\sigma_{q,t}^2 = q_0p_0F_t$. Whereas $\Delta F$ can be regarded as the rate of new inbreeding, $F_t$ is a measure of the accumulated effect of random drift over $t$ generations.

We can again consider genotype frequencies. The mean fraction of homozygotes across hypothetical lines of size $2N$ of a same base population is still given by

$$ \overline{q^2} = \overline{q}^2 + \sigma_q^2 $$

Where $\sigma_q^2$ can now be expressed in terms of the base population frequencies and the number of generations! Indeed, we have:

$$ \overline{q^2} = q_0^2 + q_0p_0F_t $$

This clearly shows the cumulating excess of homozygotes through inbreeding. Similarly, we can consider the decline in the fraction heterozygotes ($H$). Given $H_0 = 2p_0q_0$ we can see that

$$ \begin{eqnarray}
H_t &=& 2p_0q_0 - 2q_0p_0F_t \\
    &=& 2p_0q_0(1-F_t) \\
    &=& H_0(1-F_t) = H_0P_t
\end{eqnarray} $$

since the homozygote excess occurs at the expense of heterozygotes. Interestingly, this provides an interpretation for the panmictic index as the proportional decrease in the fraction of heterozygotes relative to the base population.

This is I think the core theory of inbreeding and basic population genetics in finite populations which one should always be able to reproduce on demand.
