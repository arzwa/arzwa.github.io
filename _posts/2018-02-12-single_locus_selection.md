---
layout: post
title:  "Selection on a single diallelic locus"
date:   2018-02-24 16:04:50 +0100
categories: notes popgen evolution
published: true
---

We consider the effect of selection on the allele frequency at a single locus. We restrict ourselves to the diallelic case. We denote the relative reproductive advantage of one genotype over another by the relative fitness $w$. These relative fitness derive from differential selection of genotypes, governed by some selection coefficient $s$. The basic formula for the allele frequency after one generation of selection is the following:

$$ p' = \frac{p^2 w_{AA} + p(1-p) w_{Aa}}{p^2 w_{AA} + 2p(1-p)w_{Aa} + (1-p)^2w_{aa}} $$

with $w$ the relative fitness for the genotype in subscript and $p$ the allele frequency of $A$. $p'$ denotes the allele frequency of $A$ in the next generation. The denominator ensures scaling to the novel population size obtained after selection. we further refer to it as the average fitness of the population $\bar{w}$. It is also useful to express the mean relative fitness of the allele $A$ as:

$$ \bar{w_{A}} = p w_{AA} + (1-p) w_{Aa} $$

This can be interpreted as the average relative fitness of an individual, given that it received one $A$ allele from one parent (note the analogy with average effect in quantitative genetics).

We can find the change in allele frequency due to selection as:

$$ \begin{eqnarray}
\Delta p &=& p' - p  \\
         &=& \frac{p^2 w_{AA} + p(1-p) w_{Aa} - p\bar{w}}{\bar{w}}
         &=& \frac{p(\bar{w_A} - \bar{w})}{\bar{w}}
\end{eqnarray}$$

We can use the above formulae to determine for example the change in allele frequency of a fully recessive allele $A$ with frequency $p$ and selection coefficient $s$. The average fitness of the allele is given by

$$ \begin{eqnarray}
\bar{w_A} &=& p(1+s) + (1-p) \\
          &=& 1 + sp
\end{eqnarray}$$

The average fitness of the population by

$$  \begin{eqnarray}
\bar{w} &=& p^2(1+s) + 2p(1-p) + (1-p)^2 \\
        &=& 1 + sp^2 \\
\end{eqnarray}$$

The new allele frequency after one generation with selection is than given by

$$  \begin{eqnarray}
p' &=& \frac{p\bar{w_A}}{\bar{w}} \\
   &=& \frac{p(1+sp^2)}{1 + sp^2} \\
\end{eqnarray}$$

Resulting in a per generation change in allele frequency of

$$ \Delta p = \frac{sp^2(1-p)}{1+sp^2} $$

The expressions for different forms of selection (dominance, recessiveness, additivity, *etc.*) can be readily worked out in this way.
