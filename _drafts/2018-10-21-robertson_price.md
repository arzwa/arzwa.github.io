---
title: The Robertson-Price identity
---

Consider the change in population mean $S = \mu_s - \mu$ due to selection on a
character $z$, where $\mu_s$ denotes the population mean after selection but
*before* reproduction.  Let $p(z)$ and $p_s(z)$ denote the frequency of
phenotype $z$ before and after selection respectively and $W(z)$ the abolute
fitness associated with the trait value $z$.  By definition of the mean we have 

$$ \mu_s = \int_z z p_s(z) dz $$

We need to find an expression for $p_s(z)$ to get further. $p_s(z)$ is of
course equal to $W(z)p(z)/\int_zW(z)p(z)dz$, where the denominator ensures that
$\int_z p_s(z) dz = \bar{W} = 1$ to form a proper density. We can define the
relative fitness as $w(z) = W(z) / \bar{W}$ so that we can write

$$ p_s(z) = w(z) p(z) $$ 

Which is a very intuitive expression relating the phenotypic distribution
before and after selection. We can now rewrite the expression for $mu_s$ as 

$$ \mu_s = \int_z z w(z) p(z) dz = \mathrm{E}(z w(z)) $$

Note that by definition $\mathrm{E}(w(z)) = \bar{w} = 1$ so that we have 

$$\begin{eqnarray}
 S &=& \mu_s - \mu \\
   &=& \mathrm{E}(z w(z)) - \mathrm{E}(z) \\
   &=& \mathrm{E}(z w(z)) - \mathrm{E}(z) \mathrm{E}(w(z)) \\
   &=& \sigma(z, w(z)) 
\end{eqnarray} $$

that is, *the selection differential* $S$ *is equal to the covariance of trait
value with relative fitness*. This result is known as the Robertson-Price
identity.

We now go to the offsping generation to investigate the effect of the selection
procedure.  We are interested in $\mu_o - \mu$, that is the change in trait
mean in the offspring generation due to selection. Considering the regression
of offspring trait values on the parent trait values with coefficient $\beta$,
we have that

$$ \Delta\mu = \mu_o - \mu = \beta(\mu_s - \mu) = \beta S $$

which is known as the Breeder's equation.


 
