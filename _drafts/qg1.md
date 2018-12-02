---
---

Consider a strictly additive diploid population following the following model:

```
bb                   Bb                    BB
 |--------------------|---------------------|
 0                    a                    2a
```

The population is in Hardy-Weinberg equilibrium with allele frequencies $B$:
$p$ and $b$: $(1-p)$. The expected genotypic value for the population is

$$ E(G) = 2a \times p^2 + a \times 2p(1-p) + 0 \times (1-p)^2 = 2ap $$

To calculate the variance $\sigma_A^2 = E(G^2) - E(G)^2$ we alse need 
$E(G^2)$ which is simply

$$ E(G^2) = 4a^2 \times p^2 + a^2 \times 2p(1-p) = 2pa^2(2p + (1 - p)) = 2pa^2(p+1)$$

Combining these we get

$$ \sigma_A^2 = 2pa^2(p+1) - 4p^2a^2 = 2pa^2(p + 1 - 2p) = 2a^2p(1-p) $$

Note that $2p(1-p)$ is simply the heterozygosity in the population. As expected,
the additive genetic variance (which is simply the genetic variance in this
additive model of course) increases with increased variance in allele
frequencies, that is increased heterozygosity. 

This is readily extended to the multilocus bi-allelic case.
