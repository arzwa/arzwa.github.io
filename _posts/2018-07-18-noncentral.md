---
layout: post
title:  Non-central Chi-squared
categories: notes
---


The distribution of the square of a normal but non standard-normal random variable $X \sim \mathcal{N}(\mu, sigma^2)$
is given by a scaled non-central $\chi^2$ distribution. That is the $X$ = $\sigma^2 Z$ with $Z$ a noncentral $\chi^2$
distribution with scale parameter $\lambda = (\frac{mu}{\sigma})^2$.

```python
import numpy as np
import seaborn as sns

mu, si = 0.3, 0.2
z = np.random.normal(mu, si, 10000)**2
z_ = si**2 * np.random.noncentral_chisquare(1, mu**2/si**2, 10000)
ax = sns.kdeplot(z_, shade=True, linestyle="--", color="forestgreen", label="Non-central $\chi^2$")
ax = sns.kdeplot(z, color="k", linestyle="-", label="Square of normal")
sns.despine()
```

<center><img src="{{ "/assets/img/noncentralchisq.png" | absolute_url }} " width='400px'></center>

Consider stabilizing selection on a quantitative trait, where phenotype $z$ is
related to fitness $w$ as follows:

$$ w = e^{-\frac{1}{2} \beta(z-\theta)^2} $$

We are interested in finding the expected value of $w$, given the random
variable $z$, which follows a Normal distribution with known mean
$\overline{z}$ and variance $\sigma^2$. Obviously we know that $z - \theta$
also follows a Normal distribution. Call this r.v. $X$ with mean $\mu$. Now we
want to find the distribution of $X^2$. This should have something to do with a
$\chi^2$ distribution.

It turns out that $X^2$ has the same distribution as the r.v. $\sigma^2 Y$
where $Y$ follows a non-central $\chi^2$ distribution with one d.f.
$\chi^2_{k=1}((\frac{\mu}{\sigma})^2)$.

Now we can find $E(w) = E(e^{-\frac{1}{2} \beta \sigma^2 y})$ with $y = (z -
\theta)^2$ by using the moment generating function $e^{tY}$ where we $t =
-\frac{1}{2} \beta \sigma^2$ of the non-central $\chi^2$ distribution.
