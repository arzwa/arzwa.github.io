---
layout: post
---

Consider stabilzing selection on a quantitative trait, where phenotype $z$ is
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
