
@def hascode = true

# A tour of Monte Carlo methods for sampling probability densities

'Monte Carlo' is an umbrella term for numerical methods that rely on
random numbers to obtain some desired mathematical object or result.

One area where Monte Carlo methods are used a lot is in Bayesian
statistics, where the main goal is to obtain expectations of arbitrary
functions with respect to some complicated posterior probability (target) density
$\pi(\cdot)$. If the exact form of the target is available, a convenient way
to approximate the desired expectations (of a function $f$, say) is to draw
a random sample of $n$ points and compute the arithmetic average

$$\mathbb{E}_\pi[f(X)] = \int f(x)\pi(x)dx \approx \frac{1}{n} \sum_{i=1}^n f(x_i)$$

This simple Monte Carlo method conveys an important idea: we can efficiently
approximate expectations under the target distribution given a random sample
from that distribution.

However, the problem is that often we can not obtain a random sample from the
target distribution. For instance, in all but the simplest Bayesian problems,
the posterior can only be evaluated up to a normalizing constant that is an
intractable high-dimensional integral.

In these notes I explore some Monte Carlo methods for sampling from such
unnormalized densities. The rules of the game are simple: (1) we have a density
function $\pi(\cdot)$ taking values in $\mathbb{R}^d$, (2) we can evaluate $\pi$
up to a normalizing constant, but that's all and (3) we wish to obtain a
random sample of $n$ realizations $x = (x_1, \dots, x_n) \in \mathbb{R}^d$ such
that $X \sim \pi$. Armed with such a sample, we can approximate expectations
of functions of $X$ with respect to $\pi$

```julia
using Plots, Random
Random.seed!(19081994)
```

## Importance sampling

## Rejection sampling

## Markov chain Monte Carlo

Markov chain Monte Carlo (MCMC) is a very powerful and generic approach
for sampling from complicated distributions. A vast literature exists on the
subject, and many variants and implementations of the general methodology
have been developed.

It is useful to partition the class of MCMC algorithms in two, based on
whether they use information from the differential structure of the target
or not.

## Random-walk Metropolis-Hastings MCMC

## Hamiltonian Monte Carlo

------------------------------------------------------------------------------

Suppose we have some density $\pi(x)$ for a random element $X$ taking values
$x \in \mathbb{R}^d$, and further suppose we can evaluate it up to a
normalizing constant

```julia
π(x) = sum(x .^2) <= 1. && all(x .> 0.) ? x[1]^3*x[2] : zero(x[1])
```

Note that `π(x)` coded up above is the *unnormalized version*. We assume we
can evaluate the unnormalized density, but have no further analytical tools
at our disposal (e.g. we don't know how to integrate the density). The only
thing we can do are queries like:

```julia
π([0.9, 0.1])
```

but that's it. Our goal is to obtain expectations of arbitrary functions $f$
of $X$ under the density $\pi(x)$. For most random elements $X$ this entails
that we wish to evaluate integrals of the form:

$$\mathbb{E}_\pi[f(x)] = \int f(x)\pi(x)dx $$

For a low-dimensional problem like the one we have here, a first obvious thing
to do is evaluate the density on a grid in $\mathbb{R}^2$. We use this approach
first to visualize the density:

```julia
mn, mx, step = 0, 1, 0.01
plot(surface(mn:step:mx, mn:step:mx, (x,y)->π([x,y]), grid=false),
    contourf(mn:step:mx, mn:step:mx, (x,y)->π([x,y])), size=(800,300))
savefig("assets/unnormalized-density.svg")
```

![Visualization of the unnormalized density on a fine grid.](assets/unnormalized-density.svg)

But this approach can also be used to obtain integrals, in which case it is
referred to as 'quadrature'. For instance we can compute the normalizing
constant as follows

```julia
Z = 0.
for x₁=mn:step:mx, x₂=mn:step:mx
    global Z += π([x₁, x₂])*step^2
end
1/Z
```

The full normalized density is actually $\pi(x_1, x_2) = 24 x_1^3 x_2$
for $x_1^2 + x_2^2 \le 1$ and $x_1, x_2 > 0$, so we see (surprise!) that numerical
integration using quadrature works. We can use quadrature in this way to obtain
expectations of any function of the random vector $X$. For instance the mean
of the first coordinate $\mathbb{E}_\pi[X_1] = \int x_1 \pi(x) dx$

```julia
E, Z = 0., 0.
for x₁=mn:step:mx, x₂=mn:step:mx
    p = π([x₁, x₂])*step^2
    global Z += p
    global E += x₁ * p
end
E/Z
```

## Some more 1D intuition

Assume we have the following unnormalized density, which we can only evaluate

```julia
π(x) = 0.3*exp(-(0.5-x)^2) + 0.7*exp(-(3-x)^2)
```

It looks like this

```julia
plot(x->π(x), xlim=(-2.5,6))
```

a functio of interest to get the expectation of is

```julia
f(x) = log(abs(x)+1)
```

To compute an expectation, we can perform numerical integration

```julia
mn, mx, step = -2.5, 6, 0.2
p = plot(x->π(x), xlim=(-2.5,6), color=:black, grid=false)
Z, E = 0., 0.
for x=mn:step:mx
    px = π(x)
    bar!([x], [px], bar_width=step, fillalpha=0, legend=false)
    global Z += step*px
    global E += step*f(x)*px
end
display(p)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

