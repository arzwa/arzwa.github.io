---
title: Evolutionary analysis - phylogenetic tree inference using maximum likelihood and Bayesian inference
author: Arthur Zwaenepoel
layout: post
---

# Software requirements

- IQ-TREE 1.6 ([http://www.iqtree.org/](http://www.iqtree.org/))
- MrBayes 3.2.6 ([http://nbisweden.github.io/MrBayes/](http://nbisweden.github.io/MrBayes/download.html))
- RevBayes 1.0 ([https://revbayes.github.io/](https://revbayes.github.io/software))

# Notation

- $P(A)$: probability of event $A$
- $P(A\|B)$: probability of event $A$ conditional on event $B$
- $p(X)$: **probability mass/density function** (pmf/pdf) of the **discrete random variable** $X$ (note normally a pmf is denoted using $p$ and a pdf with $f$, but here we will lump the two together)
-  $\mathcal{L}(\theta\|D,\mathcal{M})$ **likelihood** of **parameter** (vector) $\theta$ given observed **data set** $D$ and the **model** $\mathcal{M}$

# Maximum Likelihood

## Statistical principles

- $\mathcal{L}(\theta\|D,\mathcal{M}) = P(D\|\theta,\mathcal{M})$
- Discrete case: e.g. coin tossing
- Continuous case: e.g. mean and standard deviation of a Normal distribution

## Phylogenetic likelihood

- The rate matrix $Q$, and its relationship to the exponential distribution
- Markov models and the transition probability matrix $P = e^{Qt}$
- Felsenstein's pruning algorithm

## Heuristic tree search under the maximum likelihood principle

- IQ-TREE

## Model selection

- LRT, AIC, BIC (IQ-TREE)

# Bayesian inference

## Statistical principles

- $p(\)
