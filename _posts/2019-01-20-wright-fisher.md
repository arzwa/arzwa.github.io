---
title: jaynesian probability - the wright-fisher model
layout: post
categories: [popgen, jaynes]
---

This is an attempt to apply principles from Jaynes' book on population genetics/evolutionary biology. This is perhaps one of the simplest poible applications, inbreeding in a Wright-Fisher population model.

Assume we have a population of $N$ diploid individuals, which randomly mate to form a new population of $N$ diploid individuals that replaces the old one. We follow one locus which starts with $2N$ different alleles. Random mating entails that every new individual has two alleles for the focal locus that are randomly drawn with replcement from the parental population. This process is repeated for an indefinite number of dicrete generations. Every generation we pick two genes from the population (irrespective of the individual, *i.e.* we treat the population as an urn from which we draw with replacement). Consider the following propositions:

$$ A \equiv \mathrm{pick\ exactly\ the\ same\ gene\ copy\ two\ times} $$

$$ B \equiv \mathrm{the\ two\ genes\ are\ identical\ by\ descent,\ i.e.\ they\ are\ the\ same\ allele} $$

Assume futher that we ssummarize all our knowledge about the world and the problem $X$ (the prior). We are now interested in the probability of picking two genes identical by descent in generation $t$, that is $P(B_t\|X)$.

First, we note that, by the principle of indifference (*i.e.* we do not discriminate between individual gene copies) and since we are drawing with replacement $P(A_t\|X) = 1/2N$ for all $t$. Of course $P(\bar{A_t} \| X) = (1 - 1/2N)$. We have the proposition

$$ B_1 = A_1 + \bar{A_1}B_0 $$

saying, that the proposition "the two picked genes in generation 1 are identical by descent" is the same as the logical disjunction of "the two picked genes are the same gene picked twice" and "the two picked genes are different gene copies that were already identical by descent in generation 0". Now, we define $B_0$ as false, since we do not perform the sapling in generation 0, we only start picking our pairs of genes in generation 1. 

As a result we have

$$ P(B_1|X) = P(A_1 + \bar{A_1}B_0|X) = P(A_1|X) + P(\bar{A_1}|B_0X) P(B_0|X) = P(A_1|X) = \frac{1}{2N} $$

For generation 2, we get 

$$ P(B_2|X) = P(A_2 + \bar{A_2}B_1|X) = P(A_2|X) + P(\bar{A_2}|B_1X)P(B_1|X) = \frac{1}{2N} + \bigg(1 - \frac{1}{2N}\bigg)P(B_1|X) $$

in general we have the recursion

$$ P(B_t|X) = \frac{1}{2N} + \bigg(1 - \frac{1}{2N}\bigg)P(B_{t-1}|X) $$

where we call $P(B_t\|X)$ the inbreeding coefficient in generation $t$.

$$ 1 - P(B_t| X) = \bigg(1-\frac{1}{2N}\bigg)(1 - P(B_{t-1}|X)) $$

which gives

$$ 1- P(B_t|X) = \bigg(1 - \frac{1}{2N}\bigg)^tP(B_0|X) $$

So it seemed to have worked.
