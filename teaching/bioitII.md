@def title = "ml"
@def author = "Arthur Zwaenepoel"
@def hasmath = true
@def hascode = true
\toc

# Phylogenetic inference using Maximum Likelihood

Below, a rough sketch of the theory underlying tree inference by means of
maximum likelihood is given. This is not a complete nor systematic exposition,
and only serves the purpose of providing some intuition for the methods,
without delving too much in the technical details.  

## Substitution models, again, but viewed slightly differently

Last time, we saw how to infer phylogenetic trees from a multiple sequence
alignment using Maximum Parsimony (MP) and distance matrix based methods, using
the PHYLIP programs of phylogenetics pioneer Joe Felsenstein. Nowadays, with
more powerful computers, these methods are rarely used in practice[^mp], and most
phylogenetic trees you will encounter in the literature will be inferred either
using methods that employ the principle of Maximum Likelihood (ML) or Bayesian
inference. Both approaches are related in that they treat the problem of
inferring a phylogenetic tree as a problem of **statistical inference**. While
we cannot discuss the mathematical and statistical details of these methods in
detail here, a minimal exposition of the underlying theory and principles will
help to acquire some intuition on what these methods do, and how you can
properly employ them to compute a phylogenetic tree of your species or
gene family of interest.

The **substitution models** (for example, the Jukes-Cantor model) we used in
distance based methods are actually probabilistic models[^math] of how we think
the evolution of molecular sequences works. They are defined by a **rate matrix** $Q$, which represents the relative rates of substitution between character states. For instance, the rate matrix for the Jukes-Cantor model is given by

$$Q_{JC} = \begin{bmatrix}
    - & \lambda_{GA} & \lambda_{CA} & \lambda_{TA} \\
    \lambda_{AG} & - & \lambda_{CG} & \lambda_{TG} \\
    \lambda_{AC} & \lambda_{GC} & - & \lambda_{TC} \\
    \lambda_{AT} & \lambda_{GT} & \lambda_{CT} & -
\end{bmatrix} = \begin{bmatrix}
    - & \lambda & \lambda & \lambda \\
    \lambda & - & \lambda & \lambda \\
    \lambda & \lambda & - & \lambda \\
    \lambda & \lambda & \lambda & -
\end{bmatrix} $$

Where e.g. $\lambda_{AT}$ is the relative rate of substitution from `A` to `T`.
As you now, in the Jukes-Cantor model all relative rates are equal, that is,
the rate of substitution from e.g. `A` to `T` is the same as the rate from `A`
to `C`, but also `C` to `G`, `T` to `C`, etc. In distance-based methods, these
models allowed us to convert an observed number of substitutions $p_o$ into an
**expected number of substitutions** $\mathbb{E}[p]$ that actually occurred
given our modeling assumptions. This quantity, $\mathbb{E}[p]$ is what we
called the *distance* $d$.

In ML methods we will *use the exact same substitution models* in a slightly
different way: instead of using them to estimate the expected number of
substitutions based on the observed number of substitutions (which is equivalent
to saying 'correcting for multiple substitutions'), we will use them to determine
the *probability that a particular character substitutes to another character
over some time $t$*. That is, we can employ substitution models to answer
questions of the form:

>What is the probability that a nucleotide `A` substitutes into a `T` over 5
>million years under the Jukes-Cantor model with a substitution rate of 0.1
>substitution per million years?

In mathematical notation, this is a so-called **transition probability** of the
associated substitution model. It is the probability that over the time $t$ we
are considering, the `A` evolved into a `T`, but considering all possible ways
this particular evolutionary event might have happened (*i.e.* $\texttt{A}\rightarrow
\texttt{T}$ or $\texttt{A}\rightarrow \texttt{G} \rightarrow \texttt{T}$ or
$\texttt{A}\rightarrow \texttt{T} \rightarrow \texttt{A} \rightarrow
\texttt{T}$, etc.) For the Jukes & Cantor model, this transition probability
has a simple form[^complexmodels]:

$$ P_{AT}(\Delta t) = P\big[X(t+\Delta t) = \texttt{T} | X(t) = \texttt{A} ; Q_{JC}\big] = \frac{1}{4} - \frac{1}{4}e^{-4\lambda \Delta t} $$

where $\Delta t$ is the time interval over which the site evolves (in the
example, $\Delta t = 5$) and $\lambda$ is the substitution rate, which we assumed to
be $0.1$. Note that $\lambda\Delta t = d$ is the distance (this is just physics: rate
$\times$ time = distance). If we compute this we find that the probability to
find a `T` after 5 million years at this site which used to be an `A` is about
0.216.

If we now assume that different sites in the alignment **evolve independently**
we can answer the following kind of questions:

>What is the probability that a sequence `ATTTGATGAAACGTGC` evolves into a sequence
>`AGTCCACGAAAAATGC` over 10 million years under the Jukes-Cantor model with a
>substitution rate of 0.1 substitution per million years?

Indeed, basic probability theory says that if we assume that the sites evolve
independently (a quite horrible assumption, but one that *really, really* makes
the mathematics much easier), the probability we are considering in the above
question is simply given by:

$$ P(\texttt{ATT...C}\rightarrow\texttt{AGT...C}) = P_{AA}(\Delta t) \times P_{TG}(\Delta t) \times \cdots P_{CC}(\Delta t) $$

Once we can answer questions like the above, we can do the same for a bunch of
sequences related by a phylogenetic tree by using a device called
"Felsenstein's pruning algorithm". That means we can also answer questions of
the sort:

>What is the probability that an ancestral sequence `AAT` evolved into `AAA`
> for human, `AAT` for chimpanzee and `ATT` for gorilla given that human and
> chimp diverged a time $t_1$ ago and the common ancestor of chimp and human diverged
> from the lineage that gave rise to gorilla a time $t_2$ ago, assuming the
> Jukes-Cantor model and a substitution rate of $r$ substitutions per time unit?

Both the mathematics that allow us to actually compute the transition
probabilities for more complicated substitution models, as well as the
algorithmic details that allow us to compute these probabilities on general
unrooted trees are slightly too hairy to consider in this workshop[^books], but
simply take away the following crucial corollary:

>For a given data set of aligned sequences ($D$), a substitution model[^ctmc] ($Q$),
>and given parameters[^theta] ($\theta$) of the substitution model, we can compute
>the probability that the data $D$ evolved over a phylogenetic tree ($T$) that
>describes the evolutionary relationships of those sequences.

Stated more concisely in the language of probability theory, this means we can compute $P(D|T,Q,\theta)$[^condp]. This quantity, viewed as a function[^lhood]
 of $T, Q$ and $\theta$, is called **the likelihood**.

$$ L(T, Q, \theta|D) = P(D|T, Q, \theta) $$

![$L(T, Q, \theta|D) = P(D|T,Q,\theta)$ in a picture. In Maximum likelihood and
Bayesian inference, the problem of finding the phylogenetic tree is treated as
a problem of statistical inference. Of central importance in both approaches to
the problem of phylogenetic inference is the probability of the data given the
model, tree and parameters $P(D|T,Q,\theta)$.](p.pdf)

## Maximum likelihood

By now you are thinking: "All fine that we can compute $P(D|T,Q,\theta)$, but
last time I checked our goals were to find $T$, and to compute this probability
thing, we actually need $T$, which is the very stuff we are trying to find?".
Well, you're quite right! Nevertheless, we can still use the fact that we know
how to compute $P(D|T,Q,\theta)$ to find a $T$ that we like.

In parsimony based tree inference, we considered an **optimality criterion**,
and used that criterion to find the tree that provides the best possible
explanation of the data. The optimality criterion was of course the parsimony
criterion:

>**Parsimony criterion:** Trees that correspond to less character changes
>provide more parsimonious explanations of the sequence data, and are
>therefore better.

Given the optimality criterion, we can start looking in tree space for the best
tree (although this sounds very simple in principle, this is actually very tough
computationally, but that's for computer scientists to worry about).

The problem with the parsimony criterion is of course that it does not account
for unobserved substitutions, and as a result underestimates the amount of
evolution  that separates two sequences. We now consider an alternative
optimality criterion we can use to find the best tree in tree space:

>**Maximum Likelihood (ML) criterion:** Trees that correspond to a higher
>likelihood provide more likely explanations of the sequence data
>(given the substitution model) and are therefore better.

This is a very general principle in statistics: we specify a probabilistic
model, and try to find the parameters of the model for which the data was most
likely to be observed (assuming the model). For example, consider we assume that
weight is normally distributed in the human population (this is our
probabilistic  model), and we gather some data $x$ by weighing $n$ individuals.
Now we'd like to estimate the location parameter ($\mu$, or the mean) of our
probabilistic model by looking for the value of $\mu$ that makes the observed
data most likely under the Normal distribution assumption with mean $\mu$. [As it
turns out, in this simple problem the Maximum likelihood estimate of the
$\mu$ of our probabilistic model is just the sample average $\bar{x} = \sum_{i=1}^n x / n$ of our observed data.]

Inferring trees under the ML criterion therefore consists of **searching tree
space to find the tree that gives the highest likelihood value**, similarly
to what we do in maximum parsimony tree inference. However, since the likelihood
is based on a model of sequence evolution (the substitution model $Q$), we
naturally account for the problem of multiple substitutions, whereas under
the parsimony criterion we did not. That is, while maximum parsimony assumes that the observed substitutions
reflect the true evolutionary history, phylogenetic inference using
substitution models (both distance and ML based) account for our uncertainty
with regard to the true evolutionary history based on modeling assumptions (the
substitution model). Importantly however, distance based methods collapse all
evolutionary signal in the sequences to a single number (the distance), whereas
ML tree inference, like parsimony, is a character-based method, employing
the information in individual sites of the alignment.

As a last remark, note that we not only can look for the tree that maximizes
the likelihood, but we can also try to find the parameters of the substitution
model that give the highest likelihood value. So for example, we can try to
find the value of $\alpha$ for the Gamma model of rates across sites that gives
the highest likelihood value. In general mathematical terms, tree inference by
ML can be concisely stated as the following optimization problem:

$$ \hat{T}, \hat{\theta} =  \text{argmax}_{T, \theta} L(T, Q, \theta|D) $$

which can be read as "estimate the tree and parameters of the substitution
model by maximizing the likelihood function with respect to the tree and the parameters". These estimates are called **ML
estimates (MLEs)**, and they have many interesting statistical properties.

## Model selection

![Some (certainly not all) DNA substitution models and their relationships. The
simplest model (Jukes-Cantor) is shown on top, while the most complex (GTR+I+G:
General Time Reversible with invariant sites and Gamma distributed rates across
sites), of which all others are special cases, is shown below. An edge represents the addition of a parameter (I: proportion of invariant sites, G: Gamma distributed rates across sites, F: unequal base frequencies, S: unequal relative substitution rates.)](models.png)

Whereas ML allows us to infer a phylogenetic tree and the parameters of the
substitution model, we still need to decide ourselves, often rather
arbitrarily, which substitution model we are using. This is quite problematic,
because there are a very large number of possible models one could choose from
(see figure 2). Of course, this is only natural, after all, when we are modeling
something we have to decide on our modeling assumptions!  Nevertheless, given a
candidate set of substitution models, one could reason that we simply extend
the ML principle, and choose the model which gives the highest likelihood.

This is however problematic, since more complicated models will always result in
a higher likelihood, so the ML principle when applied to models will always
result in the selection of the biggest, most complicated and most flexible
model. Again, this is a general issue in statistics, known as **overfitting**.
The key issue is that we actually want to know whether some complex model fits
the data *significantly* better than a more simpler model.

![An illustration of the continuum of simple to more complex models and the
problem of overfitting.  The data is represented by black dots. From left to
right: the first model is the simplest possible model, namely a constant. It
has a single parameter $c$, and the model can be written as $f(x) = c$. The
second model is a simple linear model, $f(x) = a + bx$, which has two
parameters $a$ and $b$, being respectively the intercept and slope of the
fitted line. The third is a quadratic model, which can be represented as $f(x)
= a + bx + cx^2$ and has three parameters. The last model shows a very complex
model which perfectly fits the data. The likelihood of the model given the data
will increase for increasing model complexity (i.e. from left to right in this
example, from top to bottom in Figure 2).](overfitting.png)

Luckily, since this is a common problem in statistics, there exist devices that
can be used to select models of varying complexity. The three main methods
that are used for model selection in a ML setting are:

1. The likelihood ratio test (LRT)
2. The Akaike information criterion (AIC)
3. The Bayesian information criterion (BIC)

These are all quantities that can be used to do pairwise comparisons of models,
and they can be computed from the likelihood values under the different models
that you wish to compare. For more insight in the statistical details of these
methods for model selection, please refer to wikipedia or your favourite
statistical modeling textbook.

We briefly consider the **LRT** here, which works only for nested models. Two
models are nested if one of the two is a special case of the other (i.e. models
on different rows in Figure 2). For example the Jukes & Cantor (JC) model is
nested in the K2P model since we obtain the JC model if we set the
transition/transversion ratio $\kappa$ in the K2P model to 1 (i.e. we assume
transitions and transversions are equally likely). For models that are nested in
this way, we can compute the following test statistic

$$ D = -2\mathrm{log}\Bigg( \frac{l_{1}}{l_{0}}\Bigg) = -2\big(\mathrm{log}(l_1) - \mathrm{log}(l_0)\big) \sim \chi^2_{k_1 - k_0}$$

where

- $l_1$ is the likelihood under model 1 with $k_1$ parameters
- $l_0$ is the likelihood under model 0 with $k_0$ parameters
- $k_1 > k_0$

so $l_1$ is the likelihood for the more parameter-rich model (e.g. K2P) and
$l_0$ is the likelihood for the other model (e.g. JC). If the test statistic is
bigger than the critical value (at the desired significance level) of the
corresponding $\chi^2$ distribution, then this indicates that the more complex
model provides a statistically significant better fit to the sequence data. If
not, than we may conclude that the difference in the likelihoods is not
statistically significant, indicating that the simpler model does about an
equally good job at explaining the data. The **BIC and AIC** criterion are similar devices, which can also be computed from the likelihood values and the numbers of parameters, and which can be employed for non-nested models as well.

[^mp]: This is not entirely true, since maximum parsimony and distance methods are still used *in* maximum likelihood and Bayesian inference, e.g. to quickly provide a good starting guess of the tree.

[^math]: For the mathematics aficionados, they are a specific subclass from a family of stochastic processes referred to as "Continuous time Markov Chains" (CTMC).

[^complexmodels]: Obtaining transition probabilities for more complicated substitution models requires numerical evaluation of matrix exponentials of the form $e^{Qt}$, where $Q$ is the substitution rate matrix and $t$ is time.

[^books]: Some good books on the mathematical and algorithmic details of phylogenetic inference (using any method) are Ziheng Yang's "Computational Molecular evolution: A statistical approach" and Felsenstein's "Inferring phylogenies".

[^ctmc]: A mathematically tractable one of course, which are mostly time-reversible continuous Time Markov Models (CTMC).

[^theta]: The set of parameters $\theta$ consists of stuff like the substitution rate (in all models), the ratio of the probability that a substitution is a transition over the probability that a sustitution is a transversion (often called $\kappa$ for e.g. the K80 model), the $\alpha$ parameter in the Gamma distribution model of substitution rates across sites etc.

[^condp]: Read this as "the probability of the data given the tree, the substitution model, and the parameters of the model".

[^lhood]: This change of view, where we look at $P(D|T,Q,\theta)$ as a function of $T, Q$ and $\theta$ instead of $D$ entails some subtleties. For instance $L(T,Q,\theta|D)$ is no longer a probability density, and therefore does not have to integrate to one over its domain.


# Exercises: Inferring ML trees with IQ-TREE

For the exercises, we will use the computationally *very* efficient software
IQ-TREE (Nguyen *et al.* 2015). This is a command line program, similar to the
PHYLIP programs but not interactive. Please download IQ-TREE for your operating
system at [http://www.iqtree.org/](http://www.iqtree.org/)

In the exercises below, commands for IQ-TREE will be shown in monospaced
font `like this`.

## Basic tree inference

### The Jukes-Cantor model

Let us first reconsider the 45 species 18S rRNA data set. Remember that we had
some problematic clades in our distance-based phylogenies, where we found
different topologies depending on whether we used distances computed with the
Gamma model of rate heterogeneity across sites or not and depending on the
$\alpha$ parameter we chose. Let's see what we get when using ML.

Put the `18SrRNA_45.txt` file in the directory with the IQ-TREE executable.
Fire up IQ-TREE and run the following command:

```bash
iqtree -s 18SrRNA_45.txt -m JC -pre JC
```

The `-s` option is used for specifying the input multiple sequence alignment
file, the `-m` option is used for setting the substitution model (here Jukes &
Cantor) and the `-pre` option sets a prefix for the output file names. Note that
IQ-TREE generates three output files: a `.log`, `.treefile` and `.iqtree` file.
Most interesting stuff is in the `.iqtree` file. The tree (which can be loaded
in NJplot or FigTree for example) is in the `.treefile`. Note that likelihoods
are always shown as log-likelihoods to prevent numerical inaccuracies.

- Examine the stuff IQ-TREE prints to the screen, can you figure out the different steps IQ-TREE uses to find the ML tree?
- What is the log-likelihood associated with the ML tree?
- Do you think it is a problem that the likelihood is such a small number?
- How does the phylogeny compare with the distance based phylogeny under the JC
model? Can you explain why?

### Gamma distributed rates across sites

Now let's compare this with the JC+Gamma model.

```bash
iqtree -s 18SrRNA_45.txt -m JC+G  -pre JCG
```

- What is the MLE for the $\alpha$ parameter of the Gamma distribution?
- Does the topology differ between the ML tree found with JC and JC+G?

### The bootstrap

Since different methods gave different results for some clades, it would
definitely be worthwhile to try to get an idea of how well the different clades
in the tree are supported by the data. As you saw in the course, the most
commonly used approach to evaluate support is the bootstrap. We will use the
"ultrafast bootstrap" as implemented in IQ-TREE, which is not the same as
Felsenstein's original nonparametric bootstrapping approach.  Run IQ-TREE with
1000 ultrafast bootstrap replicates:

```bash
iqtree -s 18SrRNA_45.txt -m JC   -bb 1000 -pre JC_BSV
iqtree -s 18SrRNA_45.txt -m JC+G -bb 1000 -pre JCG_BSV
```

- Are there any nodes for which you find low support?
- Are the clades which were problematic in the distance-based analyses well-supported in the ML trees?

Now let us have a look at some other substitution models.

### The K2P model

```bash
iqtree -s msa.fasta -m K2P -pre K2P
```

- What is the difference between the K2P and JC model?
- How many more parameters does the K2P model have compared to the JC model?
- What are the ML estimate(s) of the parameter(s)
- **Extra:** Which model fits the data best according to the LRT (note: the critical value of the $\chi^2$ distribution with one degree of freedom at the 0.05 significance level is 3.845)

### The F81 model

```bash
iqtree -s 18SrRNA_45.txt -m F81 -pre F81
```

- What is the difference between the F81 model and the JC model?

## Model selection

IQ-TREE implements another program called ModelFinder (Kalyaanamoorthy *et
al.* 2017) that allows to quickly (and approximately) find the best fitting
substitution model. You can run it by just omitting the model specification
with the `-m` flag from your commands.

### The 18S rRNA set again

Run IQ-TREE with model selection:

```bash
iqtree -s 18SrRNA_45.txt
```

- Which model was selected?
- Do you get a different tree?

### Florida dentist scandal

We'll consider a different data set now. In the early 90s, a
dentist was accused of infecting several of his patients with HIV during
surgical procedures. After a "low-risk" patient was diagnosed with HIV, other
patients were screened, 10 of which had HIV.

In the file `hiv-dentist.fasta` you'll find sequences from the V3 region of the *env* gene of the HIV virus from the dentist, the patients and some local controls (AIDS patients from Florida that had no relationship to the dentist whatsoever). These sequences are unaligned, so you'll first need to align them. You can download software like MUSCLE, MAFFT or PRANK for that, or alternatively, you could run any of these tools online at EBI [`https://www.ebi.ac.uk/Tools/msa/`](https://www.ebi.ac.uk/Tools/msa/).

After you have obtained an alignment, you can visualize it with some tool like `aliview`, `seaview` or `alv`, or again, you can use some of the webservices at EBI like [`Mview`](https://www.ebi.ac.uk/Tools/msa/mview/).

Now infer a tree with IQ-TREE, try `-m JC`, `-m JC+G` and ModelFinder (by omitting the `-m` option). Remember that you can include bootstrapping by setting `-bb <number of bootstrap replicates you want>`.

- Based on the phylogeny, what can you conclude about the crime case?
- Would your conclusion change when using different substitution models?
- Why do you think they chose to sequence the *env* gene of HIV for the phylogenetic analyses?
- **Extra:** Have a look at the original paper (on minerva), how did they conduct the phylogenetic analysis? Do you obtain similar conclusions?
