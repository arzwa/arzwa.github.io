<!--This file was generated, do not modify it.-->
[back](/phylocourse/)
\toc

# Model selection for maximum likelihood based phylogenetic inference

Whereas maximum likelihood (ML) allows us to infer a phylogenetic tree and the parameters of the substitution model, we still need to decide ourselves, often rather arbitrarily, which substitution model we are using. This is quite problematic, because there are a very large number of possible models one could choose from. Consider for instance the following graph representation showing some of the more commonly used models and their relationships:

![](/assets/phylocourse/img/models.png)

Of course, this problem is only natural, after all, when we -- the scientists -- are modeling something we have to decide on our modeling assumptions!  Nevertheless, given a candidate set of substitution models, one could reason that we simply extend the ML principle, and choose the model which gives the highest likelihood.

This is however problematic, since more complicated models will always result in a higher likelihood, so the ML principle when applied to models will always result in the selection of the biggest, most complicated and most flexible model. Again, this is a general issue in statistics, known as **overfitting**. The key issue is that we actually want to know whether some complex model fits the data *significantly* better than a more simpler model. As an example of overfitting consider the following hypothetical regression example:

![](/assets/phylocourse/img/overfitting.png)

The first model fits a flat line $y = a$, i.e. a constant (or intercept if you will), to the data and involves a single parameter. This model is clearly underfitting as it does not capture any trend in the data. The second model is a linear model $y = a + bx$ and involves two parameters ($a$ and $b$). The third model fits a quadratic function and involves three parameters $y = a + bx +cx^2$. The second and third model both look fairly reasonable for the amount of data we have. The fourth model is an exceedingly complex fit (it isn't even a proper function of $x$!) And is clearly overfitting the data. Now we wish to find a quantitative criterion for assessing this fitting trade-off and preferably one that involves the likelihood, because that is something we can compute in a phylogenetics context.

Luckily, since this is a common (but challenging nevertheless) problem in statistics, there exist devices that can be used to select models of varying complexity[^mlbayes]. The three main methods that are used for model selection in a ML setting are:

1. The likelihood ratio test (LRT)
2. The Akaike information criterion (AIC)
3. The Bayesian information criterion (BIC)

These all involve quantities that can be used to do pairwise comparisons of models, and they can be computed from the likelihood values under the different models that you wish to compare. For more insight in the statistical details of these methods for model selection, please refer to wikipedia or your favourite statistical modeling textbook.

## Model selection for substitution models: nested models

We first consider the special case of **nested models**. Two models are nested if one of the two is a special case of the other (i.e. models on different rows in Figure 2). For example the Jukes & Cantor (JC) model is nested in the K2P model since we obtain the JC model if we set the transition/transversion ratio $\kappa$ in the K2P model to 1 (i.e. we assume transitions and transversions are equally likely). When we do model selection between the JC and K2P model, what we are in effect asking is

>Does a model with an extra parameter $\kappa$ for the transition /  transversion ratio does a significantly better[^allmodels] job than a model without this parameter at explaining the data?

For models that are nested in this way, we can use a **hypothesis testing approach**. Secifically we can compute the following test statistic

$$ D = 2\mathrm{log}\Bigg( \frac{l_{1}}{l_{0}}\Bigg) = 2\big(\mathrm{log}(l_1) - \mathrm{log}(l_0)\big) \sim \chi^2_{k_1 - k_0}$$

where

- $l_1$ is the likelihood under model 1 with $k_1$ parameters
- $l_0$ is the likelihood under model 0 with $k_0$ parameters
- $k_1 > k_0$

Which involves a likelihood ratio $l_1/l_0$ annd hence is called the **likelihood ratio test (LRT)** statistic.

Specifically $l_1$ is the likelihood for the more parameter-rich model (e.g. K2P) and $l_0$ is the likelihood for the other model (e.g. JC). The likelihood ratio measures how much the more complex model fits the data better than the simpler model. If the test statistic is bigger than the critical value (at the desired significance level) of the corresponding $\chi^2$ distribution, then this indicates that the more complex model provides a *statistically significant better fit* to the sequence data. If not, than we may conclude that the difference in the likelihoods is not statistically significant, indicating that the simpler model does about an equally good job at explaining the data.[^nhst]

### Example

I randomly simulated a sequence alignment from the JC model. Using the IQ-Tree program (see [the exercises on ML inference](../mliqtree)) I get a maximum log-likelihood value of -2740.380 for the JC model. For the K2P model, the likelihood value obtained is -2738.879.

```julia:ex1
lrt(ℓ_1, ℓ_0) = 2(ℓ_1 - ℓ_0)
lrt(-2738.879,-2740.380)
```

The LRT is about 3.0. The critical value at the 0.05 level for the $\chi^2$ distribution with one degree of freedom is

```julia:ex2
using Distributions
quantile(Chisq(1), 0.95)
```

So the LRT is below the critical value. In other words, the LRT suggests that the K2P model does not provide a significantly better fit to this random data set. Note that at the 0.1 significance level (where the critical value of the $\chi^2_1$ distribution is about 2.7), the LRT would suggest the K2P model to provide a better fit, even though the data was simulated according to JC. This should convince you that care should be taken not to take these statistical tests too serious.

>Consider your ML inferences for the 18SrRNA data set from the previous exercise set. Did the K2P model provide a significantly better fit to the sequence data when using the LRT?

## Model selection for substitution models: the general case

When models are not nested, other tools (derived from information theory) can be used to perform model selection, albeit not in a null hypothesis significance testing kind of approach. These involve **information criteria** which are designed to indicate the *relative* gain/loss of statistical fit when deciding between different models. We won't go in any depth here, but we note down the formulae for two such information criteria:

1. The **Akaike information criterion (AIC)** for a given model is given by $$AIC = 2k - 2l_{max}$$ where $k$ is the number of estimated parameters in the model and $l_{max}$ is the maximum log-likelihood value obtained for the model. **The lower the AIC value, the lower the estimated loss of information when representing the true data generating process by the model** (which is obviously a caricature of reality). In other words, the ower the AIC value, the better the fit. Given the AIC values of two models $AIC_0$ and $AIC_1$, where $AIC_0 > AIC_1$, we can compute the probability that model 0 provides at least as good a fit as model 1 can be estimated as $e^{(AIC_1 - AIC_0)/2}$. We see that the AIC involves a term of $2k$ penalizing the likelihood value for the number of parameters in the model.

2. The **Bayesian information criterion (BIC)** (which has not much to do with Bayesian inference or Bayesian statistics in general) is very similar and can be computed as $$BIC = \log(n)k -2l_{max}$$ where $n$ is the number of data points). It cannot exactly be interpreted in the same way, but still we are interested in the model that minimizes the BIC value.

> Compute the likelihood values for the JC, K2P, F81, HKY85 and GTR models for the 18S rRNA data using a fixed tree topology (use `iqtree -t <treefile> -s <sequence-file> -n 1 -m <model>`, which will only optimize the parameters and branch lengths for the initial topology in `treefile`)[^treefix]. Compute the AIC values for each model and rank them. Which model provides the best fit according to the AIC?

> IQ-TREE provides automated model selection. Just run `iqtree -s <sequence-file>` without the `-m` flag and it will automatically search for the best fitting model. Try to locate the model selection steps in the IQ-TREE output for a data set of choice, can you understand how IQ-TREE heuristically tries to find the best fitting model, **before doing the actual full maximum likelihood optimization**? Think a moment about this, does this make sense to you? What could be the pitfalls of such an approach?

## Partition models

Model selection does not have to be restricted to the substitution model. We may wish to use different substitution models for different sites of the alignment (for instance the codon positions), or we may wish to analyze a concatenated alignment of multiple genes, but might want to have a different substitution model for each gene. These are called partitions of the data. Partitioning the data in such ways also increases the complexity of the statistical model and may lead again to overfitting, so here we may also wish to perform model selection. We can do this using the same techniques as outlined above.

--------------------
[^mlbayes]: The notes in this page will be exclusively concerned with maximum likelihood based statistical inference. The model selection problem can be tackled in a very different (and in my opinion more natural and elegant) approach when one adopts a Bayesian framework for statistical inference.

[^allmodels]: Please do note that *better* does [not necessarily mean *good*](https://en.wikipedia.org/wiki/All_models_are_wrong).

[^nhst]: It is advised not to take the resulting $p$-values and hypothesis tests too seriously, as they hinge on a whole bunch of assumptions and the reasonableness of the assumed models that we compare. These hypothesis tests should merely be used as a guide in choosing a model.

[^treefix]: There is also the `--tree-fix` option to do this, but it didn't work in the version of IQ-TREE I was using...

