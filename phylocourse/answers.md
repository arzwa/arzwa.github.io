\toc

# Distance methods

## Distance matrix

>Question: does this distance matrix make sense to you? Can you spot clades of more closely related species already?

It does make sense at first sight. We see a couple of clear clusters. Clearly tetrapod
animals cluster strongly together with a low pairwise distance (Homo and Xenopus).
We see higher plants also forming a tight cluster (Zamia, Zea, Glycine and Oryza),
and the same holds for fungi and green algae. Palmaria (a red alga) seems to be a
bit of an outlier.

>Question: consider the portion of the matrix corresponding to the animals (human, Xenopus, Artemia and Anemonia), write down the associated part of the distance matrix. Confirm that the values are proper distances. Make a sketch of what you think the phylogeny will look like for these four species based on a quick look at the matrix.

|            | Human | *Xenopus* | *Artemia* | *Anemonia* |
| ------     | ----- | -------   | -------   | --------   |
| Human      | 0     | 0.049     | 0.20      | 0.22       |
| *Xenopus*  | 0.049 | 0         | 0.19      | 0.21       |
| *Artemia*  | 0.20  | 0.19      | 0         | 0.23       |
| *Anemonia* | 0.22  | 0.21      | 0.23      | 0          |

Clearly the values are distances, as the distance of $a$ to $b$ $d(a,b)$ is (1)
symmetric ($d(a,b) = d(b,a)$), (2) non-negative $d(a,b) \le 0$ (with equality
only when $a = b$) and (3) the distances satisfy the triangle inequality that
$d(a,c) \le d(a, c) + d(b, c)$. As a sketch of the phylogeny, I would group
Human and *Xenopus* first, for they are the closest to each other. Then I would
group the resulting clade with *Artemia*, since both human and *Xenopus* are
closer to *Artemia* than to *Anemonia*, and since *Artemia* is closer to both
Human and *Xenopus* then to *Anemonia*. Lastly I would group the resulting clade
with *Anemonia* as sister lineage.

## Clustering methods

>Question: how can you see in one glance that this is an ultrametric tree? 

The distance from tip to root is identical for all tips (leaves) of the tree.

Note that the evolutionary distance represents the estimated ‘amount’ of
evolution, i.e. the expected number of substitutions per site, so the feature
that all root-to-tip distances are equal amounts to the assumption that the
total ‘amount of evolution’ (i.e. expected number of substitutions) in the
evolutionary history of each taxon is the same. Since the actual time (in
millions of years say) to the root is the same for all tips of the tree, this
is equivalent to the assumption that the average substitution rate is equal for
all root-to-tip paths in the tree.

(Recall that an evolutionary distance is a product of the form evolutionary rate
$\times$ time. Keeping dimensions in the back of your head often helps: distance
[substitutions/site] = rate [substitutions/site/year] $\times$ time [year])

>Question: based on your knowledge of the tree of life, can you identify where
>the phylogeny is (very likely) wrong?

The most eye-catching violation of our understanding of evolution of life on earth
involves the placement of Anemonia (sea anemone, a cnidarian animal) as sister
to the Fungi, with the resulting clade in turn forming a monophyletic group with
plants that is sister to the other animals. Our current understanding is that animals
and Fungi are however both within the Opisthokonts, which is a sister lineage (in
the present taxon sampling) to the Archaeplastida (containing plants and algae).

>Question: Clustering methods like UPGMA and WPGMA produce rooted trees. Did
>the clustering algorithm identify the correct root?

No, it did not. The correct root split for this set of species would be the split of
Archaeplastida and Opisthokonta (so algae+plants vs. fungi+animals).

>Question: How do you interpret the branch lengths, what is the associated length 'unit'?

The branch lengths are in distance units which correpond the expected number of
substitutions per site. Distances reflect a statistical estimate of the
‘amount’ of molecular evolution (see also above).  In the trees constructed
using clustering methods, the branch lengths are derived from the pairwise
distance estimates using some rule (which differs for instnce for UPGMA and
WPGMA) to combine the $n \times m$ different pairwise estimates between two
clades of $n$ and $m$ taxa respectively into two branch lengths. In both UPGMA
and WPGMA the rule used to combine these pairwise estimates amounts to making
an assumption of constant rates of evolution. When interpreting the branch
lengths we should take all this into account. Using different clustering
methods will result in different internal branch lengths (but not the branch
lengths leading to tip nodes), while being derived from the same pairwise
distance matrix.

## Neighbor-joining

>Question: How can you see at one glance that this is not an ultrametric tree? What does this mean in terms of assumptions on the substitution rate?

A rooted non-ultrametric tree would be easy to spot because the distances from
root to tips would not be equal. However we have an unrooted tree (see fig. 1)
in this case. It is not possible to root this tree such that the root-to-tip
distances are identical across the tree however, so it cannot be ultrametric.
An easy way to see this is that two leaves that have a direct common ancestor
have different branch lengths (e.g. Zea and Oryza or Homo and Xenopus), which
can never be the case in an utrametric tree (see the trees for the clustering
methods if you are not convinced).

>Question: Neighbor-joining infers an unrooted tree, can you see how
>the fact that the tree is unrooted is represented in FigTree? Where
>should you root this tree (based on your knowledge of the tree of life)?
>Select the branch where you think the root should be (click on it) and
>hit the reroot button in figtree to root the tree.

An unrooted tree is usually represented by a polytomous root (i.e. a trifurcating root
node, see fig. 2). It is important to not interpret this as a true polytomy in the tree
(i.e. where an ancestral lineage would have split into tree distinct lineages at the
same time).

The current consensus view is that the Archaeplastida and Ophistokonta both form
monophyletic groups, so we should put the root in the branch connecting the last
common ancestor of the Ophistokonta (e.g. last common ancestor of *Homo* and
*Saccharomyces*) and the last common ancestor of Archaeplastida (e.g. last common
ancestor of *Palmaria* and *Zea*). See fig. 3. Note that when you root the tree in
FigTree, FigTree assumes you want to put the root in the middle of the branch
(i.e. the total distance between the two clades stemming from the root is divided
into two equal length daughter branches).

>Question: Does the NJ tree make more sense than the WPGMA tree? What do you think is causing this?

Clearly it does. Since the pairwise distances are the same for the NJ tree and
the WPGMA tree, the difference is due to the tree construction method not the
substitution model. Presumably the assumption of equal amounts of evolution
(constant rates of evolution) across the phylogeny is false. Some branches in the NJ
tree are much further from the root (or some other common ancestral node) than
others, indicating unequal rates of evolution across the tree. This is of course not
surprising, given that this tree describes the relationships between vastly different
organisms, with very different biologies that may affect key parameters of the
molecular evolutionary process.

## Gamma distances

>Question: interpret the Gamma distribution as a model for rate heterogeneity
>using the plot above. To what assumption on the rates across sites does a large
>value of $\alpha$ correspond? To what assumption does a small value of $\alpha$
>correspond?

The plot shows the probability density function for the distribution of
relative substitution rates across sites. As a model of rate heterogeneity
across sites, we can interpret the plot for instance as follows: the
probability that a particular site has a relative substitution rate (relative
to the other sites in the alignment that is) in the interval $[0.2, 0.8]$ (for
instance) is equal to the area under the curve between x-coordinates 0.2 and
0.8 (i.e. the definite integral $\int_{0.2}^{0.8} f(x)dx$ of the gamma density
function $f(x)$).

A large value of $\alpha$ corresponds to a distribution which is quite narrow
around one, so there is not much variability in rates across sites. A small
value of $\alpha$ indicates that there is a strongly asymmetric distribution of
rates across sites, with many sites having low relative rates, and many sites
having rather high relative rates.

>Question: What (if anything) is changing? What is the default $\alpha$ used in FastME? Experiment with the $\alpha$ parameter of the Gamma distribution by using for instance -g0.5 in the FastME command. What happens?

The branch lengths change considerably, both relative to each other as well as
in overall length. Assuming highly asymmetric distributions of the substitution
rate across sites (e.g. $\alpha$ = 0.2, fig. 4) leads to higher distance
estimates overall and more outspoken differences between different branches.
This is due to the combination of both allowing different substitution rates
across branches and across sites, leading to multiplicative effects (in
contrast with the clustering methods, where we saw not much change in the
relative branch lengths when using $\Gamma$ heterogeneity). Using more
symmetric distributions (e.g. $\alpha$ = 10 fig. 5) leads to a result that is
indistinguishable from the result without $\Gamma$ distances. This is because
the different rates across sites are all very close to the mean rate, so the
model is well-approximated by a model with a single substitution rate for all
sites (i.e. the normal JC model).

# Maximum likelihood phylogeny inference

## Analysis under the Jukes & Cantor model

> Examine the stuff IQ-TREE prints to the screen, can you figure out
> (roughly) the different steps IQ-TREE uses to find the ML tree?

Here I’ll highlight some parts of the output. A first informative bit is

```md
Alignment has 45 sequences with 1626 columns, 1110 distinct patterns
934 parsimony-informative, 232 singleton sites, 460 constant sites
```

This shows some features of the data. We have 1626 columns, i.e. data points
(recall that we assume independent evolution across sites, so each column of an
alignment is a data point from a statistical point of view). There are 1110
unique columns (or ‘site patterns’ as they tend to be called in phylogenetics).
Note that the likelihood will be the same for identical columns in the
alignment, so the program will not compute the log-likelihood for 1626 sites
separately, but run the basic likelihood algorithm for the 1110 unique patterns
everytime the likelihood (i.e.  $P(\text{data}|\text{tree,parameters})$) is
calculated, saving time. 934 sites are ‘parsimony informative’ which means that
for those sites, some topologies would have a higher parsimony score than
others.

Next IQ-TREE performs a compositional heterogeneity test. This is a crude
statistical test to assess whether nucleotide composition is similar or not
across sequences.  For a ‘tree of life’ scale data set like the one here,
you’ll see the sequences have quite different base compositions. More info can
be found [here](http://www.iqtree.org/doc/Frequently-Asked-Questions#what-is-the-purpose-of-composition-test).

Then the IQ-TREE algorithm starts.
- First a maximum parsimony tree is inferred to use as an initial tree for the
  ML search.
- Next, initial model parameters are estimated by doing ML estimation on the
  initial parsimony tree.
- Then pairwise distances are estimated using the model with the initial
  parameter estimates.  
- Then a Neighbor-joining tree is inferred based on the distance matrix (using
  the BIONJ algorithm), and the loglikelihood is computed for the resulting
  (BIO)NJ tree. This tree is our next ‘initial’ tree.  
  
IQ-TREE will then generate a candidate tree set based on the initial tree,
compute 1the likelihoods, and then proceed with the best scoring trees.  For
these trees it will then try to perform random topological rearrangements etc.
to search for the best tree (the one with the highest likelihood).  After that
a smaller candidate set is obtained, and IQ-TREE will do a more thorough
optimization of the model parameters, branch lengths and topology to arrive at
a final tree estimate. For the final tree, it will then do a final maximum
likelihood estimation step of the model parameters.

> Open the `.iqtree` file. What is the log-likelihood associated with the
> ML tree?

Check the line with: 

```md
BEST SCORE FOUND : -31599.974
```

> Do you think it is a problem that the likelihood is such an
> inconceivably small number?

No, this is perfectly expected. Recall that the likelihood L is defined as the
probability of the data given the tree topology, model parameters and branch
lengths, i.e. 

$$\ell(X|\Psi,b,\theta) = P(X|\Psi,b,\theta)$$

where $X$ denotes the data, $\Psi$ the topology, $b$ the branch lengths (i.e.
$\Psi$ and $b$ together define the *phylogenetic tree*) and $\theta$ denotes
the parameters of the substitution model.  The likelihood of $e^{−31600}$ is
thus the probability that if you randomly simulate the evolution of a 1626 site
long-sequence from the Jukes-Cantor with the inferred ML parameters along the
inferred ML tree, you will obtain the exact data set you have used as input. Of
course this probability really is ridiculously small. (Back of the envelope:
note that the log-likelihood per site is about -19, which you can compare to
the log-likelihood of observing a random 45 character string (one nucleotide
for each species) for a four-letter alphabet under a uniform distribution,
which is $\log(1/4^{45}) \approx −62$).

> How many parameters did we have to estimate?

Check the line:

```md
Number of free parameters (#branches + #model parameters): 87
```

## Gamma distributed rates

> - What is the MLE for the $\alpha$ parameter of the Gamma distribution?
> - How many rate classes does IQ-TREE use by default? 
> - What are the relative rates for each rate class for the estimated Gamma
>   distribution? 
> - Is the implied Gamma distribution more asymmetric or less asymmetric then an
>   Exponential distribution?

```md
Model of rate heterogeneity: Gamma with 4 categories
Gamma shape alpha: 0.46

 Category  Relative_rate  Proportion
  1         0.02623        0.25
  2         0.225          0.25
  3         0.7888         0.25
  4         2.96           0.25
Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category.
```

The ML estimate is $\alpha = 0.46$. Note that we approximate the Gamma
distribution of rates across sites by fitting a discretized Gamma model with 4
classes. The relative rates associated with the different categories are shown
in the .iqtree output. (see the plots and explanation in the exercises for more
explanation on the discretization). Clearly, IQ-TREE uses four rate classes by 
default. The relative rates are listed in the table. The distribution is more
asymmetric than an Exponential distribution ($\alpha < 1$, see plots on the
exercises page for distance methods).

## Other common substitution models

>What is the difference between the K2P and JC model? 

The K2P model allows for unequal substitution rates for transitions (pyrimidine
to pyrimidine, purine to purine) vs.  tranversions (pyrimidine to purine and
vice versa). 

> How many more parameters does the K2P model have compared to the JC model? 

The K2P has one additional parameter, the ratio of the substitution rate for
transitions to the substitution rate for transversions. This is confirmed by
the IQ-TREE output: `Number of free parameters (#branches + #model parameters): 88`

> What are the ML estimate(s) of the(se) parameter(s)?

We find in the `.iqtree` file:

```md
Rate parameter R:

  A-C: 1.0000
  A-G: 2.3152
  A-T: 1.0000
  C-G: 1.0000
  C-T: 2.3152
  G-T: 1.0000
```

IQ-TREE sets the rate of transversions to 1 and estimates the ratio of
transition rate to the transversion rate (this parameter is often denoted
$\kappa$). Here we estimated $\hat{\kappa} = 2.3152$, so given that a
substitution occurs, it is 2.3152 as likely that it will be a transition compared to a transversion.

>Extra: Consider the F81, HKY and GTR models. Look again at the .iqtree output
>files. How do these models relate to the JC and K2P model? How do they relate
>to each other? 

Briefly: the F81 model allows for unequal equilibrium frequencies (K2P and JC
assume the equilibrium frequencies -- i.e. the base frequencies when the
evolutionary process assumed by the model is run for a long time -- to be all
0.25) but assumes equal substitution rates for all substitution types. The HKY
model allows for both unequal base frequencies and a different rate for
transitions compared to transversions, so it can be seen as a combination of
F81 and K2P (if you set the equilibrium frequencies in HKY to 0.25 you get K2P,
if you set $\kappa$ to one, you get F81). The GTR model allows for unequal
equilibrium base frequences and allows all types of substitions to occur at
different rates (but assumes symmetry, i.e. the rate of `A -> T` is the same as
`T -> A`). It has all other substitution models we discussed as a special case.


## The bootstrap

> What does the `-wbtl` option do? 

This option writes the inferred trees for the bootstrap replicates to a file
with extension `.ufboot`.

> Interpret how a bootstrap support value is obtained in terms of the `.ufboot`
> file 

In the `.ufboot` file each line is a tree topology obtained for a bootstrapped
alignment (see the course on what that means). A bootstrap support value for a
particular clade will be the proportion of bootstrap trees (which are written
in that file) which contain that clade. 

Note that I say clade and not subtree.  The support value should be interpreted
as the support for a bipartition (or *split*), not a subtree. To make this
clear consider the tree for the 18SrRNA data set with 20 species:

![](/assets/phylocourse/mliqtree/20species-boot.svg)

Take for instance the bootstrap support value of 85 for the clade of green
plants. What this indicates is that in the 1000 bootstrap replicates we found
that 85% of the *unrooted* trees contain a clade which consists of the green
plants, or in other words, 85% of the trees contain a bipartition (i.e. we can
cut the unrooted tree in two trees) such that one of the resulting subtrees
consists of the green plants and the other of all the rest. Importantly, this
does not mean that 85% of the bootstrap replicates contain the exact *subtree*
we see in this figure, the branching patterns within the subtrees defined by
the bipartition may differ.
