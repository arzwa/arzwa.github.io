

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
