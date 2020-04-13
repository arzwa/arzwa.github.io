# [back](/phylocourse/)
#
# \toc

# # Maximum likelihood based inference of phylogenies

# Maximum likelihood is a powerful approach in statistics in general, and phylogenetics in particular. Besides the fairly nice statistical properties of ML estimators, ML as a methodology has the advantage that it has widespread acceptance philosophically (the frequentist school of statistics praises the consistency etc. of ML estimators, while the Bayesians often regard ML a special case of Bayesian inference with some (hidden) prior assumptions). Another advantage is that ML estimates are invariant to transformations of the parameters, so we don't have to worry about what the 'natural' scale is for our parameters (e.g. a branch length can not be negative, so should we log transform it then?). In this section, we'll first explore a bit how ML works in the phylogenetics context, and do some exercises with a modern software for ML phylogeny inference.

# ## ML on trees: substitution models (again)

# [Recall](../submod) that the principle of maximum likelihood states that

# > For a particular model with parameters $\theta$, our best estimate for the parameters based on the data $D$ is given by the parameter values $\hat{\theta}$ for which the probability of observing the data $D$ is maximized.

# In mathematical notation $$\hat{\theta} = \mathrm{argmax}_{\theta}\ p(D|\theta)$$ which you can read as "the estimate $\hat{\theta}$ is found by maximizing $p(D|\theta)$ with respect to the parameters $\theta$. Incidentally $p(D|\theta)$ is called the **the likelihood**, as it indicates how probable it is to observe data for some value of $\theta$. When the likelihood is considered in this way (i.e. as a function of $\theta$ for *fixed* $D$), it is often writtin as the likelihood function $L(\theta;D)$ or $L(\theta|D)$ to stress the fact that $\theta$ is a variable but $D$ is not (it is fixed by our observation of it, i.e. the data).

# We [already saw a bit](../submod) of how the ML principle can be used to estimate the evolutionary distance between two sequences. The important take home message here is that

# > For a given substitution model, the likelihood of a multiple sequence alignment is fully specified by the tree topology $\Psi$ that describes the relationship between the sequences, the branch lengths $b$, and the parameters of the substitution model $\theta$.

# What this means is that the quantity $p(\Psi, b, \theta |D)$ is a well defined thing. And *if we were abe to compute this quantity* we could try to adopt the ML principle *to estimate* $\Psi, b$ and $\theta$.

# We know how to do this for the case our alignment consists of two sequences. In that case, there is no tree topology, and what we're left with is the problem of estimating the distances (and parameters of the substitution model), which we already did [here](../submod) and [here](../distance). So the next question is how to compute the likelihood of a larer alignment, where the sequences are related by a tree.

# ## Computing the likelihood on a tree

# Luckily, there exists an efficient algorithm to do just that: **Felsenstein's pruning algorithm**.[^pgm] TO DO.

# ## Searching tree space

# So we know how to compute the likelihood given a tree $\Psi$, branch lengths $b$ and parameters $\theta$. For numerical parameters (i.e. $b$ and $\theta$), we may use standard techniques from mathematical optimization to find the ML estimates. To find the tree topology that leads to the maximal likelihood value is however an untractable problem that requires a **heuristic** optimization. The basic approach is the following:

# 1. start with an initial tree topology $\Psi$
# 2. set $l_{best}$ to the maximum likelihood value for $\Psi$
# 3. make a random change to $\Psi$ to obtain $\Psi'$
# 4. compute the maximum likelihood value $l'$ for $\Psi'$
# 5. if $l' > l_{best}$, set $\Psi$ to $\Psi'$ and set $l_{best}$ to $l'$
# 6. repeat from (3) until there is no improveent in $l_{best}$

# This is however a computationally very costly procedure when executed naively, and that's where it is vital to have efficient software such as RaxML, IQ-TREE and PhyML! 

# # Exercises: ML tree inference with IQ-TREE
# For the exercises, we will use the computationally *very* efficient software IQ-TREE (Nguyen *et al.* 2015). This is a command line program, similar to the FastME or PHYLIP programs but not interactive. Please download IQ-TREE for your operating system at [http://www.iqtree.org/](http://www.iqtree.org/).

# Once you obtain the program, I'd recommend to put the executable in some dedicated directory[^dir] along with the data files. In this section, we'll use the 18S rRNA data sets for [20 taxa](/assets/phylocourse/data/18SrRNA_20.phy) and [45 taxa](/assets/phylocourse/data/18SrRNA_45.phy) again. But feel free to follow along with any data set you like. For a note on using command line programs in Windows, see [the footnote here](../distance/#fndef:commandline).

# > Try to answer the questions shown in boxes like this one below. For some of the questions on substitution models you may want to google a bit or have a look in the IQ-TREE manual or the course notes or something similar.

# ## Basic tree inference

# ### The Jukes-Cantor model

# Let us first reconsider the 45 species 18S rRNA data set. Remember that we had some problematic clades in our distance-based phylogenies, where we found different topologies depending on whether we used distances computed with the Gamma model of rate heterogeneity across sites or not and depending on the $\alpha$ parameter we chose. Let's see what we get when using ML.

# Put the `18SrRNA_45.phy` file in the directory with the IQ-TREE executable. Fire up IQ-TREE and run the following command:

# ```bash
# iqtree -s 18SrRNA_45.phy -m JC -pre JC
# ```

# The `-s` option is used for specifying the input multiple sequence alignment file, the `-m` option is used for setting the substitution model (here Jukes & Cantor) and the `-pre` option sets a prefix for the output file names. Note that IQ-TREE generates three output files: a `.log`, `.treefile` and `.iqtree` file. Most interesting stuff is in the `.iqtree` file. The tree (which can be loaded FigTree or [icytree](https://icytree.org/) for example) is in the `.treefile`. Note that likelihoods are always shown as *log*-likelihoods to prevent numerical inaccuracies.

# > - Examine the stuff IQ-TREE prints to the screen, can you figure out the different steps IQ-TREE uses to find the ML tree?
# > - What is the log-likelihood associated with the ML tree?
# > - Do you think it is a problem that the likelihood is such an absurdly small number?
# > - How does the phylogeny compare with the distance based phylogeny under the JC model? Can you explain why?

# ### Gamma distributed rates across sites

# Now let's compare this with the JC+Gamma model.

# ```bash
# iqtree -s 18SrRNA_45.phy -m JC+G  -pre JCG
# ```

# > - What is the MLE for the $\alpha$ parameter of the Gamma distribution?
# > - Does the topology differ between the ML tree found with JC and JC+G?

# ## The bootstrap

# Since different methods gave different results for some clades, it would definitely be worthwhile to try to get an idea of how well the different clades in the tree are supported by the data. As you saw in the course, the most commonly used approach to evaluate the statistical support of a phylogenetic tree in ML or distance based phylogenetics is the bootstrap. We will use the "ultrafast bootstrap" as implemented in IQ-TREE (which is not exactly the same as Felsenstein's original nonparametric bootstrapping approach[^ultrafast]).  Run IQ-TREE with 1000 ultrafast bootstrap replicates:

# ```bash
# iqtree -s 18SrRNA_45.phy -m JC   -bb 1000 -pre JC_BSV
# iqtree -s 18SrRNA_45.phy -m JC+G -bb 1000 -pre JCG_BSV
# ```

# > - Are there any nodes for which you find low support?
# > - Are the clades which were problematic in the distance-based analyses well-supported in the ML trees?

# Now let us have a look at some other substitution models.

# ## The K2P model

# ```bash
# iqtree -s 18SrRNA_45.phy -m K2P -pre K2P
# ```

# > - What is the difference between the K2P and JC model?
# > - How many more parameters does the K2P model have compared to the JC model?
# > - What are the ML estimate(s) of the(se) parameter(s)

# ## The F81 model

# ```bash
# iqtree -s 18SrRNA_45.phy -m F81 -pre F81
# ```

# >- What is the difference between the F81 model and the JC model?

# Now we've performed phylogeny inference under a bunch of different substitution models, but how can we know which odel we should use? This problem is for instance similar to a well-studied issue in typical statistics courses: we fitted a linear regression with a bunch of parameters, how can we know which parameters are meaningful to include? This is a problem of **model selection**, which will be discussed in the [next section](../modsel).

# ----------------------------------------------------------------------
# [^pgm]: For the machine learning/AI aficionados, this algorithm is essentially the same as what is reffered to by 'variable elimination' in the probabilistic graphical modeling literature.

# [^ultrafast]: The ultrafast bootstrap of Minh *et al.* (2013) is for instance not only much faster, but also less biased. It is well known that Felsenstein's nonparametric bootstrap tends to be conservative (i.e. bootstrap support values tend to be underestimated relatively to the true statistical support of a clade), the ultrafast bootstrap is less so, and can be more easily interpreted as they tend to more closely resemble probabilities.

# [^dir]: I sometimes notice that Windows users that are not particular computophiles are not familiar with the word *directory*, but generally employ the term *folder*. This seems to be related to some [philosophy that originated within Microsoft Windows](9https://en.wikipedia.org/wiki/Directory_%28computing%29#Folder_metaphor) where "There is a difference between a directory, which is a file system concept, and the graphical user interface metaphor that is used to represent it (a folder)."
