# [back](/phylocourse/)
#
# \toc

# # Maximum likelihood based inference of phylogenies

# Maximum likelihood is a powerful approach in statistics in general, and
# phylogenetics in particular. Besides the fairly nice statistical properties
# of ML estimators, ML as a methodology has the advantage that it has
# widespread acceptance philosophically (the frequentist school of statistics
# praises the consistency etc. of ML estimators, while the Bayesians often
# regard ML a special case of Bayesian inference with some (hidden) prior
# assumptions). Another advantage is that ML estimates are invariant to
# transformations of the parameters, so we don't have to worry about what the
# 'natural' scale is for our parameters (e.g. a branch length can not be
# negative, so should we log transform it then?). In this section, we'll first
# explore a bit how ML works in the phylogenetics context, and do some
# exercises with a modern software for ML phylogeny inference.

# ## From distances to ML on trees

# [Recall](../submod) that the principle of maximum likelihood states that

# > For a particular model with parameters $\theta$, our best estimate for the
# > parameters based on the data $D$ is given by the parameter values
# > $\hat{\theta}$ for which the probability of observing the data $D$ is
# > maximized.

# In mathematical notation 
#
# $$\hat{\theta} = \mathrm{argmax}_{\theta} p(D|\theta)$$ 
#
# which you can read as "the estimate $\hat{\theta}$ is found by maximizing
# $p(D|\theta)$ with respect to the parameters $\theta$", where $p$ denotes a
# probability [mass](https://en.wikipedia.org/wiki/Probability_mass_function)
# or [density](https://en.wikipedia.org/wiki/Probability_density_function)
# function.  Incidentally $p(D|\theta)$ is called **the likelihood**, as it
# indicates how probable it is to observe the data $D$ for some value of
# $\theta$. When the likelihood is considered in this way (i.e. as a function
# of $\theta$ for *fixed* $D$), it is often written as **the likelihood
# function** $\ell(\theta;D)$ or $\ell(\theta|D)$ to stress the fact that
# $\theta$ is a variable but $D$ is not (it is fixed by our observation of it,
# i.e. the data).

# We [already saw a bit](../submod) of how the ML principle can be used to
# estimate the evolutionary distance between two sequences. The basic principle
# can however be generalized to multiple sequences related by a tree topology.
# Instead of estimating pairwise distances and model parameters by maximum
# likelihood and then combining these distances in a tree using some *ad hoc*
# clustering approach, we can estimate *directly* the tree topology, branch
# lengths (distances) and model parameters. In other words given
#
# 1. The parameters of the substitution model $\theta$.
# 2. The tree topology $\Psi$ that describes the relationship between the sequences
# 3. The branch lengths $b$
#
# we can compute the likelihood $p(D|\Psi,b,\theta) = \ell(\Psi, b, \theta|D)$. 
# Somewhat more graphically:
#
# ![](/assets/phylocourse/img/l.png)

# **What this means is** that the quantity $p(D|\Psi,b,\theta)$ is a well
# defined thing. And *if we were indeed able to compute this quantity* (i.e. if
# there is an *algorithm* to compute it) we could try to adopt the ML principle
# *to estimate* $\Psi, b$ and $\theta$ by searching for those topologies
# $\Psi$, branch lengths $b$ and parameter values $\theta$ that make $p(D|\Psi,
# b, \theta)$ maximal for the given data $D$. Luckily, there exists an
# algorithm to do just that: **Felsenstein's pruning algorithm**.[^pgm] I won't
# delve into that further here, but it effectively boils down to a trick which
# enables us to sum the probabilities over all possible evolutionary histories
# (all possible ancestral sequences) in an efficient way.

# ## Searching tree space

# So we know how to compute the likelihood given a tree $\Psi$, branch lengths
# $b$ and parameters $\theta$. Now we still have to maximize it with respect to
# the topologies and parameters of interest. This is a very tough
# **optimization poblem**. For numerical parameters on a fixed topology (i.e.
# $b$ and $\theta$), we may use standard techniques from mathematical
# optimization to find the ML estimates. To find the tree topology that leads
# to the maximal likelihood value is however an untractable problem that
# requires a **heuristic** optimization. (Note that because of issues with
# [underflow](https://en.wikipedia.org/wiki/Arithmetic_underflow), one usually
# works with log-likelihood values instead of likelihood values, which is
# usually denoted by $L$). The basic approach for heuristic tree search is the
# following:

# 1. start with an initial tree topology $\Psi$
# 2. set $L_{best}$ to the maximum log-likelihood value for $\Psi$, i.e. $L_{best}
#    \leftarrow \underset{b,\theta}{\mathrm{max}}\ L(\Psi,b,\theta|D)$
# 3. make a random change to $\Psi$ to obtain $\Psi'$
# 4. compute the maximum log-likelihood value $L'$ for $\Psi'$
# 5. if $L' > L_{best}$, set $\Psi \leftarrow \Psi'$ and set $L_{best} \leftarrow L'$
# 6. repeat from (3) until there is no improvement in $L_{best}$

# This is however a computationally very costly procedure when executed naively
# (note for instance that doing this as written above would involve a separate
# optimization of the parameters $\theta$ and branch lengths $b$ for each
# topology we try out, so essentially a complete optimization problem in its
# own right at each iteration!), and that's where it is vital to have efficient
# software implementing clever heuristics such as RaxML, IQ-TREE and PhyML!

# # Exercises: ML tree inference with IQ-TREE
# For the exercises, we will use the computationally *very* efficient software
# IQ-TREE (Nguyen *et al.* 2015). This is a command line program, similar to
# the FastME or PHYLIP programs but not interactive. Please download IQ-TREE
# for your operating system at
# [http://www.iqtree.org/](http://www.iqtree.org/). To install on windows, you
# can follow the same [guidelines](/phylocourse/install-windows) as for FastME.

# Once you obtain the program, I'd recommend to put the executable in some
# dedicated directory[^dir] along with the data files. In this section, we'll
# use the 18S rRNA data sets for [20
# taxa](/assets/phylocourse/data/18SrRNA_20.phy) and [45
# taxa](/assets/phylocourse/data/18SrRNA_45.phy) again. But feel free to follow
# along with any data set you like. For a note on using command line programs
# in Windows, see [the footnote here](../distance/#fndef:commandline).

# > Try to answer the questions shown in boxes like this one below. For some of
# > the questions on substitution models you may want to google a bit or have a
# > look in the IQ-TREE manual or the course notes or something similar.

# ## Basic tree inference

# ### The Jukes-Cantor model

# Let us first reconsider the [45 species 18S rRNA data
# set](/assets/phylocourse/data/18SrRNA_45.phy). Remember that we had some
# problematic clades in our distance-based phylogenies, where we found
# different topologies depending on whether we used distances computed with the
# Gamma model of rate heterogeneity across sites or not and depending on the
# $\alpha$ parameter we chose. Let's see what we get when using ML.

# Put the `18SrRNA_45.phy` file in the directory with the IQ-TREE executable.
# Fire up IQ-TREE and run the following command:

# ```bash
# iqtree -s 18SrRNA_45.phy -m JC -pre JC
# ```

# The `-s` option is used for specifying the input multiple sequence alignment
# file, the `-m` option is used for setting the substitution model (here Jukes
# & Cantor) and the `-pre` option sets a prefix for the output file names. Note
# that IQ-TREE generates three output files: a `.log`, `.treefile` and
# `.iqtree` file. **Most interesting stuff is in the `.iqtree` file**. The tree
# (which can be loaded FigTree or [icytree](https://icytree.org/) for example)
# is in the `.treefile`. 

# > - Examine the stuff IQ-TREE prints to the screen, can you figure out (roughly) the different steps IQ-TREE uses to find the ML tree?
# > - Open the `.iqtree` file. What is the log-likelihood associated with the ML tree?
# > - Do you think it is a problem that the likelihood is such an inconceivably small number?
# > - How many parameters did we have to estimate?

# ### Gamma distributed rates across sites

# Now let's compare this with the JC+Gamma model. Note that we cannot use the
# complete Gamma model as in the distance based methods (for computational
# reasons), and that we use a **discrete approximation** to the Gamma distribution,
# by chopping up the Gamma distribution into $K$ pieces of equal total
# probability, which we can think of as $K$ rate classes. If for instance
# $K=3$, it means we cut the Gamma distribution in five pieces, each of which
# accounts for 33% of the total probability distribution, like in the plot
# below
#
using Distributions, Plots, StatsPlots # hide
α = 1.0
K = 3
d = Gamma(α, 1/α)  # the Gamma distribution object
q = quantile(d, 1/K:1/K:1)  # the discretization points (defining the classes)
m = quantile(d, 1/2K:1/K:1-1/2K)  # the medians in each class
@show round.(m, digits=2)
plot(d, grid=false, legend=false, size=(500,200), xlim=(0,5), xlabel="relative rate", ylabel="probability density", guidefont=8)
vline!(q, color=:black)
vline!(m, linestyle=:dot)

savefig("_assets/phylocourse/mliqtree/gamma.svg") # hide

# ![](/assets/phylocourse/mliqtree/gamma.svg)

# In this plot the black lines divide the Gamma distribution with $\alpha = 1$
# into $K=3$ pieces of equal probability, the green dotted lines mark the
# median value in each piece. This median value will be the representative of
# our three rate classes: so we have (in this example) a slow rate class with
# relative substitution rate 0.18, a medium rate class with relative rate 0.69
# and a fast rate class with relative rate 1.79.  So what we effectively do is
# approximate the continuous gamma model, where the relative substitution rate
# for each site is distributed according to a Gamma distribution, by a model
# where each site comes from one of $K$ rate classes, each with probability
# $1/K$, where the rate classes are determined by a discretization of the Gamma
# distribution. 

# ```bash
# iqtree -s 18SrRNA_45.phy -m JC+G  -pre JCG
# ```
#
# Again, look at the `.iqtree` output file

# > - What is the MLE for the $\alpha$ parameter of the Gamma distribution?
# > - How many rate classes does IQ-TREE use by default?
# > - What are the relative rates for each rate class for the estimated Gamma distribution?
# > - Is the implied Gamma distribution more asymmetric or less asymmetric then an Exponential distribution?
# > - Does the tree topology differ between the ML tree found with JC and JC+G?

# ## Some other common substitution models

# Now let us have a look at some other substitution models, consider for
# instance Kimura's two-parameter (K2P) model. Run

# ```bash
# iqtree -s 18SrRNA_45.phy -m K2P -pre K2P
# ```

# and take a look at the `.iqtree` file 
#
# > - What is the difference between the K2P and JC model?
# > - How many more parameters does the K2P model have compared to the JC model?
# > - What are the ML estimate(s) of the(se) parameter(s)
# > - **Extra:** Consider the F81, HKY and GTR models. Look again at the `.iqtree`
# >   output files. How do these models relate to the JC and K2P model? How do they
# >   relate to each other?

# Now we've performed phylogeny inference under a bunch of different
# substitution models, but how can we know which model we should use? This
# problem is for instance similar to a well-studied issue in typical statistics
# courses: we fitted a linear regression with a bunch of parameters, how can we
# know which parameters are meaningful to include? This is a problem of **model
# selection**, which is the subject of the [next section](../modsel).

# ## The bootstrap

# Since different methods gave different results for some clades, it would
# definitely be worthwhile to try to get an idea of how well the different
# clades in the tree are supported by the data. As you saw in the course, the
# most commonly used approach to evaluate the statistical support of a
# phylogenetic tree in ML or distance based phylogenetics is by
# 'bootstrapping'. We will use the "ultrafast bootstrap" as implemented in
# IQ-TREE (which is not exactly the same as Felsenstein's original
# nonparametric bootstrapping approach[^ultrafast]).  Run IQ-TREE with 1000
# ultrafast bootstrap replicates:

# ```bash
# iqtree -s 18SrRNA_45.phy -m JC+G -bb 1000 -wbtl -pre JCG_BSV
# ```

# Open the `.treefile` in FigTree or [icytree](https://icytree.org/), and
# display the bootstrap support values along the internal nodes (In FigTree 
# select the `Node labels` checkbox and find the right field to display, in
# IcyTree you need to check `Style > Internal Node Text > Label`).
#
# > - Are there any nodes for which you find low support?
# > - Are the clades which were problematic in the distance-based analyses well-supported in the ML trees?
# > - What does the `-wbtl` option do?
# > - Interpret how a bootstrap support value is obtained in terms of the `.ufboot` file
# > - Extra: Download the [DensiTree](https://www.cs.auckland.ac.nz/~remco/DensiTree/) program
# >   and open the ufboot file with it (Download the `DensiTree.jar` file, put it in some 
# >   directory, open a terminal/command prompt and go to the directory with the jar file, 
# >   then run `java -jar DensiTree.jar <path to your ufboot file here>`. What is displayed?
# >   Interpret what you see.


# ----------------------------------------------------------------------
# [^pgm]: For the machine learning/AI aficionados, this algorithm is essentially the same as what is referred to by ['variable elimination'](https://en.wikipedia.org/wiki/Variable_elimination) in the probabilistic graphical modeling literature.

# [^dir]: I sometimes notice that Windows users that are not particular computophiles are not familiar with the word *directory*, but generally employ the term *folder*. This seems to be related to some [philosophy that originated within Microsoft Windows](9https://en.wikipedia.org/wiki/Directory_%28computing%29#Folder_metaphor) where "There is a difference between a directory, which is a file system concept, and the graphical user interface metaphor that is used to represent it (a folder)."

# [^ultrafast]: The ultrafast bootstrap of Minh *et al.* (2013) is for instance not only much faster, but also less biased. It is well known that Felsenstein's nonparametric bootstrap tends to be conservative (i.e. bootstrap support values tend to be underestimated relatively to the true statistical support of a clade), the ultrafast bootstrap is less so, and can be more easily interpreted as they tend to more closely resemble probabilities.

