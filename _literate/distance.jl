# [back](/phylocourse/)
# \toc
# # Distance-based phylogenetic inference

# >**Questions** and **Exercises** will be marked in blocks like this, try to
# >formulate a brief answer to them.

# >**Note** that there are some blocks with `julia` code included below. This
# >is just for illustrating things in a way such that those who are interested
# >can see how these illustrations (numerical or visual) were actually generated. 
# >Don't let them scare you.

# Distance-based phylogenetic methods proceed by computing a **pairwise
# distance matrix** for an input alignment based on a [substitution
# model](../submod), i.e. an assumed **model of evolution** for molecular
# sequences. Usually these estimated distances are the maximum likelihood
# estimates (MLE) of the pairwise distance under the assumed substitution
# model. For instance, consider two sequences, for an observed proportion of
# different sites $x$, the maximum likelihood estimate of the distance between
# the two sequences under the Jukes & Cantor substitution model is
#
# $$\hat{d} = -\frac{3}{4} \log \Big(1 - \frac{4}{3}x\Big)$$
# 
# So for instance for the alignment

seqa = "ATCGGGCTAGC"
seqb = "TTCGGCTTACC";

# The proportion of different sites is
differences = [seqa[i] != seqb[i] for i=1:length(seqa)]
@show differences
num_differences = sum(differences)
x = num_differences/length(seqa)

# The estimated distance is
d = -0.75 * log(1. - 4x/3)

# This distance is correctly interpreted as **an estimate of the expected
# number of substitutions per site** separating the two sequences. Note that
# the distance is a product of two other evolutionary variables: the
# **substitution rate** and the **divergence time**. If we have $n$ sequences,
# we can apply the distance formula to obtain an $n \times n$  **distance
# matrix**. Note however that such simple distance formulae are generally not
# available for more complicated substitution models. For more details, see the
# [notes on substitution models](../submod).

# > **Question (optional)**. The proportion of different sites in a pair of sequence $x$
# > is a **summary statistic** of the data, we reduce the data to a single number
# > $x$. For the Jukes - Cantor model the likelihood $P(\text{sequence pair}|\mathrm{model})$ 
# > only depends on the summary statistic $x$. Why will
# > such a summary statistic not work for more complicated substitution models,
# > for instance where transitions and transversions have different rates?
# > If we cannot use a formula based on some summary statistic, how can we
# > estimate a distance?

# After computing the distance matrix, a tree is inferred by using some sort of
# clustering algorithm or least-squares estimation step. Distance-based
# phylogenetic inference is thus essentially a two-step process:

# 1. Compute pairwise distances
# 2. Infer the tree, assuming the distances

# This two-step procedure is both the strength and weakness of distance-based
# methods. By collapsing the sequence data with $n$ sequences into a single $n
# \times n$ matrix we dramatically reduce the data, making tree inference
# computationally very fast. However, by reducing the rich sequence data to
# such a matrix we risk throwing away a lot of potentially interesting
# evolutionary information. Indeed, if we have $n$ sequences, there are $4^n$
# nucleotide patterns for any site when we compare all $n$ sequences, which we
# reduce to $n(n-1)/2$ when using methods based on pairwise distances. This is
# different from full maximum-likelihood tree inference methods, where all
# possible higher-order (higher than pairwise that is) character patterns
# contribute to phylogeny inference. One can think of the approach of phylogeny
# estimation from pairwise distances as follows: each distance is an estimate
# of the branch length between two sequences in the phylogeny, and we seek the
# phylogeny that is most compatible with these independently obtained branch
# length estimates. In other words, we estimate $n(n-1)/2$ two-species
# phylogenies and seek to somehow infer the phylogeny for all $n$ sequences
# from those two-species phylogenies.

# Another issue is that when using distance-based methods the estimated
# distances from step 1 *are treated as observed data* in step 2. However,
# distances are themselves estimates, associated with some uncertainty (a
# distance estimate has a variance for instance), and this uncertainty in the
# estimated distances is neglected in step 2. This is a common source of issues
# in many areas of bioinformatics, i.e. treating estimates as data, as if known
# without error (the multiple sequence alignments we use to infer phylogenies
# provide a good example!). 

# Both of these issues are solved in ML and Bayesian phylogenetic inference,
# however at the price of a strongly increased computational cost. Because of
# their speed and relatively good performance, distance methods are still quite
# often used, even in this age of highly performant computers. Some popular ML
# tree inference programs (for instance IQ-TREE, see the ML tutorial) start
# their tree search algorithm from a distance-based tree for instance. Also
# popular packages for phylogenomic and comparative genomic inference such as
# [OrthoFinder](https://github.com/davidemms/OrthoFinder) use distance-based
# phylogenetic methods in their workflow.

# # Software and data

# 1. Download the **FastME** software at
#    [http://www.atgc-montpellier.fr/fastme](http://www.atgc-montpellier.fr/fastme/binaries.php)
#    (on the bottom of the page, click the download button). In the downloaded
#    folder you will find a binaries directory, identify the binary (executable)
#    for your operating system and put it in some folder of your
#    convenience.[^fastmeonline] For some additional guidance with installing
#    FastME on Windows, have a look [here](/phylocourse/install-windows).

# 2. To view trees, I recommend the
#    [FigTree](https://github.com/rambaut/figtree/releases/tag/v1.4.4) tool.
#    FigTree requires Java, but that should be available on most machines.
#    Download the executable for your operating system from the link above (`.zip`
#    file for Windows users, `.tgz` for *nix users, I guess the `.dmg` file is
#    something for Mac OSX users(?))

# 3. We will use two small 'tree of life' data sets. One [18SrRNA data set with
#    20 taxa](/assets/phylocourse/data/18SrRNA_20.phy) and another [18SrRNA data
#    set with 45 taxa](/assets/phylocourse/data/18SrRNA_45.phy). These are in
#    PHYLIP format[^ancientdata], fasta formatted files can be obtained as well
#    [here](/assets/phylocourse/data/18SrRNA_20.fasta) and
#    [here](/assets/phylocourse/data/18SrRNA_45.fasta).
#
# It is always best to have a high-level look at the data before delving into analyses.
# You can view the alignments with a tool like [alv](https://github.com/arvestad/alv) 
# (personal favorite) or an online tool like [MView](https://www.ebi.ac.uk/Tools/msa/mview/).
# Here's a screenshot of the first 80 columns of the 20 taxa data set from MView.
#
# ![](/assets/phylocourse/distance/18SrRNA-aln.png)
#

# # Distance-based phyogenetics with FastME
#
# ## Computing distance matrices using FastME

# While you could easily implement a little program to generate a distance
# matrix under the Jukes and Cantor model (see the formulae for the distance
# in the section on [substitution models](../submod)), this is less
# straightforward when employing more complicated substitution models. FastME
# is probably the fastest implementation of distance-based phylogenetics
# methods using general substitution models available today. It is a software
# tool that can both be used to estimate distance matrices and infer trees using
# **Neighbor-Joining**, **least-squares**, **minimum evolution** and related
# distance matrix based methods.

# There are two ways to run FastME. There is an interactive mode (inherited
# from the influential [PHYLIP software
# package](http://evolution.genetics.washington.edu/phylip.html)) and a command
# line mode. Here I will use the command line mode, but feel free to use the
# interactive mode.  First get a look at the help message
# ```
# fastme -help
# ```

# If you have read the course notes from prof. Van de Peer, most of the options
# listed there (but not all, no worries) should make some sense to you. To
# infer a distance matrix with the Jukes & Cantor model you can run something
# like this from the command line[^commandline] :

# ```
# fastme -i 18SrRNA_20.phy -O 18SrRNA_20_JCmatrix.txt -dJC69
# ```

# Try to run this command and see whether you manage to make it work. The above
# command should generate a distance matrix. Take a look at the file
# `18SrRNA_20.phy_fastme_stat.txt`, is everything as expected?

# Let's make a heatmap to visualize the distance matrix:

using Plots
function readmatrix(file)  # little function to read the FastME distance matrix
    lines = [split(l) for l in readlines(file)[2:end] if l != ""]
    matrix = hcat([map(x->parse(Float64, x), l[2:end]) for l in lines]...)
    names = [l[1] for l in lines]
    matrix, names, length(names)
end

matrix, taxa, ntaxa = readmatrix("_assets/phylocourse/distance/18SrRNA_20_JCmatrix.txt")
heatmap(matrix, 
        yticks=(1:ntaxa, taxa), 
        xticks=(1:ntaxa, taxa), 
        xrotation=45, size=(700,650))
savefig("_assets/phylocourse/distance/hm1.svg") # hide

# ![](/assets/phylocourse/distance/hm1.svg)

# Something similar can be done in any other programming language. For instance in `R` 
# the following should work:
# ```R
# X = read.table("./18SrRNA_20_JCmatrix.txt", skip=1)
# taxa = c("taxon", as.character(X$V1))  # set column names
# names(X) = taxa
# X = as.matrix(X[,2:ncol(X)])  # remove the row names, and convert to a numeric matrix
# heatmap(X)
# ```

# >**Question**: does this distance matrix make sense to you? Can you spot
# > clades of more closely related species already?
#
# > **Question**: consider the portion of the matrix corresponding to the
# > animals  (human, *Xenopus*, *Artemia* and *Anemonia*), write down the 
# > associated part of the distance matrix. Confirm that the values are 
# > [proper distances](https://en.wikipedia.org/wiki/Metric_(mathematics)).
# > make a sketch of what you think the phylogeny will look like for these
# > four species based on a quick look at the matrix.

# ## Clustering methods

# **Clustering methods** in distance-based phylogenetics are not different from
# hierarchical clustering methods used in other applications, such as
# unsupervised machine learning for instance. Both the UPGMA and WPGMA
# clustering methods are specific cases of *average linkage clustering* (with
# different formula's used to compute the distance between two already existing
# clusters). Clustering methods can be found in many packages for scientific
# computing for most programming languages (here I'll use the `Clustering.jl`
# package for the julia progamming language, note that average linkage
# clustering as it is usually implemented is identical to what is called the
# WPGMA method in phylogenetics (see p.16 & 17 in prof. Van de Peer's course
# notes).

using Clustering, StatsPlots
hcl = hclust(matrix, linkage=:average)
plot(
    plot(hcl, xticks=false),
    heatmap(matrix[hcl.order,hcl.order], colorbar=false,
        yticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xrotation=45),
    layout=grid(2, 1, heights=[0.2,0.8]), size=(600,750))

savefig("_assets/phylocourse/distance/wpgma.svg") # hide

# ![](/assets/phylocourse/distance/wpgma.svg)

# **Note**, in `R` you can use the following to perform WPGMA clustering (see
# note above on how to load the distance matrix in R):
#
# ```R
# D = dist(X)  # convert to a 'distance matrix' object for hclust
# hc = hclust(D, method='average')  # average linkage clustering = UPGMA by default in R
# tree = as.dendrogram(hc)
# plot(tree)
# ```

# >**Question**: how can you see in one glance that this is an ultrametric tree?

# >**Question**: based on your knowledge of the tree of life, can you identify
# > where the phylogeny is (very likely) wrong?

# >**Question**: Clustering methods like UPGMA and WPGMA produce rooted trees.
# >Did the clustering algorithm identify the correct root?

# >**Question**: How do you interpret the branch lengths, what is the
# >associated length 'unit'?

# ## Neighbor-joining

# Neighbor-joining (NJ) is another method for distance-based phylogenetics. It
# is also a clustering method, but one that does not produce ultrametric trees.
# NJ and its variants are implemented in FastME. To run tree inference with NJ
# for an input alignment or distance matrix, run `fastme` with the `-o
# <output_file>` and `-m NJ` options. For instance

# ```
# fastme -i 18SrRNA_20.phy -o 18SrRNA_20_JC.nwk -dJC69 -m NJ
# ```

# Open the NJ tree in the `.nwk` file generated by FastME in FigTree (should be
# straightforward).

# >**Question**: How can you see at one glance that this is not an ultrametric
# >tree? What does this mean in terms of assumptions on the substitution rate?

# >**Question**: Neighbor-joining infers an unrooted tree, can you see how the
# >fact that the tree is unrooted is represented in FigTree? Where should you
# >root this tree (based on your knowledge of the tree of life)? Select the
# >branch where you think the root should be (click on it) and hit the `reroot`
# >button in figtree to root the tree.

# >**Question**: Does the NJ tree make more sense than the WPGMA tree? What do
# >you think is causing this?

# ## Rate heterogeneity across sites: 'Gamma distances'

# The assumption that all sites in the alignment evolve at the same rate is
# often problematic (just think of the active site of an enzyme vs. some region
# with a purely structural function for instance). Accounting for heterogenity
# in substitution rates across sites can be done by assuming some distribution
# over *relative rates* across the alignment. In distance estimation, people
# often use a Gamma distribution to model rate differences across sites. 
#
# The Gamma distribution is fixed to a mean value of 1, so that the average
# relative rate across the alignment is 1 (this is of course exactly what we
# mean by modeling *relative* rates). After fixing the mean to 1, there is
# a single free parameter left to specify the shape of the Gamma distribution,
# denoted $\alpha$. Here is a plot of the Gamma distribution for different
# values of $\alpha$

using Distributions
p = plot(title="The Gamma distribution with mean 1")
for α in [0.1, 0.25, 0.5, 1.0, 5.0, 10., 100.]
    plot!(p, Gamma(α, 1/α), label="\\alpha = $α", xlim=(0,5), ylim=(0,5))
end
savefig(p, "_assets/phylocourse/distance/gamma.svg") # hide

# ![](/assets/phylocourse/distance/gamma.svg)
#
# >**Question:** interpret the Gamma distribution as a model for rate
# > heterogeneity using the plot above. To what assumption on the rates
# > across sites does a large value of $\alpha$ correspond? To what
# > assumption does a small value of $\alpha$ correspond?

# We can infer a tree using $\Gamma$ distances nd the JC model using FastME:

# ```
# fastme -i 18SrRNA_20.phy -o 18SrRNA_20_JC.nwk -dJC69 -g -m NJ
# ```

# >**Question**: What (if anything) is changing? What is the default $\alpha$ 
# >used in FastME? Experiment with the $\alpha$
# >parameter of the Gamma distribution by using for instance `-g0.5` in the
# >FastME command. What happens?

#
# >**Question (extra)**. Why, in fact, assume a Gamma distribution for rates across sites?
# >Could we use other distributions (not in FastME, but in general)? Give a couple
# >suggestions and think about the pros and cons.

# >**Optional**: Explore some of the other substitution models available in FastME.

# # Extra: implementing Neighbor-Joining

# Implementing the neighbor-joining algorithm is fairly easy (at least if we
# don't care *too* much about efficiency). The code below is a fairly minimal
# implementation of the NJ algorithm implemented in such a way as to show the
# algorithm in action (generating the tree directly in Newick format on the
# go):

## Main NJ algorithm
using Printf  # for pretty printing

function neighbor_joining(matrix, taxa)
    clades = copy(taxa)
    nodes = collect(1:length(taxa))
    n = length(nodes)
    while length(nodes) > 1
        ## get the next neighbors to join
        (i, j, a, b, di, dj), new_dist = get_neighbors_to_join(matrix, nodes)
        ## join the chosen nodes to a new clade
        push!(clades, "($(clades[i]):$di,$(clades[j]):$dj)")
        ## update the nodes that are still left to join
        nodes[a] = n + 1
        deleteat!(nodes, b)
        ## update the matrix with the new node
        matrix = [[matrix ; new_dist[1:end-1]'] new_dist]
        ## increment internal node counter
        n += 1
        ## this will print out information so we can see the algorithm in action
        @printf "joining nodes %2d and %2d, creating internal node %2d\n" i j n
        @printf "subtree below new node %2d:\n\t%s\n" n clades[end]
        @printf "%d nodes left to join:\n\t" length(nodes)
        println(nodes, "\n", "_"^80)
    end
    clades[end]
end

## Function for a single 'join' operation
function get_neighbors_to_join(matrix, nodes)
    r = length(nodes)
    minindex = (0, 0, 0, 0, Inf)
    for a=1:r, b=a+1:r
        i = nodes[a]
        j = nodes[b]
        ## This is the neighbor joining optimality criterion, the two nodes that
        ## lead to the lowest value of `x` below are chosen in this iteration
        ## of the NJ algorithm to join and for an internal node of the tree.
        x = (r-2)*matrix[i,j] - sum([matrix[i,k] + matrix[j,k] for k in nodes])
        minindex = x < minindex[end] ? (i, j, a, b, x) : minindex
    end
    i = minindex[1]
    j = minindex[2]
    di, dj, new_distances = get_nj_distance(matrix, nodes, i, j, r)
    return (i, j, minindex[3], minindex[4], di, dj), new_distances
end

## This is the formula to compute the branch lengths leading to 
## nodes i and j, which will be joined, as well as the distances
## between the new internal node and all other nodes left to join 
function get_nj_distance(matrix, nodes, i, j, r)
    a = sum([matrix[i,k] for k in nodes])
    b = sum([matrix[j,k] for k in nodes])
    if r != 2  # we are not joining the two last nodes (root)
        di = 0.5*matrix[i,j] + 1.0/(2r - 4)*(a - b)
        dj = 0.5*matrix[i,j] + 1.0/(2r - 4)*(b - a)
    else  # we are joining the two last nodes (generating the root)
        di = matrix[i,j]
        dj = 0.
    end
    new_distances = 0.5 .* (matrix[i,:] .- di .+ matrix[j,:] .- dj)
    return di, dj, [new_distances ; 0.]
end

# Then, given a distance matrix, we can use the code like this
matrix, taxa, ntaxa = readmatrix("_assets/phylocourse/distance/18SrRNA_20_JCGamma_matrix.txt")
neighbor_joining(matrix, taxa)

# You can check this against FastME's NJ implementation, it should be correct. We have everything in place to start building our own extremely redundant phylogenetics library.

# >**Exercise**: For the diehards, try to understand the code and perhaps reimplement it in your programming language of choice.

# ------------------------------------------------------------------------------

# [^fastmeonline]: Note that you can also run FastME online [here](http://www.atgc-montpellier.fr/fastme/).

# [^ancientdata]: These data sets are from rather ancient days in phylogenetics, and are stored in the so-called PHYLIP format. In the original [PHYLIP format](http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles), taxon names were restricted to 10 characters, hence the weird truncated names for the sequences.

# [^commandline]: For those unfamiliar with the command line, it will probably be easiest to put the FastME executable for your operating system (for windows this is the `fastme.exe` file) together with the data files in a separate new folder. On windows, you can then from within your file explorer application do `Shift+right click` and click on `open command window here` or `open PowerShell window here`. Then you should be able to run the `fastme.exe` executable, e.g. `fastme -i 18SrRNA_20.phy -O 18SrRNA_20_JCmatrix.txt -d JC` provided the data file `18SrRNA_20.phy` is in the same directory as the `fastme.exe` executable. For Mac and linux users I would recommend a similar approach, make a directory where you put the executable and the data files, and open a terminal in that directory (Linux users, you know how to do this, MacOS users, I can't help you, but google is your friend). See also some notes [here](/phylocourse/install-windows)

# using Literate #src
# Literate.markdown(@__FILE__, "phylocourse", documenter=false) #src
