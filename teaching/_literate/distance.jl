# [back](/teaching/)
# \toc
# # Distance-based phylogenetic inference

# >**Questions** will be marked in blocks like this, try to formulate a brief answer to them.

# Distance-based phylogenetic methods proceed by computing a **pairwise distance matrix** for an input alignment based on a [substitution model](../submod). Usually the estimated distance is the maximum likelihood estimate (MLE) of the pairwise distance under the assumed substitution model. After computing the distance matrix, a tree is inferred by using some sort of clustering algorithm or least-squares estimation step. Distance-based phylogenetic inference is thus essentially a two-step process:

# 1. Compute distances
# 2. Infer the tree, assuming the distances

# This two step procedure is both the strength and weakness of distance-based methods. By collapsing the sequence data with $n$ sequences and $m$ sites in a single $n \times n$ matrix it dramatically reduces the data, making tree inference coputationally very fast. However, by reducing the rich sequence data to a matrix of numbers it throws away a lot of potentially interesting evolutionary information. Another issue is that when using distance-based methods the estimated distances from step 1 *are treated as observed data* in step 2. However, distances are themselves estimates, associated with some uncertainty (a distance estimate has a variance for instance), and this uncertainty in the distances is neglected in step 2. Both of these issues are solved in ML and Bayesian phylogenetic inference, however at the price of a strongly increased computational cost.

# Because of their speed, distance methods are still quite often used. Some ML tree inference programs (for instance IQ-TREE, see the ML tutorial) start their tree search algorithm from a distance-based tree for instance. Also popular packages for phylogenomic and comparative genomic inference such as [OrthoFinder](https://github.com/davidemms/OrthoFinder) use distance-based phylogenetic methods in their workflow.

# ## Software and data

# 1. Download the **FastME** software at [http://www.atgc-montpellier.fr/fastme](http://www.atgc-montpellier.fr/fastme/binaries.php) (on the bottom of the page, click the download button). In the downloaded folder you will find a binaries directory, identify the binary (executable) for your operating system and put it in some folder of your convenience.[^fastmeonline]

# 2. To view trees, I recommend the [FigTree](https://github.com/rambaut/figtree/releases/tag/v1.4.4) tool. FigTree requires Java, but that should be available on most machines. Download the executable for your operating system from the link above (`.zip` file for Windows users, `.tgz` for *nix users, I guess the `.dmg` file is something for Mc OSX users(?))

# 3. We will use two small 'tree of life' data sets. One [18SrRNA data set with 20 taxa](/assets/teaching/data/18SrRNA_20.phy) and another [18SrRNA data set with 45 taxa](/assets/teaching/data/18SrRNA_45.phy).

# ## FastME

# FastME is probably the fastest implementation of distance-based phylogenetics methods available today. It is a software tool that can both be used to compute distance matrices and infer trees using **Neighbor-Joining**, **least-squares**, **minimum evolution** and related distance matrix based methods.

# There are two ways to run FastME. There is an interactive mode (inherited from the influential [PHYLIP software package](http://evolution.genetics.washington.edu/phylip.html)) and a command line mode. Personally, I find the command line mode much more convenient. First get a look at the help message
# ```
# fastme -help
# ```

# If you have read the course notes from prof. Van de Peer, most of the options listed there (but not all, no worries) should make some sense to you. To infer a distance matrix with the Jukes & Cantor model you can run something like this from the command line[^commandline]:
# ```
# fastme -i 18SrRNA_20.phy -O 18SrRNA_20_JCmatrix.txt -dJC69
# ```
# Try to run this command and see whether you manage to make it work. The above command should generate a distance matrix. Let's make a heatmap to visualize it:

using Plots, DataFrames
function readmatrix(file)  # little function to read the FastME distance matrix
    lines = [split(l) for l in readlines(file)[2:end] if l != ""]
    DataFrame((;[Symbol(l[1])=>map(x->parse(Float64, x), l[2:end]) for l in lines]...))
end

df = readmatrix(joinpath(@__DIR__, "_assets/teaching/distance/18SrRNA_20_JCmatrix.txt"))
matrix = Matrix(df)
ntaxa = size(matrix)[1]
taxa = string.(names(df))
heatmap(matrix, yticks=(1:ntaxa, taxa), xticks=(1:ntaxa, taxa),
    xrotation=45, size=(700,650))
savefig("_assets/teaching/distance/hm1.svg") # hide

# ![](/assets/teaching/distance/hm1.svg)

# >**Question**: does this distance matrix make sense to you? Can you spot clades of more closely related species already?

# ## Clustering methods

# **Clustering methods** in distance-based phylogenetics are not different from hierarchical clustering methods used in other applications, such as unsupervised machine learning for instance. Both the UPGMA and WPGMA clustering methods are specific cases of *average linkage clustering* (with different formula's used to compute the distance between two already existing clusters). Clustering methods can befound in many packages for scientific computing for most programming languages (here I'll use the `Clustering.jl` package for the julia progamming language, note that average linkage clustering as it is usually implemented is identical to what is called the WPGMA method in phylogenetics).

using Clustering, StatsPlots
hcl = hclust(matrix, linkage=:average)
plot(
    plot(hcl, xticks=false),
    heatmap(matrix[hcl.order,hcl.order], colorbar=false,
        yticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xrotation=45),
    layout=grid(2, 1, heights=[0.2,0.8]), size=(600,750))

savefig("_assets/teaching/distance/wpgma.svg") # hide

# ![](/assets/teaching/distance/wpgma.svg)

# A pretty figure in my humble opinion.

# >**Question**: how can you see in one glance that this is an ultrametric tree?

# >**Question**: based on your knowledge of the tree of life, can you identify where the phylogeny is (very likely) wrong?

# >**Question**: Clustering methods like UPGMA and WPGMA produce rooted trees. Did the clustering algorithm identify the correct root?

# >**Question**: How do you interpret the branch lengths, what is the associated length 'unit'?

# Now let's infer a tree using $\Gamma$ **distributed rates across sites** (i.e. $\Gamma$ of 'Gamma' distances), still using the Jukes-Cantor substitution model.

# ```
# fastme -i 18SrRNA_20.phy -O 18SrRNA_20_JCGamma_matrix.txt -dJC69 -g
# ```

matrix = Matrix(readmatrix(joinpath(@__DIR__, "_assets/teaching/distance/18SrRNA_20_JCGamma_matrix.txt")))
hcl = hclust(matrix, linkage=:average)

plot(
    plot(hcl, xticks=false),
    heatmap(matrix[hcl.order,hcl.order], colorbar=false,
        yticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xrotation=45),
    layout=grid(2, 1, heights=[0.2,0.8]), size=(600,750))

savefig("_assets/teaching/distance/wpgma2.svg") # hide

# ![](/assets/teaching/distance/wpgma2.svg)

# >**Question**: did the topology change? Did the branch lengths change? Why?

# ## Neighbor-joining

# Neighbor-joining (NJ) is another method for distance-based phylogenetics. It is also a clustering method, but one that does not produce ultrametric trees. The default output of FastME includes a tree inferred using NJ (in the `.nwk` file). To run tree inference with NJ for an input alignent or distance matrix, run `fastme` with the `-o` option. For instance

# ```
# fastme -i 18SrRNA_20.phy -o 18SrRNA_20_JC.nwk -dJC69
# ```

# Open the NJ tree in the `.nwk` file generated by FastME in FigTree (should be straightforward).

# >**Question**: How can you see at one glance that this is not an ultrametric tree? What does this mean in terms of assumptions on the substitution rate?

# >**Question**: Neighbor-joining infers an unrooted tree, can you see how the fact that the tree is unrooted is represented in FigTree? Where should you root this tree (based on your knowledge of the tree of life)? Select the branch where you think the root should be (click on it) and hit the `reroot` button in figtree to root the tree.

# >**Question**: Does the NJ tree make more sense than the WPGMA tree? What do you think is causing this?

# Now infer a tree using $\Gamma$ distances

# ```
# fastme -i 18SrRNA_20.phy -o 18SrRNA_20_JC.nwk -dJC69 -g
# ```

# >**Question**: What (if anything) is changing? Experiment with the $\alpha$ parameter of the Gamma distribution by using for instance `-g0.5` in the FastME command. What happens? Below you can see a graph of the Gamma distribution for different values of $\alpha$ to help you interpret the results.

using Distributions
p = plot(title="The Gamma distribution with mean 1")
for α in [0.1, 0.25, 0.5, 1.0, 5.0, 10., 100.]
    plot!(p, Gamma(α, 1/α), label="\\alpha = $α", xlim=(0,5), ylim=(0,5))
end
savefig(p, "_assets/teaching/distance/gamma.svg") # hide

# ![](/assets/teaching/distance/gamma.svg)

# ------------------------------------------------------------------------------

# [^fastmeonline]: Note that you can also run FastME online [here](http://www.atgc-montpellier.fr/fastme/).

# [^commandline]: For those unfamiliar with the command line, it will probably be easiest to put the FastME executable for your operating system (for windows this is the `fastme.exe` file) together with the data files in a separate new folder. On windows, you can then from within your file explorer application do `Shift+right click` and click on `open command window here` or `open PowerShell window here`. Then you should be able to run the `fastme.exe` executable, e.g. `fastme -i 18SrRNA_20.phy -O 18SrRNA_20_JCmatrix.txt -d JC` provided the data file `18SrRNA_20.phy` is in the same directory as the `fastme.exe` executable. For Mac and linux users I would recommend a similar approach, make a directory where you put the executable and the data files, and open a terminal in that directory (Linux users, you know how to do this, MacOS users, I can't help you, but google i your friend).