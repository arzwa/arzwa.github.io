# \toc
# [back](/teaching/)
# # Distance-based phylogenetic inference

# ## Software requirements

# Download the **FastME** software at [http://www.atgc-montpellier.fr/fastme](http://www.atgc-montpellier.fr/fastme/binaries.php) (on the bottom of the page, click the download button). In the downloaded folder you will find a binaries directory, identify the binary (executable) for your operating system and put it in some folder of your convenience.

# ## FastME

# FastME is probably the fastest implementation of distance-based phylogenetics methods available today. It is a software tool that can both be used to compute distance matrices and infer trees using Neighbor-Joining, least-squares and related distance matrix based methods.

# There are two ways to run FastME. There is an interactive mode (inherited from the influential [PHYLIP software package](http://evolution.genetics.washington.edu/phylip.html)) and a command line mode. Personally, I find the command line mode much more convenient. First get a look at the help message
# ```
# fastme -help
# ```

# If you have read the course notes from prof. Van de Peer, most of the options listed there (but not all, no worries) should make some sense to you.

# To infer a distance matrix with the Jukes & Cantor model
# ```
# fastme -i data.phy -O JC_dist.txt -d JC69
# ```

using Plots, DataFrames
function readmatrix(file)  # little function to read the FastME distance matrix
    lines = [split(l) for l in readlines(file)[2:end] if l != ""]
    DataFrame((;[Symbol(l[1])=>map(x->parse(Float64, x), l[2:end]) for l in lines]...))
end

df = readmatrix(joinpath(@__DIR__, "_assets/teaching/distance/JC_dist.txt"))
matrix = Matrix(df)
ntaxa = size(matrix)[1]
heatmap(matrix, yticks=(1:ntaxa, string.(names(df))), size=(700,650))
savefig("_assets/teaching/distance/hm1.svg") # hide

# ![](/assets/teaching/distance/hm1.svg)

# ## Clustering methods

# **Clustering methods** in distance-based phylogenetics are not different from hierarchical clustering methods used in other applications, such as unsupervised machine learning for instance. Both the UPGMA and WPGMA clustering methods are specific cases of *average linkage clustering* (with different formula's used to compute the distance between two already existing clusters). Clustering methods can befound in many packages for scientific computing for most programming languages

using Clustering, StatsPlots
hcl = hclust(matrix, linkage=:average)
plot(
    plot(hcl, xticks=false),
    heatmap(matrix[hcl.order,hcl.order], colorbar=false,
        yticks=(1:ntaxa, string.([names(df)[i] for i in hcl.order]))),
    layout=grid(2, 1, heights=[0.2,0.8]), size=(600,750))

savefig("_assets/teaching/distance/wpgma.svg") # hide

# ![](/assets/teaching/distance/wpgma.svg)

# A pretty figure in my humble opinion.

# ## Neighbor-joining
