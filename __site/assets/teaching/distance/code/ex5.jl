# This file was generated, do not modify it. # hide
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