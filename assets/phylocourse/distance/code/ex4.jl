# This file was generated, do not modify it. # hide
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