# This file was generated, do not modify it.

seqa = "ATCGGGCTAGC"
seqb = "TTCGGCTTACC";

differences = [seqa[i] != seqb[i] for i=1:length(seqa)]
@show differences
num_differences = sum(differences)
x = num_differences/length(seqa)

d = -0.75 * log(1. - 4x/3)

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

matrix, taxa, ntaxa = readmatrix("_assets/phylocourse/distance/18SrRNA_20_JCGamma_matrix.txt")
hcl = hclust(matrix, linkage=:average)

plot(
    plot(hcl, xticks=false),
    heatmap(matrix[hcl.order,hcl.order], colorbar=false,
        yticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xticks=(1:ntaxa, [taxa[i] for i in hcl.order]),
        xrotation=45),
    layout=grid(2, 1, heights=[0.2,0.8]), size=(600,750))

savefig("_assets/phylocourse/distance/wpgma2.svg") # hide

using Distributions
p = plot(title="The Gamma distribution with mean 1")
for α in [0.1, 0.25, 0.5, 1.0, 5.0, 10., 100.]
    plot!(p, Gamma(α, 1/α), label="\\alpha = $α", xlim=(0,5), ylim=(0,5))
end
savefig(p, "_assets/phylocourse/distance/gamma.svg") # hide

# Main NJ algorithm
using Printf  # for pretty printing

function neighbor_joining(matrix, taxa)
    clades = copy(taxa)
    nodes = collect(1:length(taxa))
    n = length(nodes)
    while length(nodes) > 1
        # get the next neighbors to join
        (i, j, a, b, di, dj), new_dist = get_neighbors_to_join(matrix, nodes)
        # join the chosen nodes to a new clade
        push!(clades, "($(clades[i]):$di,$(clades[j]):$dj)")
        # update the nodes that are still left to join
        nodes[a] = n + 1
        deleteat!(nodes, b)
        # update the matrix with the new node
        matrix = [[matrix ; new_dist[1:end-1]'] new_dist]
        # increment internal node counter
        n += 1
        # this will print out information so we can see the algorithm in action
        @printf "joining nodes %2d and %2d, creating internal node %2d\n" i j n
        @printf "subtree below new node %2d:\n\t%s\n" n clades[end]
        @printf "%d nodes left to join:\n\t" length(nodes)
        println(nodes, "\n", "_"^80)
    end
    clades[end]
end

# Function for a single 'join' operation
function get_neighbors_to_join(matrix, nodes)
    r = length(nodes)
    minindex = (0, 0, 0, 0, Inf)
    for a=1:r, b=a+1:r
        i = nodes[a]
        j = nodes[b]
        # This is the neighbor joining optimality criterion, the two nodes that
        # lead to the lowest value of `x` below are chosen in this iteration
        # of the NJ algorithm to join and for an internal node of the tree.
        x = (r-2)*matrix[i,j] - sum([matrix[i,k] + matrix[j,k] for k in nodes])
        minindex = x < minindex[end] ? (i, j, a, b, x) : minindex
    end
    i = minindex[1]
    j = minindex[2]
    di, dj, new_distances = get_nj_distance(matrix, nodes, i, j, r)
    return (i, j, minindex[3], minindex[4], di, dj), new_distances
end

# This is the formula to compute the branch lengths leading to
# nodes i and j, which will be joined, as well as the distances
# between the new internal node and all other nodes left to join
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

matrix, taxa, ntaxa = readmatrix("_assets/phylocourse/distance/18SrRNA_20_JCGamma_matrix.txt")
neighbor_joining(matrix, taxa)

