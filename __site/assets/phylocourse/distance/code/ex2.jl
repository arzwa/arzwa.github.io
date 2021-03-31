# This file was generated, do not modify it. # hide
differences = [seqa[i] != seqb[i] for i=1:length(seqa)]
@show differences
num_differences = sum(differences)
x = num_differences/length(seqa)