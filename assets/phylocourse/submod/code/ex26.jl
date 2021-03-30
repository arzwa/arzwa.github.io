# This file was generated, do not modify it. # hide
distance_JC(p) = -0.75 * log(1. - 4p/3)
p = mapreduce(!=, +, seqa, seqb)/length(seqa)
distance_JC(p)