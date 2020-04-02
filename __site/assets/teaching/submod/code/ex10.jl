# This file was generated, do not modify it. # hide
Pn = Pmatrix(0.05)^10  # get the 10-step transition probabilities
x = translate(original)
y = translate(evolved)
site_probabilities = [Pn[j,i] for (i,j) in zip(x,y)]
sequence_probability = prod(site_probabilities)