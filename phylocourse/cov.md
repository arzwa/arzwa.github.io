@def title = "SARS-CoV-2 origins"
@def author = "Arthur Zwaenepoel"
@def hasmath = true
@def hascode = true

[back](/phylocourse/)

# SARS-CoV-2 origins

In this exercise we will take a look at the origins of SARS-CoV-2 -- the virus
causing Covid-19 -- by means of phylogenetic analysis. As is well-known,
SARS-CoV-2 is a [zoonotic](https://en.wikipedia.org/wiki/Zoonosis) pathogen,
meaning it jumped from some other species, which may act as a
[reservoir](https://en.wikipedia.org/wiki/Natural_reservoir), to a human host.
Many important viral pathogens are zoonotic (e.g. SARS, MERS, Ebola, HIV, ...)
and there is an obvious interest in identifying potentially zoonotic lineages
and their reservoirs.

SARS-CoV-2 is a member of the Sarbecovirus clade within Coronaviridae. It was
early on found to be relatively closely related to a sarbecovirus from
Horseshoe bats, RaTG13, and later also viral isolates from pangolins were
suggested as potential close relatives. Sarbecoviruses are, unlike many other
RNA viruses, undergoing quite a lot of recombination, which makes phylogenetic
inference challenging. [Boni *et
al.*](https://www.nature.com/articles/s41564-020-0771-4) identified
non-recombining regions (NRR) and conducted phylogenetic analyses using those. 

We will use the data from Boni *et al.*, which can be obtained from [this
Github repository](https://github.com/plemey/SARSCoV2origins). 

1. Infer a tree using maximum likelihood inference for the 'sarbecovirus' data
   set (either their NRR1, NRR2 or NRA3 data sets). Which model assumptions
   will you use? 
2. Assess the tree you found, what is your hypothesis on the zoonotic origins
   of SARS-CoV-2?
3. Why is recombination a problem for standard phylogenetic inference? Which
   assumption(s) can it possibly violate?
4. Why haven't we cared much about recombination in previous exercises? Should
   we have?
