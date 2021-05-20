@def title = "SARS-CoV-2 origins"
@def author = "Arthur Zwaenepoel"
@def hasmath = true
@def hascode = true

[back](/phylocourse/)

# SARS-CoV-2 origins

In this exercise we will take a look at the origins of SARS-CoV-2 -- the virus
causing Covid-19 -- by means of phylogenetic analysis. It is widely believed to
be of zoonotic origin like other major viral outbreaks (SARS, MERS, Ebola, HIV), 
although a lab origin is as far as I know not disproved[^1]. 

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
   of SARS-CoV-2? Could this disprove the hypothesis that the pandemic derives
   from a lab strain?
3. Why is recombination a problem for standard phylogenetic inference? Which
   assumption(s) can it possibly violate?
4. Why haven't we cared much about recombination in previous exercises? Should
   we have?

[^1]: I believe it is currently not possible to rule out either hypothesis, although much of the scientific establishment does seem to rule out the lab-origin theory. See for instance [this article](https://thebulletin.org/2021/05/the-origin-of-covid-did-people-or-nature-open-pandoras-box-at-wuhan/) for some of the reasons why the matter is complicated (politically, scientifically, biologically).
