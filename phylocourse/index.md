# Molecular phylogenetics exercises

These are some extra notes, code examples and exercises on phylogenetic
inference. Some of the material here is more directed to those interested in
bioinformatics and statistics, other material is directed to a general
biological audience.

1. [Substitution models and distances](submod)
2. [Distance-based phylogeny inference](distance)
3. [Phylogeny inference using maximum likelihood](mliqtree)
4. [Model selection for ML phylogeny inference](modsel)
5. [Bayesian phylogenetic inference](bayes)
6. [Florida dentist scandal](hiv)
7. [SARS-CoV-2 origins](cov)

Most phylogenetics software tools are distributed as **command-line
applications**.  Some extra guidance on how to install and run command-line
based programs on Microsoft Windows can be found [here](install-windows).
(Linux users probably need no help with that, and since I don't have a mac PC
at my disposal, the apple people will have to google a bit).

**Code examples** are in the julia programming language, which is a great
language for scientific computing. Even if you are unfamiliar with Julia, most
code examples should be easy to read and readily translated to your favorite
programming language, which would be a good exercise.
If you want to follow along in julia, download the latest `julia` version for
your operating system from [https://julialang.org/](https://julialang.org/). If
you encounter `using <...>` statements, this imports external julia packages.
To install a julia package, open the REPL (by typing `julia` at the command
line for instance), type `]` (now you've entered the package manager `pkg`) and
type `add package`.

