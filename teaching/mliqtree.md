@def hascode = true

[back](/teaching/)

\toc

For the exercises, we will use the computationally *very* efficient software IQ-TREE (Nguyen *et al.* 2015). This is a command line program, similar to the FastME or PHYLIP programs but not interactive. Please download IQ-TREE for your operating system at [http://www.iqtree.org/](http://www.iqtree.org/).

Once you obtain the program, I'd recommend to put the executable in some dedicated directory[^dir] along with the data files. In this section, we'll use the 18S rRNA data sets for [20 taxa](/assets/teaching/data/18SrRNA_20.phy) and [45 taxa](/assets/teaching/data/18SrRNA_45.phy) again. But feel free to follow along with any data set you like. For a note on using command line programs in Windows, see [the footnote here](../distance/#fndef:commandline).

> Try to answer the questions shown in boxes like this one below. For some of the questions on substitution models you may want to google a bit or have a look in the IQ-TREE manual or the course notes or something similar.

## Basic tree inference

### The Jukes-Cantor model

Let us first reconsider the 45 species 18S rRNA data set. Remember that we had some problematic clades in our distance-based phylogenies, where we found different topologies depending on whether we used distances computed with the Gamma model of rate heterogeneity across sites or not and depending on the $\alpha$ parameter we chose. Let's see what we get when using ML.

Put the `18SrRNA_45.txt` file in the directory with the IQ-TREE executable. Fire up IQ-TREE and run the following command:

```bash
iqtree -s 18SrRNA_45.txt -m JC -pre JC
```

The `-s` option is used for specifying the input multiple sequence alignment file, the `-m` option is used for setting the substitution model (here Jukes & Cantor) and the `-pre` option sets a prefix for the output file names. Note that IQ-TREE generates three output files: a `.log`, `.treefile` and `.iqtree` file. Most interesting stuff is in the `.iqtree` file. The tree (which can be loaded FigTree or [icytree](https://icytree.org/) for example) is in the `.treefile`. Note that likelihoods are always shown as *log*-likelihoods to prevent numerical inaccuracies.

> - Examine the stuff IQ-TREE prints to the screen, can you figure out the different steps IQ-TREE uses to find the ML tree?
> - What is the log-likelihood associated with the ML tree?
> - Do you think it is a problem that the likelihood is such an absurdly small number?
> - How does the phylogeny compare with the distance based phylogeny under the JC model? Can you explain why?

### Gamma distributed rates across sites

Now let's compare this with the JC+Gamma model.

```bash
iqtree -s 18SrRNA_45.txt -m JC+G  -pre JCG
```

> - What is the MLE for the $\alpha$ parameter of the Gamma distribution?
> - Does the topology differ between the ML tree found with JC and JC+G?

### The bootstrap

Since different methods gave different results for some clades, it would definitely be worthwhile to try to get an idea of how well the different clades in the tree are supported by the data. As you saw in the course, the most commonly used approach to evaluate the statistical support of a phylogenetic tree in ML or distance based phylogenetics is the bootstrap. We will use the "ultrafast bootstrap" as implemented in IQ-TREE (which is not exactly the same as Felsenstein's original nonparametric bootstrapping approach[^ultrafast]).  Run IQ-TREE with 1000 ultrafast bootstrap replicates:

```bash
iqtree -s 18SrRNA_45.txt -m JC   -bb 1000 -pre JC_BSV
iqtree -s 18SrRNA_45.txt -m JC+G -bb 1000 -pre JCG_BSV
```

> - Are there any nodes for which you find low support?
> - Are the clades which were problematic in the distance-based analyses well-supported in the ML trees?

Now let us have a look at some other substitution models.

### The K2P model

```bash
iqtree -s msa.fasta -m K2P -pre K2P
```

> - What is the difference between the K2P and JC model?
> - How many more parameters does the K2P model have compared to the JC model?
> - What are the ML estimate(s) of the parameter(s)
> - Which model fits the data best according to the LRT (note: the critical value of the $\chi^2$ distribution with one degree of freedom at the 0.05 significance level is 3.845)

### The F81 model

```bash
iqtree -s 18SrRNA_45.txt -m F81 -pre F81
```

- What is the difference between the F81 model and the JC model?

## Model selection

IQ-TREE implements another program called ModelFinder (Kalyaanamoorthy *et al.* 2017) that allows to quickly (and approximately) find the best fitting substitution model. You can run it by just omitting the model specification with the `-m` flag from your commands.

### The 18S rRNA set again

Run IQ-TREE with model selection:

```bash
iqtree -s 18SrRNA_45.txt
```

> - Which model was selected? What are the features of this model compared to the models we used above?
> - Do you get a different tree?

[^ultrafast]: The ultrafast bootstrap of Minh *et al.* (2013) is for instance not only much faster, but also less biased. It is well known that Felsenstein's nonparametric bootstrap tends to be conservative (i.e. bootstrap support values tend to be underestimated relatively to the true statistical support of a clade), the ultrafast bootstrap is less so, and can be more easily interpreted as they tend to more closely resemble probabilities.

[^dir]: I sometimes notice that Windows users that are not particular computophiles are not familiar with the word *directory*, but generally employ the term *folder*. This seems to be related to some [philosophy that originated within Microsoft Windows](9https://en.wikipedia.org/wiki/Directory_%28computing%29#Folder_metaphor) where "There is a difference between a directory, which is a file system concept, and the graphical user interface metaphor that is used to represent it (a folder)."
