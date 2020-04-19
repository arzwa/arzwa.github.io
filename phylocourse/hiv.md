@def title = "dentist hiv"
@def author = "Arthur Zwaenepoel"
@def hasmath = true
@def hascode = true

[back](/phylocourse/)

# Florida dentist scandal

We'll consider a different data set now. In the early 90s, a [dentist](https://en.wikipedia.org/wiki/David_J._Acer) was accused of infecting several of his patients with HIV during surgical procedures. After a "low-risk" patient was diagnosed with HIV, other patients were screened, a coupe of which had HIV.

![](/assets/phylocourse/img/tabloid.jpg)

In the file `hiv-dentist.fasta` you'll find sequences from the V3 region of the *env* gene of the HIV virus from the dentist, the patients and some local controls (AIDS patients from Florida that had no relationship to the dentist whatsoever). These sequences are unaligned, so you'll first need to align them. You can download software like MUSCLE, MAFFT or PRANK for that, or alternatively, you could run any of these tools online at EBI [`https://www.ebi.ac.uk/Tools/msa/`](https://www.ebi.ac.uk/Tools/msa/).

After you have obtained an alignment, you can visualize it with some tool like `aliview`, `seaview` or `alv`, or again, you can use some of the webservices at EBI like [`Mview`](https://www.ebi.ac.uk/Tools/msa/mview/).

Now infer a tree with IQ-TREE, try `-m JC`, `-m JC+G` and ModelFinder (by omitting the `-m` option). Remember that you can include bootstrapping by setting `-bb <number of bootstrap replicates you want>`.

>- Based on the phylogeny, what can you conclude about the crime case?
>- Would your conclusion change when using different substitution models?
>- Why do you think they chose to sequence the *env* gene of HIV for the phylogenetic analyses?
>- Have a look at the [original paper](https://science.sciencemag.org/content/256/5060/1165), how did they conduct the phylogenetic analysis? Do you obtain similar conclusions?
