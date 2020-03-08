@def title = "Ks distributions with wgd"
@def hascode = true
@def hasmath = true


# Computing $\ks$ distributions using `wgd`

Already back in 2017, [I wrote a python library](https://github.com/arzwa/wgd)
for computing so-called $\ks$[^ks] distributions, which can serve to detect
ancient whole-genome duplications from genomic data. The goal was to provide an
easy and robust pipeline to do exactly that, and it seems some people are in
fact using it. It is not very actively maintained (and I wish I'd find the time
to finish a complete rewrite, register the package and rewrite documentation)
but it does the job nicely and efficiently. Here I'll show how I generally use
the library in my research.

First, installation. I have used wgd with Python3.5+ (3.5, 3.6, 3.7 and 3.8).
I generally recommend to use a virtual environment for wgd, that way you
won't have dependency issues or clashes. The following should suffice to 
install the wgd package in a Ubuntu-ish Linux environment.

```bash
sudo apt-get install python3-pip
pip3 install virtualenv --user
virtualenv venv -p python3
source venv/bin/activate
git clone https://github.com/arzwa/wgd.git
cd wgd
pip3 install .
```

If all went well, you should try running `wgd`

```md
(venv)  $ wgd                                                                                                  
Usage: wgd [OPTIONS] COMMAND [ARGS]...

  Welcome to the wgd command line interface!

                         _______
                         \  ___ `'.
         _     _ .--./)   ' |--.\  \
   /\    \\   ///.''\\    | |    \  '
   `\\  //\\ //| |  | |   | |     |  '
     \`//  \'/  \`-' /    | |     |  |
      \|   |/   /("'`     | |     ' .'
       '        \ '---.   | |___.' /'
                 /'""'.\ /_______.'/
                ||     ||\_______|/
                \'. __//
                 `'---'
  
  wgd  Copyright (C) 2018 Arthur Zwaenepoel
  This program comes with ABSOLUTELY NO WARRANTY;
  This is free software, and you are welcome to redistribute it
  under certain conditions;

  Contact: arzwa@psb.vib-ugent.be

Options:
  -v, --verbosity [info|debug]  Verbosity level, default = info.
  -l, --logfile TEXT            File to write logs to (optional)
  --version                     Print version number
  -h, --help                    Show this message and exit.

Commands:
  dmd  All-vs.-all diamond blastp + MCL clustering.
  kde  Fit a KDE to a Ks distribution.
  ksd  Ks distribution construction.
  mcl  All-vs.-all blastp + MCL clustering.
  mix  Mixture modeling of Ks distributions.
  pre  Check and optionally rename CDS files Example usage (renaming) wgd...
  syn  Co-linearity analyses.
  viz  Plot histograms/densities (interactively).
  wf1  Standard workflow whole paranome Ks.
  wf2  Standard workflow one-vs-one ortholog Ks.
```

Some third-party tools are required to run the example below. Specifically
we will use `diamond`, `mcl`, `mafft`, `codeml` and `fasttree`.

## The data

The basic input data is a bunch of CDS sequences. Here I'll use the CDS data
available for the bladderwort (*Utricularia gibba*). To download the data 
from PLAZA:

```bash
wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/cds.selected_transcript.ugi.fasta.gz
gunzip *gz
mv cds.*.ugi.fasta ugi.fasta
```

![](https://newfs.s3.amazonaws.com/taxon-images-1000s1000/Lentibulariaceae/utricularia-gibba-fl-dcameron.jpg)

There is a little tool in `wgd` to quickly check the input data. Since $\ks$ is
a codon-model based evolutionary distance, *it is absolutely crucial that the
input data are proper codon-sequences*. If you're paranoid about this being 
the case, you can use `wgd pre` to partition a fasta file in everything that
is nicely translatable from start to stop codon and all the rest.

```bash
wgd pre ugi.fasta
```

This will spit a lot to the screen:

```
$ wgd pre ugi.fasta 
2020-03-07 19:42:14: INFO	(0) Checking ugi.fasta
2020-03-07 19:42:14: ERROR	Translation error (First codon 'CCG' is not a start codon) in sequence UGI.ctg09795.24681.1
2020-03-07 19:42:14: ERROR	Translation error (First codon 'AAA' is not a start codon) in sequence UGI.ctg09872.24683.1
2020-03-07 19:42:14: ERROR	Translation error (Sequence length 638 is not a multiple of three) in sequence UGI.ctg09892.24686.1
2020-03-07 19:42:14: ERROR	Translation error (First codon 'TAG' is not a start codon) in sequence UGI.ctg09909.24689.1
[...]
2020-03-07 19:42:21: INFO	24145/25930 (93.12%) sequences are perfect CDS (in ugi.fasta.pre.good)
2020-03-07 19:42:21: INFO	1785/25930 (6.88%) sequences are not perfect CDS (in ugi.fasta.pre.bad)
```

About 93% of the input data is textbook CDS sequence. The rest has some issue,
and is written to the `.bad` file. Of course, not all CDS sequences start with
a start codon, and if you trust your data, you should definitely not throw
those away! (you should, however, not trust any sequence with a length that
is no multiple of three!). Here I'll continue working with the 'proper' CDS
data, now in the file `ugi.fasta.pre.good`.

**Note:** `wgd pre` can also be used to rename your sequences. If you have the kind
of fasta files with very long and awkward headers, I would recommend this, as it
will make files obtained later clearer. See the `--rename` and `--prefix` options.

## Obtaining the paranome

The first step is to obtain the **paranome**, or the collection of all paralogous
genes. This is esentially a big graph, where the nodes are all genes in the
genome, and edges represent homology (paralogy) relationships. Of course, in
the absence of different genomes to compare to, homology is no well-defined
concept, as we may as well assume all genes trace back to a common ancestor 
and are thus homologous. The goal of paranome inference is of course not to lump
everything together, but simply to *cluster the genome in reasonably fine-grained 
paralogous gene families*. Ideally, we'd like to obtain paralogous gene families
with a most recent common ancestor (MRCA) that is within the time frame where 
$\ks$ can be reliably estimated, providing us a clue as to what reasonably 
fine-grained could mean.

In general, common gene family clustering methods using Markov graph clustering 
(MCL) work well for this task. In wgd an approach based on all-vs.-all protein
similarity searches and MCL clustering is implemented. Both `blastp`, and the 
much faster `diamond` are supported. The following will run `diamond` + `mcl`
to obtain the paranome:

```
wgd dmd ugi.fasta.pre.good  
```

If one wishes to use the full data (not only the textbook-CDS sequences), the
`--nostrictcds` and `--ignorestop` options can be used. Other parameters of
interest are the $e$-value threshold used to construct the sequence-similarity
graph and the inflation factor for MCL, governing the coarseness of the
inferred clusters. To see all options for `wgd dmd` run `wgd dmd --help`.

If we started from an empty directory, by now we should have obtained the 
following (using the non-default `tree` command in Linux):

```bash
$ tree   
.
├── ugi.fasta
├── ugi.fasta.pre.bad
├── ugi.fasta.pre.good
└── wgd_dmd
    └── ugi.fasta.mcl
```

The paranome consists of how many families?

```
$ wc -l wgd_dmd/ugi.fasta.pre.good.mcl
3086 wgd_dmd/ugi.fasta.mcl
```

## The paranome $\ks$ distribution

To compute the $\ks$ distribution, we'll perform multiple sequence alignment
with MAFFT, maximum likelihood estimation of codon substitution model
parameters using `codeml` and do approximate phylogenetic tree inference using
FastTree. 

Since this takes quite a bit of time (codeml performs a rather expensive 
numerical optimization step), one would usually run this on a computing 
cluster on a couple of CPUs. This data set is still feasible to process 
on a laptop with 4 cores if we filter out some big families.

If you have the required tools installed, you should be able to run
the following:

```
wgd ksd ./wgd_dmd/ugi.fasta.mcl ./ugi.fasta.pre.good -mp 1000
```

(where the `mp 1000` option filters out all families with more than a thousand
pairs of sequences, futrher decrease this number if you want to quickly test
the whole thing, since with `mp 1000` it will still take about two hours). This
will run the analysis in parallel on four CPU cores. It starts with the biggest
families (which take more time). We can monitor the analysis in the terminal,
with output that looks like this.

```
$  wgd ksd ./wgd_dmd/ugi.fasta.mcl ./ugi.fasta.pre.good  
[...]
2020-03-07 19:01:43: INFO	Performing analysis on gene family GF_003081
2020-03-07 19:01:44: INFO	Performing analysis on gene family GF_003082
2020-03-07 19:01:44: INFO	Performing analysis on gene family GF_003083
2020-03-07 19:01:44: INFO	Performing analysis on gene family GF_003084
2020-03-07 19:01:44: INFO	Performing analysis on gene family GF_003085
2020-03-07 19:01:44: INFO	Performing analysis on gene family GF_003086
2020-03-07 19:01:45: INFO	Analysis done
2020-03-07 19:01:45: INFO	Making results data frame
2020-03-07 19:02:17: INFO	Removing tmp directory
2020-03-07 19:02:18: INFO	Computing weights, outlier cut-off at Ks > 5
2020-03-07 19:02:19: INFO	Generating plots
2020-03-07 19:02:19: INFO	Will plot **node-weighted** histograms
2020-03-07 19:02:20: INFO	Done
```

The generated output is in the `wgd_ksd` directory

```
$ tree wgd_ksd
wgd_ksd
├── ugi.fasta.pre.good.ks.svg
└── ugi.fasta.pre.good.ks.tsv
```

The main output is the `.tsv` file, which contains all the results computed
in the `wgd ksd` pipeline:

```
$ head wgd_ksd/ugi.fasta.pre.good.ks.tsv 
	AlignmentCoverage	AlignmentIdentity	AlignmentLength	AlignmentLengthStripped	Distance	Family	Ka	Ks	Node	Omega	Paralog1	Paralog2	WeightOutliersIncluded	WeightOutliersExcluded
UGI.Scf00428.21433.1__UGI.Scf00767.18891.1	0.20142	0.71127	2115.0	426.0	0.51238	GF_002963	0.3227	0.6327	2.0	0.51UGI.Scf00428.21433.1	UGI.Scf00767.18891.1	1.0	1.0
UGI.Scf00223.11926.1__UGI.Scf00720.18573.1	0.31967	0.74359	732.0	234.0	0.48336	GF_002829	0.2658	0.5949	2.0	0.4467	UGI.Scf00223.11926.1	UGI.Scf00720.18573.1	1.0	1.0
UGI.Scf00013.1922.1__UGI.Scf00108.8260.1	0.84677	0.87302	372.0	315.0	0.12595	GF_002055	0.0797	0.3366	2.0	0.2369	UGI.Scf00013.1922.1	UGI.Scf00108.8260.1	1.0	1.0
UGI.Scf00015.2084.1__UGI.Scf00046.4873.1	0.30266	0.78667	1239.0	375.0	0.32056	GF_000243	0.2379	1.5895	10.0	0.1497	UGI.Scf00015.2084.1	UGI.Scf00046.4873.1	1.0	1.0
UGI.Scf00015.2084.1__UGI.Scf00051.5189.1	0.24455	0.56766	1239.0	303.0	0.99481	GF_000243	0.6446	118.4869	11.00.0054	UGI.Scf00051.5189.1	UGI.Scf00015.2084.1	0.5	0.0
UGI.Scf00046.4873.1__UGI.Scf00051.5189.1	0.23002	0.5193	1239.0	285.0	1.07317	GF_000243	4.2678	71.3586	11.0	0.0598	UGI.Scf00051.5189.1	UGI.Scf00046.4873.1	0.5	0.0
UGI.Scf00051.5189.1__UGI.Scf00090.7423.1	0.23729	0.53061	1239.0	294.0	1.09756	GF_000243	13.2138	38.8342	12.0	0.3403	UGI.Scf00090.7423.1	UGI.Scf00051.5189.1	0.33333	0.0
UGI.Scf00015.2084.1__UGI.Scf00090.7423.1	0.3293	0.625	1239.0	408.0	0.61336	GF_000243	16.5071	17.2326	12.0	0.9579	UGI.Scf00090.7423.1	UGI.Scf00015.2084.1	0.33333	0.0
UGI.Scf00046.4873.1__UGI.Scf00090.7423.1	0.33898	0.62381	1239.0	420.0	0.69171	GF_000243	3.9448	34.5095	12.0	0.1143	UGI.Scf00090.7423.1	UGI.Scf00046.4873.1	0.33333	0.0
```

The default plot outputted by `wgd` (in `wgd_ksd/*.svg`) is rather ugly, but
gives a clear overview of the parameter estimates for all gene duplication
events, being $\ks, K_\mathrm{A}$ and $\omega$ (i.e. the sysnonymous distance,
nonsynonyous distance and nonsynonymous to synonymous substitution rate ratio
respectively).

![](../assets/ugi-ks.svg)

From this distribution it's quite clear there has been some stuff going on 
with this genome.  

## Computing one-to-one ortholog $\ks$ distributions

People often also estimate one-to-one ortholog $\ks$ distributions to obtain an
idea of the relative timing of WGDs and species divergences, divergence times
or relative substitution rates. Let us compare *Utricularia* to a close relative
*Erythranthe*. First get and clean the data

```bash
wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/cds.selected_transcript.egut.fasta.gz
gunzip *gz
mv cds.*.egut.fasta egu.fasta
wgd pre egu.fasta
```

To obtain one-to-one orthologs for between species comparisons we can again 
use `diamond` with `wgd dmd`:

```bash
wgd dmd egu.fasta.pre.good ugi.fasta.pre.good
```

```
$ wgd dmd egu.fasta.pre.good ugi.fasta.pre.good
2020-03-07 20:06:51: WARNING	dir wgd_dmd exists!
2020-03-07 20:06:58: WARNING	dir wgd_dmd exists!
2020-03-07 20:07:04: INFO	Multiple CDS files: will compute RBH orthologs
2020-03-07 20:07:04: INFO	egu.fasta.pre.good vs. ugi.fasta.pre.good
```

The ouput is, as expected, a bunch of gene pairs. These are *reciprocal best hits*
(RBHs). In other words, two genes $A$ and $B$  are considered one-to-one orthologs 
if $B$ has the highest sequence similarity to $A$ among the total set of genes
(excluding $A$ itself of course) and vice versa. Have a look at the file to 
convince yourself:

```
$ head wgd_dmd/egu.fasta.pre.good_ugi.fasta.pre.good.rbh
UGI.Scf00140.9494.1	Migut.L01925.1
UGI.Scf00873.19408.1	Migut.D02258.1
UGI.Scf01221.23016.1	Migut.K00468.1
UGI.Scf00487.21516.1	Migut.F00288.1
UGI.Scf00033.3844.1	Migut.G01067.1
UGI.Scf00611.21794.1	Migut.N02364.1
UGI.Scf00611.21793.1	Migut.N02367.1
UGI.Scf00131.9250.1	Migut.F01756.1
UGI.Scf00437.15873.1	Migut.F00666.1
UGI.Scf00015.2103.1	Migut.O00089.1
```

To compute the one-to-one ortholog $\ks$ distribution, we can again 
use `wgd ksd` intuitively. For the sake of comutational feasibility 
I'll take a random subset first

```
shuf wgd_dmd/egu.fasta.pre.good_ugi.fasta.pre.good.rbh wgd_dmd/egu-ugi.2500.rbh
```

```
$ wgd ksd wgd_dmd/egu-ugi.2500.rbh egu.fasta.pre.good ugi.fasta.pre.good
2020-03-07 20:10:02: INFO	
2020-03-07 20:10:02: INFO	codeml found
2020-03-07 20:10:02: INFO	
2020-03-07 20:10:02: INFO	
2020-03-07 20:10:02: INFO	Translating CDS file
Invalid codon CTN in UGI.ctg34642.25621.1                                                       
Invalid codon TCN in UGI.Scf00004.766.1
100% (50833 of 50833) |######################################################| Elapsed Time: 0:00:14 Time:  0:00:14
2020-03-08 12:40:40: WARNING	There were 56 warnings during translation
2020-03-08 12:40:40: INFO	Started whole paranome Ks analysis
2020-03-08 12:40:40: WARNING	Filtered out the 0 largest gene families because n*(n-1)/2 > `max_pairwise`
2020-03-08 12:40:40: WARNING	If you want to analyse these large families anyhow, please raise the `max_pairwise` parameter. 
2020-03-08 12:40:40: INFO	Started analysis in parallel (n_threads = 4)
2020-03-08 12:40:40: INFO	Performing analysis on gene family GF_000001
2020-03-08 12:40:40: INFO	Performing analysis on gene family GF_000002
[...]
2020-03-08 12:51:18: INFO	Performing analysis on gene family GF_002499
2020-03-08 12:51:18: INFO	Performing analysis on gene family GF_002500
2020-03-08 12:51:19: INFO	Analysis done
2020-03-08 12:51:19: INFO	Making results data frame
2020-03-08 12:51:38: INFO	Removing tmp directory
2020-03-08 12:51:38: INFO	Computing weights, outlier cut-off at Ks > 5
2020-03-08 12:51:38: INFO	Generating plots
2020-03-08 12:51:38: INFO	Will plot **node-weighted** histograms
2020-03-08 12:51:39: INFO	Done
```

Again the default plots are pretty ugly, but show what we need to see.

![](../assets/ugi-egu-ks.svg)

## Comparing multiple distributions

Say we want to make a nice plot combining the two distributions above. 
One can use the utilities iplemented in `wgd viz`, however, if you're
somewhat familiar with Python, R, julia, or any other programming 
language with nice plotting libraries I would defintely recommend 
generating your own plots. Here I show an example of how I would go on 
visualizing the two distributions computed above comparatively using 
Python.

I use the following packages

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
```

Then load the data

```python
egu_ugi = pd.read_csv("./wgd_ksd/egu.fasta.pre.good_ugi.fasta.pre.good.ks.tsv", sep='\t', index_col=0)
ugi = pd.read_csv("./wgd_ksd/ugi.fasta.pre.good.ks.tsv", sep='\t', index_col=0)
```

Now I filter the $\ks$ data and perform node averaging. When computing
paranome $\ks$ distributions `wgd ksd` infers a phylogenetic tree for
each paralogous famiy. Using these trees, we can weigh or average 
$\ks$ estimates for ancestral duplication events in families to account
for the fact that we have multiple redundant $\ks$ estimates pertaining
to the same divergence (duplication) event. I first filter the data
frame to only retain $\ks$ estimates in a reasonable range (here I 
take $0.01 < \hat{\ks} 5$), then I will average $\ks$ estimates that
correspond to the same duplication node in each paralogous family.
For the one-to-one ortholog distribution the latter is not necessary,
and I'll just apply the same $0.01 < \hat{\ks} < 5$ filter:

```python
lower, upper = 0.01, 5
ugi = ugi[ugi["Ks"] < upper]
ugi = ugi[ugi["Ks"] > lower]
ugi_ks = ugi.groupby(["Family", "Node"])["Ks"].apply(np.mean)
egu_ugi_ks = egu_ugi[egu_ugi["Ks"] < upper]["Ks"]
```

Now we can go on plotting:

```python
fig, ax = plt.subplots(1, 1, figsize=(6,4))
ax.hist(ugi_ks, bins=np.linspace(0.0, upper, 50), 
        color="teal", alpha=0.2, 
        label="$U.\ gibba$");
ax.hist(egu_ugi_ks, bins=np.linspace(0.0, upper, 50), 
        color="black", alpha=0.2, 
        label="$U.\ gibba$ vs. $E.\ guttata$");
ax.set_xlabel("$\hat{K_\mathrm{S}}$", fontsize=14)
ax.set_ylabel("# of duplications/ortholog pairs", fontsize=14)
ax.set_xlim(0, upper)
ax.legend(frameon=False);
sns.despine(offset=5)
```

![](../assets/py-ks.svg)


-------------------
 
[^ks]: Here I'll use $\ks$ and $K_\mathrm{A}$ instead of the more commonly used $d\mathrm{S}$ and $d\mathrm{N}$. This is not due to personal preference (in fact, I would prefer following Ziheng Yang's notation), but simply because this is customary in the WGD related literature.

[^sunfish]: I have no reason for examining this species (except for its awesomeness), and have no idea whether there is something interesting to be found (although it should show the remnants of a teleost-specific WGD event).
