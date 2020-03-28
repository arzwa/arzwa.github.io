@def hascode = true
@def title = "SARS-CoV-2 phylogenetics"

\toc

In this exercise, we will have a look at the phylogenetic position of the human
coronavirus SARS-CoV-2 that is currently causing the outbreak of COVID-19. The
goal of this exercise is to go through a complete workflow for phylogenetic 
analysis. 

**Disclaimer:** I am in no way experienced with viral phylogenomics. I know
nothing about the available data, little about viral evolutionary
relationships, and barely anything about viral genomes. 

## Obtaining data

### Viral genomes

As armchair biologists, we tend not to go look for our data out there in the
wild (especially during a pandemic), but download it from some database.  Here
I describe how I got the data set we will be analyzing later. This is not
really part of the practical, but rather then providing you a data set
rightaway, it can be informative to see the typical annoying bits of
bioinformatics that necessarily come first. 

To get the data, I first did a search in the NCBI 'Nucleotide' database
using the following search term:

```
("Alphacoronavirus"[Organism] OR "Betacoronavirus"[Organism] OR 
 "Gammacoronavirus"[Organism]) AND (viruses[filter] AND refseq[filter])
```

I restrict the sequence data to the RefSeq database, which is a collection of
non-redundant and well-annotated sequences. As I am no expert in viral
phylogenomics, I considered this the best idea, since it would not be feasible
for me to assess data quality and redundancy efficiently. This provided me with
44 corona virus genome sequences.

I sent the results to a file, using the 'Complete record' option and the
GenBank format. Note that the reference sequence for nCov-2019 is the Wuhan
isolate apparently, with the following [genbank
entry](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=genbank&to=29903).

Now we could work with the entire genomes, aligning them and trying to do tree 
inference for the full viral genome alignments. Alternatively, we could analyze
the different viral proteins separetely. Besides being computationally less
involved, this has also the advantage of being somewhat clearer to work with.
I wrote the following little script using Biopython to obtain the coding DNA
sequences (CDS) and associated amino-acid (AA) sequences for each viral genome:

```py
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys

if len(sys.argv) != 3:
	print("usage: {} [genbank file] [prefix]".format(sys.argv[0]))
    exit()

path_nt = sys.argv[2] + "-nt"
path_aa = sys.argv[2] + "-aa"
os.mkdir(path_nt)
os.mkdir(path_aa)
for rec in SeqIO.parse(sys.argv[1], "genbank"):
    if rec.features:
        gid = rec.annotations["organism"].replace(" ", "_").replace("/", "_")
        recs = [SeqRecord(
            seq=ft.location.extract(rec).seq, 
            id=ft.qualifiers['protein_id'][0],
            description="") for ft in rec.features if ft.type == "CDS"]
        SeqIO.write(recs, os.path.join(path_nt, gid + ".fasta"), "fasta")
        trrec = [r.translate(id=r.id, description=r.description) for r in recs if len(r.seq)%3==0]
        SeqIO.write(trrec, os.path.join(path_aa, gid + ".fasta"), "fasta")
```

The script above can be run as follows:

```
python3 getdata.py sequence.gb genomes
```

Which gives us two directories, one with all CDS sequences for each genome and
one with all associated amino-acid sequences.

```
tree genomes-aa

genomes-aa
├── Bat_coronavirus_1A.fasta
├── Bat_coronavirus_CDPHE15_USA_2006.fasta
├── Bat_Hp-betacoronavirus_Zhejiang2013.fasta
├── Beluga_whale_coronavirus_SW1.fasta
├── Betacoronavirus_England_1.fasta
├── Betacoronavirus_Erinaceus_VMC_DEU_2012.fasta
├── Betacoronavirus_HKU24.fasta
├── Bovine_coronavirus.fasta
├── BtMr-AlphaCoV_SAX2011.fasta
├── BtNv-AlphaCoV_SC2013.fasta
├── BtRf-AlphaCoV_HuB2013.fasta
├── BtRf-AlphaCoV_YN2012.fasta
├── Camel_alphacoronavirus.fasta
├── Coronavirus_AcCoV-JC34.fasta
├── Feline_infectious_peritonitis_virus.fasta
├── Ferret_coronavirus.fasta
├── Human_coronavirus_229E.fasta
├── Human_coronavirus_HKU1.fasta
├── Human_coronavirus_NL63.fasta
├── Human_coronavirus_OC43.fasta
├── Infectious_bronchitis_virus.fasta
├── Lucheng_Rn_rat_coronavirus.fasta
├── Middle_East_respiratory_syndrome-related_coronavirus.fasta
├── Miniopterus_bat_coronavirus_HKU8.fasta
├── Mink_coronavirus_strain_WD1127.fasta
├── Murine_hepatitis_virus.fasta
├── Murine_hepatitis_virus_strain_JHM.fasta
├── NL63-related_bat_coronavirus.fasta
├── Pipistrellus_bat_coronavirus_HKU5.fasta
├── Porcine_epidemic_diarrhea_virus.fasta
├── Rabbit_coronavirus_HKU14.fasta
├── Rat_coronavirus_Parker.fasta
├── Rhinolophus_bat_coronavirus_HKU2.fasta
├── Rousettus_bat_coronavirus.fasta
├── Rousettus_bat_coronavirus_HKU10.fasta
├── Rousettus_bat_coronavirus_HKU9.fasta
├── Scotophilus_bat_coronavirus_512.fasta
├── Severe_acute_respiratory_syndrome_coronavirus_2.fasta
├── Severe_acute_respiratory_syndrome-related_coronavirus.fasta
├── Swine_enteric_coronavirus.fasta
├── Transmissible_gastroenteritis_virus.fasta
├── Turkey_coronavirus.fasta
├── Tylonycteris_bat_coronavirus_HKU4.fasta
└── Wencheng_Sm_shrew_coronavirus.fasta
```

The genomes seem to have between 5 and 14 genes, as can be quickly inspected 
using `grep -c '>' genomes-aa/* | cut -d':' -f2 | sort | uniq`.

### Gene families

With these viral genomes nicely in their respective files, we do not have our
gene families of interest yet. Since these are RefSeq records, we could
probably get gene families (e.g. the S-protein) for each genome directly by
looking at the gene names, but this is tricky (what if the names are not
consistent for instance?). I'd rather take a more general approach. We'll make
a blast database for our genomes, and fish out sequences using a bait gene of
interest (for instance the S-protein from SARS-CoV-2).

Well, actually, I'll use [diamond](https://github.com/bbuchfink/diamond)
instead of blast, which is way faster. First I make a diamond database:

```
cat genomes-aa/* > aa.fa
diamond makedb --in aa.fa --db aa
```

Now we fish for the S-protein (spike glycoprotein). This is the sequence of the
Wuhan isolate reference genome of SARS-Cov-2 (which you can easily parse out
the genbank entry linked to above):

``` 
>spike-SARS-CoV-2
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS
NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV
NNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLE
GKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT
LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETK
CTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISN
CVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIAD
YNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC
NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN
FNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITP
GTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSY
ECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI
SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE
VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDC
LGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAM
QMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN
TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRA
SANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPA
ICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDP
LQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL
QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDD
SEPVLKGVKLHYT
```

Run diamond using the above protein as query:

```
diamond blastp --db aa.dmnd -q spike-sarscov2.fasta > spike.hits
```

Now we obtain the relevant sequences using the following script
(I saved it as `gethits.py`)

```
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} <hits> <fastadir>".format(sys.argv[0]))
        exit()
    from Bio import SeqIO
    import pandas as pd
    import os
    import logging
    df = pd.read_csv(sys.argv[1], sep="\t", header=None)
    seqs = {}
    for f in os.listdir(sys.argv[2]):
        for seq in SeqIO.parse(os.path.join(sys.argv[2], f), "fasta"):
            if seq.id in list(df[1]):
                if seq[-1] == "*": seq = seq[0:-1]
                newid = "{}__{}".format(f.split(".")[0], seq.id)
                if newid in seqs: logging.warning("duplicate!")
                seqs[newid] = seq        
    if len(seqs) < len(df.index): 
        logging.warning("Found {}/{} sequences".format(
			len(seqs), len(df.index)))
    for (k,v) in seqs.items():
        print(">{}".format(k))
        print(v.seq)

```

Use it as follows: 

```
python3 gethits.py ./spike.hits ./genomes-aa > spike.aa.fasta
```

and to get the CDS sequences:

```
python3 gethits.py ./spike.hits ./genomes-nt > spike.nt.fasta
```

Clearly, you can use this approach to get other genes from the viral genome
database we assembled. Now we can do phylogenetics using the `spike.aa.fasta`
and/or `spike.nt.fasta` files.
 
## Multiple sequence alignment

In order to do phylogenetic inference, we need to have a hypothesis of
homology, that is to say, we need to know which characters in a bunch of
sequences trace back to a common ancestral character. Sadly, this is not 
given, and we have to infer this from the data. This is of course the 
problem of multiple sequence alignment (MSA).

In practice, inferring a MSA couldn't be easier:

```
mafft spike.aa.fasta > spike.aa.msa 
```

And we can have a look at this alignment using for instance the
[`alv`](https://github.com/arvestad/alv) program. Here's a little piece 
of the output of `alv spike.aa.msa`:

```
Human_coronavirus_NL63__YP_003767.1                                  --IISVVFVVLLSLLVFCCLSTGC-CGCCNCLTSSMRGCCDC-----
Coronavirus_AcCoV-JC34__YP_009380521.1                               --AIFLVLVLFSFMLLWCCCATGC-CGCCGMLGSACNGCCTK-----
Murine_hepatitis_virus__NP_045300.1                                  --LIGLAGVAVCVLLFFICCCTGCGSCCFK----KCGNCCDE-----
Rat_coronavirus_Parker__YP_003029848.1                               --LIGLAGVAVCVLLFFICCCTGCGSCCFK----KCGNCCDE-----
Infectious_bronchitis_virus__NP_040831.1                             --AIAFATIIFILILGWVFFMTGC-CGC----------CCGCFGIMP
Pipistrellus_bat_coronavirus_HKU5__YP_001039962.1                    GFIAGLVALALCVFFIL--CCTGCGTSCLGKL--KCNRCCDS-----
Severe_acute_respiratory_syndrome_coronavirus_2__YP_009724390.1      GFIAGLIAIVMVTIMLC--CMTSC-CSCLKGC-CSCGSCC-K-----
Rousettus_bat_coronavirus__YP_009273005.1                            GMIAGLVGLVLAVVLLL--CMTNC-CSCARGV-CSCKSC--A-----
Betacoronavirus_Erinaceus_VMC_DEU_2012__YP_009513010.1               GFIAGLVALALCLFFIL--CCTGCGTSCLGKI--NCTRCCDK-----
Severe_acute_respiratory_syndrome-related_coronavirus__NP_828851.1   GFIAGLIAIVMVTILLC--CMTSC-CSCLKGA-CSCGSCC-K-----
Human_coronavirus_HKU1__YP_173238.1                                  --LISFSFIIFLVLLFFICCCTGCGSACFS----KCHNCCDE-----
Lucheng_Rn_rat_coronavirus__YP_009336484.1                           --AIFLVLVLFSFLMLWCCCATGC-CGCCGLCGAACNGCCTK-----
Betacoronavirus_England_1__YP_007188579.1                            GFIAGLVALALCVFFIL--CCTGCGTNCMGKL--KCNRCCDR-----
Bovine_coronavirus__NP_150077.1                                      --LIGFAGVAMLVLLFFICCCTGCGTSCFK----KCGGCCDD-----
Bat_Hp-betacoronavirus_Zhejiang2013__YP_009072440.1                  GFIAGLIAIVMATIMLC--CMTSC-CSCFKGL-CACKRCCDS-----
Rabbit_coronavirus_HKU14__YP_005454245.1                             --LIGLAGVAVLVLLFFICCCTGCGTSCFK----KCGGCCDD-----
Middle_East_respiratory_syndrome-related_coronavirus__YP_009047204.1 GFIAGLVALALCVFFIL--CCTGCGTNCMGKL--KCNRCCDR-----
Human_coronavirus_OC43__YP_009555241.1                               --LICLAGVAMLVLLFFICCCTGCGTSCFK----KCGGCCDD-----
Turkey_coronavirus__YP_001941166.1                                   --AIGFACLIFILILGWVFFMTGC-CGC----------CCGCFGIIP
Rousettus_bat_coronavirus_HKU9__YP_001039971.1                       AMIAGIVGLVLAVIMLM--CMTNC-CSCFKGM-CDCRRCCGS-----
Tylonycteris_bat_coronavirus_HKU4__YP_001039953.1                    GFIAGLVALLLCVFFLL--CCTGCGTSCLGKM--KCKNCCDS-----
Human_coronavirus_229E__NP_073551.1                                  --CISVVLIFVVSMLLLCCCSTGC-CGFFSCFASSIRGCCE------
Murine_hepatitis_virus_strain_JHM__YP_209233.1                       --LIGLAGVAVCVLLFFICCCTGCGSCCFR----KCGSCCDE-----
BtMr-AlphaCoV_SAX2011__YP_009199609.1                                --IIVLVLILFTCLMLFCCCSTGC-CGIFSCMASSCGACCDI-----
Betacoronavirus_HKU24__YP_009113025.1                                --LIGLAGVAVLVLLFFVCCCTGCGSSCFK----KCGSCCDD-----
 ```

This however hides the crucial fact that *multiple sequence alignment is one
of the most important open problems in statistical phylogenetics*. Since every
inference is associated with uncertainty, ideally we'd like to obtain an idea
of the uncertainty in our inferred MSA, and take this uncertainty into account
in our phylogenetic analyses. In a Bayesian framework this is possible by
jointly inferring the alignment and phylogeny under a model of sequence
evolution that includes insertions and deletions (the TKF or Thorne - Kishino -
Felsenstein model for instance), which is what the program
[BAli-Phy](http://bali-phy.org/) does. However this is computationally very,
very expensive.

In the following we will rely on our MAFFT-based MSAs, but it is important to keep 
in mind that in all the following statistical phylogenetic analyses, we *are making
the terrible assumption that the MSA is known without error.*

## Phylogenetic inference using Maximum likelihood
