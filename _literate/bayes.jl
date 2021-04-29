# [back](/phylocourse/)
# \toc
#
# ## Bayesian phylogenetic inference
#
# A theoretical introduction will be given in-class. 
#
# ## MrBayes exercise
#
# Download and install
# [MrBayes](https://github.com/NBISweden/MrBayes/releases/tag/v3.2.7). Some
# more information can be found
# [here](http://nbisweden.github.io/MrBayes/download.html). Windows users
# should use the precompiled executables in the `*-WIN.zip` file. Macos users
# should be able to install it via Homebrew. Linux users should have no
# difficulty building from source as indicated in the
# [README](https://github.com/NBISweden/MrBayes).
#
# The [manual for
# MrBayes](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf)
# provides very helpful information.
#
# You will need a file in the `nexus` format to use MrBayes. Tip, if you
# have python installed, then the Biopython library provides convenient
# alignment format conversion tools:
#
# ```
# from Bio import AlignIO
# from Bio.Alphabet import DNAAlphabet
# AlignIO.convert("inputfile.fasta", "fasta", "outputfile.nex", "nexus", alphabet=DNAAlphabet())
# ```
#
# Anyway, [here](/assets/phylocourse/data/hiv-dentist.nex) is a nexus file for
# the HIV-dentist data set.
#
# Open up the MrBayes prompt by typing `mb` at the command line. You should see
# something like the following
#
# ```
# $ mb
#
#
#                            MrBayes 3.2.7a x86_64
#
#                      (Bayesian Analysis of Phylogeny)
#
#              Distributed under the GNU General Public License
#
#
#               Type "help" or "help <command>" for information
#                     on the commands that are available.
#
#                   Type "about" for authorship and general
#                       information about the program.
#
#
# MrBayes > 
# ```
#
# Now load the data file by typing `execute hiv-dentist.nex`. You can specify a
# substitution model using `lset`. For instance `lset nst=2 rates=gamma
# ngammacat=4` this will set the nucleotide substitution model (`nst`) to model
# 2 (which is the K2P model), assume Gamma distributed rates across sites and
# will approximate the Gamma distribution by four discrete rate categories.
# 
# >**Question:** To do Bayesian inference we will need to make prior
# >assumptions under the form of prior distributions for the parameters of our
# >model. In the basic model as we have set it up here (`K2P+G4`) which model
# >parameters are there? Type `help prset` and have a look at the default prior
# >settings that are currently selected. Most of the prior settings listed here 
# >are not relevant for the current model, which prior settings are relevant for
# >the currently specified model? You can check this by typing `showmodel`. 
#
# The default prior settings in MrBayes are well-thought of and designed to do
# a good job in many cases. They are more or less so-called 'uninformative' or
# 'vague' priors, so that the results are not, in general, very sensitive to
# them. Nevertheless, this **cannot be an excuse** to not be aware of the prior
# settings and parameters! If one does Bayesian inference, one does commit
# oneself to the full deal.
# 
# As an example of how we may wish to adjust the priors, consider the branch
# lengths prior. 
# 
# >**Question:** What is the default branch length prior distribution? Read about
# >it in the documentation available in `help prset`. Do you understand how this
# >how this prior is defined? Change the prior to something more simple, such
# >as an exponential distribution for each branch. Choose the parameter of the 
# >exponential prior based on your Maximum likelihood trees inferred with IQ-TREE.
#
# Now let's run the MCMC algorithm. Type `mcmc ngen=100000 samplefreq=100
# printfreq=1000 diagnfreq=10000`. 
#
# >**Question:** By default MrBayes will run two MCMC chains, why do we run a
# >second chain, what could it tell us that a single chain can't? What does the
# >standard deviation of split frequencies mean?
#
# If the standard deviation of split frequencies is low enough (below 0.01, or
# below 0.05 if you're inpatient) you can stop the run. After stopping the run
# you can type `sump` to summarize the posterior distribution for model
# parameters and `sumt` to summarize the posterior distribution for the trees.
# 
# >**Question:** Interpret the output of `sumt` printed to the screen, what would be
# >the posterior mean estimate for the transition/transversion ratio (`kappa`)?
# >What is its 95% uncertainty interval (and what does that mean?). What is the
# >ESS? Have we run the chain long enough?
#
# >**Question:** Now look at the `sumt` output, what do the clade credibility
# >values mean?
#
# To visually explore the samples from the posterior distributions, you can use
# the programs [`Tracer`](https://github.com/beast-dev/tracer/releases/tag/v1.6)
# and [`FigTree`](https://github.com/rambaut/figtree/releases/tag/v1.4.4). You
# can load the `*.t` output file(s) in FigTree and browse through the sampled
# trees. You can load the `*.p` output file(s) in Tracer to check the posterior 
# distribution (you may need to delete the first line from the `*.p` output file
# from MrBayes to load it in Tracer).
#
# >**Question:** Look at the trace plots for the different parameters. Did we
# >manage to reach the stationary distribution of the MCMC chain?


