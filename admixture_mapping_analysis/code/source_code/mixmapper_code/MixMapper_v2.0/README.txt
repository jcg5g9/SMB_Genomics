###################################
###                             ###
###      MixMapper,  v2.0       ###
###                             ###
###  Mark Lipson and Po-Ru Loh  ###
###        May 29, 2014         ###
###                             ###
###################################



Table of Contents
-----------------
  1. Overview
  2. Installation
  3. Usage
  4. Change Log
  5. License
  6. Contact


==== 1. Overview ====

The MixMapper software performs admixture inference from allele
frequency moment statistics as described in:

  Original MixMapper paper:

    Lipson M, Loh P-R, Levin A, Reich D, Patterson N, and Berger B.
    Efficient Moment-Based Inference of Population Admixture
    Parameters and Sources of Gene Flow.
    Molecular Biology and Evolution, 2013
    http://mbe.oxfordjournals.org/content/30/8/1788.full

  MixMapper 2.0:

    Lipson M, Loh P-R, Patterson N, Moorjani P, Ko Y-C, Stoneking M,
    Berger B, and Reich D.
    Reconstructing Austronesian population history in Island Southeast
    Asia.
    bioRxiv, 2014 (under revision)
    http://www.biorxiv.org/content/early/2014/05/27/005603

MixMapper can be thought of as a generalization of the qpgraph
software of Patterson et al. (Genetics, 2012), which takes as input
genotype data, along with a proposed arrangement of admixed and
unadmixed populations, and returns branch lengths and mixture
fractions that produce the best fit to allele frequency moment
statistics measured on the data.  MixMapper, by contrast, performs the
fitting in two stages, first constructing an unadmixed scaffold tree
via neighbor-joining and then automatically optimizing the placement
of admixed populations onto this initial tree.  Thus, no topological
relationships among populations need to be specified in advance.

MixMapper is also similar in spirit to the independently developed
TreeMix method of Pickrell et al. (PLoS Genetics, 2012).  Like
MixMapper, TreeMix builds admixture trees from second moments of
allele frequency divergences, although it does so via a composite
likelihood maximization approach made tractable with a multivariate
normal approximation.  Procedurally, TreeMix is structured in a
"top-down" fashion, whereby a full set of populations is initially fit
as an unadmixed tree, and gene flow edges are added sequentially to
account for the greatest errors in the fit.  This format makes TreeMix
well-suited to handling very large trees: the entire fitting process
is automated and can include arbitrarily many admixture events
simultaneously.  In contrast, MixMapper is designed as an interactive
tool to maximize flexibility and precision with a "bottom-up"
approach, beginning with a carefully screened unadmixed scaffold tree
to which admixed populations are added with best-fitting parameter
solutions.


==== 2. Installation ====

Only two standalone C++ source files need compiling, which can be done
using the following commands in a Linux environment:

g++ -O2 compute_moment_stats.cpp -o compute_moment_stats
g++ -O2 compute_most_additive_trees.cpp -o compute_most_additive_trees

The MATLAB code for mixture fitting requires a MATLAB installation
with the Bioinformatics Toolbox and Optimization Toolbox.  If the
Parallel Computing Toolbox is also installed, fitting of bootstrap
replicates can be performed in parallel.  We have tested the code on
MATLAB versions 7.14 (R2012a) and 8.0 (R2012b).

If you wish to run MixMapper but do not have access to MATLAB, we also
provide a compiled executable that can be run with an install of the
(free) MATLAB Compiler Runtime. (The MCR version is more cumbersome to
use, however.)  This package is available as a separate download at:

  http://groups.csail.mit.edu/cb/mixmapper/


==== 3. Usage ====

The MixMapper package contains two programs written in C++ that
compute allele frequency moment statistics and assist the user in
selecting a suitable scaffold tree for mixture fitting; the mixture
fitting functionality is implemented in MATLAB.  These utilities are
explained below, and an example of the full workflow is provided in
the example/ subdirectory of this package (see example/readme.txt).

---- A. compute_moment_stats ----

Computes the allele frequency moment statistics needed for MixMapper
analysis.

Input command line arguments:

- arg1 = .ind file (data in EIGENSTRAT format)
- arg2 = .snp file
- arg3 = .geno file
- arg4 = output prefix
    Output files have the specified prefix and the extensions below:
      .pops.txt, .f2_boots.txt, .h_boots.txt, .f2.tab, .neg_f3.txt
- arg5 = resample individuals in each population? (y/n)
    Choosing 'y' takes into account variability in the choice of
    samples used to represent a population.
- arg6 = number of bootstrap replicates
- arg7 = number of SNP blocks (e.g., 50)
    Bootstrapping is performed over blocks of SNPs to account for LD.

Output files:

- .pops.txt
  List of N population names.

- .f2_boots.txt
  All pairwise f2 values (N x N) computed on R replicates.  The file
  contains N x N lines, looping through populations in the same order
  as in the .pops.txt file.  Each line contains R space-separated
  numbers, one per replicate.

- .h_boots.txt
  All heterozygosity values (N) computed on R replicates.  The file
  contains N lines in the same order as the .pops.txt file.  Each line
  contains R space-separated numbers, one per replicate.

- .f2.tab:
  Tab-delimited matrix of f2 values computed without resampling (i.e.,
  on full data set); contains row and column headers.

- .neg_f3.txt:  
  Table of population triples (C,A,B) with negative values of
  f3(C;A,B) indicating admixture in population C.  Standard errors and
  z-scores are estimated from bootstrap replicates.

Output to stdout:

Filtered list of populations that the f3 test does not detect as
admixed.  This is a reasonable starting list from which to select a
scaffold subset.

Output to stderr:

Progress of the computation and miscellaneous info about the data set
for verification.


---- B. compute_most_additive_trees ----

Computes neighbor-joining trees on subsets of a user-specified set of
putatively unadmixed populations.  Outputs population subsets with
induced tree distances that are closest (i.e., have smallest maximum
error) to the measured f2 distances.  These population subsets are
good candidates for use as scaffold trees in MixMapper.

Algorithmically, subsets are generated and analyzed in increasing
order of size.  To limit the search space, for each size n, only a
user-specified maximum number of best subsets are considered; the
subsets analyzed at the next iteration (n+1) are produced by
augmenting each of the current best subsets with each possible
population.

Input command line arguments:

- arg1 = .f2.tab file (from output of compute_moment_stats)
- arg2 = max number of subsets of each size to analyze (e.g., 10000)
- arg3 = file containing list of populations to choose from
- (optional) arg4 = file containing pops required to be in the tree

Output to stdout:

Lists of population subsets satisfying the given constraints admitting
neighbor-joining trees with the smallest deviations from measured f2
distances.


---- C. MixMapper (MATLAB function) ----

Creates a scaffold tree using neighbor joining and finds the best-fit
placement of one or two additional populations as two- or three-way
admixtures.  Combines results from bootstrap replicates to obtain
confidence intervals on parameters, which are output to the console in
a summary table.

Note that the MATLAB code for mixture fitting requires a MATLAB
installation with the Bioinformatics Toolbox and Optimization
Toolbox. If the Parallel Computing Toolbox is also installed, fitting
of bootstrap replicates can be performed in parallel. We have tested
the code on MATLAB versions 7.14 (R2012a) and 8.0 (R2012b).

For full documentation, please see the MixMapper.m file, in which
MATLAB-style documentation appears at the top.


==== 4. Change Log ====

Version 2.0 (May 29, 2014):

- Implemented new features for improved three-way admixture inference:
  - More accurate 3-way mixture fraction inference
  - Option for comparing goodness of fit of 2-way vs. 3-way models
  These extensions are described in the second reference above.

- Fixed minor bug that caused a crash in rare cases (reroot_dist < 0).

- Fixed error in example data (h_boots file, which had been generated
  using outdated code).

Version 1.02 (August 9, 2013):

- Compiled MixMapper into an executable now available for use with the
  MATLAB Compiler Runtime.

- Fixed bug in compute_most_additive_trees.cpp that caused crashes in
  some cases when not requiring any populations to be in the tree or
  when analyzing greater than 32 putatively unadmixed populations.
  The program now supports up to 64 unadmixed populations as intended.

- Fixed minor bug in MixMapper.m that caused error when
  options.branch_sets was not set.

- Updated compute_moment_stats.cpp to ignore any individual with
  "ignore" or "Ignore" in the population name.

Version 1.01 (April 11, 2013):

- Fixed bug in compute_moment_stats.cpp heterozygosity computation
  (missing factor of 2).  This bug caused errors in the conversion of
  branch lengths to drift units; it did not affect the results when
  displayed in the default f2 units and did not affect the fitting
  procedure in either case.

- Fixed minor bug affecting root placement in rare cases.

- Oriented MixMapper output so that when displaying a mixture fit
  between Branch 1 and Branch 2, the major component of ancestry
  always comes from Branch 1 (i.e., the mixture fraction alpha > 0.5).

- Added option to name branches according to the sets of populations
  they split the tree into (as an alternative to the current "trace"
  naming system) and added output fields set1, set2, set3.  Enable the
  'branch_sets' option to display output using this nomenclature.

- Added support for input in "packed geno" format (2 bits/genotype).


==== 5. License ====

This software is licensed for academic and non-profit use only.


==== 6. Contact ====

If you have comments or questions about this software, please visit
our website for additional documentation and our contact info:

  http://groups.csail.mit.edu/cb/mixmapper/

Future updates will also be made available at the above link.
