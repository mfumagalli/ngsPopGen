
# ngsPopGen

Several tools to perform population genetic analyses from NGS data:
 * ` ngsFst`  - Quantificate population genetic differentiation
 * ` ngsCovar`  - Population structure via PCA (principal components analysis)
 * ` ngs2dSFS`  - Estimate 2D-SFS from posterior probabilities of sample allele frequencies
 * ` ngsStat`  - Estimates number of segregating sites, expected average heterozygosity, and number of fixed differences and Dxy (if 2 populations provided).

IMPORTANT NOTE i): 

In all analisis involving 2 populations, input data must refer to the exact same sites. If they differ (e.g. because of different filtering), you must first get the overlapping subset of sites for both populations.

To achieve this, you can follow instructions given in the tutorial ([here](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)).

A quick trick to check that everything went fine is to retrieve the dimension of each file (using "ls -l"), then divide this number by 8 and then by the double of the number of individuals plus 1. You should get the final number of sites.


IMPORTANT NOTE ii): 

The use of folded data (spectrum or sample allele frequencies probabilities) is no longer supported.
In case you do not have a reliable ancestral information, please use your reference sequence to polarise your data and follow all the steps as documented here.
However, do not attempt to make any inference based on the resulting unfolded reference/non-reference based - site frequency spectrum.


IMPORTANT NOTE iii):

It may be practical to perform a non-stringent SNP calling before running the following analyses, in order to reduce the computational time and data dimensions.
Moreover, this will reduce noise due to monomorphic sites, especially when the species' polymorphic rate is very low.


### Installation

    % git clone https://github.com/mfumagalli/ngsPopGen.git

To install these tools just run:

    % cd ngsPopGen
    % make

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

---

## ngsFST

Program to estimate FST from NGS data. It computes expected genetic variance components and estimate per-site FST from those using methods-of-moments estimator. See Fumagalli et al. Genetics 2013 for more details.
In input it receives sample allele frequencies likelihoods for each population and a 2D-SFS as a prior.

The output is a tab-separated text file. Each row represents a site. Columns are ordered as: A, AB, f, FST, Pvar; where A is the expectation of genetic variance between populations, AB is the expectation of the total genetic variance, f is the correcting factor for the ratio of expectations, FST is the per-site FST value, Pvar is the probability for the site of being variable.

### Usage
#### Examples:

 * using a 2D-SFS as a prior, estimated using ngs2dSFS:

#

    % ./ngsFST -postfiles pop1.saf pop2.saf -priorfile spectrum2D.txt -nind 20 20 -nsites 100000 -outfile pops.fst -verbose 0

#### Parameters:

    -postfiles: files with sample allele frequencies likelihoods for each population
    -priorfile: 2D-SFS to be used as a prior; you can use ngs2dSfs with parameter -relative set to 1
    -nind: number of individuals for each population
    -nsites: total number of sites; in case of a site subset this is the upper limit
    -firstbase: in case of a site subset, this is the lower limit
    -islog: boolean, are postfiles in -log? kept for compatibility but always set it to 1
    -outfile: name of the output file
    -block_size: number of sites in each chunk (for memory reasons, increase it if you can use more RAM)
    -verbose: level of verbosity

---

## ngsCovar

Program to compute the expected correlation matrix between individuals from genotype posterior probabilities. It receives as input genotype posterior probabilities. It can receive in input also posterior probabilities of sample allele frequencies for computing the probability of each site to be variant.

### Usage

#### Examples:

* not calling genotypes and weighting by each site's probability of being variable (recommended if no SNP calling is performed and depth is extremely low, otherwise please use other methods below):

#

    % ./ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 0 -norm 0 -sfsfile pop.sfs.ml.norm
    
* not calling genotypes but with SNP calling (preferred way under most circumstances):

#

    % ./ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 0 -minmaf 0.05
    
* calling genotypes (this is kept for compatibility but should not be used unless you have high-depth data, > 20X):

#

    % ./ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 1


#### Parameters:

    -probfile: file with genotype posterior probabilities
    -sfsfile: file with sample allele frequency posterior probabilities
    -nind: number of individuals
    -nsites: total number of sites; in case of a site subset this is the upper limit
    -offset: in case of a site subset, this is the lower limit
    -genoquality: text file with 'nsites' lines stating whether to use (1) or ignore (0) the site
    -norm: 0 [no normalization, recommended if no SNP calling is performed]
           1 [matrix is normalized by p(1-p) as in Patterson et al (2006)]
           2 [normalized by 2p(1-p)]
    -minmaf: ignore sites below this threshold of minor allele frequency
    -call: call genotypes based on the maximum posterior probability
    -islog: boolean, are postfiles in -log?
    -outfile: name of output file
    -block_size: number of sites in each chunk (for memory reasons)
    -verbose: level of verbosity

---

## ngs2dSFS

Program to estimate 2D-SFS from posterior probabilities of sample allele frequencies. Output file reports the occurrence of sites at distinct joint allele frequencies. This spectrum is a (2N1+1)x(2N2+2) matrix with N1 and N2 number of individuals at the two populations. Please note that cells are zero-based ordered. As an example, value reported in the cell [4,3] represents the frequency of sites with allele frequency 3 and 2 at population 1 and 2 respectively.

#### Example:

    % ./ngs2dSFS -postfiles pop.saf pop.saf -outfile spectrum.txt -relative 1 -nind 20 20 -nsites 100000

#### Parameters:

    -postfiles: file with sample allele frequency posterior probabilities (or likelihoods) for each population
    -nind: number of individuals
    -nsites: total number of sites; in case of a site subset this is the upper limit
    -offset: in case of a site subset, this is the lower limit
    -outfile: name of output file
    -maxlike: compute the MLE as the sum across sites' joint allele frequency (1, preferred) or as the sum of the products of likelihoods (0)
    -relative: boolean, whether input are absolute counts of sites with a specific joint allele frequency (0) or relative frequencies (1)
    -block_size: memory efficiency, number of sites for each chunk

ANGSD can compute a ML estimate of the 2D-SFS which should be preferred when many sites are available. However, ANGSD output file should be transformed (from log to un-log and from space-separated to tab-separated) before being used in ngsFST.


---

## ngsStat

Program to compute estimates of the number of segregating sites, the expected average heterozygosity, and the number of fixed differences (if 2 populations data is provided). It receives as input sample allele frequency posterior probabilities (from ANGSD) from 1 or 2 populations. Output is a text file with columns: start, end, segregating sites (pop 1), heterozygosity (pop 1), segregating sites (pop 2), heterozygosity (pop 2), fixed differences, dxy.

#### Example:

* 2 populations, sliding windows of 100 sites (latter recommended only if no missing data is present, however in most cases you will have some missing sites so this command should not be used):

#

    ./ngsStat -npop 2 -postfiles pop1.saf pop2.saf -nsites 1000 -iswin 1 -nind 10 10 -outfile pops.stat -verbose 0 -block_size 100
    
* 1 populations, sliding windows of 100 sites (latter recommended only if no missing data is present, so again in many cases this should not be used):

# 

    ./ngsStat -npop 1 -postfiles pop1.saf -nsites 1000 -iswin 1 -nind 10 -outfile pops.stat -block_size 100
    
* 1 population, values estimated at each site (recommended in case of missing data, and then the computation of values in sliding windows will be performed using the R script provided):

#

    ./ngsStat -npop 1 -postfiles pop1.saf -nsites 1000 -iswin 0 -nind 10 -outfile pops.stat

#### Parameters:

    -npop: number of populations (should be the first one to be specified)
    -postfiles: file with sample allele frequency posterior probabilities for each population
    -nind: number of individuals
    -nsites: total number of sites; in case of a site subset this is the upper limit
    -firstbase: in case of a site subset, this is the lower limit
    -islog: boolean, are postfiles in -log? keep 1
    -iswin: if set to 1, chuncks are considered non-overlapping sliding-windows
    -outfile: name of output file
    -block_size: number of sites in each chunk (for memory reasons)
    -verbose: level of verbosity


Further examples can be found [here](https://github.com/mfumagalli/ngsPopGen/tree/master/examples).




