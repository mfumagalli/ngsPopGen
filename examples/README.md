# Posterior probabilities of genotypes and sample allele frequencies
We use ANGSD to compute genotype posterior probabilities:

    $ANGSD/angsd -sim1 $SIM_DATA/testA.glf.gz -nInd 24 -doGeno 32 -doPost 1 -doMaf 2 -doSaf 1 -out testA -doMajorMinor 1

which will generate these files:

* `testA.arg`: parameters used;
* `testA.mafs.gz`: estimates of minor allele frequencies;
* `testA.geno.gz`: genotype posterior probabilities as doubles in binary format;
* `testA.saf`: log-likelihood ratios of sample allele frequencies

In case we prefer to weight each site rather than calling SNPs, we need to calculate posterior probabilities of sample allele frequencies:

    $ANGSD/misc/emOptim2 testA.saf 48 -nSites 10000 > testA.saf.ml
    $ANGSD/angsd -sim1 $SIM_DATA/testA.glf.gz -nInd 24 -doSaf 1 -pest testA.saf.ml -out testA.rf

At this point we can look at the estimated and true pooled site frequency spectrum:

    Rscript --vanilla --slave -e 'barplot(rbind(as.numeric(scan("../../ngsSim/examples/testA.frq", what="char")), exp(as.numeric(scan("testA.saf.ml", what="char")))), beside=T, legend=c("True","Estimated"))'

# Principal Component Analysis (PCA)
## Covariance matrix
We use ngsCovar to estimate a covariance matrix between pairs of individuals. This can be achieved in different ways. For low coverage sequencing data, we recommend to use `-norm 0` option which disables normalization proposed in [Patterson et al. (2006)](http://www.ncbi.nlm.nih.gov/pubmed/17194218).
The first way is to compute an approximation of the posterior of the covariance matrix, by weighting each site by its probability of being variable, as proposed in [Fumagalli et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23979584):

    ../ngsCovar -probfile testA.geno -outfile testA.covar1 -nind 24 -nsites 10000 -call 0 -sfsfile testA.rf.saf -norm 0

Alternatively, one can remove non-variable or low-frequency sites with the option `-minmaf`:

    ../ngsCovar -probfile testA.geno -outfile testA.covar2 -nind 24 -nsites 10000 -call 0 -minmaf 0.05

Finally, in case of high sequencing depth, one can call genotypes as the genotype with the highest posterior probability:

    ../ngsCovar -probfile testA.geno -outfile testA.covar3 -nind 24 -nsites 10000 -call 1 -minmaf 0.05

These commands will produce text files with a symmetric covariance matrix NxN (for N individuals).

## PCA plot

From the covariance matrix we can use a simple R script to perform an eigenvalue decomposition and plot the PCA results. Please make sure the use updated versions of required R packages. First, let's create a dummy plink cluster file.

    Rscript --vanilla --slave -e 'write.table(cbind(seq(1,24),rep(1,24),c(rep("A",10),rep("B",8),rep("C",6))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="testA.clst", quote=F)'

Assuming we want to plot the first 2 PCA components from the `testA.covar1` file, we can use the `plotPCA.R` script provided:

    Rscript --vanilla --slave $SCRIPTS/plotPCA.R -i testA.covar1 -c 1-2 -a testA.clst -o testA.pca.SAF.pdf

Please note that you need 'ggplot2' and 'optparse' R libraries installed. This script will output the explained genetic variance for each component and save as output the PCA plot.

# Summary Statistics

We use ngsStat to compute expectations of some basic population genetics statistics of the data, specifically the number of segregating sites, the expected heterozygosity, and the number of fixed differences between populations. We also want to compute these quantities in non-overlapping sliding windows of 100 sites each, by using the options `-iswin` and `-block_size`:

    ../ngsStat -npop 2 -postfiles testA1.rf.saf testA2.rf.saf -nsites 10000 -iswin 1 -nind 10 8 -outfile testA.stat -isfold 0 -islog 0 -block_size 100

The ouput has, for each window, start, end, number of variable sites in pop 1, expected heterozygosity in pop 1, number of variable sites in pop 2, expected heterozygosity in pop 2, number of fixed differences between populations. Values can be plotted using the script `plotSS.R` provided, either for both populations or only one:

    Rscript --vanilla --slave $SCRIPTS/plotSS.R -i testA.stat -o testA.stat.pdf -n pop1-pop2
    Rscript --vanilla --slave $SCRIPTS/plotSS.R -i testA.stat -o testA.stat.pop1.pdf -n pop1

An important note (especially for FST calculation) is that '-postfiles' of the 2 populations must correspond to the same exact sites. If they differ, you can use `GetSubSfs` program from [ngsUtils](https://github.com/mfumagalli/ngsUtils) to get an overlapping dataset.

# Population Differentiation
## Joint-SFS
To illustrate the estimation of _Fst_ between two populations, we will use the same simulated dataset as before. In case of low coverage sequencing, an improvement in the estimation accuracy can be achieved by using the joint-SFS (2D-SFS) as a prior for the sample allele frequency posterior distributions, as shown in [Fumagalli et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23979584):

    ../ngs2dSFS -postfiles testA1.rf.saf testA2.rf.saf -outfile testA.joint.spec -relative 1 -nind 10 8 -nsites 10000 -maxlike 1

Alternatively, and as proposed in the original formulation of the method, `testA1.saf` and `testA2.saf` can be used directly with `-islog 1`. For medium to high coverage data and in case of a large number of sites, another solution would be to set `-maxlike 0`.
To plot the joint SFS on a .pdf file, you can use the `plot2dSFS.R` provided script:

    Rscript --vanilla --slave $SCRIPTS/plot2dSFS.R testA.joint.spec testA.joint.spec.pdf pop1 pop2

If needed, this estimated 2D-SFS can be converted, using the script `convert.2Dsfs.to.dadi.R` provided, into `dadi` format for demographic inferences ([Gutenkunst et al. 2009)](http://www.ncbi.nlm.nih.gov/pubmed/19851460).

## _Fst_
Then we calculate a method-of-moments estimator of _Fst_, at each site, with the following command:

    ../ngsFST -postfiles testA1.saf testA2.saf -priorfile testA.joint.spec -nind 10 8 -nsites 10000 -outfile testA.fst -islog 1

Note that here we use `.saf` instead of `.rf.saf`, since we are using `.joint.spec` as a prior. In case the latter was not available (e.g. folded data or inbreeding, see below for these examples), we can run these equivalent commands:

    ../ngsFST -postfiles testA1.saf testA2.saf -priorfiles testA1.saf.ml testA2.saf.ml -nind 10 8 -nsites 10000 -outfile testA.fst2 -islog 1
    ../ngsFST -postfiles testA1.rf.saf testA2.rf.saf -nind 10 8 -nsites 10000 -outfile testA.fst3 -islog 0

We can calculate and print the overall _Fst_, as well as plot _Fst_ in sliding windows, using the provided `plotFST.R` script:

    Rscript --vanilla --slave $SCRIPTS/plotFST.R -i testA.fst -o testA.fst.pdf -w 100 -s 50
