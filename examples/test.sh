SCRIPTS=../scripts
SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd





##### Genotypes' and sample allele frequencies' posterior probabilities
$ANGSD/angsd -sim1 $SIM_DATA/testA.glf.gz -nInd 24 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 1 -out testA
$ANGSD/misc/emOptim2 testA.saf 48 -nSites 10000 > testA.saf.ml
$ANGSD/angsd -sim1 $SIM_DATA/testA.glf.gz -nInd 24 -doSaf 1 -pest testA.saf.ml -out testA.rf

# Estimated and true pooled site frequency spectrum
#Rscript --vanilla --slave -e 'barplot(rbind(as.numeric(scan("../../ngsSim/examples/testA.frq", what="char")), exp(as.numeric(scan("testA.saf.ml", what="char")))), beside=T, legend=c("True","Estimated"))'





##### PCA
# Get covariance matrix
gunzip -f testA.geno.gz
../ngsCovar -probfile testA.geno -outfile testA.covar1 -nind 24 -nsites 10000 -call 0 -sfsfile testA.rf.saf -norm 0
../ngsCovar -probfile testA.geno -outfile testA.covar2 -nind 24 -nsites 10000 -call 0 -minmaf 0.05
../ngsCovar -probfile testA.geno -outfile testA.covar3 -nind 24 -nsites 10000 -call 1 -minmaf 0.05

# Plot results
#Rscript --vanilla --slave -e 'write.table(cbind(seq(1,24),rep(1,24),c(rep("A",10),rep("B",8),rep("C",6))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="testA.clst", quote=F)'
#Rscript --vanilla --slave $SCRIPTS/plotPCA.R -i testA.covar1 -c 1-2 -a testA.clst -o testA.pca.SAF.pdf
#Rscript --vanilla --slave $SCRIPTS/plotPCA.R -i testA.covar2 -c 1-2 -a testA.clst -o testA.pca.MAF.pdf
#Rscript --vanilla --slave $SCRIPTS/plotPCA.R -i testA.covar3 -c 1-2 -a testA.clst -o testA.pca.MAFcall.pdf





##### Statistics
# Pop 1
$ANGSD/angsd -sim1 $SIM_DATA/testA1.glf.gz -nInd 10 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 1 -out testA1
$ANGSD/misc/emOptim2 testA1.saf 20 -nSites 10000 > testA1.saf.ml
$ANGSD/angsd -sim1 $SIM_DATA/testA1.glf.gz -nInd 10 -doSaf 1 -pest testA1.saf.ml -out testA1.rf
# Pop 2
$ANGSD/angsd -sim1 $SIM_DATA/testA2.glf.gz -nInd 8 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 1 -out testA2
$ANGSD/misc/emOptim2 testA2.saf 16 -nSites 10000 > testA2.saf.ml
$ANGSD/angsd -sim1 $SIM_DATA/testA2.glf.gz -nInd 8 -doSaf 1 -pest testA2.saf.ml -out testA2.rf

# Get stats
../ngsStat -npop 2 -postfiles testA1.rf.saf testA2.rf.saf -nsites 10000 -iswin 1 -nind 10 8 -outfile testA.stat -isfold 0 -islog 0 -block_size 100

# Plot results
#Rscript --vanilla --slave $SCRIPTS/plotSS.R -i testA.stat -o testA.stat.pdf -n pop1-pop2
#Rscript --vanilla --slave $SCRIPTS/plotSS.R -i testA.stat -o testA.stat.pop1.pdf -n pop1





##### Fst
# 2D-SFS
../ngs2dSFS -postfiles testA1.rf.saf testA2.rf.saf -outfile testA.joint.spec -relative 1 -nind 10 8 -nsites 10000 -maxlike 1
#Rscript --vanilla --slave $SCRIPTS/plot2dSFS.R testA.joint.spec testA.joint.spec.pdf pop1 pop2

# Estimate Fst
../ngsFST -postfiles testA1.saf testA2.saf -priorfile testA.joint.spec -nind 10 8 -nsites 10000 -outfile testA.fst -islog 1
../ngsFST -postfiles testA1.saf testA2.saf -priorfiles testA1.saf.ml testA2.saf.ml -nind 10 8 -nsites 10000 -outfile testA.fst2 -islog 1
../ngsFST -postfiles testA1.rf.saf testA2.rf.saf -nind 10 8 -nsites 10000 -outfile testA.fst3 -islog 0

# Plot
#Rscript --vanilla --slave $SCRIPTS/plotFST.R -i testA.fst -o testA.fst.pdf -w 100 -s 50





##### Check MD5
rm -f *.arg
md5sum testA* | sort -k 2,2 > /tmp/test.md5
if diff /tmp/test.md5 test.md5 > /dev/null
then
    echo "ngsPopGen: All tests OK!"
else
    echo "ngsPopGen: test(s) failed!"
fi
