
# R scripts to plot results from ngsPopGen program

## Plot FST values

     Rscript plotFST.R -i infile.fst -o outfile.eps -w win -s step      

where `infile.fst` is the output of `ngsFST`, `win` is the window size and `step` is the step for the sliding windows. Please note that this scripts assumes no missing data.

## Plot summary statistics

     Rscript plotSS.R -i infile.stat -o outfile.eps -n name1-name2 # if 2 pops     
     Rscript plotSS.R -i infile.stat -o outfile.eps -n name1 # if only 1 pop     

where `infile.stat` is the output of `ngsStat`.

## Plot PCA

     Rscript plotPCA.R -i infile.covar -c component1-component2 -a annotation.file -o outfile.eps     

where `infile.covar` is the ouput of `ngsCovar`, `component1-component2` are the two components that user wants to plot, `annotation.file` is a [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) cluster file.

## Plot 2D-SFS

     Rscript plot2dSFS.R infile.spec outfile.eps pop1 pop2    

where `infile.spec` is the output of `ngs2dSFS`, and `pop1` and `pop2` are the populations labels.




