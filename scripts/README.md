
# R scripts to plot results from ngsPopGen program

## Plot FST values

     Rscript plotFST.R -i infile.fst -o outfile -p positions.txt -w win -s step -th min_prob_variable      

where:
 `infile.fst` is the output of `ngsFST`, 
 `positions.txt` is a single column new-line separated numbers of genomic position (on the same chromosome) for corresponding values of per-site FST,
 `win` is the window size,
 `step` is the step for the sliding windows,
 'th` is the minimum threshold on the probability of being variable for a site to be kept in the analysis.

This will generate 2 files, a plot `outfile.eps` and a text file `outfile.txt` with all values for each window.

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




