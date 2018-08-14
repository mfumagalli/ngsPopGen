
print(date())

args <- commandArgs(T)

if (length(args)==0) {
	stop("No arguments given.
	Aim: Estimate the 2DSFS by sampling from .saf files.
	Usage: Rscript fast2DSFS.R input1.saf.idx.gz input1.saf.idx.gz nr_chroms1 nr_chroms2 prob_sampling_site 
	Info: Final sfs is in last argument, nr_chroms is the number of chromosomes, prob_sampling_site is the chance of sampling the current site for the estimation. This program take unfolded .saf.idx.gz data.\n")
}

fin1 <- args[1] # saf.idx file, this must be unfolded
fin2 <- args[2]
nchroms1 <- as.numeric(args[3]) # how many chroms
nchroms2 <- as.numeric(args[4]) 
psam <- as.numeric(args[5])
fout <- args[6]
rm(args)

nl <- 1000000

selpaste <- function(arr) {
	paste(arr[1:2], sep="", collapse=":")
}

expsam <- function(arr) {
	probs <- exp(as.numeric(arr[c(-1,-2)]))
	sample(x=seq(1,length(probs)), size=1, prob=probs)
}

# first file
pos <- saf <- c()

cat("Reading first file...\n")
con <- file(fin1, "r")
while (TRUE) {
	mline <- readLines(con, n=nl)
	if (length(mline)==0) {
      		break
    	}
	ss <- strsplit(mline, split="\t")
	pos <- c(pos, sapply(ss, selpaste))
	saf <- c(saf, sapply(ss, expsam))
	cat(pos[length(pos)],"\t")
}
close(con)

sfs <- matrix(0, nrow=(nchroms1+1), ncol=(nchroms2+1))

# second file
cat("\nReading second file...\n")
con <- file(fin2, "r")
while (TRUE) {
        mline <- readLines(con, n=nl)
        if (length(mline)==0) {
                break
        }
	ss <- strsplit(mline, split="\t")
	pos2 <- sapply(ss, selpaste)

	ind <- match(pos2, pos)

	take <- (!is.na(ind))*(sample(c(0,1),length(pos2),T,c(1-psam,psam)))

	saf2 <- sapply(ss, expsam)[which(take==1)]
	
	saf1 <- saf[ind[which(take==1)]]

	for (j in 1:length(saf1)) sfs[saf1[j], saf2[j]] <- sfs[saf1[j], saf2[j]] + 1

	cat(pos2[length(pos2)], "\t")

}
close(con)

final <- array(t(sfs))

cat("\nTotal sites:", sum(final))
cat("\nProportion (on pop1):", sum(final)/length(pos))

cat("\nNon-normalised SFS:\n")
cat(final)

final <- final/sum(final)
cat("\nNormalised SFS:\n")
cat(final)

cat(final, file=fout)

cat("\n")
print(date())








