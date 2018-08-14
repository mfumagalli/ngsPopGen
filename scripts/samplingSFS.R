
args <- commandArgs(T)

if (length(args)==0) {
	stop("No arguments given.\nAim: Estimate the SFS by sampling from .saf files. Fold it at the end if requested.\nUsage: Rscript fastSFS.R input.saf.idx.gz nr_chroms fold_it? plot.pdf\nInfo: Final sfs is in stdout and fold_it? is boolean (0 or 1) and nr_chroms is the number of chromosomes. This program take unfolded .saf.idx.gz data.\n")
}

fin <- args[1] # saf.idx file, this must be unfolded
nchroms <- as.numeric(args[2]) # how many chroms
tofold <- as.numeric(args[3]) # fold it at the end?
fout <- args[4] # plot in pdf
rm(args)

# output is in stdout

saf <- rep(0, nchroms+1)

con <- file(fin, "r")

while (TRUE) {
	line <- readLines(con, n=1)
	if (length(line)==0) {
      		break
    	}
	probs <- exp(as.numeric(strsplit(line, split="\t")[[1]][c(-1,-2)]))
	probs <- probs/sum(probs) # perhaps not necessary
	if (length(probs)==length(saf)) {
		sampled <- sample(x=seq(1,length(probs)), size=1, prob=probs)
		saf[sampled] <- saf[sampled]+1
	} else {
		break("Error! Quit!")
	}
}

close(con)

fold <- function(spec) {
	nt <- length(spec)
	ns <- (length(spec)-1)/2
	ff <- rep(NA, ns+1)
	for (i in 1:(length(ff))) {
		ff[i] <- spec[i]+spec[nt-i+1]
	}
	ff[length(ff)] <- spec[length(ff)]
	#print(sum(spec)==sum(ff))
	ff
}

# fold the final sfs if required
sfs <- saf
if (tofold) sfs <- fold(saf)

pdf(file=fout)

if (tofold) {
	pvar <- sfs[1]/sum(sfs)
	barplot(sfs[-1], names.arg=1:(length(sfs)-1), sub=paste("variability: ", pvar))
} else {
	pvar <- (sfs[1]+sfs[length(sfs)])/sum(sfs)
        barplot(sfs[2:(length(sfs)-1)], names.arg=1:(length(sfs)-2), sub=paste("variability: ", pvar))
}

invisible(dev.off())

cat(sfs,"\n")



















