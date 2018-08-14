
args <- commandArgs(T)

if (length(args)==0) {
        stop("No arguments given.
        Aim: Estimate NCD statistics (Bittarello) by sampling from .saf files.
        Usage: Rscript NCD.R input.saf.idx.gz nr_chroms target_freq 
	Output is in stdout: nr_sites \t NCD
        Info: this program take unfolded .saf.idx.gz data.\n")
}

fin <- args[1] # saf.idx file, this must be unfolded
nchroms <- as.numeric(args[2]) # how many chroms
tf <- as.numeric(args[3])
rm(args)

nl <- 1000000

selpaste <- function(arr) {
        paste(arr[1:2], sep="", collapse=":")
}

expsam <- function(arr) {
        probs <- exp(as.numeric(arr[c(-1,-2)]))
        sample(x=seq(1,length(probs)), size=1, prob=probs)
}

aaf <-c()
con <- file(fin, "r")
while (TRUE) {
        mline <- readLines(con, n=nl)
        if (length(mline)==0) {
                break
        }
        ss <- strsplit(mline, split="\t")
        aaf <- c(aaf, sapply(ss, expsam)/nchroms)
}
close(con)

write(length(aaf), stderr())

write(sqrt(mean((aaf-tf)^2)), stdout())



