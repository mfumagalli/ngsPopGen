
args <- commandArgs(T)

if (length(args)==0) {
        stop("No arguments given.
        Aim: Estimate NCD statistics (Bittarello) by sampling from .saf files.
        Usage: Rscript NCD.R input.saf.idx.gz nr_chroms target_freq
	Output is in stdout: nr_sites1 \t NCD1 \t nr_sites2 \t NCD2
        Info: this program take unfolded .saf.idx.gz data.\n")
}

fin <- args[1] # saf.idx file, this must be unfolded
nchroms <- as.numeric(args[2]) # how many chroms
tf <- as.numeric(args[3]) # target balanced frequency
type <- as.numeric(args[4]) # NCD1 NCD2
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
        aaf <- c(aaf, (sapply(ss, expsam)-1))
}
close(con)

# fixed differences
fd <- which(aaf==nchroms)

# snps
snps <- which(aaf>0 & aaf<nchroms)

# maf
aaf <- aaf/nchroms
aaf[which(aaf>0.5)] <- 1 - aaf[which(aaf>0.5)]

# NCD1
ind1 <- snps
ncd1 <- sqrt(mean((aaf[ind1]-tf)^2))

# NCD2
ind2 <- sort(unique(snps,fd))
ncd2 <- sqrt(mean((aaf[ind2]-tf)^2))

cat(length(ind1),"\t",ncd1,"\t",length(ind2),"\t",ncd2,"\n")



