
## Script provided by Dean Ousby and adapted

library(plot3D)

args <- commandArgs(T)
fin <- args[1]
pops <- unlist(strsplit(args[2],split="-")) # names
nchroms <- as.numeric(unlist(strsplit(args[3], split="-")))
filtmono <- as.numeric(args[4])
rm(args)

val <- scan(paste(fin, sep=""), quiet=T)
sfs <- matrix(val, nrow=nchroms[1]+1, ncol=nchroms[2]+1, byrow=T)
rm(val)

# mask non-variant sites
if (filtmono) {
	sfs[1,1] <- 0
	sfs[nrow(sfs),ncol(sfs)] <- 0
}

fout <- paste(fin,".pdf", sep="", collapse="")

pdf(file=fout)

hist3D(x = seq(0,1,length.out = nrow(sfs)), y = seq(0,1,length.out=ncol(sfs)), sfs, cex.lab=1.2, xlab=pops[1], ylab=pops[2], zlab="Frequency", main=paste("2D-SFS"),pin=c(10,0), cex.main=1.4, zlim=c(0,max(sfs)))

invisible(dev.off())



