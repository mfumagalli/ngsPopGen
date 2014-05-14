
# Usage: Rscript -i infile.fst -o outfile.eps -w win -s step

library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
make_option(c('-p','--pos_file'), action='store', type='character', default=NULL, help='Input position file'),
make_option(c('-w','--window'), action='store', type='character', default=1, help='Window length'),
make_option(c('-s','--step'), action='store', type='character', default=1, help='Step size'),
make_option(c('-t','--th'), action='store', type='character', default=0, help='Minimu probability of being variable'),
make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# Read input file
values <- read.table(opt$in_file, stringsAsFact=F);
ind <- which(values[,5]>=as.numeric(opt$th))
pos <- as.numeric(readLines(opt$pos_file))
if (length(pos)!=nrow(values)) stop("Dimensions of fst values and positions must match. Terminate.\n");
values <- values[ind,]
pos <- pos[ind]
cat("After removing sites, now there are",nrow(values),"sites going from",min(pos),"to",max(pos),"\n")
cat("Overall FST:",sum(values[,1])/sum(values[,2]),"\n");

# Windows
len=max(pos)
win=as.numeric(opt$window);
step=as.numeric(opt$step);
start=seq(min(pos), len, step);
end=start+win-1;
wpos=round(start+(win/2)); # position of the window in the plot (center)
fst=c(); 
for (i in 1:length(start)) {
	ipos=which(pos>=start[i] & pos<=end[i])
	fst[i]=sum(values[ipos,1])/sum(values[ipos,2]);
}

# Data
df=data.frame(cbind(Pop=rep(1,length(wpos)), Pos=wpos, Value=fst));
df[,2:3]=sapply(df[,2:3], as.character)
df[,2:3]=sapply(df[,2:3], as.numeric)
write.table(df, file=paste(opt$out_file,".txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=F)

# Plot
title = expression(F[ST]);
ggplot(data=df, aes(x=Pos, y=Value)) + geom_line() + ggtitle(title);
ggsave(paste(opt$out_file,".eps",sep="",collapse=""))
unlink("Rplots.pdf", force=TRUE)


