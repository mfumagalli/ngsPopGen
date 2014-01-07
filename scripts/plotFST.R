
# Usage: Rscript -i infile.fst -o outfile.eps -w win -s step

library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
make_option(c('-w','--window'), action='store', type='character', default=1, help='Window length'),
make_option(c('-s','--step'), action='store', type='character', default=1, help='Step size'),
make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# Read input file
values <- read.table(opt$in_file, stringsAsFact=F);
cat("Overall FST:",sum(values[,1])/sum(values[,2]),"\n");

# Windows
len=nrow(values);
win=as.numeric(opt$window);
step=as.numeric(opt$step);
start=seq(1, len-win, step);
end=start+win; # it requires that all windows have the same nr of sites, if less the window (usually at the end) is discarded
pos=start+win;
fst=c(); for (i in 1:length(start)) fst[i]=sum(values[start[i]:end[i],1])/sum(values[start[i]:end[i],2]);

# Data
df=data.frame(cbind(Pop=rep(1,length(pos)), Pos=pos, Value=fst));
df[,2:3]=sapply(df[,2:3], as.character)
df[,2:3]=sapply(df[,2:3], as.numeric)

# Plot
title = expression(F[ST]);
ggplot(data=df, aes(x=Pos, y=Value)) + geom_line() + ggtitle(title);
ggsave(opt$out_file)
unlink("Rplots.pdf", force=TRUE)


