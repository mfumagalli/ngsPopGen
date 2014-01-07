
# Usage: Rscript infile.spec outfile.eps pop1 pop2

# Read parameters
args <- commandArgs(trailingOnly = TRUE);
infile <- args[1];
outfile <- args[2];
pop1 <- args[3];
pop2 <- args[4];
rm(args);

# Read input file
values <- as.matrix(read.table(infile, stringsAsFact=F));
n1=nrow(values)-1;
n2=ncol(values)-1;

# Plot
pdf(outfile);
image(x=seq(0,n1), y=seq(0,n2), z=-log10(values), xlab=pop1, ylab=pop2, main="Joint SFS")
dev.off();

