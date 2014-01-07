
# Usage: Rscript -i infile.stat -o outfile.eps -n name1-name2
# Usage: Rscript -i infile.stat -o outfile.eps -n name1 # if only 1 pop

library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
                    make_option(c('-n','--names'), action='store', type='character', default=1-2, help='Name(s) of population(s)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# How many pops
pops <- as.character(strsplit(opt$names, "-", fixed=TRUE)[[1]]);
npop=length(pops);

# Read input file
values <- read.table(opt$in_file, stringsAsFact=F);
pos=as.numeric(values[,1]+(values[,2]-values[,1])/2);

# Plot

if (npop==2) {
  title <- "";

  # Data
  df = data.frame(cbind( Pop=c(rep(pops[1],length(pos)),rep(pops[2],length(pos))), Pos=pos, Segr.sites=c(values[,3], values[,5]), Exp.heterozygosity=c(values[,4],values[,6]), Fixed.differences=(rep(values[,7],2)), Pops=c(rep(paste(pops[1],"\n",pops[2]),length(pos)*2)) ) );
  df[,2:5] = sapply(df[,2:5], as.character)
  df[,2:5] = sapply(df[,2:5], as.numeric)

  p1 = ggplot(data=df, aes(x=Pos, y=Segr.sites, color=Pop)) + geom_line() + ggtitle(title)
  p2 = ggplot(data=df, aes(x=Pos, y=Exp.heterozygosity, color=Pop)) + geom_line() + ggtitle(title)
  p3 = ggplot(data=df, aes(x=Pos, y=Fixed.differences, color=Pops)) + geom_line() + ggtitle(title)

  pdf(opt$out_file);
  multiplot(p1, p2, p3, ncol=1)
  null <- dev.off();
}

if (npop==1) {
  df = data.frame(cbind( Pop=rep(pops[1],length(pos)), Pos=pos, Segr.sites=values[,3], Exp.heterozygosity=values[,4]));
  df[,2:4] = sapply(df[,2:4], as.character)
  df[,2:4] = sapply(df[,2:4], as.numeric)

  title <- "";

  p1 = ggplot(data=df, aes(x=Pos, y=Segr.sites, color=Pop)) + geom_line() + ggtitle(title)
  p2 = ggplot(data=df, aes(x=Pos, y=Exp.heterozygosity, color=Pop)) + geom_line() + ggtitle(title)

  pdf(opt$out_file);
  multiplot(p1, p2, ncol=1)
  null <- dev.off();
}


