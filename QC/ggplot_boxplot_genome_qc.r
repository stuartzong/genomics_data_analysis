library(ggplot2)
library(plyr)
library(reshape2)


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
  library(grid)
  
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

setwd('/home/szong/projects/QC/genome/HIV_Cervical')
file <- '/home/szong/projects/QC/genome/HIV_Cervical/new_bamstats_summary.csv'
df <- read.table(file, header=TRUE, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df
#factor(status)
colnames(df)

png <- paste0(file, ".png")
#ggplot(df, aes(factor("coverage", "duplicate.rate", "alignment.rate"), coverage)) +
p1 <- ggplot(df, aes(factor(status), coverage)) +
  geom_boxplot(aes(fill = factor(status))) +
  geom_jitter(width = 0.1, height = 0.1) 
  

p2 <- ggplot(df, aes(factor(status), alignment_rate)) +
  geom_boxplot(aes(fill = factor(status))) +
  geom_jitter(width = 0.1, height = 0.1) 


p3 <- ggplot(df, aes(factor(status), duplicate_rate)) +
  geom_boxplot(aes(fill = factor(status))) +
  geom_jitter(width = 0.1, height = 0.1) 

category <- colnames(df)
p4 <- ggplot(df, aes(patient)) +
        geom_point(aes(y=coverage), size=2.5, shape=21, fill="white") +
        geom_point(aes(y=alignment_rate), size=2.5, shape=21, fill="white") +
        geom_point(aes(y=duplicate_rate), size=2.5, shape=21, fill="white")

# bring your data to long format as needed by ggplot
keeps <- c("patient", "coverage", "alignment_rate", "duplicate_rate")
df.keeps <- df[keeps]

df.melt <- melt(df.keeps, id="patient")
p5 <- ggplot(data=df.melt,
       aes(x=patient, y=value, colour=variable)) +
       geom_point(aes(shape = variable), size=4, fill="white") +
       
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=6, color="black"),
    axis.text.y=element_text(face="bold"),
    axis.title.x=element_text(face="bold"))


multiplot(p1,p2,p3,p5, cols=2)






