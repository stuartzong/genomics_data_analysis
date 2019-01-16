library(ggplot2)
library(plyr)
library(reshape2)

my.path <- getwd()
files <- list.files(my.path, pattern = "bamstats_summary_transposed.csv.tmp$")
files
files.len <- length(files)

for (file in files){
  print(file)
  png <- paste0(file, ".png")
  #pdf(pdf, height=8, width=8)

#read in data
# column 1 is the row name
data <- file
df2 <- read.table(data, header=TRUE, row.names=1, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df2
print(df2)

df2$Category  <- row.names(df2)

# bring your data to long format as needed by ggplot
df2.molten <- melt(df2, value.name="QCmatrix", variable.name="patient", na.rm=TRUE)
df2.molten.sorted <- arrange(df2.molten, QCmatrix , Category) 

print("long format of the data")
print(df2.molten.sorted)
mut <- ggplot(data=df2.molten.sorted, aes(x=patient, y=QCmatrix, group = Category, colour = Category)) +
    geom_line() +
    geom_point( size=2.5, shape=21, fill="white")+
    ggtitle("HIV Cervical Genome QC Matrix")+
    #coord_cartesian(xlim = c(0, 100)+
    ## theme(axis.title.x = element_text(size = rel(3.8)))
    theme(
          axis.text.x = element_text(angle = 90, hjust = 1, size=4, color="black"),
          axis.text.y=element_text(face="bold"),
          axis.title.x=element_text(face="bold"))

ggsave(mut,file=png,width=10,height=8)
## dev.off()

}


