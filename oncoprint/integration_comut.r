# how to make a Co-mutation (comut) plot with ggplot2
#library("ggplot2")
.libPaths('/gsc/software/linux-x86_64-centos6/R-3.3.2/lib64/R/library')
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
#not sure why genom-points plotted for size=0

# RNA
# file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/69_patients/transcriptome/new_integration.txt"
file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/82_patients/transcriptome/transcriptome_integration_summary.csv.added.formated.comut"
#dna
# DNA_file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/69_patients/genome/DNA_integration.csv.formated"
DNA_file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/82_patients/genome/genome_integration_summary.csv.added.formated.comut"
df <- read.table(file,
                 header=TRUE,
                 sep="\t",
                 quote="",
                 as.is=TRUE,
                 stringsAsFactors=FALSE)
df$patient

df2 <- read.table(DNA_file,
                 header=TRUE,
                 sep="\t",
                 quote="",
                 as.is=TRUE,
                 stringsAsFactors=FALSE)
#status <- df$status

# store color blind pallette
my.colors <- c("purple",
               "#E69F00",
               "#56B4E9",
               "#009E73",
               "#F0E442",
               "#0072B2",
               "#D55E00",
               "#CC79A7",
               "#E41A1C",
               "#377EB8" )

# all categories
categories <-c("Adeno",
               "Squamous",
               "Positive",
               "Negative",
               "Multiple",
               "NON_SYNONYMOUS_CODING",
               "SPLICE_SITE_ACCEPTOR",
               "STOP_GAINED",
               "integrated",
               "unintegrated")
#my.colors <- brewer.pal(5,"Set1")

#Create a custom color scale for each category
names(my.colors) <- categories
my.colors


# now for a Comut plot with ggplot2
# this tell ggplot do not order the patient, since it is an ordered factor already
df$patient <- factor(df$patient, levels = unique(df$patient))
df$ids <- factor(df$ids, levels = unique(df$ids))

# this is to manually order the legend, need to understand better
x <- factor(df$status)
x = factor(x,levels(x)[c(1,5, 4,3, 2,6)])

## assign refactored x to the data frame
df[,"status"] = x 



mut <- ggplot(df, aes(x=patient, y=ids, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=status)) +
#scale_fill_brewer(palette = "Set1", na.value="Grey90") +
scale_fill_manual(values=my.colors, na.value="Grey90") +
#scale_colour_manual("", values = c("aaa" = "red", "bbb" = "black", "ccc" = "blue", "ddd" = "green"),labels = c('Company D','Company A','Company C','Company B')) +
#col.scale +
#geom_point(data=df2, aes(x=patient, y=ids ), size=2*df2$status)+
geom_point(data=df2, aes(x=patient, y=ids ), size=df2$status)+
xlab("") +
ggtitle("HIV Cervical HPV integration") +
theme(
legend.key = element_rect(fill='NA'),
#legend.key.size = unit(0.4, 'cm'),
legend.title = element_blank(),
legend.position="bottom",
legend.text = element_text(size=8, face="bold"),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank(),
#axis.text.x=element_blank(),
axis.text.x=element_text(angle = 90, hjust = 1,colour="black"),
axis.text.y=element_text(colour="Black"),
axis.title.x=element_text(face="bold"),
axis.title.y=element_blank(),
panel.grid.major.x=element_blank(),
panel.grid.major.y=element_blank(),
panel.grid.minor.x=element_blank(),
panel.grid.minor.y=element_blank(),
panel.background=element_blank()
))
mut
print("lllllllllllllllllllll") 
ggsave(mut,file="integration_Comutplot.png",width=10,height=8)

#getwd() #get work directory


