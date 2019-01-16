library("ggplot2")
library(ggplot2)
library(plyr)
library(reshape2)


# make a Co-mutation (comut) plot with ggplot2
## also add in the history and other pathology subtypes

file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/new_run/comut_data.txt"
df <- read.table(file,
                 header=TRUE,
                 sep="\t",
                 quote="",
                 as.is=TRUE,
                 stringsAsFactors=FALSE)
df
df$patient
df$clinic_or_mutations

# specifiy color for each categorial variables
# so that make multiple pot have the same color legend
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

## use brewer pal if no need to manually assign colors
#my.colors <- brewer.pal(5,"Set1")

## name colors using categories
names(my.colors) <- categories

## names(my.colors)


my.colors
#my.level <- levels(factor(df$clinic_or_mutations))
#my.level[c(1,7,5,3,4,8,6,2)]

# this is to manually order the legend, need to understand better
x <- factor(df$clinic_or_mutations)
x = factor(x,levels(x)[c(1,7,5,3,4,8,6,2)])

## assign refactored x to the data frame
df[,"clinic_or_mutations"] = x 
df$patient <- factor(df$patient,
                     levels=unique(df$patient))

df$clinic_or_genes <- factor(df$clinic_or_genes,
                             levels=unique(df$clinic_or_genes))

mut <- ggplot(df, aes(x=patient, y=clinic_or_genes, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=clinic_or_mutations)) +
       # use brewer to automatically assign colors 
       #scale_fill_brewer(palette = "Set1", na.value="Grey90") +
       # use predifined named color to assign colors
       scale_fill_manual(values=my.colors, na.value="Grey90") +
xlab("") +
ggtitle("HIV Cervical Mutations") +
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

ggsave(mut,file="Histology_Comutplot.png",width=10,height=8)
