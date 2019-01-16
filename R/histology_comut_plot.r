# how to make a Co-mutation (comut) plot with ggplot2
library("ggplot2")
library(ggplot2)
library(plyr)
library(reshape2)

data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/new_run/comut_data.tmp"
df <- read.table(data, header=TRUE, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df
df$patient
df$clinic_or_mutations
df$patient <- factor(df$patient, levels = unique(df$patient))
df$clinic_or_genes <- factor(df$clinic_or_genes, levels = unique(df$clinic_or_genes))
mut <- ggplot(df, aes(x=patient, y=clinic_or_genes, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=clinic_or_mutations)) +
scale_fill_brewer(palette = "Set1", na.value="Grey90") +
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
