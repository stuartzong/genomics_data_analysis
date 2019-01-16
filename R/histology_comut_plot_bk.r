# how to make a Co-mutation (comut) plot with ggplot2
library("ggplot2")
library(ggplot2)
library(plyr)
library(reshape2)

data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/new_run/comut_data.tmp"
## data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/old_run/high_moderate_SNV_summary_with_normal.txt.filtered.somatic.final.txt.comutplot"
df <- read.table(data, header=TRUE, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df
df$patient
df$status
df$patient <- factor(df$patient, levels = unique(df$patient))
df$ids <- factor(df$ids, levels = unique(df$ids))
mut <- ggplot(df, aes(x=patient, y=ids, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=status)) +
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
