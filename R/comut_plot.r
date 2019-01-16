# how to make a Co-mutation (comut) plot with ggplot2
library("ggplot2")
library(ggplot2)
library(plyr)
library(reshape2)

data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/old_run/high_moderate_SNV_summary_with_normal.txt.filtered.somatic.final.txt.filtered"
#data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/old_run/high_moderate_SNV_summary_with_normal.txt.filtered.somatic.final.txt.comutplot"
#df2 <- read.table(data, header=TRUE, row.names=1, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df <- read.table(data, header=TRUE, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df

#df2$Category  <- row.names(df2)
#df2
#print("aaaaaaaaaaaaaa")
## bring your data to long format as needed by ggplot
#df2.molten <- melt(df2, value.name="num_genes", variable.name="status", na.rm=FALSE)
#df2.molten
#print ("xxxxxxxxx")
#
#df2.molten.sorted <- arrange(df2.molten, status, Category) 
#df2.molten.sorted
#
# 
# now for a Comut plot with ggplot2
# this tell ggplot do not order the patient, since it is an ordered factor already
df$sample <- factor(df$sample, levels = df$sample)

mut <- ggplot(df, aes(x=sample, y=gene, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=mutations)) +
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
axis.text.x=element_blank(),
axis.text.y=element_text(colour="Black"),
axis.title.x=element_text(face="bold"),
axis.title.y=element_blank(),
panel.grid.major.x=element_blank(),
panel.grid.major.y=element_blank(),
panel.grid.minor.x=element_blank(),
panel.grid.minor.y=element_blank(),
panel.background=element_blank()
))
 
ggsave(mut,file="2015-5-24-ExampleComutplot.png",width=10,height=8)
