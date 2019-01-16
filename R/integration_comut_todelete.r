# how to make a Co-mutation (comut) plot with ggplot2
library("ggplot2")
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)



#data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/NCI_HIV_Cervial_BBT_results_03212016_KM-1.csv.filtered"
#data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/new_integration_2.txt"
#data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/old_run/new_mutations_2.txt"
#df2 <- read.table(data, header=TRUE, row.names=1, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
data <- "/projects/trans_scratch/validations/workspace/szong/Cervical/variant_bwamem/new_run/comut_data.txt"
df <- read.table(data, header=TRUE, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df$patient
status <- df$status

# store color blind pallette
my.colors <- c("purple", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E41A1C", "#377EB8" )

# all categories
categories <-c("Adeno", "Squamous", "Positive", "Negative", "Multiple", "NON_SYNONYMOUS_CODING","SPLICE_SITE_ACCEPTOR","STOP_GAINED", "integrated", "unintegrated") 
#my.colors <- brewer.pal(5,"Set1")

#Create a custom color scale for each category
names(my.colors) <- categories
print("xxxxxxxxx88888888")
names(my.colors)
print("my.colors are:")
my.colors

#col.scale <- scale_colour_manual(name ="categories", values = my.colors, na.value="Grey90")
#col.scale
print("---------------------")

#manul assign brew color
#my.cols <- brewer.pal(6, "Blues")
#my.cols <- brewer.pal(status)
#my.cols
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
df$patient <- factor(df$patient, levels = unique(df$patient))
df$ids <- factor(df$ids, levels = unique(df$ids))
print("yyyyyyyyyyyyyyyyyyy")
mut <- ggplot(df, aes(x=patient, y=ids, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=status)) +
#scale_fill_brewer(palette = "Set1", na.value="Grey90") +
scale_fill_manual(values=my.colors, na.value="Grey90") +
#scale_colour_manual("", values = c("aaa" = "red", "bbb" = "black", "ccc" = "blue", "ddd" = "green"),labels = c('Company D','Company A','Company C','Company B')) +
#col.scale +
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
print("lllllllllllllllllllll") 
ggsave(mut,file="2015-5-24-ExampleComutplot.png",width=10,height=8)
