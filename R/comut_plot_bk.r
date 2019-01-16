# how to make a Co-mutation (comut) plot with ggplot2
library("ggplot2")
# library('biomaRt')
source("~/Downloads/biocLite.R")
# biocLite("biomaRt")
library("biomaRt") 

# first - lets get random gene names from biomaRt
#listMarts() # will show all available databases
listMarts(host="www.ensembl.org") 
print("aaaaaaa")
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
print("bbbbbbbb")
#listDatasets(mart) # will list all available datasets
mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
head(listAttributes(mart), n = 50)
 
# Extract information from biomart for uniprot gene names
results <- getBM(attributes = c("uniprot_genename"), mart = mart)
results

 
#random sample of 40 genes
RandomGenes <- sample(results$uniprot_genename,40)
#RandomGenes
 
# Now we will make an example data frame for a Co-mutation plot
 
# This data is generally obtained from MutSig output (or Genome MuSiC, etc)
 
# with 40 genes of interest in 100 subjects
df <- expand.grid(gene=1:40, subj=1:100)
df$Mutation <- as.factor(sample(c(rep("Missense",1500),
rep("Nonsense",500),
rep("Frame Shift",1000),
rep("Indel", 250),
rep("Splice Site", 749), 4000)))
df$Mutation[sample(4000,3700)] <- NA
 
df$gene <- RandomGenes
df 
# now for a Comut plot with ggplot2
mut <- ggplot(df, aes(x=subj, y=gene, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=Mutation)) +
scale_fill_brewer(palette = "Set1", na.value="Grey90") +
xlab("Subject") +
ggtitle("Example Comut plot") +
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
