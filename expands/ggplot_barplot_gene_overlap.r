library(ggplot2)
library(plyr)
library(reshape2)

#barplot iwth multiple groups
#Data derived from ToothGrowth data sets are used. 
#ToothGrowth describes the effect of Vitamin C on tooth growth in Guinea pigs. 
#Three dose levels of Vitamin C (0.5, 1, and 2 mg) with each of 
#two delivery methods [orange juice (OJ) or ascorbic acid (VC)] are used :


#read in data
# tell column 1 is the row names
data <- "clonal_evolution.txt"
df2 <- read.table(data, header=TRUE, row.names=1, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df2


df2$Category  <- row.names(df2)

# bring your data to long format as needed by ggplot

df2.molten <- melt(df2, value.name="num_genes", variable.name="status", na.rm=TRUE)


df2.molten.sorted <- arrange(df2.molten, status, Category) 
df2.molten.sorted

df_cumsum <- ddply(df2.molten.sorted, "status",
                   transform, 
                   label_ypos=cumsum(num_genes) - 0.5*num_genes)


# Create the barplot
ggplot(data=df_cumsum, aes(x=status, y=num_genes, fill=Category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=num_genes), vjust=1.6, 
            color="black", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

