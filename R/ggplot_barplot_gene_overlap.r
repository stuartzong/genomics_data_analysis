library(ggplot2)
library(plyr)
library(reshape2)

#barplot with multiple groups, 3 types of tissue and 7 categories of somatic mutations (SNVs)
#type	diagnosis	relapse 	postmortem
#unique2Diagnosis	27	0	0
#Diagnosis+PM	35	0	35
#unique2Relapse	0	49	0
#unique2PM	0	0	80
#sharedByAll	213	213	213
#Relapse+PM	0	19	19
#Diagnosis+Relapse	16	16	0

#2 groups
#type	NB333R_hypoxic	NB333L_hypoxic
#NB333R	25	0
#NB333L	0	48
#shared	14	14


#pass input file name from command line
#args <- commandArgs(trailingOnly = TRUE)
#data <- args[1]


# list all input files
#my.path <- "/projects/trans_scratch/validations/workspace/szong/David_Kaplan/SV"
my.path <- getwd()
files <- list.files(my.path, pattern = "128153295_somatic_SNVs.csv$")
files
files.len <- length(files)

for (file in files){
  print(file)
  pdf <- paste0(file, ".pdf")
  pdf(pdf, height=8, width=8)



#read in data
# column 1 is the row name
data <- file
df2 <- read.table(data, header=TRUE, row.names=1, sep="\t", quote="", as.is=TRUE, stringsAsFactors=FALSE)
df2
#                  diagnosis relapse postmortem
#unique2Diagnosis         27       0          0
#Diagnosis+PM             35       0         35
#unique2Relapse            0      49          0
#unique2PM                 0       0         80
#sharedByAll             213     213        213
#Relapse+PM                0      19         19
#Diagnosis+Relapse        16      16          0

df2$Category  <- row.names(df2)

# bring your data to long format as needed by ggplot
df2.molten <- melt(df2, value.name="num_genes", variable.name="status", na.rm=TRUE)
df2.molten.sorted <- arrange(df2.molten, status, Category) 
df2.molten.sorted

#Using Category as id variables
#            Category     status num_genes
#1       Diagnosis+PM  diagnosis        35
#2  Diagnosis+Relapse  diagnosis        16
#3         Relapse+PM  diagnosis         0
#4        sharedByAll  diagnosis       213
#5   unique2Diagnosis  diagnosis        27
#6          unique2PM  diagnosis         0
#7     unique2Relapse  diagnosis         0
#8       Diagnosis+PM    relapse         0
#9  Diagnosis+Relapse    relapse        16
#10        Relapse+PM    relapse        19
#11       sharedByAll    relapse       213
#12  unique2Diagnosis    relapse         0
#13         unique2PM    relapse         0
#14    unique2Relapse    relapse        49
#15      Diagnosis+PM postmortem        35
#16 Diagnosis+Relapse postmortem         0
#17        Relapse+PM postmortem        19
#18       sharedByAll postmortem       213
#19  unique2Diagnosis postmortem         0
#20         unique2PM postmortem        80
#21    unique2Relapse postmortem         0



# figure out cumulative y position for bar labels (number of genes)
df_cumsum <- ddply(df2.molten.sorted, "status",transform, label_ypos=cumsum(num_genes) - 0.35*num_genes)
df_cumsum
#            Category     status num_genes label_ypos
#1       Diagnosis+PM  diagnosis        35       17.5
#2  Diagnosis+Relapse  diagnosis        16       43.0
#3         Relapse+PM  diagnosis         0       51.0
#4        sharedByAll  diagnosis       213      157.5
#5   unique2Diagnosis  diagnosis        27      277.5
#6          unique2PM  diagnosis         0      291.0
#7     unique2Relapse  diagnosis         0      291.0
#8       Diagnosis+PM    relapse         0        0.0
#9  Diagnosis+Relapse    relapse        16        8.0
#10        Relapse+PM    relapse        19       25.5
#11       sharedByAll    relapse       213      141.5
#12  unique2Diagnosis    relapse         0      248.0
#13         unique2PM    relapse         0      248.0
#14    unique2Relapse    relapse        49      272.5
#15      Diagnosis+PM postmortem        35       17.5
#16 Diagnosis+Relapse postmortem         0       35.0
#17        Relapse+PM postmortem        19       44.5
#18       sharedByAll postmortem       213      160.5
#19  unique2Diagnosis postmortem         0      267.0
#20         unique2PM postmortem        80      307.0
#21    unique2Relapse postmortem         0      347.0




# Create the barplot
# if not print, pdf files will be empty, weired
# supply a data frame to geom_text so that 0 will not labelled for the bar height

print(
  ggplot(data=df_cumsum, aes(x=status, y=num_genes, fill=Category)) +
  geom_bar(stat="identity")+
  geom_text(data=subset(df_cumsum, num_genes>0), aes(y=label_ypos, label=num_genes), vjust=0.6, color="black", size=3.5)+
  xlab("Neuroblastoma Cell Lines") +
  ylab("Number of Events") +
  ggtitle(pdf) +
  scale_fill_brewer(palette="Paired")#+
  #theme(axis.title.x = element_text(size = rel(3.8), angle = 00))
  #theme_minimal()
  ) 

dev.off()

}

# use the following command line command to merge multiple PDF files
#pdftk *RNA*.pdf cat output all_RNA.pdf

#convert pdf to png #convert -verbose -density 150 -trim NB333L_normoxic_vs_hypoxic_somatic_SNVs.csv.pdf -quality 100 -sharpen 0x1.0 aa.png

