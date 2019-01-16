# /gsc/software/linux-x86_64-centos6/R-3.1.1/bin/Rscript use this R to run it locally
# load the ggplot2 and grid packages
library(ggplot2)
library(grid)
require(directlabels)
require(scales)
library(stringr)
pdf(file.path("bbb.pdf"), height=8, width=10)
for (file in list.files(path=".", pattern = "_cf_plot.csv$")) {
    #file <- "TARGET-21-PARZIA_cf_plot.csv"
    data <- read.table(file, header=TRUE, sep="\t", quote="", as.is=TRUE)
    #order data frame
    data <- data[with(data, order(refractory,primary)), ]
    outfile <- paste(file, "pdf", sep= ".")
    title <- paste(strsplit(file, "_")[[1]][1], "cellular frequency", sep= " ")
    print(title)
    #pdf(file.path("./", outfile), height=8, width=10)
    primary_cf <- data$primary
    refractory_cf <-data$refractory
    non_coding <- subset(data, impact != "MODERATE" & impact != "HIGH")
    non_coding_x <- non_coding$primary
    non_coding_y <- non_coding$refractory
    non_silent <- subset(data, impact == "MODERATE" | impact == "HIGH")
    #qplot(primary_cf, refractory_cf, data=data)
    #snv <- str_split_fixed(data$mutation, "_", 2)[,1]
    snv <- non_silent$gene
    non_silent_x <- non_silent$primary
    non_silent_y <- non_silent$refractory
    num_snvs <- length(snv)
     
    qplot <- qplot(x=non_silent_x, y=non_silent_y, data=non_silent,xlab="Primary cellular frequency", ylab="Refractory cellular frequency",main=title,color=snv,size=I(6))+
             scale_colour_discrete(guide=FALSE)
    final_plot <- direct.label(qplot)+
                  geom_point(data=non_coding, aes(x=non_coding_x, y=non_coding_y), colour="black", size=2)
    print(final_plot)

}
dev.off()
