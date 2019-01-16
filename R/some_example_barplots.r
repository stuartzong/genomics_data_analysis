# Fitting Labels
par(las=2) # make label text perpendicular to axis
#set the default device
par(mar= c(5,4,4,2))
#pdf(file.path(".", "tc.pdf"), height=10, width=10)
#png(file.path(".", "tc.png"))
par(mar=c(7,5,4,2)) # increase plot margin. bottom,left,top,right
#par(mar=c(7.5,7,4,2)) # increase y-axis margin.

#current work directory
#setwd("/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/variant")
getwd()


## example 1: stacked barplot
#file <-"tumour_content.txt"
file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/bbt/HPV_read_level_results_ordered.txt"
##file <-"/home/szong/projects/ALL/knownRecurrentFusions.tsv.unsplit"
data <- read.table(file,
                   header=TRUE,
                   sep="\t",
                   quote="",
                   as.is=TRUE)
patient.names <- data$patient
genome_reads <- data$genome_reads
rna_reads <- data$transcriptome_reads
both_reads <- rbind(rna_reads, genome_reads*10)

barplot(both_reads,
        main="HIV Cervical HPV BBT read support",
        ylim=c(0,0.01),
        xlab="", ylab="adjusted normalized percentage",
        col=c("darkblue", "red"),
        legend=rownames(both_reads),
        args.legend = list(x = "topleft"), # specifiy legend position
        names.arg=patient.names, # x axis tick label
        cex.axis=1, # for y axis tick font
        cex.names=0.5)# for x axis tick label font
       
#reset parameters
par(mgp = c(3, 1, 0))
par(mar= c(5,4,4,2)) # default margin
dev.off()
quit()



## example 2, grouped barplot
par(mgp = c(2.6, 0.5, 0))
file <- "/projects/trans_scratch/validations/workspace/szong/David_Kaplan/probing/fusions_probing_af_by_aligned_reads.csv"
##file <-"/home/szong/projects/ALL/knownRecurrentFusions.tsv.unsplit"
par(mar=c(7,18,4,2)) # increase plot margin. bottom,left,top,right
#png(file.path(".", "af.png"))
data <- read.table(file,
                   header=TRUE,
                   sep="\t",
                   quote="",
                   as.is=TRUE)
data <- data[ order(-data[,3], data[,5]), ]


postmortem_af <- data$postmortem_af
diagnosis_af <- data$diagnosis_af
relapse_af <- data$relapse_af
fusion_names <- data$fusion
all_af <-rbind(postmortem_af, diagnosis_af, relapse_af)



mp <- barplot(all_af, 
              axes=TRUE,
              col=c('red', 'darkblue', 'green'),
              border=NA,
              cex.axis=1.2,
              cex.lab=1.3,
              #xlim=c(0, 0.000001),
              #ylim=c(0,ylim),
              ylab="", 
              xlab="fusion allele frequency %",
              xaxs="i",
              yaxs="i",
              horiz=TRUE, 
              names.arg=fusion_names, 
              legend=rownames(all_af),
              args.legend = list(x = "topright"), # specifiy legend position
              cex.names=0.6,
              las=1, 
              beside=TRUE)


#reset parameters
par(mgp = c(3, 1, 0))
par(mar= c(5,4,4,2))


dev.off()




# example 3: line plot, three point data

par(mgp = c(2.6, 0.5, 0))
file <- "/projects/trans_scratch/validations/workspace/szong/David_Kaplan/variants/run3/high_moderate_SNV_summary_no_normal.txt.128153295.somatic.filtered.3d.RSD_0.1_AF_0.1.csv"
##file <-"/home/szong/projects/ALL/knownRecurrentFusions.tsv.unsplit"
par(mar=c(7,4,4,2)) # increase plot margin. bottom,left,top,right
#png(file.path(".", "af.png"))
data <- read.table(file,
                   header=TRUE,
                   sep="\t",
                   quote="",
                   as.is=TRUE)

postmortem_af <- data$postmortem_af
diagnosis_af <- data$diagnosis_af
relapse_af <- data$relapse_af
variant_names <- data$variant

all_af <-cbind(diagnosis_af, relapse_af, postmortem_af)

plot(0,0, type="n", xlab="status",xlim=c(0,1), ylim=c(0,1), ylab="af" )

num_variants <- length(variant_names)

for (i in 1:num_variants){
  lines(c(0, 0.5,1), all_af[i,])

 }

#reset parameters
par(mgp = c(3, 1, 0))
par(mar= c(5,4,4,2))
    
  dev.off()
  


# use the following command line command to merge multiple PDF files
#pdftk *RNA*.pdf cat output all_RNA.pdf

#convert pdf to png #convert -verbose -density 150 -trim NB333L_normoxic_vs_hypoxic_somatic_SNVs.csv.pdf -quality 100 -sharpen 0x1.0 aa.png


dev.off()





#legend(x= 80, y=40, legend = c("patients", "uniq variants"), cex=1.8, fill = cols, bty="n")
#title("Number of patients and uniq variants per gene",line=0)
legend("right", legend = c("primary", "refractory"), cex=1.3, fill = cols, bty="n")
# Draw the bar values above the bars
#text(mp, height, labels = format(height, 4), pos = 3, cex = .75)
text(text.coordinate,mp, height, pos = 4, cex = .75)
#text(x, bp, signif(x,2), pos=4)


##data$num
#pdf(file.path("/home/szong/projects/NCI-AML/ALL/", "recurrentFusionRplot.pdf"), height=8, width=16)
#pdf(file.path("/home/szong/projects/NCI-AML/ALL/", "recurrentFusionRplot.pdf"))
gene <- data$Gene
# as a function applied over a list of filenames, c concatenate two coloumns, paste join them with , as separator
#fusion.names <- lapply(data$genes, function(x) paste(unlist(strsplit(x, ","))[c(1,2)],collapse=","))
gene.names <- data$patient
num_pat <- data$primary_tc
num_var <- data$refractory_tc
str(num_pat)
#cols <- c("blue", "red")[data$status == "known"]

cols <- c("deepskyblue2", "pink")



#cols <- c(data$status)
#barplot(frequency, main="fusion occurrence", col=cols,border=NA,xlim=c(0,30),xlab="number of libraries", horiz=TRUE, names.arg=fusion.names, cex.names=0.8)
#barplot(t(data[,4:5]),main="gene occurrence",names.arg = fusion.names,cex.names = 0.4,col = cols,border=NA,ylim=c(0,80),xlab="genes",ylab="number of libraries", horiz=FALSE,)

###barplot(t(data[,4:5]), main="fusion occurrence", axes=TRUE,col=cols,border=NA,ylim=c(0,30),ylab="number of libraries", horiz=FALSE, names.arg=fusion.names, cex.axis=1.2,cex.main=1.6,cex.lab=1.2,cex.names=0.85)
#below is for plotting without splitting primary and relapse




# create a two row matrix with x and y
height <- rbind((num_pat), num_var )
text.coordinate <-rbind(num_pat+0.2, num_var)
# Use height and set 'beside = TRUE' to get pairs
# save the bar midpoints in 'mp'
# Set the bar pair labels to A:D
#mp <- barplot(height, beside = TRUE, ylim = c(0, 100), names.arg = gene.names)


#par(mgp) control distance of axis label, tick label, tick mark symbol from axes, default is c(3,1,0)
# xaxs="i" set the distance between plot and x axis to 0
par(mgp = c(2.6, 0.5, 0))
#las specifiy the style of axis labels. (0=parallel, 1=all horizontal, 2=all perpendicular to axis, 3=all vertical)las=3 vertically x axis ticks, las=2 horizontal
#mp <- barplot(height, axes=TRUE,col=cols,border=NA,cex.axis=1.2,main="", cex.main=1.3,cex.lab=1.3,xlim=c(0, 50),ylim=c(0,160),ylab="", xlab="Occurrence",horiz=TRUE, names.arg=gene.names, cex.names=1,las=1, beside=TRUE)

#determine ylim
ylim <- 3.1*length(height)/2
mp <- barplot(height, axes=TRUE,col=cols,border=NA,cex.axis=1.2,cex.lab=1.3,xlim=c(50, 100),ylim=c(0,ylim),ylab="", xlab="Tumour Content %",xaxs="i",yaxs="i",horiz=TRUE, names.arg=gene.names, cex.names=0.9,las=1, beside=TRUE)
#legend(x= 80, y=40, legend = c("patients", "uniq variants"), cex=1.8, fill = cols, bty="n")
#title("Number of patients and uniq variants per gene",line=0)
legend("right", legend = c("primary", "refractory"), cex=1.3, fill = cols, bty="n")
# Draw the bar values above the bars
#text(mp, height, labels = format(height, 4), pos = 3, cex = .75)
text(text.coordinate,mp, height, pos = 4, cex = .75)
#text(x, bp, signif(x,2), pos=4)



#barplot(counts, main="Number of patients", axes=TRUE,col=cols,border=NA,cex.axis=1.2,cex.main=1.6,cex.lab=1.2,ylim=c(0,100),ylab="number of patients", horiz=FALSE, names.arg=gene.names, cex.names=0.85,las=3, beside=TRUE)


#legend("topright", legend = c("known ALL fusions", "possible novel ALL fusions"), cex=1.2,fill = c("seagreen4", "mediumpurple3"),bty="n")
#add a box around the plotting area
###box()

###lines(c(9.7,9.7),c(0,100),lwd=1, col="gray")
###text(0, 29, "known ALL fusions", cex=1,pos = 4)
###text(10, 29, "potential novel ALL fusions", cex=1,pos = 4)
###legend("topright", legend = c("primary", "relapse"), fill = c("deepskyblue2", "pink"),bty="n")
#barplot(frequency, main="gene occurrence", axes=TRUE,col=cols,border=NA,ylim=c(0,80),ylab="number of libraries", horiz=FALSE, names.arg=fusion.names,cex.axis=1.2,cex.main=1.6,cex.lab=1.2, cex.names=0.85,las=3)

#reset parameters
par(mgp = c(3, 1, 0))
par(mar= c(5,4,4,2))


dev.off()


