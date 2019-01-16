library(RColorBrewer)

#read in chromosome sizes
chr.file <- "/home/szong/projects/resource/chrominfo.txt.test"                                                                                                              
chr.df <- read.table(chr.file, 
                     header=TRUE, 
                     sep="\t", 
                     quote="", 
                     as.is=TRUE, 
                     stringsAsFactors=FALSE)
chr.df
row.names(chr.df)

# plot barplot for the chromosomes based on size
bp <- barplot(chr.df$size, 
              xlab="Chromosome 23=chrX, 24=chrY", 
              ylim=c(0, 250000000), 
              axis.lty=0, 
              names.arg=row.names(chr.df), 
              border=NA, 
              col="grey88")

#mark1: DNA unitegrated virus
#mark1.file <- "/home/szong/projects/resource/fake_marks.txt" 
#mark1s.df <- read.table(mark1.file, 
#                        header=TRUE, 
#                        sep="\t", 
#                        quote="", 
#                        as.is=TRUE, 
#                        stringsAsFactors=FALSE)
#mark1s.df
#type1 <- mark1s.df$Type



# mark2: DNA integrated virus
mark2.file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/all_integ.tmp" 
mark2s.df <- read.table(mark2.file, 
                        header=TRUE, 
                        sep="\t", 
                        quote="", 
                        as.is=TRUE, 
                        stringsAsFactors=FALSE)
mark2s.df
type2 <- mark2s.df$Type




#mark3: RNA unitegrated virus
#mark3.file <- "/home/szong/projects/resource/fake_marks.txt" 
#mark3s.df <- read.table(mark3.file, 
#                        header=TRUE, 
#                        sep="\t", 
#                        quote="", 
#                        as.is=TRUE, 
#                        stringsAsFactors=FALSE)
#mark3s.df
#type3 <- mark3s.df$Type



# mark4: RNA integrated virus
mark4.file <- "/projects/trans_scratch/validations/workspace/szong/Cervical/integration/RNA_all_integ.tmp" 
mark4s.df <- read.table(mark4.file, 
                        header=TRUE, 
                        sep="\t", 
                        quote="", 
                        as.is=TRUE, 
                        stringsAsFactors=FALSE)
mark4s.df
type4 <- mark4s.df$Type


# mark5: chromosome centermere
mark5.file <- "/home/szong/projects/resource/centromere.pos.txt.test" 
mark5s.df <- read.table(mark5.file, 
                        header=TRUE, 
                        sep="\t", 
                        quote="", 
                        as.is=TRUE, 
                        stringsAsFactors=FALSE)
mark5s.df



# merge vectors to get the full list of types
#all.types <- unique(as.vector(rbind(type1, type2, type3, type4)))
all.types <- sort(unique(as.vector(rbind(type2, type4))))
num.colors <- length(all.types)

#use color brewer to generate a list of colors
#my.colors <- brewer.pal(num.colors,"Set1")

# use this to get more than 9 colors
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

my.colors <- getPalette(num.colors)
# use all.types factors to name the color
names(my.colors) <- factor(all.types)
names(my.colors)
my.colors

# this is how you access those colors
my.colors[mark2s.df$Type]
my.colors[mark1s.df$Type]


# add DNA unintegrated marks
#with(mark1s.df,
#     points(
#       bp[Chromosome,]-0.5, Position,
#       #bp[Chromosome,]+0, Position,
#       col=my.colors[mark1s.df$Type],
#       pch=5,
#       lwd=2, 
#       lend=2
#     )
#)


# add DNA integrated marks
with(mark2s.df,
     segments(
       bp[Chromosome,]-0.5, Position,
       bp[Chromosome,]+0, Position,
       col=my.colors[mark2s.df$Type],
       lwd=2, 
       lend=1
     )
)


# add RNA unintegrated marks
#with(mark3s.df,
#     points(
#       bp[Chromosome,]-0.5, Position,
#       #bp[Chromosome,]+0, Position,
#       col=my.colors[mark3s.df$Type],
#       pch=3,
#       lwd=2, 
#       lend=2
#     )
#)


# add RNA integrated marks
with(mark4s.df,
     segments(
       bp[Chromosome,]+0, Position,
       bp[Chromosome,]+0.5, Position,
       col=my.colors[mark4s.df$Type],
       lwd=2, 
       lend=1,
      
     )
)


# add chromosome centermere positions
with(mark5s.df,
     points(
       bp[Chromosome,]+0, Position,
       #bp[Chromosome,]+0.5, Position,
       col='black',
       pch=4,
       cex=2,
       lend=1,
       
     )
)
title("HIV Cervical Human papillomavirus Integration Sites(69 genome + transcriptome)", cex=1)
legend('topright', all.types, lty=1, col=my.colors, bty='n', cex=.75)
legend('top', 'centermere', pch=4, col='black', bty='n', cex=.75)
#legend('topleft', all.types, pch=5, col=my.colors, bty='n', cex=.75)

# get the type for all global variables 
eapply(.GlobalEnv,typeof)

