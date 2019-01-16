my.path <- "/Users/stuartzong/GSC/alleleFrequency/"
data.path <- "/Users/stuartzong/GSC/alleleFrequency/"
pdf.path <- file.path(data.path,"pdfs")
pdf.path

dir.create(pdf.path, showWarnings = TRUE)
setwd(pdf.path)

file.list <- list.files(my.path, pattern = "^TARGET")  
file.list
file.list.len <- length(file.list)
max.value.to.plot <- 2

chrom.info.df <- read.table("/Users/stuartzong/GSC/alleleFrequency/hg19.chromInfo.simple.txt", row.names=1, header=FALSE, sep="\t", quote="", as.is=TRUE)
colnames(chrom.info.df) <- "length"
head(chrom.info.df)
rownames(chrom.info.df) <- gsub("chr","", rownames(chrom.info.df))
head(chrom.info.df)
chrom.info.df["1",]

# as a function applied over a list of filenames
sample.names <- lapply(file.list, function(x) unlist( strsplit( unlist(strsplit(x, "-"))[3], "_") )[1] )

write(unlist(sample.names),file="sample.names.txt")
write(unlist(file.list),file="file.names.txt")

seq( 1, file.list.len*max.value.to.plot, 2 )
length(seq( 1, file.list.len*max.value.to.plot, 2 ))
length(sample.names)
sample.name.list <- unlist(sample.names)
sample.name.list

for (chr in 1:25){
  pdf.chr.name <- paste0("chr",chr,".pdf")
  pdf(file.path(pdf.path, pdf.chr.name), height=8, width=8)
  chr.length.Mb <- chrom.info.df[as.character(chr),]/1000000
  par(mar=c(4.5,6,1.5,1) + 0.1) 
  
  plot(0,0,col="white", xlim=c(0,chr.length.Mb), ylim=c(0,file.list.len*max.value.to.plot), xlab="position (Mb)", ylab="", cex.lab=1.5, main=paste0("chr",chr), yaxt='n',xaxt='n')
  axis(1,cex.axis=0.7)
  axis(2, cex.axis=0.6,at=seq( 1, file.list.len*max.value.to.plot, 2 ), labels=sample.name.list, las=1, cex.lab=0.5)
  #plot(0,0,col="white", xlim=c(0,chr.length.Mb), ylim=c(0,file.list.len*max.value.to.plot), cex.lab=1.2, xlab="position (Mb)", ylab="AF", main=paste0("chr",chr))
  # vertical lines at start and end of chromosome
  lines( c(0,0),c(0,file.list.len*max.value.to.plot), col="gray")
  lines( c(chr.length.Mb,chr.length.Mb),c(0,file.list.len*max.value.to.plot), col="gray")
  # line at the top
  lines(c(0,chr.length.Mb),c(file.list.len*max.value.to.plot,file.list.len*max.value.to.plot),col="gray")
  
  for (f in 1:file.list.len){
    current.file <- file.list[f]
    current.file.split <- strsplit(current.file,'[-_]')
    current.file.name <- current.file.split[[1]][3]
    current.path.file <- file.path( data.path, current.file )
    data <- read.table(current.path.file, header=FALSE, sep="\t", quote="", as.is=TRUE)
    df <- data.frame(data[,c(1,13,19,20)])
    colnames(df) <- c("location","PriAF","MalAF","Pvalue")
    #print (df[5,])
    locations.raw <- as.character(data[,1])
    locations.chrs <- apply( as.matrix(locations.raw), 1, function(x) unlist( strsplit( x, "\\_" ) )[1] )
    locations.coords <- as.numeric(apply( as.matrix(locations.raw), 1, function(x) unlist( strsplit( x, "\\_" ) )[2] ))
    processed.df <- data.frame( "chr"=as.character(locations.chrs), "coord"=locations.coords, "PriAF"=df$PriAF,"MalAF"=df$MalAF,"Pvalue"=df$Pvalue )
    #replace chr X Y MT with numbers
    levels(processed.df$"chr")[levels(processed.df$"chr")=='MT'] <- '25'
    levels(processed.df$"chr")[levels(processed.df$"chr")=='Y'] <- '24'
    levels(processed.df$"chr")[levels(processed.df$"chr")=='X'] <- '23'
    
    table(processed.df$"chr")
    
    
    processed.ordered.df <- processed.df[ order(processed.df[,1],processed.df[,2]), ]
    processed.df.chr1 <- processed.df[processed.df$"chr"==chr, ]
    print (af.Pvalue <- processed.df.chr1$"Pvalue")
    y.baseline <- f*max.value.to.plot-1
    
    
    #draw baseline 
    lines(c(0,chr.length.Mb),c(y.baseline,y.baseline),lwd=0.5, col="lightgray")
    
    # a separator line
    lines(c(0,chr.length.Mb),c(y.baseline-1,y.baseline-1),lwd=0.5, col="gray")
    

    if (nrow(processed.df.chr1)>0){
      
      
    for (i in 1:nrow(processed.df.chr1)){
     
      if (processed.df.chr1[i,5]<=0.05){
      lines(c( processed.df.chr1[i,2]/1000000, processed.df.chr1[i,2]/1000000 ),c( (y.baseline-processed.df.chr1[i,4]), y.baseline ), type="l",col="blue")
      lines(c( processed.df.chr1[i,2]/1000000, processed.df.chr1[i,2]/1000000 ),c( y.baseline, (y.baseline+processed.df.chr1[i,3]) ), type="l",col="red")
      }
      #else{
       # lines(c( processed.df.chr1[i,2]/1000000, processed.df.chr1[i,2]/1000000 ),c( (y.baseline-processed.df.chr1[i,4]), y.baseline ), type="l",col="red")
        #lines(c( processed.df.chr1[i,2]/1000000, processed.df.chr1[i,2]/1000000 ),c( y.baseline, (y.baseline+processed.df.chr1[i,3]) ), type="l",col="red")
        
        
      #}
      
      }
    #text(25,y.baseline, current.file.name, cex=0.5, pos=2, offset=0.25)
    
    }
    else{ 
      #lines(c(0,chr.length.Mb),c(y.baseline,y.baseline),lwd=0.5, col="lightgray")
      #text(25,y.baseline, current.file.name, cex=0.1, pos=2, offset=0.25)
      
    }
  }
dev.off()
  
}
    
    
  
