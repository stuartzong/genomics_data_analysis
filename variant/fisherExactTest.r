#options(digits=5)
wkdir <- getwd()
for (file in list.files(path=wkdir, pattern = ".combined.adjusted$")) {
#file="A44821_A44813~TARGET-30-PAUFSR.combined"
cat("Current working file: ", file)
cat("\n")
data <- read.table(file=paste(wkdir, "/", file, sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE, 
            colClasses = c("character"))
pvals <- matrix(0,nrow=nrow(data),ncol=1)
for(i in 1:nrow(data))
{
	altP <- round(as.numeric(data[i,7]))
	refP <- round(as.numeric(data[i,6]))
	totP <- round(as.numeric(data[i,5]))

	altR <- round(as.numeric(data[i,16]))
	refR <- round(as.numeric(data[i,15]))
	totR <- round(as.numeric(data[i,14]))

	if(totR!=0 && totP!=0){
		contingencyTable1 <- matrix(c(altP, refP, altR, refR),nrow=2,ncol=2)
		fisherTest1 <- fisher.test(contingencyTable1)
                pvals[i,1]=sprintf("%.3f",fisherTest1$p.value)
	}else{
		pvals[i,1]="."
	}
}
newData<-cbind(data,pvals)
write.table(format(newData,nsmall=3,scientific=FALSE), file=paste(wkdir, "/", file, ".fisher", sep=""),sep="\t",quote=FALSE)
}


