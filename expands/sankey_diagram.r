

require(rCharts)

#require(json)
wkdir <- getwd()
setwd(wkdir)
pdf(file.path(wkdir, "sankey.pdf"), height=100, width=160)

file <-"128153295.sankey" 
df.sps <- read.table(file, header=TRUE, sep="\t", quote="", as.is=TRUE, stringsAsFactors = FALSE)
workingdata=df.sps
colnames(workingdata)=c('Sector','source','target','value')
mut.num <- nrow(workingdata)

sankeyPlot=function(df){
sankeyPlot <- rCharts$new()
sankeyPlot$setLib('/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/expands/test/rCharts_d3_sankey-gh-pages/')
sankeyPlot$setTemplate(script = "/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/expands/test/rCharts_d3_sankey-gh-pages/layouts/chart.html")

sankeyPlot$set(
  data = df,
  nodeWidth = 15,
  nodePadding = 10,
  layout = 32,
  width = 750,
  height = 12*mut.num,
  labelFormat = ".1%"
)
sankeyPlot


}



agg.pri.sps.2.mutations=aggregate(value ~ primary_sp + mutation, df.sps, sum)
colnames(agg.pri.sps.2.mutations)=c('source','target','value')

#We can now generate a single data file combing all source and target data
mutation.2.refractory.sps=subset(workingdata,select=c('source','target','value'))
all.source.target.df=rbind(mutation.2.refractory.sps,agg.pri.sps.2.mutations)


 

sankeyPlot(all.source.target.df)


dev.off()
