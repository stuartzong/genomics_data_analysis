#simulate some data
#runif: random uniform distribution data, min to max
#gl: generate factor levels, 
#?rep > rep(1:4, each = 2, times = 3)
#[1] 1 1 2 2 3 3 4 4 1 1 2 2 3 3 4 4 1 1 2 2 3 3 4 4
dat<-data.frame(X=runif(100,-2,2),T1=gl(n=4,k=25,labels=c("Small","Medium","Large","Big")),Site=rep(c("Site1","Site2"),time=50))
mm<-model.matrix(~Site+X*T1,dat)
betas<-runif(9,-2,2)
dat$Y<-rnorm(100,mm%*%betas,1)
summary(dat)
head(dat)

#color by treatment
#select the colors that will be used
library(RColorBrewer)
#all palette available from RColorBrewer
display.brewer.all()
#we will select the first 4 colors in the Set1 palette
cols<-brewer.pal(n=4,name="Set1")
#cols contain the names of four different colors
#create a color vector corresponding to levels in the T1 variable in dat
cols_t1<-cols[dat$T1]
#plot
plot(Y~X,dat,col=cols_t1,pch=16)

#change plot symbols
#symbols available: http://www.endmemo.com/program/R/pchsymbols.php
pch_site<-c(16,18)[factor(dat$Site)]
#the argument that control the plotting symbols is pch
plot(Y~X,dat,col=cols_t1,pch=pch_site)


#add legend
#bty: legend box type: default o=with a box, n=no box
plot(Y~X,dat,col=cols_t1,pch=pch_site)
legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="o",ncol=2,cex=0.7,pt.cex=0.7)

#legend outside the plot
plot(Y~X,dat,col=cols_t1,pch=pch_site)
legend(x=-1,y=11,legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7,xpd=TRUE)


#generate a new data frame with ordered X values
new_X<-expand.grid(X=seq(-2,2,length=10),T1=c("Small","Medium","Large","Big"),Site=c("Site1","Site2"))
#the model
m<-lm(Y~Site+X*T1,dat)
#get the predicted Y values
pred<-predict(m,new_X)
#plot
xs<-seq(-2,2,length=10)
plot(Y~X,dat,col=cols_t1,pch=pch_site)
lines(xs,pred[1:10],col=cols[1],lty=1,lwd=3)
lines(xs,pred[11:20],col=cols[2],lty=1,lwd=3)
lines(xs,pred[21:30],col=cols[3],lty=1,lwd=3)
lines(xs,pred[31:40],col=cols[4],lty=1,lwd=3)
lines(xs,pred[41:50],col=cols[1],lty=2,lwd=3)
lines(xs,pred[51:60],col=cols[2],lty=2,lwd=3)
lines(xs,pred[61:70],col=cols[3],lty=2,lwd=3)
lines(xs,pred[71:80],col=cols[4],lty=2,lwd=3)
legend(x=-1,y=13,legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),lwd=1,lty=rep(c(1,2),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7,xpd=TRUE)





#axis title etc
#some data
x<-1:100
y<-runif(100,-2,2)
#a usual plot with per default settings
plot(x,y)
#changing the axis title is pretty straightforward
plot(x,y,xlab="Index",ylab="Uniform draws")


#change the sizes of the axis labels and axis title
op<-par(no.readonly=TRUE) #this is done to save the default settings 
par(cex.lab=1.5,cex.axis=1.3)
plot(x,y,xlab="Index",ylab="Uniform draws")


#if we want big axis titles and labels we need to set more space for them
par(mar=c(6,6,3,3),cex.axis=1.5,cex.lab=2)
plot(x,y,xlab="Index",ylab="Uniform draws")



#some data
x<-1:100
y<-runif(100,-2,2)
#a usual plot with per default settings
plot(x,y)
#changing the axis title is pretty straightforward
plot(x,y,xlab="Index",ylab="Uniform draws")

Here is the plot:
  axis1

The settings of the plot are usually controlled by the par function (see ?par for the many possible arguments), once the arguments are set in par they apply to all subsequent plots. Some arguments in par (for example cex.axis) can also be set in other plot functions like axis or text. When these arguments are set in these other functions they will then apply only to the current plot. One can then control if he/she wants all plots to be affected by the change or only the current one.

#change the sizes of the axis labels and axis title
op<-par(no.readonly=TRUE) #this is done to save the default settings 
par(cex.lab=1.5,cex.axis=1.3)
plot(x,y,xlab="Index",ylab="Uniform draws")
#if we want big axis titles and labels we need to set more space for them
par(mar=c(6,6,3,3),cex.axis=1.5,cex.lab=2)
plot(x,y,xlab="Index",ylab="Uniform draws")

Here is the plot:
  axis2

A handy function to gain deeper control into the axis is the axis function which can control among other things at which values the tick marks are drawn, what axis labels to put under the tick marks, the line type and width of the axis line, the width of the tick marks, the color of the tick marks and axis line.

#we can further control the axis using the axis function
par(op) #re-set the plot to the default settings
plot(x,y,xaxt="n") #this prevent the plotting of x-axis labels 
axis(side=1,at=c(5,50,100)) #force the tick marks to be drawn at x=5, 50 and 100
#one can also specify the labels at the tick marks
plot(x,y,yaxt="n")
axis(side=2,at=c(-2,0,2),labels=c("Small","Medium","Big"))
#the axis function also control the axis line and tick marks
plot(x,y)
axis(side=3,at=c(5,25,75),lwd=4,lwd.ticks=2,col.ticks="red")
#some time you may want to remove the box around the plot and only show the axis lines
plot(x,y,bty="n",xaxt="n",yaxt="n")
axis(side=1,at=seq(0,100,20),lwd=3)
axis(side=2,at=seq(-2,2,2),lwd=3)


#tick marks finer control goes through par or axis
par(tcl=0.4,mgp=c(1.5,0,0)) #tcl control the length of the tick marks
#positive values will make the tick being drawn inside the plot
#negative values will make tick go outside
#mgp takes three values, the first one control how much line between plot and axis title
#the second between plot and axis labels and the third between plot and axis line
plot(x,y)
#another example using axis
par(op)
plot(x,y,xaxt="n",yaxt="n",xlab="",ylab="")
axis(side=1,at=seq(5,95,30),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0))
mtext(side=1,text="X axis",line=1.5) #cannot set the axis title with axis so need to use mtext
axis(side=2,at=seq(-2,2,2),tcl=0.3,lwd.ticks=3,col.ticks="orange",mgp=c(0,0,2))
mtext(side=2,text="Numbers taken randomly",line=2.2)





#some data
x<-1:100
y<-runif(100,-2,2)
#a usual plot with per default settings
plot(x,y)
#changing the axis title is pretty straightforward
plot(x,y,xlab="Index",ylab="Uniform draws")

#change the sizes of the axis labels and axis title
op<-par(no.readonly=TRUE) #this is done to save the default settings 
par(cex.lab=1.5,cex.axis=1.3)
plot(x,y,xlab="Index",ylab="Uniform draws")
#if we want big axis titles and labels we need to set more space for them
par(mar=c(6,6,3,3),cex.axis=1.5,cex.lab=2)
plot(x,y,xlab="Index",ylab="Uniform draws")


#we can further control the axis using the axis function
par(op) #re-set the plot to the default settings
plot(x,y,xaxt="n") #this prevent the plotting of x-axis labels 
axis(side=1,at=c(5,50,100)) #force the tick marks to be drawn at x=5, 50 and 100
#one can also specify the labels at the tick marks
plot(x,y,yaxt="n")
axis(side=2,at=c(-2,0,2),labels=c("Small","Medium","Big"))
#the axis function also control the axis line and tick marks
plot(x,y)
axis(side=3,at=c(5,25,75),lwd=4,lwd.ticks=2,col.ticks="red")
#some time you may want to remove the box around the plot and only show the axis lines
plot(x,y,bty="n",xaxt="n",yaxt="n")
axis(side=1,at=seq(0,100,20),lwd=3)
axis(side=2,at=seq(-2,2,2),lwd=3)

#tick marks finer control goes through par or axis
par(tcl=0.4,mgp=c(1.5,0,0)) #tcl control the length of the tick marks
#positive values will make the tick being drawn inside the plot
#negative values will make tick go outside
#mgp takes three values, the first one control how much line between plot and axis title
#the second between plot and axis labels and the third between plot and axis line
plot(x,y)
#another example using axis
par(op)
plot(x,y,xaxt="n",yaxt="n",xlab="",ylab="")
axis(side=1,at=seq(5,95,30),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0))
mtext(side=1,text="X axis",line=1.5) #cannot set the axis title with axis so need to use mtext
axis(side=2,at=seq(-2,2,2),tcl=0.3,lwd.ticks=3,col.ticks="orange",mgp=c(0,0,2))
mtext(side=2,text="Numbers taken randomly",line=2.2)


  #understanding the lines
plot(1:10,1:10,xlab="",ylab="",xaxt="n",yaxt="n")
for(i in 0:4){
  mtext(side=1,text=paste0("Line ",i),line=i, las=1)
}
for(i in 0:3){
  mtext(side=2,text=paste0("Line ",i),line=i, las=0)
}
#of course this can be changed with the mar argument in par
par(mar=c(7,2,2,2))
plot(1:10,1:10,xlab="",ylab="",xaxt="n",yaxt="n")
for(i in 0:6){
  mtext(side=1,text=paste0("Line ",i),line=i, las=1)
}
for(i in 0:1){
  mtext(side=2,text=paste0("Line ",i),line=i, las=0)
}



#but we can add some
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
for(side in 1:4){
  inner<-round(par()$mar[side],0)-1
  for(line in 0:inner){
    mtext(text=paste0("Inner line ",line),side=side,line=line, las=1)
  }
  outer<-round(par()$oma[side],0)-1
  for(line in 0:inner){
    mtext(text=paste0("Outer line ",line),side=side,line=line,las=1,outer=TRUE)
  }
}





  
  #a plot has inner and outer margins
  #by default there is no outer margins
  par()$oma

[1] 0 0 0 0

#but we can add some
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
for(side in 1:4){
  inner<-round(par()$mar[side],0)-1
  for(line in 0:inner){
    mtext(text=paste0("Inner line ",line),side=side,line=line)
  }
  outer<-round(par()$oma[side],0)-1
  for(line in 0:inner){
    mtext(text=paste0("Outer line ",line),side=side,line=line,outer=TRUE)
  }
}


  
  #Outer margins are useful in various context
  #when axis label is long and one does not want to shrink plot area
  par(op)
#example
par(cex.lab=1.7)
plot(1,1,ylab="A very very long axis titlenthat need special care",xlab="",type="n")
#one option would be to increase inner margin size
par(mar=c(5,7,4,2))
plot(1,1,ylab="A very very long axis titlenthat need special care",xlab="",type="n")
#sometime this is not desirable so one may plot the axis text outside of the plotting area
par(op)
par(oma=c(0,4,0,0))
plot(1,1,ylab="",xlab="",type="n")
mtext(text="A very very long axis titlenthat need special care",side=2,line=0,outer=TRUE,cex=1.7,las=0)



#this is particularly useful when having a plot with multiple panels and similar axis labels
par(op)
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,2))

plot(1,1,ylab="",xlab="",type="n")
plot(1,1,ylab="",xlab="",type="n")
plot(1,1,ylab="",xlab="",type="n")
plot(1,1,ylab="",xlab="",type="n")

mtext(text="A common x-axis label",side=1,line=0,outer=TRUE,las=1)
mtext(text="A common y-axis label",side=2,line=0,outer=TRUE, las=0)
dev.off()

set.seed(20160228)
#outer margins can also be used for plotting legend in them
x<-runif(10)
y<-runif(10)
cols<-rep(c("red","green","orange","yellow","black"),each=2)

par(op)
par(oma=c(2,2,0,4),mar=c(3,3,2,0),mfrow=c(2,2),pch=16)
#par(oma=c(3,3,0,4),mar=c(3,3,2,0),mfrow=c(2,2),pch=16)
for(i in 1:4){
  plot(x,y,col=cols,ylab="",xlab="")
}

mtext(text="A common x-axis label",side=1,line=0,outer=TRUE, las=1)
mtext(text="A common y-axis label",side=2,line=0,outer=TRUE, las=0)

legend(x=1,y=1.7,legend=LETTERS[1:5],col=unique(cols),pch=16,bty="n",xpd=NA)
dev.off()