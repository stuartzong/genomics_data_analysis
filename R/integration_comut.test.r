# how to make a Co-mutation (comut) plot with ggplot2
library("ggplot2")
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)



print("lllllllllllllllllllll") 
#Some test data
dat <- data.frame(x=runif(10),y=runif(10),
        grp = rep(LETTERS[1:5],each = 2),stringsAsFactors = TRUE)
dat
dat$grp
#Create a custom color scale
myColors <- brewer.pal(5,"Set1")
myColors
names(myColors) <- levels(dat$grp)
levels(dat$grp)
names(myColors)
colScale <- scale_colour_manual(name = "grp",values = myColors)
colScale
#One plot with all the data
p <- ggplot(dat,aes(x,y,colour = grp)) + geom_point()
p1 <- p + colScale

#A second plot with only four of the levels
p2 <- p %+% droplevels(subset(dat[4:10,])) + colScale
ggsave(p1,file="2015-5-24-ExampleComutplot.png",width=10,height=8)
