df <- data.frame(dose=c("D0.5", "D1", "D2"),
                 len=c(4.2, 10, 29.5))

head(df)
library(ggplot2)
# Basic barplot
p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")
p

# Horizontal bar plot
p + coord_flip()

# Change the width of bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", width=0.5)

# Change colors
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", color="blue", fill="white")

# Minimal theme + blue fill color
p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()
p

# choose which item to display
p + scale_x_discrete(limits=c("D0.5", "D2"))

# barplot with labels
# Outside bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal()

# Inside bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_minimal()



# use mtcars data

head(mtcars)

# Don't map a variable to y
ggplot(mtcars, aes(x=factor(cyl)))+
  geom_bar(stat="bin", width=0.7, fill="steelblue")+
  theme_minimal()

# Change barplot line colors by groups
p<-ggplot(df, aes(x=dose, y=len, color=dose)) +
  geom_bar(stat="identity", fill="white")
p

#it is also possible to change manually barplot line colors using the functions :  
#scale_color_manual() : to use custom colors
#scale_color_brewer() : to use color palettes from RColorBrewer package
#scale_color_grey() : to use grey color palettes


# Use custom color palettes
p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# Use brewer color palettes
p+scale_color_brewer(palette="Dark2")

# Use grey scale
p + scale_color_grey() + theme_classic()

#in the R code below, barplot fill colors are automatically controlled by the levels of dose
# Change barplot fill colors by groups
p<-ggplot(df, aes(x=dose, y=len, fill=dose)) +
  geom_bar(stat="identity")+theme_minimal()
p

#It is also possible to change manually barplot fill colors using the functions :
#scale_fill_manual() : to use custom colors
#scale_fill_brewer() : to use color palettes from RColorBrewer package
#scale_fill_grey() : to use grey color palettes

# Use custom color palettes
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# use brewer color palettes
p+scale_fill_brewer(palette="Dark2")

# Use grey scale
p + scale_fill_grey()


#use black outline color
ggplot(df, aes(x=dose, y=len, fill=dose))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()

# Change bar fill colors to blues
p <- p+scale_fill_brewer(palette="Blues")

p + theme(legend.position="top")

p + theme(legend.position="bottom")

# Remove legend
p + theme(legend.position="none")


#The function scale_x_discrete can be used to change the order of items to “2”, “0.5”, “1” :
p + scale_x_discrete(limits=c("D2", "D0.5", "D1"))



#barplot iwth multiple groups
#Data derived from ToothGrowth data sets are used. 
#ToothGrowth describes the effect of Vitamin C on tooth growth in Guinea pigs. 
#Three dose levels of Vitamin C (0.5, 1, and 2 mg) with each of 
#two delivery methods [orange juice (OJ) or ascorbic acid (VC)] are used :
df2 <- data.frame(supp=rep(c("VC", "OJ"), each=3),
                  dose=rep(c("D0.5", "D1", "D2"),2),
                  len=c(6.8, 15, 33, 4.2, 10, 29.5))

head(df2)

#A stacked barplot is created by default. 
#You can use the function position_dodge() to change this. 
#The barplot fill color is controlled by the levels of dose :
# Stacked barplot with multiple groups
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")

# Use position=position_dodge()
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", position=position_dodge())


# Change the colors manually
p <- ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))

# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")

ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=len), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

#Add labels to a stacked barplot : 3 steps are required

#Sort the data by dose and supp : the package plyr is used
#Calculate the cumulative sum of the variable len for each dose
#Create the plot
library(plyr)
# Sort by dose and supp
df_sorted <- arrange(df2, dose, supp) 
head(df_sorted)


# Calculate the cumulative sum of len for each dose
df_cumsum <- ddply(df_sorted, "dose",
                   transform, label_ypos=cumsum(len))
head(df_cumsum)

# Create the barplot
ggplot(data=df_cumsum, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()


#if you want to place the labels at the middle of bars, you have to modify the cumulative sum as follow

df_cumsum <- ddply(df_sorted, "dose",
                   transform, 
                   label_ypos=cumsum(len) - 0.5*len)



# Create the barplot
ggplot(data=df_cumsum, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()





