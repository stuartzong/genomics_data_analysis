library(ggplot2)
river<-data.frame(
    long=c(-2.816452494909265,-2.845487331898639,-2.883036393822358),
    lat=c(56.38229290416972,56.36346886284386,56.36577994637793))

samploc<-data.frame(
    site=c("Site1","Site2","Site3"),
    long=c(-2.826213585663894,-2.816519300644918,-2.868437228090127),
    lat=c(56.3649482229089,56.38166100310631,56.36716019476281))

ggplot(samploc, aes(x = long, y = lat)) + 
  geom_point() + 
  geom_text(aes(label = site), vjust = 2) + 
  geom_line(data = river, aes(y = lat))
