library(ggplot2)
mtcars$name <- row.names(mtcars)
mtcars
p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point()
p + geom_point() + 
  geom_text(data=subset(mtcars, wt > 4 | mpg > 25),
            aes(wt,mpg,label=name))
