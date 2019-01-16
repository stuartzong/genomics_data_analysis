library(plot3D)
data(iris)
head(iris)
# x, y and z coordinates
x <- sep.l <- iris$Sepal.Length
y <- pet.l <- iris$Petal.Length
z <- sep.w <- iris$Sepal.Width
scatter3D(x, y, z, clab = c("Sepal", "Width (cm)"))

df <- read.table(header=TRUE, '/projects/trans_scratch/validations/workspace/szong/David_Kaplan/variants/run3/high_moderate_SNV_summary_no_normal.txt.128153295.somatic.filtered.3d')
df
diagnosis_af <- df$diagnois_af
relapse_af <- df$relapse_af
postmortem_af <- df$postmortem
# by default, points colored by z
scatter3D(diagnosis_af, relapse_af, postmortem_af, clab = c("diagnosis_af","%"), 
          bty = "g", 
          colkey = FALSE, 
          main =" af 3d plot bty= 'g'", 
          theta = 0,
          phi = 40)


library("plot3Drgl")
plotrgl()