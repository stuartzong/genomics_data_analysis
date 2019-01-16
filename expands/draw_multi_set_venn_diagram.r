library(VennDiagram)
library(grDevices)

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}


#3 set venn diagram
#tum_single1 <- read.table('NB128R_primary_snvs.txt')
#tum_single2 <- read.table('NB153_relapse_snvs.txt')
#tum_single3 <- read.table('NB295_normal_snvs.txt')

#temp <- venn.diagram(list(NB128R_primary=tum_single1[,1],NB153_relapse=tum_single2[,1],NB295_normal=tum_single3[,1]),
#    fill = c("red", "green", "blue"), alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)
#pdf(file="PATENX.Primary.pdf")
#grid.draw(temp)
#dev.off()

#3 set venn diagram
#tum_paired1 <- read.table('NB128R_NB295_strelka_snv.txt')
#tum_paired2 <- read.table('NB153_NB295_strelka_snv.txt')
#tum_paired3 <- read.table('NB153_NB128R_strelka_snv.txt')

#temp <- venn.diagram(list(NB128R_NB295=tum_paired1[,1],B153_NB295=tum_paired2[,1],NB153_NB128R=tum_paired3[,1]),
#    fill = c("red", "green", "blue"), alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)
#pdf(file="PATENX.Primary.pdf")
#grid.draw(temp)



#3 set venn diagram

tum_paired1 <- read.table('diagnosis_variants.txt.tmp')
tum_paired2 <- read.table('relapse_variants.txt.tmp')
tum_paired3 <- read.table('postmortem_variants.txt.tmp')

temp <- venn.diagram(list(diagnosis=tum_paired1[,1],relapse=tum_paired2[,1],post_mortem=tum_paired3[,1]),
    fill = c("red", "green", "blue"), alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)
pdf(file="128153295.pdf")
grid.draw(temp)


# to figure out what are shared and what are unique to each group
xx.1 <- list(diagnosis=tum_paired1[,1],relapse=tum_paired2[,1],post_mortem=tum_paired3[,1])

# Create a list of all the combinations
combs <- 
  unlist(lapply(1:length(xx.1), 
                function(j) combn(names(xx.1), j, simplify = FALSE)),
         recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)
print(combs)


elements <- 
  lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))

n.elements <- sapply(elements, length)
print(n.elements)

print(elements)

# write list into a text file
lapply("genes_overlap", write, "128153295_genes_overlap.txt", append=FALSE, ncolumns=1000)
lapply(combs, write, "128153295_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(n.elements, write, "128153295_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(elements, write, "128153295_genes_overlap.txt", append=TRUE, ncolumns=1000)




# process 333_rl_normoxic
#2 set venn diagram

tum_paired1 <- read.table('333_right_variants.txt')
tum_paired2 <- read.table('333_left_variants.txt')

temp <- venn.diagram(list(right=tum_paired1[,1],left=tum_paired2[,1]),
    fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)
pdf(file="333_rl_normoxic.pdf")
grid.draw(temp)



xx.1 <- list(right=tum_paired1[,1],left=tum_paired2[,1])

# Create a list of all the combinations
combs <-
  unlist(lapply(1:length(xx.1),
                function(j) combn(names(xx.1), j, simplify = FALSE)),
         recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)
print(combs)
elements <-
  lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))

n.elements <- sapply(elements, length)
print(n.elements)

print(elements)

# write list into a text file
lapply("genes_overlap", write, "333_rl_normoxic_genes_overlap.txt", append=FALSE, ncolumns=1000)
lapply(combs, write, "333_rl_normoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(n.elements, write, "333_rl_normoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(elements, write, "333_rl_normoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)



#2 set venn diagram

tum_paired1 <- read.table('333_right_hypo_variants.txt')
tum_paired2 <- read.table('333_right_normoxic_variants.txt')

temp <- venn.diagram(list(right_hypoxic=tum_paired1[,1],right_normoxic=tum_paired2[,1]),
    fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)
pdf(file="333_righ_normoxic_vs_hypoxic.pdf")
grid.draw(temp)



xx.1 <- list(right_hypoxic=tum_paired1[,1],right_normoxic=tum_paired2[,1])

# Create a list of all the combinations
combs <-
  unlist(lapply(1:length(xx.1),
                function(j) combn(names(xx.1), j, simplify = FALSE)),
         recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)
print(combs)
elements <-
  lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))

n.elements <- sapply(elements, length)
print(n.elements)

print(elements)

# write list into a text file
lapply("genes_overlap", write, "333_righ_normoxic_vs_hypoxic_genes_overlap.txt", append=FALSE, ncolumns=1000)
lapply(combs, write, "333_righ_normoxic_vs_hypoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(n.elements, write, "333_righ_normoxic_vs_hypoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(elements, write, "333_righ_normoxic_vs_hypoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)


#2 set venn diagram left hpyoxic vs normoxic

tum_paired1 <- read.table('333_left_hypo_variants.txt')
tum_paired2 <- read.table('333_left_normoxic_variants.txt')

temp <- venn.diagram(list(left_hypoxic=tum_paired1[,1],left_normoxic=tum_paired2[,1]),
    fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)
pdf(file="333_left_normoxic_vs_hypoxic.pdf")
grid.draw(temp)



xx.1 <- list(left_hypoxic=tum_paired1[,1],left_normoxic=tum_paired2[,1])

# Create a list of all the combinations
combs <-
  unlist(lapply(1:length(xx.1),
                function(j) combn(names(xx.1), j, simplify = FALSE)),
         recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)
print(combs)
elements <-
  lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))

n.elements <- sapply(elements, length)
print(n.elements)

print(elements)

# write list into a text file
lapply("genes_overlap", write, "333_left_normoxic_vs_hypoxic_genes_overlap.txt", append=FALSE, ncolumns=1000)
lapply(combs, write, "333_left_normoxic_vs_hypoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(n.elements, write, "333_left_normoxic_vs_hypoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)
lapply(elements, write, "333_left_normoxic_vs_hypoxic_genes_overlap.txt", append=TRUE, ncolumns=1000)



#3 set venn diagram
#DN: generate Venn Diagram
#DN_r1_data <- read.table('/projects/CIRM_HALT/RNASeq/analysis/data/plate7/Sasan/G1_G2/DEfine/DNRecurrence.w_gene_names_r1.txt')
#DN_r2_data <- read.table('/projects/CIRM_HALT/RNASeq/analysis/data/plate7/Sasan/G1_G2/DEfine/DNRecurrence.w_gene_names_r2.txt')
#DN_DESeq_data <- read.table('/projects/CIRM_HALT/RNASeq/analysis/data/plate7/Sasan/G1_G2/DESeq/DN_G1_G2_by_DESeq_FDR0.05.txt_no_header_no_1stcol_w_gene_names')

#temp <- venn.diagram(list(DN_r1_data=DN_r1_data[,1],DN_r2_data=DN_r2_data[,1],DN_DESeq_data=DN_DESeq_data[,1]),
#    fill = c("red", "green", "blue"), alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.05)

#pdf(file="G1_G2_DN_venn.pdf")
#grid.draw(temp)
#dev.off()



#4 set venn diagram
#DN: generate Venn Diagram
#tum_single <- read.table('PASTIE.Primary.single_mpileup.txt')
#tum_paired <- read.table('PASTIE.Primary.paired_mpileup.txt')
#tum_strelka <- read.table('PASTIE.Primary.strelka.txt')
#nor_single <- read.table('PASTIE.Remission.single_mpileup.txt')

#temp <- venn.diagram(list(tum_single=tum_single[,1],nor_single=nor_single[,1],tum_paired=tum_paired[,1],tum_strelka=tum_strelka[,1]),
#    fill = c("red", "green", "blue", 'yellow'), alpha = c(0.5, 0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.195)

# # # #pdf(file="PATENX.Primary.pdf")
# # # #grid.draw(temp)
# # # #dev.off()


#4 set venn diagram
#DN: generate Venn Diagram
# tum_left_normoxic <- read.table('333_left_normoxic_snvs.txt')
# tum_right_normoxic <- read.table('333_right_normoxic_snvs.txt')
# tum_left_hypoxic <- read.table('333_left_hypoxic_snvs.txt')
# tum_right_hypoxic <- read.table('333_right_hypoxic_snvs.txt')
# 
# temp <- venn.diagram(list(left_normoxic=tum_left_normoxic[,1],right_hypoxic=tum_right_hypoxic[,1],right_normoxic=tum_right_normoxic[,1],left_hypoxic=tum_left_hypoxic[,1]),
#     fill = c("red", "green", "blue", 'yellow'), alpha = c(0.5, 0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.195)
# 
# #pdf(file="PATENX.Primary.pdf")
# grid.draw(temp)
# #dev.off()




#4 set venn diagram
#DN: generate Venn Diagram
#tum_left_normoxic <- read.table('333_left_normoxic.tmp')
#tum_right_normoxic <- read.table('333_right_normoxic.tmp')
#tum_left_hypoxic <- read.table('333_left_hypoxic.tmp')
#tum_right_hypoxic <- read.table('333_right_hypoxic.tmp')

#temp <- venn.diagram(list(left_normoxic=tum_left_normoxic[,1],right_hypoxic=tum_right_hypoxic[,1],right_normoxic=tum_right_normoxic[,1],left_hypoxic=tum_left_hypoxic[,1]),
#    fill = c("red", "green", "blue", 'yellow'), alpha = c(0.5, 0.5, 0.5, 0.5), cex = 2,cat.fontface = 4, sub.just=c(0.5,1,2),lty =2, filename = NULL, resolution = 1080, cat.dist=0.195)

#pdf(file="PATENX.Primary.pdf")
#grid.draw(temp)
#dev.off()

