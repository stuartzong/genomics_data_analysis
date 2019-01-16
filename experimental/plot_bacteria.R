#setwd("~/Desktop/TCGA_microbial_detection/MESO/WES/")
#MESO.WES.tab <- read.delim(file="/projects/trans_scratch/validations/workspace/szong/Cervical/bbt/CESC_RNAseq_BBT_counts.txt", header = TRUE, sep = "\t", quote = "\"",strip.white = T, stringsAsFactors = F)
#MESO.WES.tab <- read.delim(file="/projects/trans_scratch/validations/workspace/szong/Cervical/bbt/82_patients/bbt_raw_counts.tmp", header = TRUE, sep = "\t", quote = "\"",strip.white = T, stringsAsFactors = F)
MESO.WES.tab <- read.delim(file="/projects/trans_scratch/validations/workspace/szong/Cervical/bbt/101_patients/RNA_bbt_raw_counts_formatted.txt.new", header = TRUE, sep = "\t", quote = "\"",strip.white = T, stringsAsFactors = F)
dim(MESO.WES.tab)
#166 51
is.data.frame(MESO.WES.tab)


# separating diseased from normal for main data 
MESO.WES.tab.diseased.rowNumber <- which(substr(MESO.WES.tab[,2],0,8)=="diseased")
MESO.WES.tab.normal.rowNumber <-  which(substr(MESO.WES.tab[,2],0,6)=="normal")

length(MESO.WES.tab.normal.rowNumber)
#83
length(MESO.WES.tab.diseased.rowNumber)
#83


MESO.WES.tab_diseased <- MESO.WES.tab[MESO.WES.tab.diseased.rowNumber,]
MESO.WES.tab_normal <- MESO.WES.tab[MESO.WES.tab.normal.rowNumber,]

N_Aspergillus <- as.numeric(MESO.WES.tab_normal$Aspergillus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Bacteriodes <- as.numeric(MESO.WES.tab_normal$Bacteriodes)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Bradyrhizobium <- as.numeric(MESO.WES.tab_normal$Bradyrhizobium)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Burkholderia <- as.numeric(MESO.WES.tab_normal$Burkholderia)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Campylobacter <- as.numeric(MESO.WES.tab_normal$Campylobacter)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Candida_albican <- as.numeric(MESO.WES.tab_normal$Candida_albican)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Chlamydia <- as.numeric(MESO.WES.tab_normal$Chlamydia)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Clostridium <- as.numeric(MESO.WES.tab_normal$Clostridium)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Cryptococcus <- as.numeric(MESO.WES.tab_normal$cryptococcus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Entamoeba_histolytica <- as.numeric(MESO.WES.tab_normal$Entamoeba_histolytica)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Escherichia_coli <- as.numeric(MESO.WES.tab_normal$Escherichia_coli)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Fusobacterium_nucleatum <- as.numeric(MESO.WES.tab_normal$Fusobacterium_nucleatum)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Gordonia_polyisoprenivorans <- as.numeric(MESO.WES.tab_normal$Gordonia_polyisoprenivorans)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Helicobacter_pylori <- as.numeric(MESO.WES.tab_normal$Helicobacter_pylori)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Hepatitis_A_B_C <- as.numeric(MESO.WES.tab_normal$Hepatitis_A_B_C)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_T_lymphotropic <- as.numeric(MESO.WES.tab_normal$Human_T_lymphotropic_virus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_adenoviruses <- as.numeric(MESO.WES.tab_normal$Human_adenovirus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_1 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_1)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_2 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_2)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_3 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_3)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_4 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_4)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_5 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_5)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_6A <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_6A)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_6B <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_6B)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_7 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_7)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_herpesvirus_8 <- as.numeric(MESO.WES.tab_normal$Human_herpesvirus_8)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_immunodeficiency_virus <- as.numeric(MESO.WES.tab_normal$Human_immunodeficiency_virus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Human_papillomavirus <- as.numeric(MESO.WES.tab_normal$Human_papillomavirus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Klebsiella <- as.numeric(MESO.WES.tab_normal$Klebsiella)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Listeria <- as.numeric(MESO.WES.tab_normal$Listeria)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Mycobacterium <- as.numeric(MESO.WES.tab_normal$Mycobacterium)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Polyomaviruses <- as.numeric(MESO.WES.tab_normal$Polyomaviruses)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Propionibacterium <- as.numeric(MESO.WES.tab_normal$Propionibacterium)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Pseudomonas <- as.numeric(MESO.WES.tab_normal$Pseudomonas)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Ralstonia <- as.numeric(MESO.WES.tab_normal$Ralstonia)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Rotaviruses <- as.numeric(MESO.WES.tab_normal$Rotaviruses)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Saccharomyces_cerevisiae <- as.numeric(MESO.WES.tab_normal$Saccharomyces_cerevisiae)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Shigella <- as.numeric(MESO.WES.tab_normal$Shigella)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Sphingomonas <- as.numeric(MESO.WES.tab_normal$Sphingomonas)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Staphylococcus <- as.numeric(MESO.WES.tab_normal$Staphylococcus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Streptococcus <- as.numeric(MESO.WES.tab_normal$Streptococcus)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_Yersinia <- as.numeric(MESO.WES.tab_normal$Yersinia)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_other_bacteria <- as.numeric(MESO.WES.tab_normal$other_bacteria)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_other_viruses <- as.numeric(MESO.WES.tab_normal$other_viruses)/as.numeric(MESO.WES.tab_normal$total)*10^6
N_vectors <- as.numeric(MESO.WES.tab_normal$vectors)/as.numeric(MESO.WES.tab_normal$total)*10^6

D_Aspergillus <- as.numeric(MESO.WES.tab_diseased$Aspergillus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Bacteriodes <- as.numeric(MESO.WES.tab_diseased$Bacteriodes)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Bradyrhizobium <- as.numeric(MESO.WES.tab_diseased$Bradyrhizobium)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Burkholderia <- as.numeric(MESO.WES.tab_diseased$Burkholderia)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Campylobacter <- as.numeric(MESO.WES.tab_diseased$Campylobacter)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Candida_albican <- as.numeric(MESO.WES.tab_diseased$Candida_albican)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Chlamydia <- as.numeric(MESO.WES.tab_diseased$Chlamydia)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Clostridium <- as.numeric(MESO.WES.tab_diseased$Clostridium)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Cryptococcus <- as.numeric(MESO.WES.tab_diseased$cryptococcus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Entamoeba_histolytica <- as.numeric(MESO.WES.tab_diseased$Entamoeba_histolytica)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Escherichia_coli <- as.numeric(MESO.WES.tab_diseased$Escherichia_coli)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Fusobacterium_nucleatum <- as.numeric(MESO.WES.tab_diseased$Fusobacterium_nucleatum)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Gordonia_polyisoprenivorans <- as.numeric(MESO.WES.tab_diseased$Gordonia_polyisoprenivorans)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Helicobacter_pylori <- as.numeric(MESO.WES.tab_diseased$Helicobacter_pylori)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Hepatitis_A_B_C <- as.numeric(MESO.WES.tab_diseased$Hepatitis_A_B_C)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_T_lymphotropic <- as.numeric(MESO.WES.tab_diseased$Human_T_lymphotropic_virus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_adenoviruses <- as.numeric(MESO.WES.tab_diseased$Human_adenovirus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_1 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_1)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_2 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_2)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_3 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_3)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_4 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_4)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_5 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_5)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_6A <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_6A)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_6B <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_6B)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_7 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_7)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_herpesvirus_8 <- as.numeric(MESO.WES.tab_diseased$Human_herpesvirus_8)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_immunodeficiency_virus <- as.numeric(MESO.WES.tab_diseased$Human_immunodeficiency_virus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
#D_Human_papillomavirus <- as.numeric(MESO.WES.tab_diseased$HPV_Nov_2014)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Human_papillomavirus <- as.numeric(MESO.WES.tab_diseased$Human_papillomavirus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Klebsiella <- as.numeric(MESO.WES.tab_diseased$Klebsiella)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Listeria <- as.numeric(MESO.WES.tab_diseased$Listeria)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Mycobacterium <- as.numeric(MESO.WES.tab_diseased$Mycobacterium)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Polyomaviruses <- as.numeric(MESO.WES.tab_diseased$Polyomaviruses)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Propionibacterium <- as.numeric(MESO.WES.tab_diseased$Propionibacterium)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Pseudomonas <- as.numeric(MESO.WES.tab_diseased$Pseudomonas)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Ralstonia <- as.numeric(MESO.WES.tab_diseased$Ralstonia)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Rotaviruses <- as.numeric(MESO.WES.tab_diseased$Rotaviruses)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Saccharomyces_cerevisiae <- as.numeric(MESO.WES.tab_diseased$Saccharomyces_cerevisiae)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Shigella <- as.numeric(MESO.WES.tab_diseased$Shigella)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Sphingomonas <- as.numeric(MESO.WES.tab_diseased$Sphingomonas)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Staphylococcus <- as.numeric(MESO.WES.tab_diseased$Staphylococcus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Streptococcus <- as.numeric(MESO.WES.tab_diseased$Streptococcus)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_Yersinia <- as.numeric(MESO.WES.tab_diseased$Yersinia)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_other_bacteria <- as.numeric(MESO.WES.tab_diseased$other_bacteria)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_other_viruses <- as.numeric(MESO.WES.tab_diseased$other_viruses)/as.numeric(MESO.WES.tab_diseased$total)*10^6
D_vectors <- as.numeric(MESO.WES.tab_diseased$vectors)/as.numeric(MESO.WES.tab_diseased$total)*10^6

Dnames2 <- c('Aspergillus', 'Bacteriodes', 'Bradyrhizobium', 'Burkholderia',  
             'Campylobacter',  'Candida_albican','Chlamydia','Clostridium','Cryptococcus',
             'Entamoeba_histolytica', 'Escherichia_coli', 'Fusobacterium_nucleatum','Gordonia_polyisoprenivorans',
             'Helicobacter_pylori','Hepatitis_A_B_C','Human_T_lymphotropic','Human_adenoviruses',
             'Human_herpesvirus_1','Human_herpesvirus_2','Human_herpesvirus_3','Human_herpesvirus_4', 
             'Human_herpesvirus_5','Human_herpesvirus_6A', 'Human_herpesvirus_6B', 'Human_herpesvirus_7', 
             'Human_herpesvirus_8', 'Human_immunodeficiency_virus', 'Human_papillomavirus',  'Klebsiella', 
             'Listeria', 'Mycobacterium',  'Polyomaviruse', 'Propionibacterium', 'Pseudomonas',  
             'Ralstonia', 'Rotaviruses', 'Saccharomyces_cerevisiae', 'Shigella', 'Sphingomonas','Staphylococcus',  
             'Streptococcus', 'Yersinia', 'other_bacteria', 'other_viruses', 'vectors')


All_RPMs <- c(D_Aspergillus, D_Bacteriodes, D_Bradyrhizobium, D_Burkholderia,  
              D_Campylobacter,  D_Candida_albican,D_Chlamydia,D_Clostridium,D_Cryptococcus,
              D_Entamoeba_histolytica, D_Escherichia_coli, D_Fusobacterium_nucleatum,D_Gordonia_polyisoprenivorans,
              D_Helicobacter_pylori,D_Hepatitis_A_B_C,D_Human_T_lymphotropic,D_Human_adenoviruses,
              D_Human_herpesvirus_1,D_Human_herpesvirus_2,D_Human_herpesvirus_3,D_Human_herpesvirus_4, 
              D_Human_herpesvirus_5,D_Human_herpesvirus_6A, D_Human_herpesvirus_6B, D_Human_herpesvirus_7, 
              D_Human_herpesvirus_8, D_Human_immunodeficiency_virus, D_Human_papillomavirus,  D_Klebsiella, 
              D_Listeria, D_Mycobacterium,  D_Polyomaviruses, D_Propionibacterium, D_Pseudomonas,  
              D_Ralstonia, D_Rotaviruses, D_Saccharomyces_cerevisiae, D_Shigella, D_Sphingomonas,D_Staphylococcus,  
              D_Streptococcus, D_Yersinia, D_other_bacteria, D_other_viruses, D_vectors,
              N_Aspergillus, N_Bacteriodes, N_Bradyrhizobium, N_Burkholderia,  
              N_Campylobacter,  N_Candida_albican,N_Chlamydia,N_Clostridium,N_Cryptococcus,
              N_Entamoeba_histolytica, N_Escherichia_coli, N_Fusobacterium_nucleatum,N_Gordonia_polyisoprenivorans,
              N_Helicobacter_pylori,N_Hepatitis_A_B_C,N_Human_T_lymphotropic,N_Human_adenoviruses,
              N_Human_herpesvirus_1,N_Human_herpesvirus_2,N_Human_herpesvirus_3,N_Human_herpesvirus_4, 
              N_Human_herpesvirus_5,N_Human_herpesvirus_6A, N_Human_herpesvirus_6B, N_Human_herpesvirus_7, 
              N_Human_herpesvirus_8, N_Human_immunodeficiency_virus, N_Human_papillomavirus,  N_Klebsiella, 
              N_Listeria, N_Mycobacterium,  N_Polyomaviruses, N_Propionibacterium, N_Pseudomonas,  
              N_Ralstonia, N_Rotaviruses, N_Saccharomyces_cerevisiae, N_Shigella, N_Sphingomonas,N_Staphylococcus,  
              N_Streptococcus, N_Yersinia, N_other_bacteria, N_other_viruses, N_vectors)
max(log10(All_RPMs+1))
#colgrid=c("gray90","white")
colgrid=c("gray90","", "white", "")
pdf("MESO.WES_boxplot_new.pdf", height=8, width=7)
par(mar=c(5.1,11,4.1,2.1))
plot(1,1,xlim=c(0,max(log10(All_RPMs+1))), ylim=c(0.8,2*length(Dnames2)),col="white", xaxt="n",yaxt="n",xlab="",ylab="",main="")
for (i in seq(from=1, to=2*length(Dnames2), by=2)){
  #print(i, i%%4+1)
  #print(c(0,i,3,i+2,colgrid))
  rect(-0.05,i-0.5,max(log10(All_RPMs+1))+0.1,i+1.5, col = colgrid[i %% 4],border=NA)
}
boxplot(log10(N_Aspergillus+1),log10(D_Aspergillus+1), log10(N_Bacteriodes+1), log10(D_Bacteriodes+1), log10(N_Bradyrhizobium+1),log10(D_Bradyrhizobium+1),log10(N_Burkholderia+1), 
        log10(D_Burkholderia+1),log10(N_Campylobacter+1),log10(D_Campylobacter+1),log10(N_Candida_albican+1),log10(D_Candida_albican+1),
        log10(N_Chlamydia+1),log10(D_Chlamydia+1),log10(N_Clostridium+1),log10(D_Clostridium+1),log10(N_Cryptococcus+1),log10(D_Cryptococcus+1),
        log10(N_Entamoeba_histolytica+1),log10(D_Entamoeba_histolytica+1),log10(N_Escherichia_coli+1),log10(D_Escherichia_coli+1),
        log10(N_Fusobacterium_nucleatum+1),log10(D_Fusobacterium_nucleatum+1),log10(N_Gordonia_polyisoprenivorans+1),
        log10(D_Gordonia_polyisoprenivorans+1),log10(N_Helicobacter_pylori+1),log10(D_Helicobacter_pylori+1),
        log10(N_Hepatitis_A_B_C+1),log10(D_Hepatitis_A_B_C+1),
        log10(N_Human_T_lymphotropic+1),log10(D_Human_T_lymphotropic+1),
        log10(N_Human_adenoviruses+1),log10(D_Human_adenoviruses+1),log10(N_Human_herpesvirus_1+1),
        log10(D_Human_herpesvirus_1+1),log10(N_Human_herpesvirus_2+1),
        log10(D_Human_herpesvirus_2+1),log10(N_Human_herpesvirus_3+1),log10(D_Human_herpesvirus_3+1),
        log10(N_Human_herpesvirus_4+1),log10(D_Human_herpesvirus_4+1),
        log10(N_Human_herpesvirus_5+1),log10(D_Human_herpesvirus_5+1), log10(N_Human_herpesvirus_6A+1),
        log10(D_Human_herpesvirus_6A+1), log10(N_Human_herpesvirus_6B+1),
        log10(D_Human_herpesvirus_6B+1),log10(N_Human_herpesvirus_7+1),log10(D_Human_herpesvirus_7+1),
        log10(N_Human_herpesvirus_8+1),log10(D_Human_herpesvirus_8+1),
        log10(N_Human_immunodeficiency_virus+1),log10(D_Human_immunodeficiency_virus+1),
        log10(N_Human_papillomavirus+1),log10(D_Human_papillomavirus+1), 
        log10(N_Klebsiella+1),log10(D_Klebsiella+1),log10(N_Listeria+1),log10(D_Listeria+1),
        log10(N_Mycobacterium+1),log10(D_Mycobacterium+1),log10(N_Polyomaviruses+1),log10(D_Polyomaviruses+1),
        log10(N_Propionibacterium+1),log10(D_Propionibacterium+1),log10(N_Pseudomonas+1),
        log10(D_Pseudomonas+1),log10(N_Ralstonia+1),log10(D_Ralstonia+1),log10(N_Rotaviruses+1),
        log10(D_Rotaviruses+1),log10(N_Saccharomyces_cerevisiae+1),log10(D_Saccharomyces_cerevisiae+1),
        log10(N_Shigella+1),log10(D_Shigella+1),log10(N_Sphingomonas+1), 
        log10(D_Sphingomonas+1),log10(N_Staphylococcus+1),log10(D_Staphylococcus+1),log10(N_Streptococcus+1),
        log10(D_Streptococcus+1),log10(N_Yersinia+1),log10(D_Yersinia+1),
        log10(N_other_bacteria+1),log10(D_other_bacteria+1),log10(N_other_viruses+1),log10(D_other_viruses+1),log10(N_vectors+1),log10(D_vectors+1),
        col=c("lightblue","pink1"),las =1 ,horizontal = TRUE, cex.axis=0.6, main='HIV_CERVICAL microbial load', 
        xlab=expression(paste("Log"[10],"(RPM+1)")), yaxt="n",add=T,
        outpch=c(21,21), outcol=c("cornflowerblue","pink3"), outbg=c("lightblue","pink1"),outcex=0.7)
axis(2,at=seq(from=1.5, to=2*length(Dnames2), by=2), labels= Dnames2, las=2, cex.axis=0.7,mgp=c(3,0.6,0),tck=-0.02)
grid(nx=NULL, ny=NA)
legend("bottomright", c("Tumour - 101 samples","Normal - 0 samples"), cex=0.8,
       pt.bg = c("pink1","lightblue"), pt.cex = 1, pt.lwd = 0.5,
       pch=c(21,21), col=c("pink3","cornflowerblue"), bg="white")
abline(v=log10(1.2), col=1,lty=5)
text(0.2, -1, 'RPM = 0.2', cex = 0.8 )
dev.off()


# pdf("MESO_WES_boxplot_manuscript.pdf", height=10, width=8)
# par(mar=c(5.1,0.5,0.5,1))
# plot(1,1,xlim=c(0,max(log10(All_RPMs+1))), ylim=c(0.8,2*length(Dnames2)),col="white", xaxt="n",yaxt="n",xlab="",ylab="",main="")
# for (i in seq(from=1, to=2*length(Dnames2), by=2)){
#   #print(i, i%%4+1)
#   #print(c(0,i,3,i+2,colgrid))
#   rect(-0.05,i-0.5,max(log10(All_RPMs+1))+0.1,i+1.5, col = colgrid[i %% 4],border=NA)
# }
# boxplot(log10(N_Aspergillus+1),log10(D_Aspergillus+1),log10(N_Bradyrhizobium+1),log10(D_Bradyrhizobium+1),log10(N_Burkholderia+1), 
#         log10(D_Burkholderia+1),log10(N_Campylobacter+1),log10(D_Campylobacter+1),log10(N_Candida_albican+1),log10(D_Candida_albican+1),
#         log10(N_Chlamydia+1),log10(D_Chlamydia+1),log10(N_Clostridium+1),log10(D_Clostridium+1),log10(N_Cryptococcus+1),log10(D_Cryptococcus+1),
#         log10(N_Entamoeba_histolytica+1),log10(D_Entamoeba_histolytica+1),log10(N_Escherichia_coli+1),log10(D_Escherichia_coli+1),
#         log10(N_Fusobacterium_nucleatum+1),log10(D_Fusobacterium_nucleatum+1),log10(N_Gordonia_polyisoprenivorans+1),
#         log10(D_Gordonia_polyisoprenivorans+1),log10(N_Helicobacter_pylori+1),log10(D_Helicobacter_pylori+1),
#         log10(N_Hepatitis_A_B_C+1),log10(D_Hepatitis_A_B_C+1),
#         log10(N_Human_T_lymphotropic+1),log10(D_Human_T_lymphotropic+1),
#         log10(N_Human_adenoviruses+1),log10(D_Human_adenoviruses+1),log10(N_Human_herpesvirus_1+1),
#         log10(D_Human_herpesvirus_1+1),log10(N_Human_herpesvirus_2+1),
#         log10(D_Human_herpesvirus_2+1),log10(N_Human_herpesvirus_3+1),log10(D_Human_herpesvirus_3+1),
#         log10(N_Human_herpesvirus_4+1),log10(D_Human_herpesvirus_4+1),
#         log10(N_Human_herpesvirus_5+1),log10(D_Human_herpesvirus_5+1), log10(N_Human_herpesvirus_6A+1),
#         log10(D_Human_herpesvirus_6A+1), log10(N_Human_herpesvirus_6B+1),
#         log10(D_Human_herpesvirus_6B+1),log10(N_Human_herpesvirus_7+1),log10(D_Human_herpesvirus_7+1),
#         log10(N_Human_herpesvirus_8+1),log10(D_Human_herpesvirus_8+1),
#         log10(N_Human_immunodeficiency_virus+1),log10(D_Human_immunodeficiency_virus+1),
#         log10(N_Human_papillomavirus+1),log10(D_Human_papillomavirus+1), 
#         log10(N_Klebsiella+1),log10(D_Klebsiella+1),log10(N_Listeria+1),log10(D_Listeria+1),
#         log10(N_Mycobacterium+1),log10(D_Mycobacterium+1),log10(N_Polyomaviruses+1),log10(D_Polyomaviruses+1),
#         log10(N_Propionibacterium+1),log10(D_Propionibacterium+1),log10(N_Pseudomonas+1),
#         log10(D_Pseudomonas+1),log10(N_Ralstonia+1),log10(D_Ralstonia+1),log10(N_Rotaviruses+1),
#         log10(D_Rotaviruses+1),log10(N_Saccharomyces_cerevisiae+1),log10(D_Saccharomyces_cerevisiae+1),
#         log10(N_Shigella+1),log10(D_Shigella+1),log10(N_Sphingomonas+1), 
#         log10(D_Sphingomonas+1),log10(N_Staphylococcus+1),log10(D_Staphylococcus+1),log10(N_Streptococcus+1),
#         log10(D_Streptococcus+1),log10(N_Yersinia+1),log10(D_Yersinia+1),
#         log10(N_other_bacteria+1),log10(D_other_bacteria+1),log10(N_other_viruses+1),log10(D_other_viruses+1),log10(N_vectors+1),log10(D_vectors+1),
#         col=c("lightblue","pink1"),las =1 ,horizontal = TRUE, cex.axis=0.6, main='', 
#         xlab=expression(paste("Log"[10],"(RPM+1)")), yaxt="n",add=T)
# axis(4,at=seq(from=1.5, to=2*length(Dnames2), by=2), labels= Dnames2, las=2, cex.axis=0.7)
# grid(nx=NULL, ny=NA)
# legend("bottomright", c("Normal - 11 samples","Tumour - 39 samples"),fill = c("lightblue","pink1"), cex=0.8,
#        bg="white")
# 
# abline(v=log10(1.2), col=1,lty=5)
# text(0.1, -1, 'RPM = 0.2', cex = 0.8 )
# dev.off()






