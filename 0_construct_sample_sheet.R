setwd("C:/Work/Anna_Plaksienko/Anna_methylation/Data/")

#read both samples sheets
ss1 <- read.csv("20220620_FT-CoxCov-02-DNAMetilazione.csv")
ss2 <- read.csv("20220830_CoxCoc_recuperi.csv")
#bind them 
ss_new <- rbind(ss1, ss2[c(9:13), ])
#add colnames an delete the header
colnames(ss_new) <- ss_new[8, ]
ss_new <- ss_new[-c(1:8), ]
#make sure Sentrix ID is numeric because csv can think it's a character
ss_new$Sentrix_ID <- as.numeric(ss_new$Sentrix_ID)

#delete columns that are empty
ss_new$X <- ss_new$Sample_Plate <- ss_new$Sample_Well <-
  ss_new$Pool_ID <- ss_new$Sample_Group <- NULL

#add run 1 and run 2
ss_new$Run <- c(rep("Run1", 48), rep("Run2", 5))

#add other metadata info. We use same file we had for RNA-Seq data
load("C:/Work/Anna_Plaksienko/RNA_Seq_workflow/RNA_Seq_workflow/metadata.Rdata")
metadata$ID <- as.character(metadata$ID)
#match the order of the samples in metadata to sample sheet
metadata <- metadata[order(match(metadata$ID, ss_new$Sample_Name)), ]

ss_new$Tempo <- metadata$tempo
ss_new$Stato <- metadata$stato
ss_new$Stimolo <- metadata$stimolo
ss_new$Sample_Group <- paste(metadata$stato, 
                                metadata$tempo, metadata$stimolo, sep = ":")
ss_new

#save sample sheet. Remember to add it to the data folder!
write.csv(ss_new, file = "sample_sheet_2.csv")
