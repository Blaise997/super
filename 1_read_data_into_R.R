library("ChAMP")

MyDir <- "C:/Work/Anna_Plaksienko/Anna_methylation/Data/"
MyOutputDir <- paste0(MyDir,"ChAMP_Output/")
if (!file.exists(MyOutputDir)){
  dir.create(file.path(MyOutputDir))
}
MyDataDir <- paste0(MyDir,"data")

myLoad <- champ.load(MyDir ,arraytype = "EPIC")

#make sure everything is factors in the phenodata!
{
  #make sure sample names are characters, not numbers
  myLoad$pd$Sample_Name <- as.character(myLoad$pd$Sample_Name)
  #make sure everything is factors
  myLoad$pd$Slide <- as.factor(myLoad$pd$Slide)
  levels(myLoad$pd$Slide)
  myLoad$pd$Array <- as.factor(myLoad$pd$Array)
  levels(myLoad$pd$Array)
  myLoad$pd$Run <- as.factor(myLoad$pd$Run)
  levels(myLoad$pd$Run)
  myLoad$pd$Tempo <- as.factor(myLoad$pd$Tempo)
  levels(myLoad$pd$Tempo)
  myLoad$pd$Stato <- as.factor(myLoad$pd$Stato)
  levels(myLoad$pd$Stato)
  myLoad$pd$Stimolo <- as.factor(myLoad$pd$Stimolo)
  levels(myLoad$pd$Stimolo)
  myLoad$pd$Sample_Group <- as.factor(myLoad$pd$Sample_Group)
  levels(myLoad$pd$Sample_Group)
  
}

#sav raw data object
save(myLoad, file = paste0(MyOutputDir, "MyLoad.Rdata"))
load(paste0(MyOutputDir, "MyLoad.Rdata"))

library(openxlsx)
df<-merge(myLoad$beta,myLoad$detP)
write.xlsx(df,"C:/Work/GEO/methylation_df.xlsx")
