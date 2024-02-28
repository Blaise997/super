library(ChAMP)
library(missMethyl)
library(ggplot2)
library(dplyr)

MyDir <- "C:/Work/Anna_Plaksienko/Anna_methylation/Data/"
MyOutputDir <- paste0(MyDir,"ChAMP_Output/")

load(paste0(MyOutputDir, "questionable_results_from_fitting.Rdata"))
load("C:/Work/Anna_Plaksienko/RNA_Seq_workflow/RNA_Seq_workflow/contrasts.Rdata")
load("C:/Work/Anna_Plaksienko/TranscriptomicsCronXCov/contrasts.Rdata")
setwd(MyOutputDir)
#data("probe.features.epic")
#annot=tibble::rownames_to_column(probe.features,"PROBE")
#allcpgs <- annot[,1]



# Create an empty list to store results for each contrast
go<-list()
top_results_list <- list()

for (z in colnames(contrasts_full)){
  go[[z]]=gometh(res[[z]]$PROBE,all.cpg = NULL,collection = c("GO","KEGG"),array.type = "EPIC",anno = NULL,plot.bias = F)
  
  # Store the top results for each contrast in the list
  top_results_list[[z]] <- topGSA(go[[z]], n = 10)
  #go_k=gometh(res[[z]]$PROBE,all.cpg = allcpgs,collection = "KEGG",array.type = "EPIC",anno=NULL,plot.bias = T)
  #top1=topGSA(go_k,n=10)
}

# Loop over each element in the list
for(i in 1:length(gsa)) {
  
}

ncontrasts<-dim(contrasts_t)[2]

# Loop over each element in the list
for(c in 1:length(gsa)) {
  # Convert row names to a column
  gsa[[c]]$ID <- rownames(gsa[[c]])
  gsa[[c]]$corrected_qvalue<-gsa[[c]]$FDR*ncontrasts
  
  # Filter the results
  filtered_results <- gsa[[c]][gsa[[c]]$FDR < 0.01, ]
  
  # Create a dot plot
  p <- ggplot(filtered_results, aes(x=reorder(ID, -log10(FDR)), y=-log10(FDR), size=N)) +
    geom_point(alpha=0.6) +
    scale_size(range = c(1, 8)) +  # Adjust the range as needed
    coord_flip() +  # Flip the axes
    labs(x="GO Term", y="-log10 FDR", title=paste("Enriched GO Terms for Contrast", c)) +
    theme_minimal()
  
  # Print the plot
  print(p)
}



