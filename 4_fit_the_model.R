library(limma)
library(minfi)
library(dplyr)
library(DMRcate)

MyDir <- "C:/Work/Anna_Plaksienko/Anna_methylation/Data/"
MyOutputDir <- paste0(MyDir,"ChAMP_Output/")

load(paste0(MyOutputDir, "mybetaCombat_run_slide_array.Rdata"))

Mvalue <- logit2(mybetaCombat_run_array)

#construct design matrix
#all possible combinations of state, time and stimulus
design <- model.matrix(~0 + pd$Stato:pd$Tempo:pd$Stimolo)
colnames(design)=gsub(":",".",colnames(design))
#stimulated at t0 doesn't exist so we drop those columns
design <- design[, colSums(design) != 0]
rownames(design) <- colnames(mybetaCombat_run_array)
#make more understandable column names
colnames(design) <- gsub("pd\\$Stato", "", colnames(design))
colnames(design) <- gsub("pd\\$Tempo", "", colnames(design))
colnames(design) <- gsub("pd\\$Stimolo", "", colnames(design))
colnames(design)
#print(design)
dim(design)
ngroups <- dim(design)[2]
ngroups

load("C:/Work/Anna_Plaksienko/RNA_Seq_workflow/RNA_Seq_workflow/contrasts.Rdata")

contrasts <- contrasts_full
ncontrasts <- dim(contrasts)[2]
ncontrasts
rownames(contrasts)=gsub(":",".",rownames(contrasts))
design <- design[ , rownames(contrasts)]

#computes coefficients, residual variances and standard errors
fit <- lmFit(Mvalue, design)
fit
#rownames=row.names(contrasts)
#rownames=rownames%>%gsub('[:]','_',rownames)
#rownames
contrast.matrix <- makeContrasts(INFL.0.NO-FIS.0.NO, INFL.1.NO-FIS.1.NO, INFL.48.NO-FIS.48.NO,
                                 levels=design)

#contrast.matrix=makeContrasts(row.names(contrasts),levels = design)
#converts the coefficients and standard errors to reflect the contrasts 
#rather than the original design matrix, but does not compute t-statistics or p-values.
fit2 <- contrasts.fit(fit, contrast.matrix)
#computes t-statistics and p-values from the coefficients and standard errors
fit3 <- eBayes(fit2)

res<-list()
adjPVal = 0.05
for (i in c(1:dim(contrast.matrix)[2]))
{res[[i]]<-topTable(fit3, coef=i,adjust="BH", #Benjamini,Hochberg è lo standard method
                    number = nrow(mybetaCombat_run_array))}
#a=gsub(" ","",unlist(strsplit(colnames(contrast.matrix)[i],split="-")))
#ind=match(a,colnames(fit$coefficients)) tutto questo nel ciclo


col_1<-fit$coefficients[,12]
Coeff_betaINFL.0.NO<-(2^col_1)/(1+2^col_1)

col_2<-fit$coefficients[,1]
Coeff_betaFIS.0.NO<-(2^col_2)/(1+2^col_2)

col_3<-fit$coefficients[,13]
Coeff_betaINFL.1.NO<-(2^col_3)/(1+2^col_3)

col_4<-fit$coefficients[,2]
Coeff_betaFIS.1.NO<-(2^col_4)/(1+2^col_4)

col_5<-fit$coefficients[,14]
Coeff_betaINFL.48.NO<-(2^col_5)/(1+2^col_5)

col_6<-fit$coefficients[,3]
Coeff_betaFIS.48.NO<-(2^col_6)/(1+2^col_6)


D_beta_IF.0=Coeff_betaINFL.0.NO-Coeff_betaFIS.0.NO
D_beta_IF.1=Coeff_betaINFL.1.NO-Coeff_betaFIS.1.NO
D_beta_IF.48=Coeff_betaINFL.48.NO-Coeff_betaFIS.48.NO

contrasti_B<-data.frame(D_beta_IF.0,D_beta_IF.1,D_beta_IF.48)
contrasti_B=tibble::rownames_to_column(contrasti_B,"PROBE")

threshold.padj=0.05
threshold.dbeta=0.1

c=1#primo contrasto
res <- topTable(fit3, coef = c, adjust = "BH",
                number = nrow(mybetaCombat_run_array))

A=contrasti_B[match(rownames(res),contrasti_B$PROBE),]
res1<-cbind(res,A)

res_df <- subset(res1, adj.P.Val < threshold.padj)

HyperProbes1<-rownames(res_df[res_df$D_beta_IF.0>threshold.dbeta,])
print(paste("#hypermethylated probes=",length(HyperProbes1)))

HypoProbes1<-rownames(res_df[res_df$D_beta_IF.0<(-threshold.dbeta),])#prova senza parentesi per vedere problemi freccia
print(paste("#hypomethylated probes=", length(HypoProbes1)))


c=2#secondo contrasto
res <- topTable(fit3, coef = c, adjust = "BH",
                number = nrow(mybetaCombat_run_array))

A=contrasti_B[match(rownames(res),contrasti_B$PROBE),]
res1<-cbind(res,A)

res_df <- subset(res1, adj.P.Val < threshold.padj)

HyperProbes2<-rownames(res_df[res_df$D_beta_IF.1>threshold.dbeta,])
print(paste("#hypermethylated probes=",length(HyperProbes2)))

HypoProbes2<-rownames(res_df[res_df$D_beta_IF.1<(-threshold.dbeta),])
print(paste("#hypomethylated probes=", length(HypoProbes2)))

c=3#secondo contrasto
res <- topTable(fit3, coef = c, adjust = "BH",
                number = nrow(mybetaCombat_run_array))

A=contrasti_B[match(rownames(res),contrasti_B$PROBE),]
res1<-cbind(res,A)

res_df <- subset(res1, adj.P.Val < threshold.padj)

HyperProbes3<-rownames(res_df[res_df$D_beta_IF.48>threshold.dbeta,])
print(paste("#hypermethylated probes=",length(HyperProbes3)))

HypoProbes3<-rownames(res_df[res_df$D_beta_IF.48<(-threshold.dbeta),])
print(paste("#hypomethylated probes=", length(HypoProbes3)))



myAnnotation<-cpg.annotate(object = Mvalue,datatype = "array",what = "M",
                           analysis.type = "differential",design = design,
                           contrasts = T,cont.matrix = contrast.matrix,
                           coef = colnames(contrast.matrix)[1],arraytype = "EPIC")


#fit4=contrasts.fit(fit,contrasts = contrasts[,1])
FC_thresh <- 2
FDR_thresh <- 0.05

fun_DMP <- function(c) {
  
  res <- topTable(fit3, coef = c, adjust = "BH",
                  number = nrow(mybetaCombat_run_array))
  res_df <- subset(res, adj.P.Val < 0.05)
  
  HyperProbes <- rownames(res_df[res_df$logFC > 0, ])
  print(paste("#hypermethylated probes = ", length(HyperProbes)))
  HypoProbes <- rownames(res_df[res_df$logFC < 0, ])
  print(paste("#hypomethylated probes = ", length(HypoProbes)))
  
  curr_contrast <- contrasts[ , c]
  
  delta_beta <- 0
  for (i in 1:ngroups) {
    if (curr_contrast[i] == 0) {
      delta_beta <- delta_beta
    } else {
      coeff <- fit$coefficients[, i]
      coeff <- (2^coeff)/(1 + 2^coeff)
      delta_beta <- delta_beta + coeff * curr_contrast[i]
    }
  }

  return(list(res = res, delta_beta = delta_beta))
}

#specify the path to save contrasts plots
#path <- "/home/anna/Christine_methylation/contrasts_plots"
#setwd(path)

#prepare a list where to save results
DMP <- vector(mode = "list", length = ncontrasts)
names(DMP) <- colnames(contrasts)

for (c in 1:ncontrasts) {
  DMP[[c]] <- fun_DMP(c)
}

#I get a plot like this
plot(DMP[[2]]$delta_beta, -log10(DMP[[2]]$res$adj.P.Val))


#which is not how it's supposed to look like
#so probably there is something wrong

save(list = c("DMP"),
     file = paste0(MyOutputDir, "questionable_results_from_fitting.Rdata"))
