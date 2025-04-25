
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)
library(tidyverse)
library(stats)

dirin <- "Spin_Results/";
dirout <- "Spin_Results/Perm_testing/"
dir.create(dirout)

GP_F1_act <- 0.6955
GP_F2_act <- 0.2344
GP_SA_act <- 0.4086
GP_AB_act <- 0.7967

F1_F2_act <- 0.3832
F1_SA_act <- 0.4580
F1_AB_act <- 0.6795

F2_SA_act <- 0.3224
F2_AB_act <- 0.7230

infile <- paste0(dirin, 'Corr_GP_nulls_avg.mat')
readin <- readMat(infile);
GP_F1_null <- t(readin$GP.F1)
GP_F2_null <- t(readin$GP.F2)
GP_SA_null <- t(readin$GP.SA)
infile <- paste0(dirin, 'Corr_GP_nulls_AB.mat')
readin <- readMat(infile);
GP_AB_null <- t(readin$GP.AB)

infile <- paste0(dirin, 'Corr_F1_nulls_avg.mat')
readin <- readMat(infile);
F1_F2_null <- t(readin$F1.F2)
F1_SA_null <- t(readin$F1.SA)
infile <- paste0(dirin, 'Corr_F1_nulls_AB.mat')
readin <- readMat(infile);
F1_AB_null <- t(readin$F1.AB)

infile <- paste0(dirin, 'Corr_F2_nulls_avg_SA.mat')
readin <- readMat(infile);
F2_SA_null <- t(readin$F2.SA)
infile <- paste0(dirin, 'Corr_F2_nulls_AB.mat')
readin <- readMat(infile);
F2_AB_null <- t(readin$F2.AB)

p_counter <- matrix(1,9,1)
rownames(p_counter) <- c('GP_F1','GP_F2','GP_SA','GP_AB','F1_F2','F1_SA','F1_AB','F2_SA','F2_AB')

for (fold in 1:1000) {
  if (GP_F1_null[fold] >= GP_F1_act)
  {
    p_counter[1,1] <- p_counter[1,1] + 1  
  }
  if (GP_F2_null[fold] >= GP_F2_act)
  {
    p_counter[2,1] <- p_counter[2,1] + 1  
  }
  if (GP_SA_null[fold] >= GP_SA_act)
  {
    p_counter[3,1] <- p_counter[3,1] + 1  
  }
  if (GP_AB_null[fold] >= GP_AB_act)
  {
    p_counter[4,1] <- p_counter[4,1] + 1  
  }
  if (F1_F2_null[fold] >= F1_F2_act)
  {
    p_counter[5,1] <- p_counter[5,1] + 1  
  }
  if (F1_SA_null[fold] >= F1_SA_act)
  {
    p_counter[6,1] <- p_counter[6,1] + 1  
  }
  if (F1_AB_null[fold] >= F1_AB_act)
  {
    p_counter[7,1] <- p_counter[7,1] + 1  
  }
  if (F2_SA_null[fold] >= F2_SA_act)
  {
    p_counter[8,1] <- p_counter[8,1] + 1  
  }
  if (F2_AB_null[fold] >= F2_AB_act)
  {
    p_counter[9,1] <- p_counter[9,1] + 1  
  }
}

Corr_p_values <- p_counter/1001
write.csv(Corr_p_values, file = paste0(dirout,"Corr_perm_p_vals.csv")) #save out SOW p values

Corr_p_values_FDR <- p.adjust(Corr_p_values[c('GP_F1','GP_F2','F1_F2','GP_SA','F1_SA','F2_SA'),1],"fdr") #FDR corrected p-vals
write.csv(Corr_p_values_FDR, file = paste0(dirout,"Corr_perm_p_vals_FDR.csv")) #save out p FDR values of corrs between F1, F2, GP

CorrAB_p_values_FDR <- p.adjust(Corr_p_values[c('GP_AB','F1_AB','F2_AB'),1],"fdr") #FDR corrected p-vals
write.csv(CorrAB_p_values_FDR, file = paste0(dirout,"CorrAB_perm_p_vals_FDR.csv")) #save out p FDR values of corrs between F1, F2, GP


for (corr in c('GP_F1','GP_F2','GP_SA','GP_AB','F1_F2','F1_SA','F1_AB','F2_SA','F2_AB')) {
  if (corr == 'GP_F1') {
    null_dist = as.data.frame(GP_F1_null)
    xint = GP_F1_act
  }
  else if (corr == 'GP_F2') {
    null_dist = as.data.frame(GP_F2_null)
    xint = GP_F2_act
  }
  else if (corr == 'GP_SA') {
    null_dist = as.data.frame(GP_SA_null)
    xint = GP_SA_act
  }
  else if (corr == 'GP_AB') {
    null_dist = as.data.frame(GP_AB_null)
    xint = GP_AB_act
  }
  else if (corr == 'F1_F2') {
    null_dist = as.data.frame(F1_F2_null)
    xint = F1_F2_act
  }
  else if (corr == 'F1_SA') {
    null_dist = as.data.frame(F1_SA_null)
    xint = F1_SA_act
  }
  else if (corr == 'F1_AB') {
    null_dist = as.data.frame(F1_AB_null)
    xint = F1_AB_act
  }
  else if (corr == 'F2_SA') {
    null_dist = as.data.frame(F2_SA_null)
    xint = F2_SA_act
  }
  else if (corr == 'F2_AB') {
    null_dist = as.data.frame(F2_AB_null)
    xint = F2_AB_act
  }
  colnames(null_dist) <- c('Correlations')
  ggplot(null_dist, aes(x=Correlations)) + geom_histogram()
  ggplot(null_dist, aes(x=Correlations)) + geom_histogram(binwidth=0.5) # Change the width of bins
  p<-ggplot(null_dist, aes(x=Correlations)) + geom_histogram(color="black", fill="black")
  #add null mean line, actual sum of weights line, and title
  p <- p + geom_vline(aes(xintercept=mean(x=Correlations)),
              color="gray", linetype="dashed", size=0.5) + geom_vline(aes(xintercept=xint),
              color="red", linetype="solid", size=0.5) + ggtitle(corr)
  p
  ggsave(paste0(dirout,"Correlation_perm_",corr,".tiff"))
}


