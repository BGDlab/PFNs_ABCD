
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)
library(tidyverse)
library(stats)

dirin <- "../../../Step_5/FULL/Weight_Haufe_Nulls_041024/";
dirout <- "weight_maps_HaufeAAB/sig_testing/"
dir.create(dirout)

# SIG TESTING AVG B/W A AND B
for (catg in c('PRS_1','PRS_2')) { 
  sum_of_weights_actual <- read.csv(paste0('weight_maps_HaufeAAB/',catg, '_vertex_ABS_sum_Haufe_weights.csv'))
  
  sum_of_weights_ABS_nulls <- matrix(0,59412,1000)
  p_counter_ABS <- matrix(1,59412,1)
  count <- 0
  for (dir in c('A','B','C','D','E','F','G','H','I','J')) { #cycle through 10 directories
    for (fold in 0:99) { #cycle through 100 files
      count <- count + 1 #count to determine placement of values in sum_of_weights null matrices
      for (grp in c('A','B')) {
        
        infile <- paste0(dirin, 'Weight_Haufe_Nulls_041024_All_cov_', dir, '/', catg, '_', grp, '_Weight_haufetrans_all_', fold,'.mat')
        
        # read in weights_by_vertex.mat
        readin <- readMat(infile);
        read_weights <- readin$Weights
        
        #format Haufe weights into a matrix of 17 PFNs
        all_weights <- matrix(0,59412,17)
        for (vert in c(1:59412)) {
          for (PFN in c(1:17)) {
            all_weights[vert,PFN] <- read_weights[vert+(PFN-1)*59412]
          }
        }
        
        if (grp == 'A') {
          all_weights_A <- matrix(0,59412,17)
          all_weights_A <- all_weights
        }
        if (grp == 'B') {
          all_weights_B <- matrix(0,59412,17)
          all_weights_B <- all_weights
        }
      }
      
      #average Haufe weights between groups A and B
      all_weights_avg <- matrix(0,59412,17) 
      for (vert in c(1:59412)) {
        for (PFN in c(1:17)) {
          all_weights_avg[vert,PFN] <- sum(all_weights_A[vert,PFN]+all_weights_B[vert,PFN])/2
        }
      }
  
      for (vert in c(1:59412)) { # cycle through all vertices to compare null sum of weights across networks to actual
        summed_abs_weights <- 0 #records sum of weights across 17 networks
        for (PFN in c(1:17)) {
          summed_abs_weights <- summed_abs_weights + abs(all_weights_avg[vert,PFN]) # Sum absolute values across PFNs
        }
        sum_of_weights_ABS_nulls[vert,count] <- sum_of_weights_ABS_nulls[vert,count] + summed_abs_weights # assign abs SOW to matrix
        
        if (sum_of_weights_ABS_nulls[vert,count] >= sum_of_weights_actual[vert,2])
        {
          p_counter_ABS[vert,1] <- p_counter_ABS[vert,1] + 1  
        }
      }
      print(count)
    }
  }
  
  SOW_p_values <- p_counter_ABS/1001
  write.csv(SOW_p_values, file = paste0(dirout,catg,"_vertex_sum_of_weights_p_vals.csv")) #save out SOW p values
  
  SOW_p_values_FDR <- p.adjust(SOW_p_values,"fdr") #FDR corrected p-vals
  write.csv(SOW_p_values_FDR, file = paste0(dirout,catg,"_vertex_sum_of_weights_p_vals_FDR.csv")) #save out SOW p FDR values 
  
}

