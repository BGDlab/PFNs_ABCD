
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Program Files/workbench')
library(ggplot2)
library(tidyverse)
library(stats)

dirin <- "../../Step_5/BASELINE/Haufe_Nulls_031824/";
dirout <- "sum_of_weights_HaufeAAB/sig_testing/"
dir.create(dirout)

# SIG TESTING AVG B/W A AND B
sum_of_weights_actual <- read.csv('sum_of_weights_HaufeAAB/Sum_of_weights_HaufeAAB.csv') #don't need to normalize by vert for significance testing
#(would be the same normalization for actual and null sum of weights, so can ignore)

catg <- 'General_PB1'
sum_of_weights_NET_nulls <- matrix(0,17,1000) #matrices for nulls for NET and ABS
sum_of_weights_ABS_nulls <- matrix(0,17,1000)
p_counter_NET <- matrix(1,17,1) # p-counter for each PFN starts at 1
p_counter_ABS <- matrix(1,17,1)
count <- 0
for (dir in c('A','B','C','D','E','F','G','H','I','J')) { #cycle through 10 directories
  for (fold in 0:99) { #cycle through 100 files
    count <- count + 1 #count to determine placement of values in sum_of_weights null matrices
    for (grp in c('A','B')) {
      
      infile <- paste0(dirin, 'Weight_Haufe_Nulls_031824_', dir, '/', catg, '_', grp, '_Weight_haufetrans_all_', fold,'.mat')
      
      # read in weights_by_vertex.mat
      readin <- readMat(infile);
      read_weights <- readin$PB.test
      
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

    for (PFN in c(1:17)) { # cycle through all networks to compare null sum of weights to actual
      for (vert in c(1:59412)) {
        sum_of_weights_NET_nulls[PFN,count] <- sum_of_weights_NET_nulls[PFN,count] + all_weights_avg[vert,PFN] #NET
        sum_of_weights_ABS_nulls[PFN,count] <- sum_of_weights_ABS_nulls[PFN,count] + abs(all_weights_avg[vert,PFN]) #ABS
      }
      if (abs(sum_of_weights_NET_nulls[PFN,count]) >= abs(sum_of_weights_actual[PFN,2]))
      {
        p_counter_NET[PFN,1] <- p_counter_NET[PFN,1] + 1  
      }
      if (sum_of_weights_ABS_nulls[PFN,count] >= sum_of_weights_actual[PFN,3])
      {
        p_counter_ABS[PFN,1] <- p_counter_ABS[PFN,1] + 1  
      }
    }
    print(count)
  }
}

for (looper in c(1:2)) {
  if (looper == 1) {
    p_counter <- p_counter_NET
    nulls <- sum_of_weights_NET_nulls
    name <- "NET"
  } else {
    p_counter <- p_counter_ABS
    nulls <- sum_of_weights_ABS_nulls
    name <- "ABS"
  }
  SOW_p_values <- p_counter/1001
  write.csv(SOW_p_values, file = paste0(dirout,"Sum_of_weights_",name,"_p_vals.csv")) #save out SOW p values

  SOW_p_values_FDR <- p.adjust(SOW_p_values,"fdr") #FDR corrected p-vals
  write.csv(SOW_p_values_FDR, file = paste0(dirout,"Sum_of_weights_",name,"_p_vals_FDR.csv")) #save out SOW p FDR values 

  for (PFN in c(1:17)) { #plot and save histogram for each PFN
    null_dist = as.data.frame(nulls[PFN,]) #convert matrix row into data frame
    colnames(null_dist) <- c('Sum_of_weights')
    ggplot(null_dist, aes(x=Sum_of_weights)) + geom_histogram()
    ggplot(null_dist, aes(x=Sum_of_weights)) + geom_histogram(binwidth=0.5) # Change the width of bins
    p<-ggplot(null_dist, aes(x=Sum_of_weights)) + geom_histogram(color="black", fill="black")
    #add null mean line, actual sum of weights line, and title
    p <- p + geom_vline(aes(xintercept=mean(x=Sum_of_weights)),
                color="gray", linetype="dashed", size=0.5) + geom_vline(aes(xintercept=sum_of_weights_actual[PFN,looper+1]),
                color="red", linetype="solid", size=0.5) + ggtitle(PFN)
    p
    ggsave(paste0(dirout,"/Sum_of_weights_",name,"_PFN",PFN,".tiff"))
  }
}

#SIG TESTING OF A AND B SEPARATELY
catg <- 'General_PB1'
for (grp in c('A','B')) {
  sum_of_weights_NET_nulls <- matrix(0,17,1000) #matrices for nulls for NET and ABS
  sum_of_weights_ABS_nulls <- matrix(0,17,1000)
  p_counter_NET <- matrix(1,17,1) # p-counter for each PFN starts at 1
  p_counter_ABS <- matrix(1,17,1)
  count <- 0
  for (dir in c('A','B','C','D','E','F','G','H','I','J')) { #cycle through 10 directories
    for (fold in 0:99) { #cycle through 100 files
      count <- count + 1 #count to determine placement of values in sum_of_weights null matrices
      infile <- paste0(dirin, 'Weight_Haufe_Nulls_031824_', dir, '/', catg, '_', grp, '_Weight_haufetrans_all_', fold,'.mat')
      
      # read in weights_by_vertex.mat
      readin <- readMat(infile);
      read_weights <- readin$PB.test
      
      #format Haufe weights into a matrix of 17 PFNs
      all_weights <- matrix(0,59412,17)
      for (vert in c(1:59412)) {
        for (PFN in c(1:17)) {
          all_weights[vert,PFN] <- read_weights[vert+(PFN-1)*59412]
        }
      }
      
      sum_of_weights_actual <- read.csv(paste0('sum_of_weights_HaufeAAB/Sum_of_weights_',grp,'_HaufeAAB.csv'))
      
      for (PFN in c(1:17)) { # cycle through all networks to compare null sum of weights to actual
        for (vert in c(1:59412)) {
          sum_of_weights_NET_nulls[PFN,count] <- sum_of_weights_NET_nulls[PFN,count] + all_weights[vert,PFN] #NET
          sum_of_weights_ABS_nulls[PFN,count] <- sum_of_weights_ABS_nulls[PFN,count] + abs(all_weights[vert,PFN]) #ABS
        }
        if (abs(sum_of_weights_NET_nulls[PFN,count]) >= abs(sum_of_weights_actual[PFN,2]))
        {
          p_counter_NET[PFN,1] <- p_counter_NET[PFN,1] + 1  
        }
        if (sum_of_weights_ABS_nulls[PFN,count] >= sum_of_weights_actual[PFN,3])
        {
          p_counter_ABS[PFN,1] <- p_counter_ABS[PFN,1] + 1  
        }
      }
      print(count)
    }
  }
  for (looper in c(1:2)) {
    if (looper == 1) {
      p_counter <- p_counter_NET
      nulls <- sum_of_weights_NET_nulls
      name <- "NET"
    } else {
      p_counter <- p_counter_ABS
      nulls <- sum_of_weights_ABS_nulls
      name <- "ABS"
    }
    SOW_p_values <- p_counter/1001
    write.csv(SOW_p_values, file = paste0(dirout,"Sum_of_weights_",grp,"_",name,"_p_vals.csv")) #save out SOW p values
    
    SOW_p_values_FDR <- p.adjust(SOW_p_values,"fdr") #FDR corrected p-vals
    write.csv(SOW_p_values_FDR, file = paste0(dirout,"Sum_of_weights_",grp,"_",name,"_p_vals_FDR.csv")) #save out SOW p FDR values
  }
}
