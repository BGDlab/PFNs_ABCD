
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "../../Step_5/BASELINE/Weight_Haufe_031324/";
dirout <- "weight_map_HaufeAAB/"
dir.create(dirout)

PFNs_hardparcel<-read_cifti('../hardparcel_group.dscalar.nii')

catg <- 'General_PB1'
for (grp in c('A','B')) { #write and save out A and B weight maps

  infile <- paste0(dirin, catg, '_', grp, '_Weight_haufetrans_all','.mat')
  outfile <- paste0(dirout,catg, '_', grp, '_Weight_haufetrans_all_med_wall')
  
  # read in weights_by_vertex.mat
  readin <- readMat(infile);
  read_weights <- readin$PB.test
  
  #absolute value of weights
  abs_weights <- abs(read_weights)  
  
  #sum abs weights across 17 networks
  summed_weights <- matrix(0,59412)
  for (i in c(1:59412)) {
    for (f in c(1:17)) {
      summed_weights[i] <- summed_weights[i]+abs_weights[i+(f-1)*59412]
    }
  }

  #assign summed weights to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_weights_lh[i,1] <- summed_weights[i]
  }

  #assign summed weights to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_weights_rh[i,1] <- summed_weights[i+29696]
  }

  if (grp == 'A') {
    vertex_weights_lh_A <- matrix(0,29696,1)
    vertex_weights_lh_A <- vertex_weights_lh
    vertex_weights_rh_A <- matrix(0,29716,1)
    vertex_weights_rh_A <- vertex_weights_rh 
  }
  if (grp == 'B') {
    vertex_weights_lh_B <- matrix(0,29696,1)
    vertex_weights_lh_B <- vertex_weights_lh
    vertex_weights_rh_B <- matrix(0,29716,1)
    vertex_weights_rh_B <- vertex_weights_rh 
  }

  #define map variable using hard parcel cifti file
  Weight_map <- PFNs_hardparcel

  #assign weight values into hard parcel
  Weight_map$data$cortex_left <- vertex_weights_lh
  Weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Weight_map_med <- move_from_mwall(Weight_map)

  #save out cifti files in A or B dir
  #write_cifti(Weight_map_med,outfile)
}


#AVERAGE WEIGHTS FIRST
#save out A and B weights separately
catg <- 'General_PB1'
for (grp in c('A','B')) {
  
  infile <- paste0(dirin, catg, '_', grp, '_Weight_haufetrans_all','.mat')
  
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
    concat_weights_A <- read_weights
  }
  if (grp == 'B') {
    all_weights_B <- matrix(0,59412,17)
    all_weights_B <- all_weights
    concat_weights_B <- read_weights
  }
}

#average Haufe weights between groups A and B
all_weights_avg <- matrix(0,59412,17) 
for (vert in c(1:59412)) {
  for (PFN in c(1:17)) {
    all_weights_avg[vert,PFN] <- sum(all_weights_A[vert,PFN]+all_weights_B[vert,PFN])/2
  }
}
#write.csv(all_weights_avg,file = paste0(dirout,catg,"_avg_Haufe_weights.csv"))

#average Haufe weights (all PFNs concatenated) to more easily visualize distribution of weights
concat_weight_avg <- matrix(0,1010004)
for (vert in c(1:1010004)) {
  concat_weight_avg[vert] <- sum(concat_weights_A[vert]+concat_weights_B[vert])/2
}
#write.csv(concat_weight_avg,file = paste0(dirout,catg,"_avg_Haufe_weights_concat.csv"))

#ABS VAL OF WEIGHTS
abs_avg_weights <- abs(all_weights_avg)

#SUM ABS WEIGHTS across 17 networks
summed_weights <- matrix(0,59412)
for (i in c(1:59412)) {
  for (f in c(1:17)) {
    summed_weights[i] <- summed_weights[i]+abs_avg_weights[i+(f-1)*59412]
  }
}

#assign summed weights to left hemi
vertex_weights_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weights_lh[i,1] <- summed_weights[i]
}

#assign summed weights to right hemi
vertex_weights_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weights_rh[i,1] <- summed_weights[i+29696]
}

#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel

#assign weight values into hard parcel
Weight_map$data$cortex_left <- vertex_weights_lh
Weight_map$data$cortex_right <- vertex_weights_rh

#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)

#define output file for averaged map
outfile <- paste0("./",dirout,catg,'_Weight_haufetrans_all_AVG_FIRST')

#save out cifti files
#write_cifti(Weight_map_med,outfile)


#SIG ONLY (TFCE corrected)
dirout <- "weight_map_HaufeAAB/sig_testing/"
sig_weights <- read.csv(paste0(dirout,"pfac_TFCE_p_vals.csv"))

#Only take summed weight if sig
sig_summed_weights <- matrix(0,59412)
for (i in c(1:59412)) {
  if (sig_weights[i,1] < 0.05) {
    sig_summed_weights[i] <- summed_weights[i]
  }
}
write.csv(sig_summed_weights,file = paste0(dirout,"pfac_vertex_ABS_sum_Haufe_weights_sig_only_TFCE.csv"))

#assign sig summed weights to left hemi
vertex_weights_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weights_lh[i,1] <- sig_summed_weights[i]
}

#assign sig summed weights to right hemi
vertex_weights_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weights_rh[i,1] <- sig_summed_weights[i+29696]
}

#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel

#assign weight values into hard parcel
Weight_map$data$cortex_left <- vertex_weights_lh
Weight_map$data$cortex_right <- vertex_weights_rh

#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)

#define output file for averaged map
outfile <- paste0("./",dirout,'pfac_SIG_TFCE_Weight_haufetrans_all_AVG_FIRST')

#save out cifti files
write_cifti(Weight_map_med,outfile)


#AVERAGED INDIVIDUAL PFN MAPS
for (PFN in c(1:17))
{
  #assign averaged weights for 1 of the 17 networks
  PFN_weights <- matrix(0,59412)
  for (i in c(1:59412)) {
    PFN_weights[i] <- all_weights_avg[i+(PFN-1)*59412]
  }
  
  #assign averaged weights to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_weights_lh[i,1] <- PFN_weights[i]
  }
  
  #assign averaged weights to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_weights_rh[i,1] <- PFN_weights[i+29696]
  }
  
  #define map variable using hard parcel cifti file
  Weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel
  Weight_map$data$cortex_left <- vertex_weights_lh
  Weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Weight_map_med <- move_from_mwall(Weight_map)
  
  #define output file for averaged map
  outfile <- paste0("./",dirout,catg,'_PFN',PFN,'_Weight_haufetrans_all_AVG_FIRST')
  
  #save out cifti files
  write_cifti(Weight_map_med,outfile)
}


#Max Pos/Neg Weight Maps (NO LONGER USED- now calculated in "Rel_max_maps_PFNID.R" using wb -cifti-find-clusters func)
sum_threshold <- 0 #0.0075
weight_threshold <- 0.0025 #0.003
max_pos_weight_PFNID <- matrix(0,59412)
max_neg_weight_PFNID <- matrix(0,59412)

for (vert in c(1:59412)) {
  if (summed_weights[vert] > sum_threshold) { #check if summed weight is > sum_threshold
    max_tracker_pos <- 0
    max_tracker_neg <- 0
    pos_PFNID <- 0
    neg_PFNID <- 0
    for (PFN in c(1:17)) {
      if (all_weights_avg[vert,PFN] > max_tracker_pos & all_weights_avg[vert,PFN] > weight_threshold) { #check if PFN weight is larger than max pos weight AND larger than weight_threshold
        max_tracker_pos <- all_weights_avg[vert,PFN]
        pos_PFNID <- PFN
      } else if (all_weights_avg[vert,PFN] < max_tracker_neg  & abs(all_weights_avg[vert,PFN]) > weight_threshold) { #check if PFN weight is smaller than max neg weight AND (abs) larger than weight_threshold
        max_tracker_neg <- all_weights_avg[vert,PFN]
        neg_PFNID <- PFN
      }
    }
    max_pos_weight_PFNID[vert] <- pos_PFNID
    max_neg_weight_PFNID[vert] <- neg_PFNID
  }
}

#MAX POS MAP
#assign IDs to left hemi
vertex_ID_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_ID_lh[i,1] <- max_pos_weight_PFNID[i]
}
#assign IDs to right hemi
vertex_ID_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_ID_rh[i,1] <- max_pos_weight_PFNID[i+29696]
}

#Save left + right hem to mat files, converted to gifit--> cifti in matlab "Labeling_PFNs_into_cifti.m"

# For L hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_lh <- vertex_ID_lh
writeMat("Max_Pos_L.mat", x = names_list)
# For R hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_rh <- vertex_ID_rh
writeMat("Max_Pos_R.mat", x = names_list)


#define map variable using hard parcel cifti file
Pos_ID_map <- PFNs_hardparcel

#assign weight values into hard parcel
Pos_ID_map$data$cortex_left <- vertex_ID_lh
Pos_ID_map$data$cortex_right <- vertex_ID_rh

#add back in medial wall
Pos_ID_map_med <- move_from_mwall(Pos_ID_map)

#define output file for averaged map
outfile <- paste0("./",dirout,catg,'_PFNID_Pos_no_sum_thres')

#save out cifti files
write_cifti(Pos_ID_map_med,outfile)


#MAX NEG MAP
#assign IDs to left hemi
vertex_ID_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_ID_lh[i,1] <- max_neg_weight_PFNID[i]
}
#assign IDs to right hemi
vertex_ID_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_ID_rh[i,1] <- max_neg_weight_PFNID[i+29696]
}

#Save left + right hem to mat files, converted to gifit--> cifti in matlab "Labeling_PFNs_into_cifti.m"

# For L hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_lh <- vertex_ID_lh
writeMat("Max_Neg_L.mat", x = names_list)
# For R hemi, ensure that both the list and its elements are named
names_list <- list()
names_list$vertex_ID_rh <- vertex_ID_rh
writeMat("Max_Neg_R.mat", x = names_list)


#define map variable using hard parcel cifti file
Neg_ID_map <- PFNs_hardparcel

#assign weight values into hard parcel
Neg_ID_map$data$cortex_left <- vertex_ID_lh
Neg_ID_map$data$cortex_right <- vertex_ID_rh

#add back in medial wall
Neg_ID_map_med <- move_from_mwall(Neg_ID_map)

#define output file for averaged map
outfile <- paste0("./",dirout,catg,'_PFNID_Neg_no_sum_thres')

#save out cifti files
write_cifti(Neg_ID_map_med,outfile)


#TOP 1% OF SUMMED WEIGHTS AFTER AVERAGING FIRST
summed_weights_ranked = rank(summed_weights) #rank summed weights
top_1_summed_weights <- matrix(0,59412)
for (i in c(1:59412)) {
  if (summed_weights_ranked[i] > 58818) { #add summed weight to new matrix if in the top 1% of summed weights
   top_1_summed_weights[i] <- summed_weights[i]
  }
}
write.csv(top_1_summed_weights, file = paste0(dirout,catg,"_Sum_of_weights_top_1%.csv")) #save out top 1% summed weights
  
#assign top 1% summed weights to left hemi
vertex_weights_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weights_lh[i,1] <- top_1_summed_weights[i]
}
  
#assign top 1% summed weights to right hemi
vertex_weights_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weights_rh[i,1] <- top_1_summed_weights[i+29696]
}
  
#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel
  
#assign weight values into hard parcel
Weight_map$data$cortex_left <- vertex_weights_lh
Weight_map$data$cortex_right <- vertex_weights_rh
  
#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)
  
#define output file for top 1% map
outfile <- paste0("./",dirout,catg,'_Weight_haufetrans_top_1%')
  
#save out cifti files
write_cifti(Weight_map_med,outfile)


#TOP 5% OF SUMMED WEIGHTS AFTER AVERAGING FIRST
summed_weights_ranked = rank(summed_weights) #rank summed weights
top_1_summed_weights <- matrix(0,59412)
for (i in c(1:59412)) {
  if (summed_weights_ranked[i] > 56441) { #add summed weight to new matrix if in the top 5% of summed weights
    top_1_summed_weights[i] <- summed_weights[i]
  }
}
write.csv(top_1_summed_weights, file = paste0(dirout,catg,"_Sum_of_weights_top_5%.csv")) #save out top 5% summed weights

#assign top 5% summed weights to left hemi
vertex_weights_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weights_lh[i,1] <- top_1_summed_weights[i]
}

#assign top 5% summed weights to right hemi
vertex_weights_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weights_rh[i,1] <- top_1_summed_weights[i+29696]
}

#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel

#assign weight values into hard parcel
Weight_map$data$cortex_left <- vertex_weights_lh
Weight_map$data$cortex_right <- vertex_weights_rh

#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)

#define output file for top 5% map
outfile <- paste0("./",dirout,catg,'_Weight_haufetrans_top_5%')

#save out cifti files
write_cifti(Weight_map_med,outfile)


#TOP 10% OF SUMMED WEIGHTS AFTER AVERAGING FIRST
summed_weights_ranked = rank(summed_weights) #rank summed weights
top10_summed_weights <- matrix(0,59412)
for (i in c(1:59412)) {
  if (summed_weights_ranked[i] > 53471) { #add summed weight to new matrix if in the top 10% of summed weights
    top10_summed_weights[i] <- summed_weights[i]
  }
}
write.csv(top10_summed_weights, file = paste0(dirout,catg,"_Sum_of_weights_top_10%.csv")) #save out top 10% summed weights

#assign top 10% summed weights to left hemi
vertex_weights_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_weights_lh[i,1] <- top10_summed_weights[i]
}

#assign top 10% summed weights to right hemi
vertex_weights_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_weights_rh[i,1] <- top10_summed_weights[i+29696]
}

#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel

#assign weight values into hard parcel
Weight_map$data$cortex_left <- vertex_weights_lh
Weight_map$data$cortex_right <- vertex_weights_rh

#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)

#define output file for top 10% map
outfile <- paste0("./",dirout,catg,'_Weight_haufetrans_top10%')

#save out cifti files
write_cifti(Weight_map_med,outfile)
