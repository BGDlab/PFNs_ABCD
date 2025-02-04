
clear all

% Add package to read NPY files into matlab
addpath('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Matlab_package/npy-matlab')

% Load in prediction results
folder_name='../../2_Server_all/FULL/scanner_results_010225/all_network';
file_name_pred1_samp1='PRS_1_prediction_testA.npy';
file_name_pred2_samp1='PRS_2_prediction_testA.npy';
file_name_pred1_samp2='PRS_1_prediction_testB.npy';
file_name_pred2_samp2='PRS_2_prediction_testB.npy';
pred_PRS_1_samp1 = readNPY([folder_name '/' file_name_pred1_samp1]);
pred_PRS_2_samp1 = readNPY([folder_name '/' file_name_pred2_samp1]);
pred_PRS_1_samp2 = readNPY([folder_name '/' file_name_pred1_samp2]);
pred_PRS_2_samp2 = readNPY([folder_name '/' file_name_pred2_samp2]);

% put predictions on same scale as actual
pred_PRS_1_samp1 = zscore(pred_PRS_1_samp1);
pred_PRS_2_samp1 = zscore(pred_PRS_2_samp1);
pred_PRS_1_samp2 = zscore(pred_PRS_1_samp2);
pred_PRS_2_samp2 = zscore(pred_PRS_2_samp2);

% fifth column has the matched_group
data_for_ridge=xlsread('../../1_Prepare_Data/data_for_ridge_PRS.csv')

% Load in actual PRSs and Exposome scores
actual_PRS_1 = data_for_ridge(:,3);
actual_PRS_2 = data_for_ridge(:,4);

% split data out by matched group
actual_PRS_1_samp1 = actual_PRS_1(data_for_ridge(:,5)==1);
actual_PRS_1_samp2 = actual_PRS_1(data_for_ridge(:,5)==2);
actual_PRS_2_samp1 = actual_PRS_2(data_for_ridge(:,5)==1);
actual_PRS_2_samp2 = actual_PRS_2(data_for_ridge(:,5)==2);


%% Plotting

% load red colors
color_fill_1 = [0.75, 0.31, 0.27]
color_fill_2 = [0.85, 0.57, 0.58]

% PRS_1
figure()
scatter(actual_PRS_1_samp2,pred_PRS_1_samp2,15,color_fill_2,'Marker','^','MarkerFaceColor',color_fill_2)
[r,p,ci,stats]=corrcoef(actual_PRS_1_samp2,pred_PRS_1_samp2)
PRS_1_samp1_r = r(1,2)
PRS_1_samp1_p = p(1,2)
PRS_1_samp1_CI = [ci(1,2),stats(1,2)]
hold on
scatter(actual_PRS_1_samp1,pred_PRS_1_samp1,75,color_fill_1,'Marker','.')
[r,p,ci,stats]=corrcoef(actual_PRS_1_samp1,pred_PRS_1_samp1);
PRS_1_samp2_r = r(1,2)
PRS_1_samp2_p = p(1,2)
PRS_1_samp2_CI = [ci(1,2),stats(1,2)]
hold on
linearFit = polyfit(actual_PRS_1_samp1,pred_PRS_1_samp1,1);
hline = refline(linearFit);
hline.Color=color_fill_1;
hline.LineWidth=2;
hold on
linearFit = polyfit(actual_PRS_1_samp2,pred_PRS_1_samp2,1);
hline = refline(linearFit);
hline.Color=color_fill_2;
hline.LineWidth=2;
title(["PRS-F1 Model Performance","with Scanner Covariate"],'FontSize',20)
xlabel("Actual PRS-F1 Score",'FontSize',16)
ylabel("Model-Fit PRS-F1 Score",'FontSize',16)
axis([-4 4 -4 4])

% load purple colors
color_fill_1 = [0.38, 0.33, 0.64]
color_fill_2 = [0.67, 0.64, 0.82]

% PRS_2
figure()
scatter(actual_PRS_2_samp2,pred_PRS_2_samp2,15,color_fill_2,'Marker','^','MarkerFaceColor',color_fill_2)
[r,p,ci,stats]=corrcoef(actual_PRS_2_samp2,pred_PRS_2_samp2);
PRS_2_samp1_r = r(1,2)
PRS_2_samp1_p = p(1,2)
PRS_2_samp1_CI = [ci(1,2),stats(1,2)]
hold on
scatter(actual_PRS_2_samp1,pred_PRS_2_samp1,75,color_fill_1,'Marker','.')
[r,p,ci,stats]=corrcoef(actual_PRS_2_samp1,pred_PRS_2_samp1);
PRS_2_samp2_r = r(1,2)
PRS_2_samp2_p = p(1,2)
PRS_2_samp2_CI = [ci(1,2),stats(1,2)]
hold on
linearFit = polyfit(actual_PRS_2_samp1,pred_PRS_2_samp1,1);
hline = refline(linearFit);
hline.Color=color_fill_1;
hline.LineWidth=2;
hold on
linearFit = polyfit(actual_PRS_2_samp2,pred_PRS_2_samp2,1);
hline = refline(linearFit);
hline.Color=color_fill_2;
hline.LineWidth=2;
title(["PRS-F2 Model Performance","with Scanner Covariate"],'FontSize',20)
xlabel("Actual PRS-F2 Score",'FontSize',16)
ylabel("Model-Fit PRS-F2 Score",'FontSize',16)
axis([-4 4 -4 4])

