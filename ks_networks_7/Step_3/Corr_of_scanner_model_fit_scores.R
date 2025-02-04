library(R.matlab)
library(ggplot2)

# Model-fit data for A and B, not including scanner
Pfac_mat_A <- readMat('../../Step_2/2_Server_All/BASELINE/baseline_results_031324/all_network/General_PB1_prediction_testA.mat')
Pfac_mat_B <- readMat('../../Step_2/2_Server_All/BASELINE/baseline_results_031324/all_network/General_PB1_prediction_testB.mat')
F1_mat_A <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_1_prediction_testA.mat')
F1_mat_B <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_1_prediction_testB.mat')
F2_mat_A <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_2_prediction_testA.mat')
F2_mat_B <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_2_prediction_testB.mat')

Pfac_A <- t(Pfac_mat_A$PB.test)
Pfac_B <- t(Pfac_mat_B$PB.test)
F1_A <- t(F1_mat_A$PB.test)
F1_B <- t(F1_mat_B$PB.test)
F2_A <- t(F2_mat_A$PB.test)
F2_B <- t(F2_mat_B$PB.test)

# Model-fit data for A and B, including scanner 
Pfac_mat_A_scan <- readMat('../../Step_2/2_Server_All/BASELINE/scanner_results_010225/all_network/General_PB1_prediction_testA.mat')
Pfac_mat_B_scan <- readMat('../../Step_2/2_Server_All/BASELINE/scanner_results_010225/all_network/General_PB1_prediction_testB.mat')
F1_mat_A_scan <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/scanner_results_010225/all_network/PRS_1_prediction_testA.mat')
F1_mat_B_scan <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/scanner_results_010225/all_network/PRS_1_prediction_testB.mat')
F2_mat_A_scan <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/scanner_results_010225/all_network/PRS_2_prediction_testA.mat')
F2_mat_B_scan <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/scanner_results_010225/all_network/PRS_2_prediction_testB.mat')

Pfac_A_scan <- t(Pfac_mat_A_scan$PB.test)
Pfac_B_scan <- t(Pfac_mat_B_scan$PB.test)
F1_A_scan <- t(F1_mat_A_scan$PB.test)
F1_B_scan <- t(F1_mat_B_scan$PB.test)
F2_A_scan <- t(F2_mat_A_scan$PB.test)
F2_B_scan <- t(F2_mat_B_scan$PB.test)

scanner_corrs <- matrix(0,3)
#Concatenate across subsets A and B and create data frame where x values are model-fit scores w/o scanner cov and y values are model-fit scores w/ scanner cov
Pfac_x <- c(Pfac_A,Pfac_B)
Pfac_y_scan <- c(Pfac_A_scan,Pfac_B_scan)
Pfac_data <- data.frame(Pfac_x, Pfac_y_scan)
Pfac_data$Pfac_x <-(Pfac_data$Pfac_x-mean(Pfac_data$Pfac_x))/sd(Pfac_data$Pfac_x) # z-score x vals
Pfac_data$Pfac_y_scan <-(Pfac_data$Pfac_y_scan-mean(Pfac_data$Pfac_y_scan))/sd(Pfac_data$Pfac_y_scan) # z-score y vals
scanner_corrs[1] <- cor(Pfac_data$Pfac_x, Pfac_data$Pfac_y_scan, method = 'pearson')

F1_x <- c(F1_A,F1_B)
F1_y_scan <- c(F1_A_scan,F1_B_scan)
F1_data <- data.frame(F1_x, F1_y_scan)
F1_data$F1_x <-(F1_data$F1_x-mean(F1_data$F1_x))/sd(F1_data$F1_x) # z-score x vals
F1_data$F1_y_scan <-(F1_data$F1_y_scan-mean(F1_data$F1_y_scan))/sd(F1_data$F1_y_scan) # z-score y vals
scanner_corrs[2] <- cor(F1_data$F1_x, F1_data$F1_y_scan, method = 'pearson')

F2_x <- c(F2_A,F2_B)
F2_y_scan <- c(F2_A_scan,F2_B_scan)
F2_data <- data.frame(F2_x, F2_y_scan)
F2_data$F2_x <-(F2_data$F2_x-mean(F2_data$F2_x))/sd(F2_data$F2_x) # z-score x vals
F2_data$F2_y_scan <-(F2_data$F2_y_scan-mean(F2_data$F2_y_scan))/sd(F2_data$F2_y_scan) # z-score y vals
scanner_corrs[3] <- cor(F2_data$F2_x, F2_data$F2_y_scan, method = 'pearson')

write.csv(scanner_corrs, file = "Pearson_Rs_of_scanner_model_fit_scores.csv") #save out top 1% summed weights

#P-factor plot
ggplot(Pfac_data, aes(x = Pfac_x, y = Pfac_y_scan)) +
  geom_point(color = "black", size = 1) +  # Scatter plot
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +  # Add a regression line
  labs(title = "P-factor Scanner Model Comparison",
       x = "Model-Fit Score, no SM covariate",
       y = "Model-Fit Score, SM covariate") +  # Titles for plot and axes
  scale_x_continuous(breaks = seq(-4, 5, by = 1)) +  # Set x-axis breaks
  scale_y_continuous(breaks = seq(-4, 5, by = 1)) +  # Set y-axis breaks
  theme_minimal()+  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center the plot title
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black", size = 0.2),  # Add axis lines
    axis.ticks = element_line(color = "black", size = 0.2)  # Add tick marks on axes
  )
ggsave("P-factor_scanner_corrplot.tiff", bg = "white", width = 10, height = 10, dpi = 600, units = "cm")

#F1 plot
ggplot(F1_data, aes(x = F1_x, y = F1_y_scan)) +
  geom_point(color = "black", size = 1) +  # Scatter plot
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +  # Add a regression line
  labs(title = "PRS-F1 Scanner Model Comparison",
       x = "Model-Fit Score, no SM covariate",
       y = "Model-Fit Score, SM covariate") +  # Titles for plot and axes
  scale_x_continuous(breaks = seq(-4, 4, by = 1)) +  # Set x-axis breaks
  scale_y_continuous(breaks = seq(-4, 4, by = 1)) +  # Set y-axis breaks
  theme_minimal()+  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center the plot title
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black")  # Add tick marks on axes
  )
ggsave("PRS-F1_scanner_corrplot.tiff", bg = "white", width = 10, height = 10, dpi = 600, units = "cm")

#F2 plot
ggplot(F2_data, aes(x = F2_x, y = F2_y_scan)) +
  geom_point(color = "black", size = 1) +  # Scatter plot
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +  # Add a regression line
  labs(title = "PRS-F2 Scanner Model Comparison",
       x = "Model-Fit Score, no SM covariate",
       y = "Model-Fit Score, SM covariate") +  # Titles for plot and axes
  scale_x_continuous(breaks = seq(-4, 4, by = 1)) +  # Set x-axis breaks
  scale_y_continuous(breaks = seq(-4, 4, by = 1)) +  # Set y-axis breaks
  theme_minimal()+  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center the plot title
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black")  # Add tick marks on axes
  )
ggsave("PRS-F2_scanner_corrplot.tiff", bg = "white", width = 10, height = 10, dpi = 600, units = "cm")

