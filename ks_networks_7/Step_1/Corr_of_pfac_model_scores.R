library(ggplot2)

# Model-fit data for A and B, not including scanner
Pfac_long_full <- read.csv('abcd_longitudinal_psychopathology_factors_scores_full_sample.csv')
Pfac_long_full <- Pfac_long_full[Pfac_long_full$eventname=="1",] #filter for baseline timepoint
Pfac_CS_full <- read.csv('Sun_ABCD_T1_psychopathology_factors_scores_full_sample.csv')

Pfac_long <- (Pfac_long_full$General_p)
Pfac_CS <- (Pfac_CS_full$General_p)
Pfac_data <- data.frame(Pfac_long, Pfac_CS)
Pfac_corr <- cor(Pfac_long, Pfac_CS, method = 'pearson')

#P-factor plot
ggplot(Pfac_data, aes(x = Pfac_long, y = Pfac_CS)) +
  geom_point(color = "black", size = 1) +  # Scatter plot
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +  # Add a regression line
  labs(title = "Bifactor Model Comparison\nLongitudinal vs. Cross-Sectional",
       x = "Baseline P-factor, Longitudinal Model",
       y = "Baseline P-factor, Cross-Sectional Model") +  # Titles for plot and axes
  scale_x_continuous(breaks = seq(-2, 4, by = 1)) +  # Set x-axis breaks
  scale_y_continuous(breaks = seq(-2, 4, by = 1)) +  # Set y-axis breaks
  theme_minimal()+  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center the plot title
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black", size = 0.2),  # Add axis lines
    axis.ticks = element_line(color = "black", size = 0.2)  # Add tick marks on axes
  )
ggsave("P-factor_long_vs_CS.tiff", bg = "white", width = 10, height = 10, dpi = 600, units = "cm")
