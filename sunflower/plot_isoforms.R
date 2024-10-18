# Name: plot_isoforms
# Author: EY
# Date: 10/17/2024
# Version:4.2.1
# Description: will plot isoform quantification in TPM


# Load libraries
library(ggplot2)
library(dplyr)

count_data <- read.csv("/scratch/ely67071/sunflower_inflo_dev_data_b3/STAR_output/BAM_files/extracted/output_counts_transcript_counts.csv")
count_data
# Extract treatment groups
count_data2 <- count_data %>%
  mutate(Samples = case_when(
    grepl("10D", Sample) ~ "10D",
    grepl("20D", Sample) ~ "20D",
    grepl("30D", Sample) ~ "30D",
    grepl("35D", Sample) ~ "35D",
  ))


# average the TPM values across biological replicates for each Transcript_ID and sample
average_data <- count_data2 %>%
  group_by(Transcript_ID, Samples) %>%
  summarize(Average_TPM = mean(TPM, na.rm = TRUE), .groups = 'drop')  # Calculate the average across replicates

# plot
ggplot(average_data, aes(x = Samples, y = Average_TPM, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge") +  # 'dodge' positions bars side by side
  facet_wrap(~ Transcript_ID) +  
  labs(title = "Average TPM Values by Stage for Each Transcript",
       y = "Average TPM") +
  scale_fill_manual(values = c("10D" = "blue", "20D" = "red", "30D" = "green", "35D" = "orange")) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  theme_minimal() 