library(stringr)
library(dplyr)
library(ggplot2)

# Function to perform 2D PCA with color coding for different sheets
PCA_2D <- function() {
  numsheets <- as.numeric(readline("How many sheets would you like to integrate: "))
  while(is.na(numsheets)) {
    numsheets <- as.numeric(readline("Invalid input. How many sheets would you like to integrate: "))
  }
  
  combined_df <- data.frame()
  sheet_labels <- c()
  for (i in 1:numsheets) {
    datasheet <- as.data.frame(read.csv(file.choose()))
    sheet_label <- rep(paste0("Sheet_", i), nrow(datasheet))
    combined_df <- rbind(combined_df, datasheet)
    sheet_labels <- c(sheet_labels, sheet_label)
  }
  
  combined_df$Sheet <- sheet_labels  # Add a new column indicating the sheet
  
  avbl_ftrs <- colnames(combined_df)
  print("Available Features: ")
  print(avbl_ftrs)
  
  num_slctd_ftrs <- as.numeric(readline("How many features would you like to integrate: "))
  while(is.na(num_slctd_ftrs)) {
    num_slctd_ftrs <- as.numeric(readline("Invalid input. How many features would you like to integrate: "))
  }
  
  slctd_ftrs <- c()
  for (i in 1:num_slctd_ftrs) {
    feature <- readline("Which feature(s) would you like to integrate into your PCA: ")
    feature <- str_to_lower(feature)
    while(!(feature %in% tolower(avbl_ftrs))) {
      feature <- readline("Invalid feature. Please enter a valid feature: ")
      feature <- str_to_lower(feature)
    }
    slctd_ftrs[i] <- feature
  }
  
  # Filter the dataframe with selected features
  slctd_df <- combined_df[, colnames(combined_df) %in% slctd_ftrs]
  
  # Perform PCA
  pca <- prcomp(slctd_df, scale = TRUE)
  pca_df <- as.data.frame(pca$x)
  pca_df$Sheet <- combined_df$Sheet  # Add the sheet information for color coding
  
  # Create the 2D PCA plot with color coding
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Sheet)) +
    geom_point(size = 3) +
    labs(title = "PCA",
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
  return(pca)  # Return PCA object for further analysis
}

# Function to create a scree plot
scree_plot <- function(pca) {
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
  
  # Create a data frame for plotting
  pca.var.per.df <- data.frame(
    Number = seq_along(pca.var.per),
    Variance = pca.var.per
  )
  
  # Create the scree plot
  ggplot(data = pca.var.per.df, aes(x = Number, y = Variance)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Variance), position = position_dodge(width = 0.9), vjust = -0.25) +
    labs(x = "Principal Component", y = "Percentage of Variance Explained") +
    ggtitle("Scree Plot") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

# Function to get the top 5 contributing variables to PC1
get_top_5 <- function(pca) {
  loading_scores <- pca$rotation[,1]
  gene_scores <- abs(loading_scores)
  gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
  top_5_variables <- names(gene_score_ranked[1:5])
  print("Top 5 contributing variables to PC1:")
  print(top_5_variables)
  print("Scores:")
  print(pca$rotation[top_5_variables,1])
}
