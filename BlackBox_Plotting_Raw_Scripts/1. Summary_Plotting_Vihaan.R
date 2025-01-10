library(ggplot2)
library(stringr)
library(patchwork)  # Load the patchwork package
library(knitr)

# Function to load and process tracking data
summary_load <- function(tag) {
  file_path <- file.choose()  # Prompt user to select file
  
  # Make a tag dataset to save all tag datasets
  if (exists("summary_tag_memory", envir = .GlobalEnv)) {
    summary_tag_memory[length(summary_tag_memory) + 1] <<- tag
  } else {
    summary_tag_memory <<- c(tag)
  }
  
  suppressWarnings({
    assign(paste0(tag, "_summary_data"), read.csv(file_path, row.names = 1, header = TRUE, fileEncoding = "UTF-8"), envir = .GlobalEnv)
  })
  
  metrics <<- colnames(get(paste0(tag, "_summary_data")))
  mouse_names <<- rownames(get(paste0(tag, "_summary_data")))
}

remove_periods_and_underscores <- function(text) {
  # Replace all underscores with spaces
  text <- gsub("_", " ", text)
  
  # Replace multiple periods in a row with a single period
  text <- gsub("\\.{2,}", ".", text)
  
  # Replace all periods with spaces
  text <- gsub("\\.", " ", text)
  
  # Replace multiple spaces in a row with a single space
  text <- gsub(" +", " ", text)
  
  # Trim leading and trailing spaces
  text <- trimws(text)
  
  return(text)
}

plot_summary_data <- function(tag, mouse_name, metric) {
  # Find the index of the mouse and metric
  row_index <- grep(mouse_name, rownames(get(paste0(tag, "_summary_data"))))
  col_index <- grep(metric, colnames(get(paste0(tag, "_summary_data"))))
  
  # Extract the data point for the mouse and metric
  datapt <- data.frame(Value = get(paste0(tag, "_summary_data"))[row_index, col_index])
  datapt$Treat_group <- rownames(get(paste0(tag, "_summary_data")))[row_index]
  
  # Determine colors for positive and negative values
  datapt$Color <- ifelse(datapt$Value >= 0, "positive", "negative")
  
  # Define custom colors
  colors <- c("positive" = "blue", "negative" = "red")
  
  # Create the plot
  plot <- ggplot(data = datapt, aes(x = Treat_group, y = Value, fill = Color)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = colors) + 
    geom_text(aes(label = round(Value, 2)), vjust = ifelse(datapt$Value >= 0, -0.5, 1.5), color = "black", size = 3.5) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    labs(
      title = str_to_title(paste0("Average ",remove_periods_and_underscores(metric), " Bar Plot")), 
      y = paste0(str_to_title(paste(remove_periods_and_underscores(metric), " Value"))),
      x = ""
    ) +
    theme_classic(base_size = 22) + 
    theme(legend.position = "none")
  
  return(plot)
}

summary_load_and_plot <- function() {
  plot_count <- as.integer(readline("# of Plots: "))
  combined_plot <- NULL
  
  for (i in 1:plot_count) {
    tag <- readline(paste0("Enter Plot ", i, " tag: "))
    summary_load(tag)
    
    print("Available Metrics: ")
    print(metrics)
    
    mtric <- readline(paste0("Enter Plot ", i, " Metric: "))
    
    print("Mouse Names: ")
    print(mouse_names)

    animal <- readline(paste0("Enter Plot ", i, " Mouse Name: "))
    
    current_plot <- plot_summary_data(tag = tag, mouse_name = animal, metric = mtric)
    
    if (is.null(combined_plot)) {
      combined_plot <- current_plot
    } else {
      combined_plot <- combined_plot + current_plot + plot_layout(ncol = 1)
    }
  }
  
  print(combined_plot)
}