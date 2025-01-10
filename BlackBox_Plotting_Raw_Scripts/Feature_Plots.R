###################################################
#### Plot behavioral feature time-series data ####
###################################################

# Load packages
# Can install rhdf5 by install BiocManager, then running: BiocManager::install("rhdf5")
library(rhdf5)
library(ggplot2)

features_load_and_plot <- function() {
  # Load the .h5 file (have to have .h5 file saved locally for R to have accessibility)
  bboxdataPath <<- file.choose()
  bboxdata <<- H5Fopen(bboxdataPath)
  
  # Get group name (can find by running h5ls(bboxdata))
  groupls <- h5ls(bboxdata)
  bboxgroupname <- unique(groupls$group)
  
  # Provide the group index of the group which has the data
  bboxgroupname <- bboxgroupname[2]
  
  # Get unique tag (for purpose of naming the files, isn't necessary but helps distinguish different datasets)
  bboxdsetnm <<- readline("Enter unique tag (eg baseline): ")
  
  # Get list of dataset names within the group
  behavioral_metrics <<- names(h5read(bboxdataPath,bboxgroupname))
  
  #Data extraction and plotting
  for (i in 1:length(behavioral_metrics)) {
    # Assign each object name to it's corresponding dataset
    objectpath <- paste(bboxgroupname, behavioral_metrics[i], sep="/")
    assign(paste0(behavioral_metrics[i],"_", bboxdsetnm), h5read(bboxdataPath, objectpath))
    
    # Plot objects of class array
    obj <- get(paste0(behavioral_metrics[i],"_", bboxdsetnm))
    
    # Make timelist
    objframe_count <- length(obj)
    fps <- 45
    vid_length <- objframe_count / fps
    time_increments <<- 1 / fps
    timelist <- seq(from = 1 / fps, to = vid_length, by = time_increments)
    
    # Create data frame with obj and time increments
    obj_df <- data.frame(values = obj, times = timelist[1:objframe_count])
    
    # Create a dynamic name for the data frame
    df_name <- paste0(behavioral_metrics[i], "_", bboxdsetnm, "_df")
    assign(df_name, obj_df, envir = .GlobalEnv)
    
    if (is.factor(obj) || is.logical(obj)) {
      # Factor or Logical: Convert to numeric for plotting
      obj_df$values <- as.numeric(obj_df$values)
      
      # Generate the plot
      plot_title <- paste0(behavioral_metrics[i], "_", bboxdsetnm," Time Series Plot")
      plot_name <- paste0(behavioral_metrics[i], "_", bboxdsetnm, "_plot")
      p <- ggplot(obj_df, aes(x = times, y = values)) +
        geom_tile(aes(width = time_increments, height = Inf, fill = factor(values))) +
        scale_fill_manual(values = c("1" = "blue", "0" = "red"), na.value = "transparent") +
        labs(x = "Time (s)", y = behavioral_metrics[i], fill = "State") +
        ggtitle(plot_title) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
      
      # Assign the plot to a dynamic variable name
      assign(plot_name, p, envir = .GlobalEnv)
      
      # Clean up the data frame variable
      rm(df_name)
    }
    else if (is.array(obj)) {
      # Generate the plot
      plot_title <- paste0(behavioral_metrics[i], "_", bboxdsetnm, " Time Series Plot")
      plot_name <- paste0(behavioral_metrics[i], "_", bboxdsetnm, "_plot")  # Dynamic plot name
      
      print(str(obj_df))
      p <- ggplot(obj_df, aes(x = times, y = values)) +
        geom_area() +
        labs(x = "Time (s)", y = behavioral_metrics[i]) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(hjust = 0.5))
      
      # Assign the plot to a unique variable name
      assign(plot_name, p, envir = .GlobalEnv)
      
      # Clean up the data frame
      rm(df_name)
    }
  }
  
  # Clean up unnecessary variables
  rm(objectpath)
}

## USAGE:
# mouse_1 = extract_and_plot_h5() # Load in "feature.h5" file for a given mouse and give the dataset a tag/name
# mouse_2 = extract_and_plot_h5() # Add how many ever mice (or datasets) you would like

# Access desired behavioral feature data or plot (look up from the "Environment" tab)
# Example:
# ankle_distance_cont_1_df
# ankle_distance_cont_1_plot
