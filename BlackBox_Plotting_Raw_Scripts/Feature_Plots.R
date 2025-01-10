#####################################################
#### 2. Plot behavioral feature time-series data ####
#####################################################

# Load packages
# Can install rhdf5 by install BiocManager, then running: BiocManager::install("rhdf5")
library(rhdf5)
library(ggplot2)
library(manipulate)

#### A. Enter basic functions

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



# Function to create checkbox plot and store it globally
checkbox_front_paws_luminance_plot <- function() {
  
  # Convert values to numeric
  both_front_paws_lifted_baseline_df$values <- as.numeric(both_front_paws_lifted_baseline_df$values)
  
  # Create a plot object and store it globally
  global_plot <<- NULL  # Initialize a global variable to store the plot
  
  manipulate(
    {
      # Create the base plot
      front_paws_luminance_baseline_plot <- ggplot(both_front_paws_lifted_baseline_df, aes(x = times, y = values)) +
        guides(fill = "none") +
        coord_cartesian(ylim = range(front_right_luminance_baseline_df$values)) +
        labs(title = "Front Paws Luminance Plot Baseline", x = "Time", y = "Luminance") +
        theme(plot.title = element_text(hjust = 0.5))
      
      # Add layers based on checkbox selections
      if (do_tile) {
        front_paws_luminance_baseline_plot <- front_paws_luminance_baseline_plot + 
          geom_tile(aes(width = time_increments, height = Inf, fill = factor(values))) + 
          scale_fill_manual(values = c("0" = "transparent", "1" = "transparent"), na.value = "green")
      }
      
      if (do_right_line) {
        front_paws_luminance_baseline_plot <- front_paws_luminance_baseline_plot + 
          front_right_luminance_baseline_plot
      }
      
      if (do_left_line) {
        front_paws_luminance_baseline_plot <- front_paws_luminance_baseline_plot + 
          geom_line(data = front_left_luminance_baseline_df, aes(x = times, y = values), color = "blue")
      }
      
      if (do_left_point) {
        front_paws_luminance_baseline_plot <- front_paws_luminance_baseline_plot + 
          geom_point(data = front_left_luminance_baseline_df, aes(x = times, y = values), color = "blue")
      }
      
      if (do_right_point) {
        front_paws_luminance_baseline_plot <- front_paws_luminance_baseline_plot + 
          geom_point(data = front_right_luminance_baseline_df, aes(x = times, y = values), color = "red")
      }
      
      # Store the plot globally
      global_plot <<- front_paws_luminance_baseline_plot
      
      # Print the plot
      print(global_plot)
    },
    
    # Define checkboxes for user interaction
    do_tile = checkbox("Add both paws lifted tile", initial = FALSE),
    do_right_line = checkbox("Add front right luminance line", initial = FALSE),
    do_left_line = checkbox("Add front left luminance line", initial = FALSE),
    do_right_point = checkbox("Add front right luminance point", initial = FALSE),
    do_left_point = checkbox("Add front left luminance point", initial = FALSE)
  ) 
} 




#### B. Load in data

# mouse_1 = extract_and_plot_h5() # Load in "feature.h5" file for a given mouse and give the dataset a tag/name
# mouse_2 = extract_and_plot_h5() # Add how many ever mice (or datasets) you would like

# Access desired behavioral feature data or plot (look up from the "Environment" tab)
# Example:
# ankle_distance_cont_1_df
# ankle_distance_cont_1_plot


#### C. Overlay Plots

# Define a function to create overlay plots (Need to have ran the extract_and_plot_h5 function at least twice)
# tag1 and tag2 are the unique tags you assigned your datasets within the extract and plot function (in the line bboxdsetnm <- readline("Enter unique tag (eg baseline): "))

overlay_features_plots <- function(tag1, tag2) {
  
  # Get data frame names from global environment
  all_objects <- ls(envir = .GlobalEnv)
  
  # Create a list of data frame names for data set 1
  df1_list <- all_objects[grepl(tag1, all_objects) & grepl("_df", all_objects)]
  
  # Remove specific entries from the df1 list that are not relevant
  df1_list <- df1_list[!df1_list %in% c("df1_list", paste0("fps_", tag1, "_df"), paste0("frame_count_", tag1, "_df"))]
  
  # Create a list of data frame names for data set 2
  df2_list <- all_objects[grepl(tag2, all_objects) & grepl("_df", all_objects)]
  
  # Remove specific entries from the df2 list that are not relevant
  df2_list <- df2_list[!df2_list %in% c("df2_list", paste0("fps_", tag2, "_df"), paste0("frame_count_", tag2, "_df"))]
  
  # Loop through each object name in the df1 list
  for (i in 1:length(df1_list)) {
    # Check that both data sets are not NA and exist
    if (!is.na(df1_list[i]) && !is.na(df2_list[i])) {
      
      # Retrieve the data frame object by its name
      df1_obj <- get(df1_list[i], envir = .GlobalEnv)
      df2_obj <- get(df2_list[i], envir = .GlobalEnv)
      
      # Create a title for the plot based on the object name
      plot_title <- paste0(df1_list[i], " ", tag1, " & ", tag2, " Time Series Overlay Plot")
      
      # Create a ggplot object for the data set 1
      p <- ggplot(df1_obj, aes(x = times, y = values)) +
        geom_line(color = "red") +  # Add a red line for data set 1 with some transparency
        geom_line(data = df2_obj,  # Overlay a blue line for data set 2 with some transparency
                  aes(x = times, y = values), color = "blue") +
        labs(x = "Time (s)", y = "Value") +  # Label the x-axis as 'Time (s)' and y-axis with the object name
        ggtitle(paste0(df1_list[i], " ", tag1, " vs ", tag2)) +  # Add a title to the plot
        theme(plot.title = element_text(hjust = 0.5))  # Center the plot title horizontally
      
      # Assign the ggplot object to a variable with a dynamic name
      assign(paste0(df1_list[i], "_", tag1, "_", tag2, "_overlay"), p, envir = .GlobalEnv)
    }
  }
}