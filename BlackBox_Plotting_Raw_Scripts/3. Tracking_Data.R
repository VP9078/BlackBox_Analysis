#####################################
#### 3. Tracking data processing ####
#####################################

# Load required libraries
library(ggplot2)
library(rhdf5)
library(stringr)


#### Run Functions (defined below)
# load_tracking_data() # loads in the data
# body_parts # plots body parts available for tracking

# plot_likelihood("control", "hip") # mouse centroid

#plot_speed("control", "snout")

#### 1. Load input tracking ####

pixel_to_cm <- 0.03

# Function to load and process tracking data
tracking_load <- function(){
  file_path = file.choose()  # Prompt user to select file
  
  tag <- readline("Give a unique tag: ")  # Prompt user to provide a tag
  
  # Make a tag dataset to save all tag datasets
  if (exists("tag_memory", envir = .GlobalEnv)) {
    tag_memory[length(tag_memory) + 1] <<- tag
  } else {
    tag_memory <<- c(tag)
  }
  
  
  # Open the HDF5 file and read tracking data
  assign(paste0(tag, "_tracking_data"), H5Fopen(file_path), envir = .GlobalEnv)
  tracking_data <- get(paste0(tag, "_tracking_data"))
  
  # Read coordinates data from the HDF5 file
  assign(paste0(tag, "_tracking_coordinates"), as.data.frame(h5read(tracking_data, name="/df_with_missing/table", compoundAsDataFrame = FALSE)$values_block_0), envir = .GlobalEnv)
  tracking_coordinates <- get(paste0(tag, "_tracking_coordinates"))
  
  # Extract column names
  tracking_names = h5readAttributes(tracking_data, name="/df_with_missing/table")$values_block_0_kind

  # Use stringr to extract body part names (as tracking names gives us a long, cryptic string that contains the body parts)
  pattern <- "(?<=\\nV)(.*?)(?=\\np)" # Looking for the whatever is in between nV and np (as that is where the body part names are in the patter)
  all_matches <- str_extract_all(tracking_names, pattern)[[1]]
  exclude_indices <- c(2, 3, 4)  # Indices to exclude
  
  # Filter body part names
  body_parts <<- all_matches[-exclude_indices]
  
  # Initialize row names list
  assign(paste0(tag, "_row_names"), list(), envir = .GlobalEnv)
  
  # Create row names with body parts and corresponding attributes
  assign(paste0(tag, "_row_names"),
         unlist(lapply(body_parts, function(part) {
           c(paste(part, "x", sep = " "),
             paste(part, "y", sep = " "),
             paste(part, "likelihood", sep = " "))
         }), use.names = FALSE), envir = .GlobalEnv)
  
  # Clean up temporary variables
  rm(all_matches, exclude_indices, pattern, tracking_names)
  
  # Set row names to tracking coordinates
  rownames(tracking_coordinates) <- get(paste0(tag, "_row_names"))
  assign(paste0(tag, "_tracking_coordinates"), tracking_coordinates, envir = .GlobalEnv)
  
  # Clean up
  rm(tracking_coordinates, tracking_data)
}

###################
#### Plotting ####
##################

# REQUIRES: load_tracking_data() to be used at least once

# Function to compute and plot speed for a specific body part
plot_speed <- function (tag, body_part, subset_size = "all") {
  
  # Setting the default subset size
  if (any(grepl("all", subset_size)))
    subset_size = 1:ncol(get(paste0(tag,"_tracking_coordinates")))
  
  # Initialize speed and distance vectors
  assign(paste0(tag, "_distances"), numeric())
  distances <- get(paste0(tag, "_distances"))
  
  assign(paste0(tag, "_speeds"), numeric())
  speeds <- get(paste0(tag, "_speeds"))
  
  # Set frame rate and calculate video length
  fps <- 45
  time_between_frames <- 1/fps
  vid_length <- time_between_frames * ncol(get(paste0(tag, "_tracking_coordinates")))
  
  # Get the x and y row indices for the specified body part
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- get(paste0(tag, "_tracking_coordinates"))[x_row_index,]
  x_row <- unlist(x_row)
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- get(paste0(tag, "_tracking_coordinates"))[y_row_index,]
  y_row <- unlist(y_row)
  
  # Calculate distances between consecutive points
  for (i in 2:ncol(get(paste0(tag, "_tracking_coordinates")))) {
    distances[i-1] <- sqrt((x_row[i] - x_row[i-1])^2 + (y_row[i] - y_row[i-1])^2)
  }
  distances <- unlist(distances)
  distances <- distances * pixel_to_cm
  
  # Compute speeds as distances per frame
  for (i in 1:length(distances)) {
    speeds[i] <- distances[i]/time_between_frames
  }
  assign(paste0(tag, "_",body_part, "_speed_values"), speeds, envir = .GlobalEnv)
  
  # Create a time vector for plotting
  time_vector <- seq(from = time_between_frames, to = (vid_length - time_between_frames), by = time_between_frames)
  
  # Combine time vector and speed values into a dataframe
  assign(paste0(tag, "_speed_df"), data.frame(time_vector, get(paste0(tag, "_",body_part, "_speed_values"))), envir = .GlobalEnv)
  speed_df <- get(paste0(tag, "_speed_df"))[subset_size,]
  
  # Generate speed plot
  assign(paste0(tag, "_", body_part, "_speed_plot"),
         ggplot(speed_df, aes(x = time_vector, y = (get(paste0(tag, "_",body_part, "_speed_values"))[subset_size]))) +
           geom_line(color = "blue") +
           labs(
             title = str_to_title(paste0(tag, " ", body_part, " Speed")),
             x = "Time (s)",
             y = "Speed (cm/s)"
           ) +
           theme_minimal() +
           theme(
             plot.title = element_text(face = "bold", hjust = 0.5)  # Center the title and make it bold
           )
         , envir = .GlobalEnv)
  
  # Return the speed plot
  get(paste0(tag, "_",body_part, "_speed_plot"))
}

plot_likelihood <- function(tag, body_part, subset_size = "all") {
  
  # Setting the default subset size
  if (any(grepl("all", subset_size)))
    subset_size = 1:ncol(get(paste0(tag,"_tracking_coordinates")))
  
  # Set frame rate and calculate video length
  fps <- 45
  time_between_frames <- 1/fps
  vid_length <- time_between_frames * ncol(get(paste0(tag, "_tracking_coordinates")))
  
  ## Likelihood plotting
  likelihood_row_index <- which(grepl("likelihood", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  likelihood_row <- get(paste0(tag, "_tracking_coordinates"))[likelihood_row_index,]
  likelihood_row <- unlist(likelihood_row)
  
  # Normalize time vector:
  time_vector <- seq(from = time_between_frames, to = (vid_length), by = time_between_frames)
  
  # Make likelihood vs time dataframe
  assign(paste0(tag, "_", body_part,"_likelihood_df"), data.frame(time_vector, likelihood_row), envir = .GlobalEnv)
  likelihood_df <- get(paste0(tag, "_", body_part,"_likelihood_df"))[subset_size, ]
  
  # Plot likelihood
  assign(paste0(tag, "_", body_part, "_likelihood_plot"),
         ggplot(likelihood_df, aes(x = time_vector, y = likelihood_row)) +
           geom_line(color = "blue") +
           labs(
             title = str_to_title(paste0(tag, " ", body_part, " Likelihood")),
             x = "Time (s)",
             y = "Likelihood"
           ) +
           theme_minimal() +
           theme(
             plot.title = element_text(face = "bold", hjust = 0.5)  # Center the title and make it bold
           ) +
           coord_cartesian(ylim = c(0, max(likelihood_row)))
         , envir = .GlobalEnv)
  
  # Return the likelihood plot
  get(paste0(tag, "_",body_part, "_likelihood_plot"))
}


# Plotting Trajectory function (subset_size should be an integer between 1 and the amount of frames of your video -> (ncol(get(paste0(tag,"_tracking_coordinates")))))
plot_trajectory <- function (tag, body_part, subset_size = "all", points = FALSE, grid_size = "full") {
  
  # Setting the default subset size
  if (any(grepl("all", subset_size)))
    subset_size = 1:ncol(get(paste0(tag,"_tracking_coordinates")))
  
  # Get the x and y row indices for the specified body part
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- get(paste0(tag, "_tracking_coordinates"))[x_row_index,]
  x_row <- unlist(x_row)
  x_row <- x_row * pixel_to_cm
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- get(paste0(tag, "_tracking_coordinates"))[y_row_index,]
  y_row <- unlist(y_row)
  y_row <- y_row * pixel_to_cm
  
  # Make trajectory data frame
  trajectory_df <- data.frame(x_row, y_row)[subset_size,]
  
  # Make an index column in the dataframe
  trajectory_df$index <- 1:nrow(trajectory_df)
  
  # Set grid size
  if (grid_size == "full")
    grid_size <- max(x_row)
  
  # Generate trajectory plot
  trajectory_plot <-
    ggplot(trajectory_df, aes(x = x_row, y = y_row)) +
    geom_path(aes(color = index)) +
    scale_color_gradientn(colors = c("blue", "green", "yellow", "red")) +
    labs(
      title = str_to_title(paste0(tag, " ", body_part, " trajectory")),
      x = "Box Length (cm)",
      y = "Box Width (cm)"
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, grid_size), ylim = c(0, grid_size)) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)  # Center the title and make it bold
    )
  
  if (points == TRUE)
    trajectory_plot <- trajectory_plot +
    geom_point(aes(color = index))
  
  assign(paste0(tag, "_", body_part, "_trajectory_plot"), trajectory_plot, envir = .GlobalEnv)
  
  
  # Return the trajectory plot
  get(paste0(tag, "_",body_part, "_trajectory_plot"))
}


# Plot heatmap
plot_heatmap <- function (tag, body_part, subset_size = "all", grid_size = "full") {
  
  # Setting the default subset size
  if (any(grepl("all", subset_size)))
    subset_size = 1:ncol(get(paste0(tag,"_tracking_coordinates")))
  
  # Get the x and y row indices for the specified body part
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- get(paste0(tag, "_tracking_coordinates"))[x_row_index,]
  x_row <- unlist(x_row)
  x_row <- x_row * pixel_to_cm
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- get(paste0(tag, "_tracking_coordinates"))[y_row_index,]
  y_row <- unlist(y_row)
  y_row <- y_row * pixel_to_cm
  
  # Set grid size
  if (grid_size == "full")
    grid_size <- 50
  
  # Make trajectory data frame
  trajectory_df <- as.matrix(data.frame(x_row, y_row)[subset_size,])
  
  heatmap <- ggplot(trajectory_df, aes(x = x_row, y = y_row)) +
    geom_bin_2d(bins = 50) +
    scale_fill_gradient(low = "blue", high = "red", name = "Density") +
    theme_minimal() +
    coord_cartesian(xlim = c(0, grid_size), ylim = c(0, grid_size)) +
    labs(title = "Mouse Cage Occupancy Heatmap",
         x = "Box Length (cm)",
         y = "Box Width (cm)") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)  # Center the title and make it bold
    )
  
  print(heatmap)
  
  assign(paste0(tag, "_", body_part, "_heatmap_plot"), heatmap, envir = .GlobalEnv)
}