library(stringr)
library(rhdf5)
library(writexl)

# Set working directory to directory with all tracking files

### Load Tracking Data code:
tracking_load <- function(file_path, tag){
  pixels_to_cm <<- 0.03
  grid_size <<- 30
  
  # Open the HDF5 file and read tracking data
  assign(paste0(tag, "_tracking"), H5Fopen(file_path))
  tracking_data <- get(paste0(tag, "_tracking"))
  
  # Read coordinates data from the HDF5 file
  assign(paste0(tag, "_tracking_data"), as.data.frame(h5read(tracking_data, name="/df_with_missing/table", compoundAsDataFrame = FALSE)$values_block_0), envir = .GlobalEnv)
  tracking_coordinates <- get(paste0(tag, "_tracking_data"))
  
  # Extract column names
  tracking_names = h5readAttributes(tracking_data, name="/df_with_missing/table")$values_block_0_kind
  
  # Use stringr to extract body part names (as tracking names gives us a long, cryptic string that contains the body parts)
  pattern <- "(?<=\\nV)(.*?)(?=\\np)" # Looking for the whatever is in between nV and np (as that is where the body part names are in the pattern)
  all_matches <- str_extract_all(tracking_names, pattern)[[1]]
  exclude_indices <- c(2, 3, 4)  # Indices to exclude
  
  # Filter body part names
  body_parts <- all_matches[-exclude_indices]
  
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
  assign(paste0(tag, "_tracking_data"), tracking_coordinates, envir = .GlobalEnv)
  
  # Clean up
  rm(tracking_coordinates, tracking_data)
}

get_th_ratio <- function(tag, body_part) {
  # Setting the default subset size
  subset_size = 1:ncol(get(paste0(tag,"_tracking_data")))
  
  # Get the x and y row indices for the specified body part
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- get(paste0(tag, "_tracking_data"))[x_row_index,]
  x_row <- unlist(x_row)
  x_row <- x_row * pixels_to_cm
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- get(paste0(tag, "_tracking_data"))[y_row_index,]
  y_row <- unlist(y_row)
  y_row <- y_row * pixels_to_cm
  
  # Make coordinate data frame
  coord_df <- as.data.frame(data.frame(x = x_row, y = y_row)[subset_size,])

  cg_start <- grid_size/4
  cg_end <- 3*(grid_size/4)
  
  is_in_center <- with(coord_df, x >= cg_start & x <= cg_end & y >= cg_start & y <= cg_end)
  cg_count <- sum(is_in_center)
  outside_count <- sum(!is_in_center)
  
  # Calculate thigmotaxis ratio (handle division by zero)
  th_ratio <- ifelse(outside_count > 0, cg_count / outside_count, NA)
  
  return(th_ratio)
}

get_accel_metrics <- function(tag, body_part) {
  data <- get(paste0(tag, "_tracking_data"))
  row_names <- rownames(data)
  
  # Frame rate and time
  fps <- 45
  time_between_frames <- 1 / fps
  
  # Extract x and y coordinates
  x_row <- unlist(data[grep(paste0(body_part, " x"), row_names), ])
  y_row <- unlist(data[grep(paste0(body_part, " y"), row_names), ])
  
  # Compute distances, speeds, and accelerations
  distances <- sqrt(diff(x_row)^2 + diff(y_row)^2) * pixels_to_cm
  speeds <- distances / time_between_frames
  accels <- diff(speeds) / time_between_frames
  
  # Calculate metrics
  c(mean = mean(accels), sd = sd(accels), max = max(accels))
}

get_empty_bins_count <- function(tag, body_part) {
  
  # Setting the default subset size
  subset_size <- 1:ncol(get(paste0(tag, "_tracking_data")))
  
  # Get the x and y row indices for the specified body part
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- get(paste0(tag, "_tracking_data"))[x_row_index,]
  x_row <- unlist(x_row)
  x_row <- x_row * pixels_to_cm
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- get(paste0(tag, "_tracking_data"))[y_row_index,]
  y_row <- unlist(y_row)
  y_row <- y_row * pixels_to_cm
  
  # Make coordinate data frame
  coord_df <- data.frame(x = x_row, y = y_row)[subset_size, ]
  
  # Define grid edges for 50x50 bins over 30x30 area
  x_edges <- seq(0, 30, length.out = 51)
  y_edges <- seq(0, 30, length.out = 51)
  
  # Initialize a matrix to count occurrences in each bin
  bin_counts <- matrix(0, nrow = 50, ncol = 50)
  
  # Assign data points to bins
  for (i in seq_len(nrow(coord_df))) {
    x_bin <- findInterval(coord_df$x[i], x_edges, rightmost.closed = TRUE)
    y_bin <- findInterval(coord_df$y[i], y_edges, rightmost.closed = TRUE)
    if (x_bin > 0 && x_bin <= 50 && y_bin > 0 && y_bin <= 50) {
      bin_counts[x_bin, y_bin] <- bin_counts[x_bin, y_bin] + 1
    }
  }
  
  empty_bins_count <- sum(bin_counts == 0)
  
  return(empty_bins_count)
}

get_speed_avg <- function(tag, body_part) {
  # Initialize speed and distance vectors
  assign(paste0(tag, "_distances"), numeric())
  distances <- get(paste0(tag, "_distances"))
  
  assign(paste0(tag, "_speeds"), numeric())
  speeds <- get(paste0(tag, "_speeds"))
  
  # Set frame rate and calculate video length
  fps <- 45
  time_between_frames <- 1/fps
  vid_length <- time_between_frames * ncol(get(paste0(tag, "_tracking_data")))
  
  # Get the x and y row indices for the specified body part
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- get(paste0(tag, "_tracking_data"))[x_row_index,]
  x_row <- unlist(x_row)
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- get(paste0(tag, "_tracking_data"))[y_row_index,]
  y_row <- unlist(y_row)
  
  # Calculate distances between consecutive points
  for (i in 2:ncol(get(paste0(tag, "_tracking_data")))) {
    distances[i-1] <- sqrt((x_row[i] - x_row[i-1])^2 + (y_row[i] - y_row[i-1])^2)
  }
  distances <- unlist(distances)
  distances <- distances * pixels_to_cm
  
  # Compute speeds as distances per frame
  for (i in 1:length(distances)) {
    speeds[i] <- distances[i]/time_between_frames
  }
  
  # Average the speeds together
  return(mean(speeds))
}

### Process Files
files <- list.files()

# Time points and dataframe
time_pts <- c("Baseline", "1 day post SCI", "7 days post SCI", "14 days post SCI", "21 days post SCI", "28 days post SCI", "35 days post SCI", "42 days post SCI")

accel_metrics <- data.frame(tag = character(), mean = numeric(), sd = numeric(), max = numeric(), stringsAsFactors = FALSE)
th_ratios <- c()
empty_bins_counts <- c()
speed_avgs <- c()
fthp_spd_avgs <- c()

k <- 0
for (i in seq_along(files)) {
  if ((i - 1) %% 16 == 0) k <- k + 1
  
  # Load tracking data
  file_path <- file.path(getwd(), files[i])
  mouse_num <- str_extract(files[i], "(?<=_)[^_]*(?=-tracking)")
  tag <- paste(time_pts[k], mouse_num, sep = "_")
  
  tracking_load(file_path, tag)
  
  # Compute accel metrics
  accel_metrics <- rbind(accel_metrics, c(tag, get_accel_metrics(tag, "tailbase")))
  
  # Compute thigmotaxis ratio
  th_ratios <- c(th_ratios, get_th_ratio(tag, "tailbase"))
  
  # Compute empty bins counts
  empty_bins_counts <- c(empty_bins_counts, get_empty_bins_count(tag, "tailbase"))

  # Compute speed averages
  speed_avgs <- c(speed_avgs, get_speed_avg(tag, "tailbase"))
  
  # Compute fthp speed ratio
  rfpaw_speed_avg <- get_speed_avg(tag, "rfpaw")
  lfpaw_speed_avg <- get_speed_avg(tag, "lfpaw")
  rhpaw_speed_avg <- get_speed_avg(tag, "rhpaw")
  lhpaw_speed_avg <- get_speed_avg(tag, "lhpaw")
  
  favg <- mean(rfpaw_speed_avg, lfpaw_speed_avg)
  havg <- mean(rhpaw_speed_avg, lhpaw_speed_avg)
  
  fthp_spd_avgs <- c(fthp_spd_avgs, favg/havg)
}

tracking_metrics <- accel_metrics
tracking_metrics <- cbind(tracking_metrics, th_ratios)
tracking_metrics <- cbind(tracking_metrics, empty_bins_counts)
tracking_metrics <- cbind(tracking_metrics, speed_avgs)
tracking_metrics <- cbind(tracking_metrics, fthp_spd_avgs)

colnames(tracking_metrics) <- c(
  "tag", "acceleration_mean",
  "acceleration_sd",
  "acceleration_max",
  "thigmotaxis_ratio",
  "empty_bins_counts",
  "speed_avgs",
  "fthp_spd_avgs"
  )

tracking_metrics_orig <- tracking_metrics

index_map <- c(9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121)
orig_index_map <- c(17, 33, 49, 65, 81, 97, 113, 9, 25, 41, 57, 73, 89, 105, 121)

for (i in seq_along(index_map)) {
  tracking_metrics[index_map[i]:(index_map[i] + 7), ] <- 
    tracking_metrics_orig[orig_index_map[i]:(orig_index_map[i] + 7), ]
}

write_xlsx(tracking_metrics, "tracking_metrics.xlsx", format_headers = T)
