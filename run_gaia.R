#!/usr/bin/env Rscript

# Basic logging function to use before logging package is installed
basic_log <- function(level, message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

# Function to get number of SLURM cores
get_slurm_cores <- function() {
  # Check SLURM environment variables
  cpus_per_task <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (cpus_per_task != "") {
    return(as.integer(cpus_per_task))
  }
  
  # If SLURM_CPUS_PER_TASK is not set, try SLURM_CPUS_ON_NODE
  cpus_on_node <- Sys.getenv("SLURM_CPUS_ON_NODE")
  if (cpus_on_node != "") {
    return(as.integer(cpus_on_node))
  }
  
  # If no SLURM variables are set, default to 1 core
  return(1)
}

# Function to ensure a directory exists
create_dir_if_missing <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    basic_log("INFO", sprintf("Created directory: %s", dir_path))
  }
}

# Function to ensure required packages are installed
ensure_packages <- function() {
  required_packages <- c("data.table", "logging", "future", "future.apply", "remotes")
  
  # First, handle base packages needed for installation
  for (pkg in c("utils", "tools")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      basic_log("ERROR", sprintf("Base R package %s not available!", pkg))
      stop(sprintf("Base R package %s not available!", pkg))
    }
  }
  
  # Install required packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      basic_log("INFO", sprintf("Installing package: %s", pkg))
      utils::install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  
  # Install gaia if not present
  if (!requireNamespace("gaia", quietly = TRUE)) {
    basic_log("INFO", "Installing gaia package from GitHub")
    remotes::install_github("blueraleigh/gaia")
  }
  
  # Load all required packages
  for (pkg in c(required_packages, "gaia")) {
    library(pkg, character.only = TRUE)
  }
}

# Function to check if output files already exist
check_outputs_exist <- function(tree_file, mpr_dir, locations_dir) {
  base_name <- tools::file_path_sans_ext(basename(tree_file))
  mpr_output_file <- file.path(mpr_dir, paste0(base_name, "_mpr.csv"))
  locations_output_file <- file.path(locations_dir, paste0(base_name, "_locations.csv"))
  
  # Check if both output files exist and are not empty
  mpr_exists <- file.exists(mpr_output_file) && file.size(mpr_output_file) > 0
  locations_exist <- file.exists(locations_output_file) && file.size(locations_output_file) > 0
  
  return(list(
    skip = mpr_exists && locations_exist,
    mpr_file = mpr_output_file,
    locations_file = locations_output_file
  ))
}

# Main processing function
process_tree_pair <- function(tree_file, sample_file, mpr_dir, locations_dir) {
  tryCatch({
    # Check if outputs already exist
    output_status <- check_outputs_exist(tree_file, mpr_dir, locations_dir)
    
    if (output_status$skip) {
      logging::loginfo(sprintf("Skipping existing outputs for tree: %s", basename(tree_file)))
      return(TRUE)
    }
    
    logging::loginfo(sprintf("Processing tree pair: %s and %s", tree_file, sample_file))
    
    # Read tree sequence and sample locations
    tree <- gaia::treeseq_load(tree_file)
    sample_locations <- as.matrix(data.table::fread(sample_file))
    
    # Run MPR analysis
    mpr_output <- gaia::treeseq_quadratic_mpr(tree, sample_locations, TRUE)
    
    # Save MPR output
    saveRDS(mpr_output, output_status$mpr_file)
    logging::loginfo(sprintf("Saved MPR output to: %s", output_status$mpr_file))
    
    # Run minimize and save locations
    locations_output <- gaia::treeseq_quadratic_mpr_minimize(mpr_output)
    data.table::fwrite(as.data.table(locations_output), output_status$locations_file)
    logging::loginfo(sprintf("Saved locations output to: %s", output_status$locations_file))
    
    return(TRUE)
  }, error = function(e) {
    logging::logerror(sprintf("Error processing %s: %s", tree_file, e$message))
    return(FALSE)
  })
}

# Main execution function
main <- function() {
  # Get CWD from command line args or use current working directory
  args <- commandArgs(trailingOnly = TRUE)
  CWD <- if (length(args) > 0) args[1] else getwd()
  
  # Create log directory first
  log_dir <- file.path(CWD, "logs", "processing_logs")
  create_dir_if_missing(log_dir)
  
  # Create log file with timestamp to avoid conflicts
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(log_dir, sprintf("run_gaia_sh_logs_%s.txt", timestamp))
  
  basic_log("INFO", sprintf("Starting processing with CWD: %s", CWD))
  basic_log("INFO", sprintf("Logs will be written to: %s", log_file))
  
  # Ensure required packages are installed
  ensure_packages()
  
  # Now that logging is installed, set up proper logging
  logging::basicConfig(level = "DEBUG")  # Changed to DEBUG for more verbose logging
  logging::addHandler(logging::writeToFile, file = log_file, append = TRUE)
  
  # Create output directories
  logging::loginfo("Gathering and creating necessary directories...")
  mpr_dir <- file.path(CWD, "inferred_locations", "mpr")
  locations_dir <- file.path(CWD, "inferred_locations", "locations")
  create_dir_if_missing(mpr_dir)
  create_dir_if_missing(locations_dir)
  logging::loginfo("All directories present!")
  
  # Get list of input files
  trees_dir <- file.path(CWD, "trees", "simplified")
  samples_dir <- file.path(CWD, "sample_locations")
  
  tree_files <- list.files(trees_dir, pattern = "\\.trees$", full.names = TRUE)
  
  if (length(tree_files) == 0) {
    logging::logerror("No tree files found")
    return(FALSE)
  }
  
  logging::loginfo(sprintf("Found %d tree files to process", length(tree_files)))
  
  # Create corresponding sample file paths
  sample_files <- file.path(
    samples_dir,
    paste0(tools::file_path_sans_ext(basename(tree_files)), ".csv")
  )
  
  # Verify all sample files exist
  missing_samples <- !file.exists(sample_files)
  if (any(missing_samples)) {
    logging::logerror(sprintf(
      "Missing sample files for trees: %s",
      paste(basename(tree_files[missing_samples]), collapse = ", ")
    ))
    return(FALSE)
  }
  
  # Set up parallel processing
  num_cores <- get_slurm_cores()
  num_cores <- max(1, num_cores)  # Ensure at least 1 core
  
  logging::loginfo(sprintf("Processing using %d cores", num_cores))
  logging::loginfo("Setting up parallel processing with multicore...")
  
  # Set up future to use multicore instead of multisession
  tryCatch({
    future::plan(future::multicore, workers = num_cores)
    logging::loginfo("Successfully initialized multicore parallel processing")
  }, error = function(e) {
    logging::logerror(sprintf("Failed to initialize parallel processing: %s", e$message))
    logging::loginfo("Falling back to sequential processing")
    future::plan(future::sequential)
  })
  
  # Add some debug logging
  logging::logdebug(sprintf("Memory available: %s", system("free -h", intern = TRUE)))
  logging::logdebug(sprintf("Current working directory: %s", getwd()))
  
  # Process files in parallel with progress logging
  logging::loginfo("Starting parallel processing of files...")
  results <- vector("list", length(tree_files))
  
  for (i in seq_along(tree_files)) {
    logging::loginfo(sprintf("Submitting job %d of %d", i, length(tree_files)))
    results[[i]] <- future::future({
      process_tree_pair(tree_files[i], sample_files[i], mpr_dir, locations_dir)
    })
  }
  
  # Resolve futures with progress logging
  for (i in seq_along(results)) {
    logging::loginfo(sprintf("Resolving job %d of %d", i, length(results)))
    results[[i]] <- future::value(results[[i]])
  }
  
  # Clean up parallel backend
  future::plan(future::sequential)
  logging::loginfo("Parallel processing completed")
  
  # Check results
  success_count <- sum(unlist(results))
  total_count <- length(results)
  skipped_count <- sum(sapply(tree_files, function(tf) check_outputs_exist(tf, mpr_dir, locations_dir)$skip))
  
  logging::loginfo(sprintf(
    "Processing complete. Successfully processed %d out of %d files (%d files skipped)",
    success_count, total_count, skipped_count
  ))
  
  return(TRUE)
}

# Execute main function
if (!interactive()) {
  main()
}