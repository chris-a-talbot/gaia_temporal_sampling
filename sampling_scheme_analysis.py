import os
import numpy as np
import pandas as pd
import tskit
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import mean_squared_error, mean_absolute_error
import re
import argparse
import logging
import sys
from datetime import datetime
import socket
from concurrent_log_handler import ConcurrentRotatingFileHandler


def setup_logging(cwd=None):
    """Set up logging configuration for both file and console output."""
    # Create logs directory structure
    log_dir = Path(cwd) / "logs" / "processing_logs" / "sampling_scheme_analysis_sh_logs" if cwd else Path(
        "logs/processing_logs/sampling_scheme_analysis_sh_logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    # Create a timestamp for the log file
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    hostname = socket.gethostname()
    pid = os.getpid()
    log_file = log_dir / f"analysis_{timestamp}_{hostname}_{pid}.log"

    # Create logger
    logger = logging.getLogger('sampling_scheme_analysis')
    logger.setLevel(logging.DEBUG)

    # Create handlers
    # Console handler with INFO level
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    # File handler with DEBUG level - using ConcurrentRotatingFileHandler for cluster safety
    file_handler = ConcurrentRotatingFileHandler(
        str(log_file),
        maxBytes=10 * 1024 * 1024,  # 10MB
        backupCount=5,
        encoding='utf-8'
    )
    file_handler.setLevel(logging.DEBUG)

    # Create formatters and add them to the handlers
    file_formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(processName)s-%(process)d | %(threadName)s | '
        '%(filename)s:%(lineno)d | %(message)s'
    )
    console_formatter = logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s')

    file_handler.setFormatter(file_formatter)
    console_handler.setFormatter(console_formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # Log initial information
    logger.info(f"Starting analysis on host: {hostname}")
    logger.info(f"Process ID: {pid}")
    logger.info(f"Log file location: {log_file}")

    return logger


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze location inference accuracy across different sampling schemes.')
    parser.add_argument('--cwd', type=str, default=None,
                        help='Working directory containing the required folder structure. '
                             'Expected: {CWD}/trees/simplified, {CWD}/inferred_locations/locations, '
                             'and outputs will be saved to {CWD}/accuracy_analysis')
    parser.add_argument('--ntrees', type=int, default=None,
                        help='Number of trees to randomly select for processing. '
                             'If not specified, all matched trees will be processed.')
    return parser.parse_args()


def setup_paths(cwd=None):
    """Set up directory paths based on optional CWD argument."""
    logger = logging.getLogger('sampling_scheme_analysis')

    if cwd is not None:
        base_path = Path(cwd)
        tree_dir = base_path / "trees" / "simplified"
        loc_dir = base_path / "inferred_locations" / "locations"
        output_dir = base_path / "accuracy_analysis"
        metrics_dir = output_dir / "individual_metrics"  # New directory for individual file metrics
        logger.info(f"Using provided working directory: {cwd}")
    else:
        tree_dir = Path("/home/christ/PycharmProjects/gaia_temporal_sampling/gaia_temporal_testing/trees/simplified")
        loc_dir = Path(
            "/home/christ/PycharmProjects/gaia_temporal_sampling/gaia_temporal_testing/inferred_locations/locations")
        output_dir = Path("/home/christ/PycharmProjects/gaia_temporal_sampling/gaia_temporal_testing/accuracy_analysis")
        metrics_dir = output_dir / "independent_metrics"
        logger.warning("No working directory provided, using default paths")

    # Create output directories if they don't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_dir.mkdir(parents=True, exist_ok=True)

    logger.debug(f"Tree directory: {tree_dir}")
    logger.debug(f"Location directory: {loc_dir}")
    logger.debug(f"Output directory: {output_dir}")
    logger.debug(f"Individual metrics directory: {metrics_dir}")

    return tree_dir, loc_dir, output_dir, metrics_dir


def load_and_match_files(tree_dir, loc_dir, ntrees=None):
    """Match tree files with their corresponding location files."""
    logger = logging.getLogger('sampling_scheme_analysis')

    tree_files = {f.stem: f for f in Path(tree_dir).glob("*.trees")}
    loc_files = {f.stem.replace("_locations", ""): f for f in Path(loc_dir).glob("*_locations.csv")}

    logger.info(f"Found {len(tree_files)} tree files and {len(loc_files)} location files")

    # Match files and return paired paths
    matched_files = []
    for name in set(tree_files.keys()) & set(loc_files.keys()):
        matched_files.append((tree_files[name], loc_files[name]))

    total_matches = len(matched_files)
    logger.info(f"Successfully matched {total_matches} file pairs")

    if total_matches == 0:
        logger.error("No matching file pairs found!")
        return matched_files

    # If ntrees is specified and valid, randomly select that many trees
    if ntrees is not None:
        if ntrees > total_matches:
            logger.warning(f"Requested {ntrees} trees but only {total_matches} pairs available. Using all pairs.")
        else:
            logger.info(f"Randomly selecting {ntrees} trees from {total_matches} available pairs")
            # Convert to array for random selection
            indices = np.arange(total_matches)
            selected_indices = np.random.choice(indices, size=min(ntrees, total_matches), replace=False)
            matched_files = [matched_files[i] for i in selected_indices]
            logger.info(f"Selected {len(matched_files)} trees for processing")

    return matched_files


def parse_filename(filename):
    """
    Parse tree filename to extract sampling scheme information.
    Handles both standard format: 'simplified_tree-{prefix}_{a|f}_c{c}_nc{nc}_{q}q'
    and simplified format: 'simplified_tree-{prefix}_{a|f}_c{c}'
    """
    logger = logging.getLogger('sampling_scheme_analysis')

    # Try standard format first
    standard_pattern = r"simplified_tree-(.+?)_(a|f)_c(\d+)_nc(\d+)_(\d+)q"
    match = re.match(standard_pattern, filename)

    if match:
        prefix, adj_type, c_samples, nc_samples, quarters = match.groups()
        info = {
            'prefix': prefix,
            'adjustment_type': adj_type,
            'contemporary_samples': int(c_samples),
            'noncontemporary_samples': int(nc_samples),
            'quarters': int(quarters),
            'total_samples': int(c_samples) + int(nc_samples)
        }
        logger.debug(f"Successfully parsed filename (standard format): {filename} -> {info}")
        return info

    # Try simplified format
    simplified_pattern = r"simplified_tree-(.+?)_(a|f)_c(\d+)"
    match = re.match(simplified_pattern, filename)

    if match:
        prefix, adj_type, c_samples = match.groups()
        info = {
            'prefix': prefix,
            'adjustment_type': adj_type,
            'contemporary_samples': int(c_samples),
            'noncontemporary_samples': 0,
            'quarters': 0,
            'total_samples': int(c_samples)
        }
        logger.debug(f"Successfully parsed filename (simplified format): {filename} -> {info}")
        return info

    logger.error(f"Failed to parse filename: {filename}")
    return None


def calculate_accuracy_metrics(true_locations, inferred_locations):
    """Calculate various accuracy metrics for location inference."""
    logger = logging.getLogger('sampling_scheme_analysis')

    try:
        metrics = {}
        # Calculate metrics for x and y coordinates separately
        for coord in ['x', 'y']:
            true_coord = true_locations[:, 0 if coord == 'x' else 1]
            inferred_coord = inferred_locations[:, 0 if coord == 'x' else 1]

            mse = mean_squared_error(true_coord, inferred_coord)
            rmse = np.sqrt(mse)
            mae = mean_absolute_error(true_coord, inferred_coord)

            try:
                mape = np.mean(np.abs((true_coord - inferred_coord) / true_coord)) * 100
            except RuntimeWarning:
                logger.warning(f"MAPE calculation generated warnings for {coord} coordinate")
                mape = np.nan

            metrics[f'{coord}_mse'] = mse
            metrics[f'{coord}_rmse'] = rmse
            metrics[f'{coord}_mae'] = mae
            metrics[f'{coord}_mape'] = mape

        # Calculate combined x,y metrics
        combined_mse = mean_squared_error(true_locations, inferred_locations)
        combined_rmse = np.sqrt(combined_mse)
        combined_mae = mean_absolute_error(true_locations.ravel(), inferred_locations.ravel())

        metrics.update({
            'combined_mse': combined_mse,
            'combined_rmse': combined_rmse,
            'combined_mae': combined_mae
        })

        logger.debug(f"Calculated metrics: {metrics}")
        return metrics

    except Exception as e:
        logger.error(f"Error calculating metrics: {str(e)}", exc_info=True)
        raise


def get_unprocessed_files(matched_files, metrics_dir):
    """Filter out already processed files."""
    logger = logging.getLogger('sampling_scheme_analysis')

    unprocessed = []
    for tree_path, loc_path in matched_files:
        metrics_file = metrics_dir / f"{tree_path.stem}_metrics.csv"
        if not metrics_file.exists():
            unprocessed.append((tree_path, loc_path))

    logger.info(f"Found {len(unprocessed)} unprocessed files out of {len(matched_files)} total files")
    return unprocessed


def save_individual_metrics(metrics, tree_path, metrics_dir):
    """Save metrics for an individual tree to CSV."""
    logger = logging.getLogger('sampling_scheme_analysis')

    output_path = metrics_dir / f"{tree_path.stem}_metrics.csv"
    try:
        pd.DataFrame([metrics]).to_csv(output_path, index=False)
        logger.debug(f"Saved individual metrics to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save metrics for {tree_path.stem}: {str(e)}")
        return False


def load_all_metrics(metrics_dir):
    """Load all individual metrics files from the metrics directory."""
    logger = logging.getLogger('sampling_scheme_analysis')

    all_metrics = []
    for metrics_file in metrics_dir.glob("*_metrics.csv"):
        try:
            metrics = pd.read_csv(metrics_file)
            all_metrics.append(metrics)
        except Exception as e:
            logger.error(f"Failed to load metrics from {metrics_file}: {str(e)}")

    if all_metrics:
        combined = pd.concat(all_metrics, ignore_index=True)
        logger.info(f"Loaded metrics from {len(all_metrics)} files")
        return combined
    else:
        logger.warning("No metrics files found to load")
        return pd.DataFrame()


def get_node_locations(ts):
    """Extract x,y locations for each node in a tree sequence."""
    logger = logging.getLogger('sampling_scheme_analysis')

    try:
        # Initialize array to store locations
        locations = np.full((ts.num_nodes, 2), np.nan)

        # Iterate through all nodes
        for node in ts.nodes():
            # Check if node has an associated individual
            if node.individual != -1:
                # Get individual object
                ind = ts.individual(node.individual)
                # Extract x,y coordinates
                if ind.location is not None and len(ind.location) >= 2:
                    locations[node.id] = [ind.location[0], ind.location[1]]

        logger.debug(f"Extracted locations for {np.sum(~np.isnan(locations[:, 0]))} nodes")
        return locations

    except Exception as e:
        logger.error(f"Error extracting node locations: {str(e)}", exc_info=True)
        raise


def analyze_tree_locations(tree_path, loc_path):
    """Analyze accuracy of inferred locations for a single tree."""
    logger = logging.getLogger('sampling_scheme_analysis')

    try:
        logger.info(f"Analyzing tree: {tree_path.name}")

        # Load tree and get true locations
        ts = tskit.load(str(tree_path))
        true_locations = get_node_locations(ts)

        # Load inferred locations
        inferred_locations = pd.read_csv(loc_path).values
        logger.debug(f"Loaded {len(inferred_locations)} inferred locations")

        # Get sample node IDs
        sample_nodes = ts.samples()
        logger.debug(f"Found {len(sample_nodes)} sample nodes")

        # Remove sample nodes from analysis
        non_sample_mask = ~np.isin(np.arange(len(true_locations)), sample_nodes)
        true_locations_filtered = true_locations[non_sample_mask]
        inferred_locations_filtered = inferred_locations[non_sample_mask]
        logger.debug(f"Analyzing {len(true_locations_filtered)} non-sample nodes")

        # Calculate accuracy metrics
        metrics = calculate_accuracy_metrics(true_locations_filtered, inferred_locations_filtered)

        # Parse filename info
        filename_info = parse_filename(tree_path.stem)

        # Combine both dictionaries
        if filename_info is not None:
            combined_metrics = {**metrics, **filename_info}
        else:
            combined_metrics = metrics

        logger.info(f"Completed analysis for {tree_path.name}")
        return combined_metrics

    except Exception as e:
        logger.error(f"Error analyzing tree {tree_path.name}: {str(e)}", exc_info=True)
        raise


def aggregate_results(results):
    """Aggregate results across trees with identical sampling schemes."""
    logger = logging.getLogger('sampling_scheme_analysis')

    try:
        df = pd.DataFrame(results)
        logger.debug(f"Aggregating results for {len(results)} trees")

        # Group by sampling scheme parameters
        grouping_cols = ['contemporary_samples', 'noncontemporary_samples', 'quarters', 'total_samples']
        aggregated = df.groupby(grouping_cols).agg({
            'x_mse': ['mean', 'std'],
            'y_mse': ['mean', 'std'],
            'combined_mse': ['mean', 'std'],
            'x_rmse': ['mean', 'std'],
            'y_rmse': ['mean', 'std'],
            'combined_rmse': ['mean', 'std'],
            'x_mae': ['mean', 'std'],
            'y_mae': ['mean', 'std'],
            'combined_mae': ['mean', 'std']
        }).reset_index()

        logger.info(f"Generated aggregated results for {len(aggregated)} unique sampling schemes")
        return aggregated

    except Exception as e:
        logger.error(f"Error aggregating results: {str(e)}", exc_info=True)
        raise


def rank_sampling_schemes(aggregated_results):
    """Rank sampling schemes based on combined RMSE."""
    logger = logging.getLogger('sampling_scheme_analysis')

    try:
        rankings = aggregated_results.sort_values(('combined_rmse', 'mean'))
        rankings['rank'] = range(1, len(rankings) + 1)
        logger.info(f"Generated rankings for {len(rankings)} sampling schemes")
        return rankings

    except Exception as e:
        logger.error(f"Error ranking sampling schemes: {str(e)}", exc_info=True)
        raise


def main():
    try:
        # Parse command line arguments
        args = parse_args()

        # Set up logging
        logger = setup_logging(args.cwd)
        logger.info("Starting sampling scheme analysis")

        # Set random seed for reproducibility
        np.random.seed(42)
        logger.debug("Set random seed to 42 for reproducible tree selection")

        # Set up directory paths
        tree_dir, loc_dir, output_dir, metrics_dir = setup_paths(args.cwd)

        # Get matched files
        matched_files = load_and_match_files(tree_dir, loc_dir, args.ntrees)

        # Filter out already processed files
        unprocessed_files = get_unprocessed_files(matched_files, metrics_dir)

        # Process unprocessed files
        for tree_path, loc_path in unprocessed_files:
            try:
                metrics = analyze_tree_locations(tree_path, loc_path)
                save_individual_metrics(metrics, tree_path, metrics_dir)
            except Exception as e:
                logger.error(f"Failed to analyze {tree_path.name}: {str(e)}")
                continue

        # Load all metrics for aggregation
        logger.info("Loading all metrics for aggregation")
        all_results = load_all_metrics(metrics_dir)

        if not all_results.empty:
            # Aggregate results
            logger.info("Aggregating results")
            aggregated_results = aggregate_results(all_results.to_dict('records'))

            # Rank sampling schemes
            logger.info("Ranking sampling schemes")
            rankings = rank_sampling_schemes(aggregated_results)

            # Save aggregated results
            output_files = {
                'accuracy_analysis.csv': aggregated_results,
                'sampling_scheme_rankings.csv': rankings
            }

            for filename, data in output_files.items():
                output_path = output_dir / filename
                try:
                    data.to_csv(output_path)
                    logger.info(f"Successfully saved {filename} to {output_path}")
                except Exception as e:
                    logger.error(f"Failed to save {filename}: {str(e)}")

            logger.info("Analysis completed successfully")
            return aggregated_results, rankings
        else:
            logger.warning("No results to aggregate")
            return None, None

    except Exception as e:
        logger.critical(f"Critical error in main execution: {str(e)}", exc_info=True)
        raise


if __name__ == "__main__":
    main()