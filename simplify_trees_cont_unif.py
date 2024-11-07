import tskit
import numpy as np
from datetime import datetime
import itertools
import os
import pandas as pd
import glob
import argparse
import sys
import logging
from pathlib import Path
import multiprocessing as mp
from functools import partial
import threading
from queue import Empty


# Add environment check
def check_environment():
    """Verify all required packages are available."""
    required_packages = ['tskit', 'numpy', 'pandas']
    missing_packages = []
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)

    if missing_packages:
        raise ImportError(f"Missing required packages: {', '.join(missing_packages)}")


class MPLogHandler(logging.Handler):
    """A logging handler that puts logs into a multiprocessing queue."""

    def __init__(self, queue):
        super().__init__()
        self.queue = queue

    def emit(self, record):
        self.queue.put_nowait(record)


def logger_thread(queue, log_file):
    """Separate thread to handle logging from multiple processes."""
    file_handler = logging.FileHandler(log_file)
    console_handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    while True:
        try:
            record = queue.get()
            if record is None:  # Sentinel value to stop the thread
                break
            file_handler.emit(record)
            console_handler.emit(record)
        except Empty:
            continue
        except Exception as e:
            print(f"Error in logger thread: {e}")
            break

    file_handler.close()
    console_handler.close()


def setup_parallel_logging(log_dir, job_id, task_id):
    """Configure logging for parallel processing with SLURM support."""
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'tree_processing_job{job_id}_task{task_id}_{timestamp}.log')

    # Set up queue and thread as before
    queue = mp.Queue()
    threading.Thread(target=logger_thread, args=(queue, log_file)).start()

    logger = logging.getLogger(__name__)
    logger.addHandler(MPLogHandler(queue))
    logger.setLevel(logging.INFO)

    # Add basic system information to log
    logger.info(f"Running on host: {os.uname().nodename}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Number of CPUs: {mp.cpu_count()}")

    return logger, queue


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process tree sequences with various sampling schemes.')

    parser.add_argument('--job-id', type=str, default=os.getenv('SLURM_JOB_ID', 'local'),
                        help='SLURM job ID (automatically detected if running under SLURM)')
    parser.add_argument('--task-id', type=str, default=os.getenv('SLURM_ARRAY_TASK_ID', '0'),
                        help='SLURM array task ID (automatically detected if running under SLURM)')
    parser.add_argument('--trees-dir', type=str, default='./simulation_data/trees',
                        help='Directory containing tree files')
    parser.add_argument('--logs-dir', type=str, default='./simulation_data/logs',
                        help='Directory containing log files')
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed for reproducibility')
    parser.add_argument('--contemporary-n-max', type=int, default=250,
                        help='Maximum number of contemporary samples')
    parser.add_argument('--noncontemporary-n-max', type=int, default=100,
                        help='Maximum number of non-contemporary samples')
    parser.add_argument('--ploidy', type=int, default=2,
                        help='Ploidy of the organisms')
    parser.add_argument('--tree-numbers', type=int, nargs='+',
                        help='Optional: Only process trees with these numbers')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Optional: Custom output directory for simplified trees')
    parser.add_argument('--num-processes', type=int, default=None,
                        help='Number of parallel processes to use (default: number of CPU cores)')
    parser.add_argument('--chunk-size', type=int, default=1,
                        help='Number of trees to process per chunk in parallel processing')

    return parser.parse_args()


# Naming scheme mappings remain the same
NAMING_KEY = {
    "fixed": "f",
    "adjusted": "a",
    "contemporaries": "c",
    "non_contemporaries": "nc",
    "last_quarter_percent": "025p",  # 0.25%
    "last_half_percent": "050p",     # 0.50%
    "last_three_quarters_percent": "075p",  # 0.75%
    "last_percent": "100p"           # 1.00%
}

def get_percentage_generations(total_gens):
    """Calculate generation ranges for different percentage spans."""
    spans = {}

    # Calculate the number of generations for each percentage
    one_percent_gens = max(1, int(round(total_gens * 0.01)))  # At least 1 generation
    zero_point_75_gens = max(1, int(round(total_gens * 0.0075)))
    zero_point_50_gens = max(1, int(round(total_gens * 0.005)))
    zero_point_25_gens = max(1, int(round(total_gens * 0.0025)))

    # Calculate the ranges for each span
    spans = {
        'last_quarter_percent': range(total_gens - zero_point_25_gens + 1, total_gens + 1),
        'last_half_percent': range(total_gens - zero_point_50_gens + 1, total_gens + 1),
        'last_three_quarters_percent': range(total_gens - zero_point_75_gens + 1, total_gens + 1),
        'last_percent': range(total_gens - one_percent_gens + 1, total_gens + 1)
    }

    return spans


def read_pedigree_log(log_path):
    """Read and parse the pedigree log file, removing duplicate sections if present."""
    try:
        # Read the raw file first
        df = pd.read_csv(log_path)

        # Check for duplicate headers within the data
        header_rows = df[df['cycle'] == 'cycle'].index

        if len(header_rows) > 0:
            # If we found duplicate headers, keep only the most recent section
            last_header = header_rows[-1]
            # Get the latest section (from after the last header to the end)
            df = df[last_header + 1:].reset_index(drop=True)

            # Write the cleaned data to a new file
            output_path = log_path.rsplit('.', 1)[0] + '_trimmed.' + log_path.rsplit('.', 1)[1]
            df.to_csv(output_path, index=False)
            logging.info(f"Cleaned file saved to: {output_path}")

        # Convert the pedigree_IDs column to lists of integers
        df['pedigree_IDs'] = df['pedigree_IDs'].apply(lambda x: [int(i) for i in x.split(',')])

        # Create the generation to IDs dictionary
        gen_to_ids = dict(zip(df['cycle'], df['pedigree_IDs']))
        total_gens = df['cycle'].max()

        return gen_to_ids, total_gens

    except Exception as e:
        logging.error(f"Error reading pedigree log {log_path}: {str(e)}")
        raise


def get_quarter_generations(total_gens):
    """Divide generations into quarters and return ranges for different spans."""
    quarters = np.array_split(range(total_gens + 1), 4)
    spans = {
        'last_quarter': quarters[-1],
        'last_half': np.concatenate([quarters[-2], quarters[-1]]),
        'last_three_quarters': np.concatenate([quarters[-3], quarters[-2], quarters[-1]])
    }
    return spans



def process_tree_file(tree_path, args, process_seed):
    """Process a single tree file and its corresponding log file."""
    # Set a unique random seed for this process
    np.random.seed(args.random_seed + process_seed)

    logger = logging.getLogger(__name__)
    logger.info(f"Processing tree file: {tree_path}")

    try:
        tree_name = Path(tree_path).stem
        log_path = os.path.join(args.logs_dir, f"{tree_name}.txt")

        if not os.path.exists(log_path):
            logger.error(f"No corresponding log file found for {tree_path}")
            return

        ts = tskit.load(tree_path)
        gen_to_ids, total_gens = read_pedigree_log(log_path)

        individuals_times = ts.individual_times
        contemporary_individuals_ids = np.where(individuals_times == 0)[0]
        contemporary_individuals_sample_ids = np.random.choice(
            contemporary_individuals_ids,
            size=min(len(contemporary_individuals_ids), args.contemporary_n_max),
            replace=False
        )

        percentage_spans = get_percentage_generations(total_gens)
        sampling_schemes = {}

        sampling_scheme_options = {
            'contemporary_fixed': [True, False],
            'n_non_contemporaries': [2, 10, args.noncontemporary_n_max],
            'non_contemporary_span': ['last_quarter_percent', 'last_half_percent',
                                    'last_three_quarters_percent', 'last_percent']
        }

        pedigree_to_ind = {
            individual.metadata['pedigree_id']: individual.id
            for individual in ts.individuals()
            if hasattr(individual, 'metadata') and 'pedigree_id' in individual.metadata
        }

        # Create sampling schemes
        for contemporary_fixed, n_non_contemporaries, span in itertools.product(
                sampling_scheme_options['contemporary_fixed'],
                sampling_scheme_options['n_non_contemporaries'],
                sampling_scheme_options['non_contemporary_span']
        ):
            # Create shortened scheme name
            prefix = f"{NAMING_KEY['fixed']}" if contemporary_fixed else f"{NAMING_KEY['adjusted']}"
            if contemporary_fixed:
                c_number = 250
            else:
                c_number = args.contemporary_n_max - n_non_contemporaries

            scheme_name = f"{prefix}_c{c_number}_nc{n_non_contemporaries}_{NAMING_KEY[span]}"

            # Select contemporary samples
            if contemporary_fixed:
                contemporaries = contemporary_individuals_sample_ids
            else:
                contemporaries = np.random.choice(contemporary_individuals_sample_ids,
                                                size=c_number, replace=False)

            # Get generations for the specified span
            span_generations = percentage_spans[span]

            # Pool all pedigree IDs from the relevant generations
            pedigree_pool = []
            for gen in span_generations:
                if gen in gen_to_ids:
                    pedigree_pool.extend(gen_to_ids[gen])

            # Sample non-contemporary individuals by pedigree ID
            if len(pedigree_pool) < n_non_contemporaries:
                logger.warning(f"Warning: Not enough non-contemporary individuals for scheme '{scheme_name}'")
                sampled_pedigree_ids = pedigree_pool
            else:
                sampled_pedigree_ids = np.random.choice(pedigree_pool, size=n_non_contemporaries, replace=False)

            # Get the individual IDs for our sampled pedigree IDs
            non_contemporary_individuals = [
                pedigree_to_ind[ped_id]
                for ped_id in sampled_pedigree_ids
                if ped_id in pedigree_to_ind
            ]

            # Combine samples
            sample_ids = np.concatenate([contemporaries, non_contemporary_individuals])
            sample_ids = np.unique(sample_ids)
            sampling_schemes[scheme_name] = sample_ids.tolist()

            logger.info(f"Scheme '{scheme_name}' created with {len(sample_ids)} individuals.")

        del pedigree_to_ind, gen_to_ids, total_gens, individuals_times, contemporary_individuals_ids

        # Add contemporary-only scheme
        scheme_name = "f_c250"
        sampling_schemes[scheme_name] = contemporary_individuals_sample_ids.tolist()

        # Define output directory with process-specific subdirectory to avoid conflicts
        if args.output_dir:
            output_dir = os.path.join(args.output_dir, "simplified")
        else:
            output_dir = os.path.join(os.path.dirname(tree_path), "simplified")

        os.makedirs(output_dir, exist_ok=True)

        for scheme_name, individual_ids in sampling_schemes.items():
            logger.info(f"Simplifying tree sequence for scheme '{scheme_name}'...")

            selected_node_ids = []
            for ind_id in individual_ids:
                individual = ts.individual(ind_id)
                if individual is None:
                    logger.warning(f"Warning: Individual ID {ind_id} not found in tree sequence.")
                    continue
                selected_node_ids.extend(individual.nodes)

            selected_node_ids = [node_id for node_id in selected_node_ids if node_id is not None]

            if not selected_node_ids:
                logger.warning(f"Warning: No valid node IDs found for scheme '{scheme_name}'. Skipping simplification.")
                continue

            simplified_ts = ts.simplify(samples=selected_node_ids, keep_unary_in_individuals=True)

            output_path = os.path.join(output_dir, f"simplified_{tree_name}_{scheme_name}.trees")
            simplified_ts.dump(output_path)

            logger.info(f"Saved simplified tree sequence to {output_path}")


    except Exception as e:
        logger.error(f"Error processing tree file {tree_path}: {str(e)}")
        return False


def main():
    """Main function to process tree files in parallel."""

    # Check environment first
    check_environment()

    args = parse_arguments()

    # Set up parallel logging with SLURM job information
    log_dir = os.path.join(args.logs_dir, 'processing_logs')
    logger, log_queue = setup_parallel_logging(log_dir, args.job_id, args.task_id)

    try:
        # Get tree files
        if args.tree_numbers:
            tree_files = []
            for number in args.tree_numbers:
                pattern = f"*-{number}.trees"
                matches = glob.glob(os.path.join(args.trees_dir, pattern))
                tree_files.extend(matches)
        else:
            tree_files = glob.glob(os.path.join(args.trees_dir, "*.trees"))

        if not tree_files:
            logger.warning("No tree files found matching the criteria.")
            return

        # Use SLURM_CPUS_PER_TASK if available, otherwise use specified or default
        num_processes = int(os.getenv('SLURM_CPUS_PER_TASK', args.num_processes or mp.cpu_count()))
        logger.info(f"Using {num_processes} processes to handle {len(tree_files)} tree files")

        # Create process pool and process files in parallel
        with mp.Pool(processes=num_processes) as pool:
            # Create the list of arguments for each file
            work_items = [(tree_file, args, idx) for idx, tree_file in enumerate(tree_files)]

            # Process files directly with starmap
            results = pool.starmap(
                process_tree_file,
                work_items,
                chunksize=args.chunk_size
            )

        # Report results
        successful = sum(1 for r in results if r)
        failed = len(tree_files) - successful
        logger.info(f"Processing complete. Successfully processed {successful} files, {failed} failed.")

    except Exception as e:
        logger.error(f"Fatal error in main processing: {str(e)}")
        sys.exit(1)
    finally:
        # Stop the logger thread
        log_queue.put_nowait(None)


if __name__ == "__main__":
    mp.freeze_support()  # Required for Windows compatibility
    main()
