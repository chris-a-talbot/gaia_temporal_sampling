import tskit
import numpy as np
from datetime import datetime
import itertools
import os
import pandas as pd
import glob
import logging
from pathlib import Path
import multiprocessing as mp
from functools import partial


def setup_logging():
    """Configure basic logging."""
    os.makedirs('logs', exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = f'logs/tree_processing_{timestamp}.log'

    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


# Naming scheme mappings
NAMING_KEY = {
    "fixed": "f",
    "adjusted": "a",
    "contemporaries": "c",
    "non_contemporaries": "nc",
    "last_quarter_percent": "025p",
    "last_half_percent": "050p",
    "last_three_quarters_percent": "075p",
    "last_percent": "100p"
}


def get_percentage_generations(total_gens):
    """Calculate generation ranges for different percentage spans."""
    spans = {}

    one_percent_gens = max(1, int(round(total_gens * 0.01)))
    zero_point_75_gens = max(1, int(round(total_gens * 0.0075)))
    zero_point_50_gens = max(1, int(round(total_gens * 0.005)))
    zero_point_25_gens = max(1, int(round(total_gens * 0.0025)))

    spans = {
        'last_quarter_percent': range(total_gens - zero_point_25_gens + 1, total_gens + 1),
        'last_half_percent': range(total_gens - zero_point_50_gens + 1, total_gens + 1),
        'last_three_quarters_percent': range(total_gens - zero_point_75_gens + 1, total_gens + 1),
        'last_percent': range(total_gens - one_percent_gens + 1, total_gens + 1)
    }

    return spans


def read_pedigree_log(log_path):
    """Read and parse the pedigree log file."""
    try:
        df = pd.read_csv(log_path)
        df['pedigree_IDs'] = df['pedigree_IDs'].apply(lambda x: [int(i) for i in x.split(',')])
        gen_to_ids = dict(zip(df['cycle'], df['pedigree_IDs']))
        total_gens = df['cycle'].max()
        return gen_to_ids, total_gens
    except Exception as e:
        logging.error(f"Error reading pedigree log {log_path}: {str(e)}")
        raise


def process_tree_file(tree_path, contemporary_n_max=250, noncontemporary_n_max=100,
                      random_seed=42, process_seed=0):
    """Process a single tree file and its corresponding log file."""
    logger = logging.getLogger(__name__)

    # Set unique random seed for this process
    np.random.seed(random_seed + process_seed)

    try:
        tree_name = Path(tree_path).stem
        log_path = f"./gaia_temporal_testing/logs/SLiM_logs/{tree_name}.txt"

        if not os.path.exists(log_path):
            logger.error(f"No corresponding log file found for {tree_path}")
            return False

        ts = tskit.load(tree_path)
        gen_to_ids, total_gens = read_pedigree_log(log_path)

        individuals_times = ts.individual_times
        contemporary_individuals_ids = np.where(individuals_times == 0)[0]
        contemporary_individuals_sample_ids = np.random.choice(
            contemporary_individuals_ids,
            size=min(len(contemporary_individuals_ids), contemporary_n_max),
            replace=False
        )

        percentage_spans = get_percentage_generations(total_gens)
        sampling_schemes = {}

        sampling_scheme_options = {
            'contemporary_fixed': [True, False],
            'n_non_contemporaries': [2, 10, noncontemporary_n_max],
            'non_contemporary_span': ['last_quarter_percent', 'last_half_percent',
                                      'last_three_quarters_percent', 'last_percent']
        }

        pedigree_to_ind = {
            individual.metadata['pedigree_id']: individual.id
            for individual in ts.individuals()
            if hasattr(individual, 'metadata') and 'pedigree_id' in individual.metadata
        }

        # Create sampling schemes
        for contemporary_fixed, n_non_contemporaries, span in itertools.product(*sampling_scheme_options.values()):
            prefix = f"{NAMING_KEY['fixed']}" if contemporary_fixed else f"{NAMING_KEY['adjusted']}"
            c_number = 250 if contemporary_fixed else (contemporary_n_max - n_non_contemporaries)

            scheme_name = f"{prefix}_c{c_number}_nc{n_non_contemporaries}_{NAMING_KEY[span]}"

            if contemporary_fixed:
                contemporaries = contemporary_individuals_sample_ids
            else:
                contemporaries = np.random.choice(contemporary_individuals_sample_ids,
                                                  size=c_number, replace=False)

            span_generations = percentage_spans[span]
            pedigree_pool = []
            for gen in span_generations:
                if gen in gen_to_ids:
                    pedigree_pool.extend(gen_to_ids[gen])

            if len(pedigree_pool) < n_non_contemporaries:
                logger.warning(f"Not enough non-contemporary individuals for scheme '{scheme_name}'")
                sampled_pedigree_ids = pedigree_pool
            else:
                sampled_pedigree_ids = np.random.choice(pedigree_pool, size=n_non_contemporaries, replace=False)

            non_contemporary_individuals = [
                pedigree_to_ind[ped_id]
                for ped_id in sampled_pedigree_ids
                if ped_id in pedigree_to_ind
            ]

            sample_ids = np.concatenate([contemporaries, non_contemporary_individuals])
            sample_ids = np.unique(sample_ids)
            sampling_schemes[scheme_name] = sample_ids.tolist()

            logger.info(f"Created scheme '{scheme_name}' with {len(sample_ids)} individuals")

        # Add contemporary-only scheme
        scheme_name = "f_c250"
        sampling_schemes[scheme_name] = contemporary_individuals_sample_ids.tolist()

        # Create simplified trees
        output_dir = os.path.join('output', 'simplified')
        os.makedirs(output_dir, exist_ok=True)

        for scheme_name, individual_ids in sampling_schemes.items():
            selected_node_ids = []
            for ind_id in individual_ids:
                individual = ts.individual(ind_id)
                if individual is None:
                    continue
                selected_node_ids.extend(individual.nodes)

            selected_node_ids = [node_id for node_id in selected_node_ids if node_id is not None]

            if not selected_node_ids:
                logger.warning(f"No valid node IDs found for scheme '{scheme_name}'. Skipping.")
                continue

            simplified_ts = ts.simplify(samples=selected_node_ids, keep_unary_in_individuals=True)

            output_path = os.path.join(output_dir, f"simplified_{tree_name}_{scheme_name}.trees")
            simplified_ts.dump(output_path)

            logger.info(f"Saved simplified tree sequence to {output_path}")

        return True

    except Exception as e:
        logger.error(f"Error processing tree file {tree_path}: {str(e)}")
        return False


def main():
    """Main function to process tree files using local multiprocessing."""
    logger = setup_logging()

    # Configuration
    trees_dir = "./gaia_temporal_testing/trees"  # Directory containing your tree files
    num_cores = max(1, mp.cpu_count() // 2)  # Use half of available cores

    try:
        # Get all tree files
        tree_files = glob.glob(os.path.join(trees_dir, "*.trees"))

        if not tree_files:
            logger.warning("No tree files found in the specified directory.")
            return

        logger.info(f"Found {len(tree_files)} tree files to process")
        logger.info(f"Using {num_cores} CPU cores")

        # Process files in parallel
        with mp.Pool(processes=num_cores) as pool:
            work_items = [(tree_file, 250, 100, 42, idx)
                          for idx, tree_file in enumerate(tree_files)]
            results = pool.starmap(process_tree_file, work_items)

        # Report results
        successful = sum(1 for r in results if r)
        failed = len(tree_files) - successful
        logger.info(f"Processing complete. Successfully processed {successful} files, {failed} failed.")

    except Exception as e:
        logger.error(f"Fatal error in main processing: {str(e)}")
        raise


if __name__ == "__main__":
    main()