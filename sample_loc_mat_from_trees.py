import os
import tskit
import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial


def process_tree_file(tree_file_path, output_dir):
    """
    Process a single .trees file and create corresponding CSV with sample locations.

    Args:
        tree_file_path (str): Path to the .trees file
        output_dir (str): Directory to save the output CSV
    """
    try:
        # Load the tree sequence
        ts = tskit.load(tree_file_path)

        # Create list to store sample location data
        sample_locations = []

        # Get all sample nodes
        samples = ts.samples()

        # For each sample, get its associated individual and location
        for sample_id in samples:
            # Get the individual associated with this node
            individual_id = ts.node(sample_id).individual
            if individual_id is not None:
                # Get the location for this individual
                individual = ts.individual(individual_id)
                if len(individual.location) >= 2:
                    x, y = individual.location[0], individual.location[1]
                    sample_locations.append({
                        'node_id': sample_id,
                        'x': x,
                        'y': y
                    })

        # Create DataFrame and save to CSV
        df = pd.DataFrame(sample_locations)

        # Create output filename based on input filename
        tree_name = os.path.splitext(os.path.basename(tree_file_path))[0]
        output_path = os.path.join(output_dir, f"{tree_name}.csv")

        # Save to CSV without index
        df.to_csv(output_path, index=False)
        return f"Processed {tree_name} successfully"
    except Exception as e:
        return f"Error processing {os.path.basename(tree_file_path)}: {str(e)}"


def main():
    # Input and output directories
    input_dir = "./trees/simplified"
    output_dir = "./sample_locations"

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Get list of all .trees files
    tree_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)
                  if f.endswith('.trees')]

    # Determine number of processes to use (leave one core free for system)
    num_processes = max(1, cpu_count() - 1)

    # Create partial function with fixed output_dir
    process_func = partial(process_tree_file, output_dir=output_dir)

    # Process files in parallel
    with Pool(processes=num_processes) as pool:
        # Use tqdm if available for progress bar
        try:
            from tqdm import tqdm
            results = list(tqdm(
                pool.imap(process_func, tree_files),
                total=len(tree_files),
                desc="Processing tree files"
            ))
        except ImportError:
            results = pool.map(process_func, tree_files)

    # Print results
    for result in results:
        print(result)


if __name__ == "__main__":
    main()