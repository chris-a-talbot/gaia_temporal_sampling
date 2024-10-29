# gaia temporal sampling

A set of scripts for testing the effects of various temporal sampling schemes on the accuracy of gaia ancestor location inference on the cluster.

## Workflow:

- Run the SLiM simulation up to 20 times to generate unique tree sequence output (using slim.sh and run-pareto.slim)
- Simplify the output tree sequences to fit various sampling schemes (using simplify_trees.sh/.py)
- Create sample location matrices from the simplified trees (using sample_loc_mat_from_trees.sh/.py)
- Run gaia on the simplified trees with sample locations (using run_gaia.sh/.R)
- Analyze the accuracy of gaia inference (using sampling_scheme_analysis.py)
- Visualize the results across sampling schemes (using sampling_scheme_comparison.py)
