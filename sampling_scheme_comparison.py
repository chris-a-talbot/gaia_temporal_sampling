import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_and_clean_data(filepath):
    # Read CSV, skipping the first row (sub-headers)
    df = pd.read_csv(filepath, header=0)

    # Drop the unnamed index column
    df = df.drop(df.columns[0], axis=1)

    # Drop the row with std/mean subheaders
    df = df.drop(0, axis=0)

    # Convert numeric columns to float
    numeric_columns = df.columns.drop(['quarters'])
    df[numeric_columns] = df[numeric_columns].astype(float)

    # Calculate contemporary/noncontemporary ratio
    df['noncontem_over_contem'] = df['noncontemporary_samples'] / df['contemporary_samples']

    return df


def create_focused_heatmap(df):
    # Define parameters and metrics with new labels
    sampling_params = ['contemporary_samples', 'noncontemporary_samples',
                       'total_samples', 'quarters', 'noncontem_over_contem']

    # Define the new labels for x-axis
    param_labels = {
        'contemporary_samples': '# contemporary samples',
        'noncontemporary_samples': '# noncontemporary samples',
        'total_samples': '# samples, total',
        'quarters': '% generations sampled',
        'noncontem_over_contem': '# noncontemporary / # contemporary'
    }

    accuracy_metrics = ['combined_mse', 'combined_rmse', 'combined_mae']

    # Calculate correlations between sampling parameters and accuracy metrics only
    correlation_matrix = pd.DataFrame()
    for param in sampling_params:
        correlations = []
        for metric in accuracy_metrics:
            corr = df[param].corr(df[metric])
            correlations.append(corr)
        correlation_matrix[param_labels[param]] = correlations

    correlation_matrix.index = accuracy_metrics

    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(correlation_matrix,
                annot=True,  # Show numbers
                fmt='.2f',  # Round to 2 decimal places
                cmap='coolwarm',
                center=0,  # Center colormap at 0
                vmin=-1,  # Set minimum correlation value
                vmax=1)  # Set maximum correlation value

    plt.title('Correlation between Sampling Parameters and Accuracy Metrics for Gaia Inference (All Sampling Schemes)')
    plt.xlabel('Sampling Parameters')
    plt.ylabel('Accuracy Metrics')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    return plt.gcf()


def create_filtered_heatmap(df, contemporary_filter=250):
    # Filter the dataframe
    filtered_df = df[df['contemporary_samples'] == contemporary_filter].copy()

    if filtered_df.empty:
        raise ValueError(f"No data found for contemporary_samples = {contemporary_filter}")

    # Define parameters and metrics with new labels
    sampling_params = ['noncontemporary_samples', 'total_samples', 'quarters', 'noncontem_over_contem']

    # Define the new labels for x-axis
    param_labels = {
        'noncontemporary_samples': '# noncontemporary samples',
        'total_samples': '# samples, total',
        'quarters': '% generations sampled',
        'noncontem_over_contem': '# noncontemporary / # contemporary'
    }

    accuracy_metrics = ['combined_mse', 'combined_rmse', 'combined_mae']

    # Calculate correlations between sampling parameters and accuracy metrics only
    correlation_matrix = pd.DataFrame()
    for param in sampling_params:
        correlations = []
        for metric in accuracy_metrics:
            corr = filtered_df[param].corr(filtered_df[metric])
            correlations.append(corr)
        correlation_matrix[param_labels[param]] = correlations

    correlation_matrix.index = accuracy_metrics

    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(correlation_matrix,
                annot=True,  # Show numbers
                fmt='.2f',  # Round to 2 decimal places
                cmap='coolwarm',
                center=0,  # Center colormap at 0
                vmin=-1,  # Set minimum correlation value
                vmax=1)  # Set maximum correlation value

    plt.title(f'Correlation between Sampling Parameters and Accuracy Metrics for Gaia Inference (Contemporary Samples = {contemporary_filter})')
    plt.xlabel('Sampling Parameters')
    plt.ylabel('Accuracy Metrics')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    return plt.gcf()


def create_filtered_heatmap_excl(df, total_filter=250):
    # Filter the dataframe
    filtered_df = df[df['total_samples'] == total_filter].copy()

    if filtered_df.empty:
        raise ValueError(f"No data found for total_samples = {total_filter}")

    # Define parameters and metrics with new labels
    sampling_params = ['noncontemporary_samples',
                       'quarters', 'noncontem_over_contem']

    # Define the new labels for x-axis
    param_labels = {
        'noncontemporary_samples': '# noncontemporary samples',
        'quarters': '% generations sampled',
        'noncontem_over_contem': '# noncontemporary / # contemporary'
    }

    accuracy_metrics = ['combined_mse', 'combined_rmse', 'combined_mae']

    # Calculate correlations between sampling parameters and accuracy metrics only
    correlation_matrix = pd.DataFrame()
    for param in sampling_params:
        correlations = []
        for metric in accuracy_metrics:
            corr = filtered_df[param].corr(filtered_df[metric])
            correlations.append(corr)
        correlation_matrix[param_labels[param]] = correlations

    correlation_matrix.index = accuracy_metrics

    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(correlation_matrix,
                annot=True,  # Show numbers
                fmt='.2f',  # Round to 2 decimal places
                cmap='coolwarm',
                center=0,  # Center colormap at 0
                vmin=-1,  # Set minimum correlation value
                vmax=1)  # Set maximum correlation value

    plt.title(f'Correlation between Sampling Parameters and Accuracy Metrics for Gaia Inference (Total Samples = {total_filter})')
    plt.xlabel('Sampling Parameters')
    plt.ylabel('Accuracy Metrics')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    return plt.gcf()


def main():
    # Load and process the data
    filepath = '/home/christ/PycharmProjects/gaia_temporal_sampling/gaia_temporal_testing/accuracy_analysis/accuracy_analysis.csv'
    df = load_and_clean_data(filepath)

    # Create focused heatmap
    heatmap_fig = create_focused_heatmap(df)
    heatmap_fig_2 = create_filtered_heatmap(df, 250)
    heatmap_fig_3 = create_filtered_heatmap_excl(df, 250)
    plt.show()


if __name__ == "__main__":
    main()