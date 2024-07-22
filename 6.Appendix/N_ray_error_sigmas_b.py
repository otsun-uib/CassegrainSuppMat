import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# List of file names
file_names = ['results_10K.txt', 'results_50K.txt', 'results_100K.txt', 'results_200K.txt', 'results_500K.txt', 'results_1M.txt']

# Create an empty DataFrame to store the data
all_data = pd.DataFrame()

# Read each file and append the data to the DataFrame
for file_name in file_names:
    data = pd.read_csv(file_name, delimiter=',', header=None)
    all_data = pd.concat([all_data, data], ignore_index=True)

# Identify unique values of "Number of Rays Emitted"
ray_counts = all_data[8].unique()

# Assuming column 9 (index 8) is the number of rays emitted and column 15 (index 14) is the optical efficiency
num_rays_column = 8
efficiency_column = 14

# Check if column indices exist
if num_rays_column in all_data.columns and efficiency_column in all_data.columns:
    # Extract relevant columns
    data_subset = all_data[[num_rays_column, efficiency_column]]
    data_subset.columns = ['Number of Rays Emitted', 'Optical Efficiency']

    # Group by number of rays and calculate mean and standard deviation of efficiency
    grouped_data = data_subset.groupby('Number of Rays Emitted').agg(['mean', 'std'])
    grouped_data.columns = ['Mean Optical Efficiency', 'Std Deviation']

    # Reset index for easier access
    grouped_data.reset_index(inplace=True)

    # Visualize the results
    plt.errorbar(grouped_data['Number of Rays Emitted'], grouped_data['Mean Optical Efficiency'], yerr=grouped_data['Std Deviation'], fmt='o', capsize=5)
    plt.xscale('log')
    
    # Format X-axis labels
    plt.xticks(ray_counts, [f'${int(x / 10**int(np.log10(x)))} \\times 10^{{{int(np.log10(x))}}}$' for x in ray_counts], fontsize=10)

    plt.xlabel('Number of Rays Emitted')
    plt.ylabel('Mean Optical Efficiency')
    plt.grid(True)

    # Save the figure as a PDF file
    plt.savefig('plot_errors_sigma.pdf')
    plt.show()

    # Calculate interval based on different confidence levels
    confidence_levels = [0.95, 0.98, 0.9975]  # Confidence levels for interval calculation
    z_criticals = [stats.norm.ppf(1 - (1 - cl) / 2) for cl in confidence_levels]  # Calculate Z values for given confidence levels
    for cl, z_critical in zip(confidence_levels, z_criticals):
        grouped_data[f'Error ({cl * 100:.2f}\%)'] = z_critical * grouped_data['Std Deviation'] * 100

    # Generate LaTeX code for the table
    latex_table = "\\begin{table}[htb]\n\\centering\n\\begin{tabular}{l" + "c" * len(confidence_levels) + "}\n\\hline\n"
    latex_table += "Number of Rays Emitted" + " & ".join([f"Error ({cl * 100:.2f}\\%)" for cl in confidence_levels]) + " \\\\\n\\hline\n"

    for _, row in grouped_data.iterrows():
        row_str = f"${int(row['Number of Rays Emitted'] / 10**int(np.log10(row['Number of Rays Emitted'])))} \\times 10^{{{int(np.log10(row['Number of Rays Emitted']))}}}$"
        for cl in confidence_levels:
            error_col = f'Error ({cl * 100:.2f}\\%)'
            row_str += f" & $\\pm${row[error_col]:.1f}\\%"
        row_str += " \\\\\n"
        latex_table += row_str

    latex_table += "\\hline\n\\end{tabular}\n"
    latex_table += "\\caption{Error Intervals at Different Confidence Levels}\n"
    latex_table += "\\label{table:error_intervals}\n\\end{table}"

    # Print LaTeX code for the table
    print(latex_table)
