import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from tabulate import tabulate

file_names = ['results_10K.txt', 'results_50K.txt', 'results_100K.txt', 'results_200K.txt', 'results_500K.txt', 'results_1M.txt']

def read_and_process_files(file_names):
    """
    Read and process data from multiple files.
    
    Args:
    - file_names (list): List of file names to read
    
    Returns:
    - data (list): List of tuples, each containing (num_rays, opt_efficiency)
    """
    data = []
    for file_name in file_names:
        try:
            df = pd.read_csv(file_name, header=None)
            if df.shape[1] >= 15:  # Check if there are at least 15 columns
                num_rays = df.iloc[:, 8].unique()[0]  # Unique values from the 9th column (index 8)
                opt_efficiency = df.iloc[:, 14]  # 15th column (index 14) for all rows
                data.append((num_rays, opt_efficiency))
            else:
                print(f"Error: {file_name} does not have enough columns.")
        except Exception as e:
            print(f"Error reading {file_name}: {e}")
    return data

def gaussian(x, mu, sigma):
    """
    Gaussian (normal) distribution function.

    Args:
    - x (np.ndarray): Input values
    - mu (float): Mean of the distribution
    - sigma (float): Standard deviation of the distribution

    Returns:
    - np.ndarray: Gaussian distribution evaluated at x
    """
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def plot_histograms(data):
    """
    Plot histograms and Gaussian fits for optical efficiency data.

    Args:
    - data (list): List of tuples (num_rays, opt_efficiency)
    """
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    for i, (num_rays, opt_efficiency) in enumerate(data):
        row = i // 3
        col = i % 3
        ax = axs[row, col]

        # Normalized histogram
        n, bins, patches = ax.hist(opt_efficiency, bins=13, density=True, alpha=0.6, color='blue', edgecolor='black', linewidth=1.2, label='Histogram')
        
        # Get bin centers and calculate Gaussian fit
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        x_fit = np.linspace(min(bin_centers), max(bin_centers), 1000)
        p0 = [np.mean(opt_efficiency), np.std(opt_efficiency)]
        popt, _ = curve_fit(gaussian, bin_centers, n, p0=p0)
        mu, sigma = popt
        y_fit = gaussian(x_fit, mu, sigma)
        
        # Calculate R²
        r_squared = r2_score(n, gaussian(bin_centers, mu, sigma))
        
        # Plot Gaussian fit
        ax.plot(x_fit, y_fit, 'r-', label=f'Gaussian fit (R²={r_squared:.4f})')
        
        ax.set_xlabel('Optical Efficiency')
        ax.set_ylabel('Density')
        ax.set_title(f'Number of emitted rays = ${num_rays // 10**int(np.log10(num_rays))}\\times10^{{{int(np.log10(num_rays))}}}$')
        ax.legend()
        
        # Set x-axis ticks
        x_ticks = np.linspace(min(opt_efficiency), max(opt_efficiency), 7)
        ax.set_xticks(x_ticks)
        
        # Limit to 3 decimals
        ax.set_xticklabels([f'{x:.3f}' for x in x_ticks])
        
        # Place legend at lower center
        ax.legend(loc='lower center', fontsize=10)
        
    plt.tight_layout()
    
    # Save figure to a PDF file
    plt.savefig('plot_errors.pdf')
    plt.show()

def format_latex_table(data):
    """
    Format data into a LaTeX table.

    Args:
    - data (list): List of tuples (num_rays, opt_efficiency)

    Returns:
    - str: LaTeX formatted table string
    """
    headers = ["Number of Rays", "$\mu$", "$\sigma$", "$\sigma_{data}$"]
    rows = []
    for num_rays, opt_efficiency in data:
        exponent = int(np.log10(num_rays))
        base = num_rays // 10**exponent
        num_rays_str = f'${base} \\times 10^{{{exponent}}}$'
        
        n, bins, _ = plt.hist(opt_efficiency, bins=13, density=True, alpha=0.6, color='blue', edgecolor='black', linewidth=1.2, label='Histogram')
        
        # Get bin centers and calculate Gaussian fit
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        x_fit = np.linspace(min(bin_centers), max(bin_centers), 1000)
        p0 = [np.mean(opt_efficiency), np.std(opt_efficiency)]
        popt, _ = curve_fit(gaussian, bin_centers, n, p0=p0)
        mu, sigma = popt
        y_fit = gaussian(x_fit, mu, sigma)
        
        # Calculate R²
        # r_squared = r2_score(n, gaussian(bin_centers, mu, sigma))
        
        mu_str = f'{mu:.4f}'
        sigma_str = f'{sigma:.4f}'
        # r_squared_str = f'{r_squared:.4f}'
        std_from_data_str = f'{np.std(opt_efficiency):.4f}'
        rows.append([num_rays_str, mu_str, sigma_str, std_from_data_str])
    
    plt.close()  # Close figure created by plt.hist to avoid overlapping
    
    table = tabulate(rows, headers=headers, tablefmt="plain")
    latex_table = "\\begin{table}[htb]\n\\centering\n\\begin{tabular}{lccccl}\n\\hline\n"
    latex_table += " & ".join(headers) + " \\\\\n\\hline\n"
    
    for row in rows:
        latex_table += " & ".join(row) + " \\\\\n"
    
    caption = 'Parameters from Gaussian fits of optical efficiency distributions for various numbers of emitted rays, including the standard deviation of the data values for comparison with the obtained $\sigma$ value.'
    label = 'table_errors'
    latex_table += "\\hline\n\\end{tabular}\n"
    latex_table += f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{table}}"
    return latex_table

def main():
    # Read and process data
    data = read_and_process_files(file_names)
    
    # Plot histograms and Gaussian fits
    plot_histograms(data)
    
    # Format data into LaTeX table and print
    print(format_latex_table(data))

if __name__ == "__main__":
    main()
