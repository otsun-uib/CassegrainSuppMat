import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from matplotlib.ticker import AutoMinorLocator

# Define the file names
file_T = 'results_effiopt_vs_angle_T.txt'
file_L = 'results_effiopt_vs_angle_L.txt'

# Read data from files
data_T = pd.read_csv(file_T, sep=',', header=None)
data_L = pd.read_csv(file_L, sep=',', header=None)

# Define design parameter values for two cases
case_1_values = [3000, 400, 0, 240, 19.0, 0.5]
case_2_values = [2000, 200, 0, 240, 20.0, 0.4]

# Function to filter data for a given case
def filter_case(data, case_values, plane_value):
    conditions = (data.iloc[:, 0] == case_values[0]) & \
                 (data.iloc[:, 1] == case_values[1]) & \
                 (data.iloc[:, 2] == case_values[2]) & \
                 (data.iloc[:, 3] == case_values[3]) & \
                 (data.iloc[:, 4] == case_values[4]) & \
                 (data.iloc[:, 5] == case_values[5]) & \
                 (data.iloc[:, 11] == plane_value)
    filtered_data = data[conditions]
    return filtered_data

# Filter data for each case in the transversal plane (column 11 = 0)
case_1_T = filter_case(data_T, case_1_values, 0)
case_2_T = filter_case(data_T, case_2_values, 0)

# Filter data for each case in the longitudinal plane (column 11 = 90)
case_1_L = filter_case(data_L, case_1_values, 90)
case_2_L = filter_case(data_L, case_2_values, 90)

# Get incidence angles and optical efficiencies for each case in the transversal plane
angles1_T = case_1_T.iloc[:, 12].values
efficiencies1_T = case_1_T.iloc[:, 13].values

angles2_T = case_2_T.iloc[:, 12].values
efficiencies2_T = case_2_T.iloc[:, 13].values

# Get incidence angles and optical efficiencies for each case in the longitudinal plane
angles1_L = case_1_L.iloc[:, 12].values
efficiencies1_L = case_1_L.iloc[:, 13].values

angles2_L = case_2_L.iloc[:, 12].values
efficiencies2_L = case_2_L.iloc[:, 13].values

# Sort data by incidence angle
sorted_indices_1_T = np.argsort(angles1_T)
angles1_T_sorted = angles1_T[sorted_indices_1_T]
efficiencies1_T_sorted = efficiencies1_T[sorted_indices_1_T]

sorted_indices_2_T = np.argsort(angles2_T)
angles2_T_sorted = angles2_T[sorted_indices_2_T]
efficiencies2_T_sorted = efficiencies2_T[sorted_indices_2_T]

sorted_indices_1_L = np.argsort(angles1_L)
angles1_L_sorted = angles1_L[sorted_indices_1_L]
efficiencies1_L_sorted = efficiencies1_L[sorted_indices_1_L]

sorted_indices_2_L = np.argsort(angles2_L)
angles2_L_sorted = angles2_L[sorted_indices_2_L]
efficiencies2_L_sorted = efficiencies2_L[sorted_indices_2_L]

# Create a PDF file to save the figures
with PdfPages('tracking_error.pdf') as pdf:
    # Create a single figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot the curves
    ax.plot(angles1_T_sorted, efficiencies1_T_sorted, 'bo:', label='Case 1 Transversal', markersize=8)
    ax.plot(angles2_T_sorted, efficiencies2_T_sorted, 'go:', label='Case 2 Transversal', markersize=8)
    ax.plot(angles1_L_sorted, efficiencies1_L_sorted, 'b^:', label='Case 1 Longitudinal', markersize=8)
    ax.plot(angles2_L_sorted, efficiencies2_L_sorted, 'g^:', label='Case 2 Longitudinal', markersize=8)
    
    # Set labels and title
    ax.set_xlabel('Incidence angle [degrees]', fontsize=16)
    ax.set_ylabel('Optical efficiency', fontsize=16)
    ax.legend(fontsize=14)
    ax.grid(True)
    
    # Set limits for X and Y axes
    ax.set_xlim(angles1_T_sorted.min() - 0.1, angles1_T_sorted.max() + 0.1)
    ax.set_xticks(list(np.arange(angles1_T_sorted.min(), angles1_T_sorted.max() + 0.1, 0.2)))
    ax.set_ylim(0, 0.8)
    ax.set_yticks(list(np.arange(0., 0.81, 0.1)))
    
    # Setting minor ticks on Y axis
    minor_locator_y = AutoMinorLocator(5)  # 5 minor ticks per major tick
    ax.yaxis.set_minor_locator(minor_locator_y)
    ax.tick_params(axis='y', which='minor', length=4)  # Minor ticks properties
    
    # Adjust font size of ticks on X and Y axes
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Adjust layout
    plt.tight_layout()

    # Save the figure to the PDF file
    pdf.savefig(fig)
    plt.close(fig)

# Display the plot (optional)
plt.show()
