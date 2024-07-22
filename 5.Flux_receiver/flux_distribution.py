import numpy as np
import pandas as pd
import re
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from scipy.ndimage import convolve

num_rays = 10000

# Center of the sphere
center_z = 0
sphere_radius = 240

# Function to extract only the numeric value from a line
def extract_value(line):
    return float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])

# Load S_wl
with open('source_wavelengths_Case_1.txt', 'r') as file:
    S_wl = extract_value(file.readline().strip())

# Load ray data
ray_data = pd.read_csv('Th_points_absorber_Case_1.txt', delim_whitespace=True, comment='#', header=None)

# Adjust the column name list according to the number of columns
ray_data.columns = ['Factor', 'X', 'Y', 'Z', 'Unknown1', 'Unknown2', 'Unknown3', 'Unknown4', 'Unknown5', 'Unknown6', 'Wavelength']

# Load the solar spectrum
solar_spectrum = pd.read_csv('ASTMG173-direct.txt', delim_whitespace=True, header=None)
solar_spectrum.columns = ['Wavelength', 'Power']

# Define the wavelength range
min_wavelength = solar_spectrum['Wavelength'].min()
max_wavelength = solar_spectrum['Wavelength'].max()
step = 1  # Step to define the wavelength range
wavelengths_range = np.arange(min_wavelength, max_wavelength + step, step)

# Perform cubic interpolation to adjust the data to 10 nm intervals
binned_solar_spectrum = interp1d(solar_spectrum['Wavelength'], solar_spectrum['Power'], kind='cubic')(wavelengths_range)

# Calculate the integral
integral = np.trapz(binned_solar_spectrum, wavelengths_range)

print("Integral value of the solar spectrum (every 10 nm):", integral)

# Calculate the power of each ray
def calculate_ray_power(row):
    wl = row['Wavelength']
    rounded_wl = round(wl / 10) * 10  # Round to the nearest wavelength in steps of 10 nm
    wl_left = rounded_wl - 5  # Left limit of integration
    wl_right = rounded_wl + 5  # Right limit of integration
    
    # Find the corresponding indices in the wavelength array
    try:
        index_left = np.where(wavelengths_range == wl_left)[0][0]
        index_right = np.where(wavelengths_range == wl_right)[0][0]
        
        # Calculate the integral of the solar spectrum between the limits
        integral = np.trapz(binned_solar_spectrum[index_left:index_right+1], wavelengths_range[index_left:index_right+1])
        power_per_wavelength = integral  # No need to multiply by 10 as the integral is already calculated in 10 nm intervals
    except IndexError:
        power_per_wavelength = 0  # Assign 0 power if the wavelength is not found (error handling)
    
    return (power_per_wavelength / num_rays) * row['Factor'] * S_wl

ray_data['RayPower'] = ray_data.apply(calculate_ray_power, axis=1)

print("Power absorbed in the sphere:", sum(ray_data['RayPower']))

# Convert Cartesian coordinates to spherical coordinates
ray_data['r'] = np.sqrt(ray_data['X']**2 + ray_data['Y']**2 + (ray_data['Z'] - center_z)**2)
ray_data['theta'] = np.arccos((ray_data['Z'] - center_z) / ray_data['r'])
ray_data['phi'] = np.arctan2(ray_data['Y'], ray_data['X'])

# Define the minimum theta value according to the results
theta_min = min(ray_data['theta'])

# Filter data to exclude theta values smaller than theta_min
ray_data = ray_data[ray_data['theta'] >= theta_min]

# Define the number of divisions for discretization
theta_bins = 80
phi_bins = 30

# Create a mesh to discretize the sphere
theta_edges = np.linspace(theta_min, np.pi, theta_bins + 1)
phi_edges = np.linspace(-np.pi, np.pi, phi_bins + 1)

# Initialize the heatmap
heatmap = np.zeros((theta_bins, phi_bins))

# Calculate the total power in each cell
for i in range(theta_bins):
    for j in range(phi_bins):
        in_bin = ((ray_data['theta'] >= theta_edges[i]) & (ray_data['theta'] < theta_edges[i + 1]) &
                  (ray_data['phi'] >= phi_edges[j]) & (ray_data['phi'] < phi_edges[j + 1]))
        heatmap[i, j] = ray_data[in_bin]['RayPower'].sum()

# Calculate the area of each cell
dtheta = np.diff(theta_edges)
dphi = np.diff(phi_edges)
cell_areas = np.outer(np.sin(theta_edges[:-1] + dtheta / 2) * dtheta, dphi)

# Power per unit area
heatmap /= cell_areas

# Smooth the heatmap using an average filter
kernel = np.ones((5, 5)) / 25  # 5x5 convolution kernel for a wider average
heatmap_smoothed = convolve(heatmap, kernel, mode='reflect')  # Use 'reflect' mode for edges

# Convert from radians to degrees
theta_edges_deg = np.degrees(theta_edges)
phi_edges_deg = np.degrees(phi_edges)

# Invert the theta angle scale (Y-axis labels)
theta_edges_deg_reversed = np.flip(theta_edges_deg)
heatmap_reversed = np.flip(heatmap_smoothed, axis=0)

# Power per unit area (convert from W/m^2 to kW/m^2)
heatmap_kW = heatmap_reversed / 1000

# Maximum power per unit area
print("Maximum power per unit area:", heatmap_kW.max())

# Create the heatmap with a custom color scale using PowerNorm
norm = mcolors.PowerNorm(gamma=0.8)  # Adjust gamma for greater sensitivity to low values

plt.figure(figsize=(10, 8))
heatmap_plot = plt.pcolormesh(phi_edges_deg, theta_edges_deg_reversed, heatmap_kW, shading='auto', cmap='inferno', norm=norm)

# Create an inverted color bar with adjusted pad
cbar = plt.colorbar(heatmap_plot, location='left', pad=0.15)
cbar.set_label('Power per unit area [kW/m$^2$]', fontsize=18)
# Set the Y-axis color bar ticks every 3
cbar.set_ticks(np.arange(0, np.max(heatmap_kW) + 1, 3))
cbar.ax.invert_yaxis()  # Invert the color bar

# Labels and title with increased font size
plt.xlabel('Phi [degrees]', fontsize=18)
plt.ylabel('Theta [degrees]', fontsize=16)
plt.gca().invert_yaxis()  # Invert the Y-axis

# Manually set the tick locations and labels for the X-axis
plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], [-180, -135, -90, -45, 0, 45, 90, 135, 180])

# Increase the size of the numbers on the axes
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Save the plot to a PDF file
plt.savefig('heat_map_Case_1.pdf', format='pdf')

# Show the plot
plt.show()

# Convert the heatmap from kW/m² to W/m²
heatmap_W = heatmap_kW * 1000

# Calculate the areas of the cells on the sphere
dtheta = np.diff(theta_edges)
dphi = np.diff(phi_edges)
cell_areas = np.outer(np.sin(theta_edges[:-1] + dtheta / 2) * dtheta, dphi)

# Invert the cell areas array to match the heatmap
cell_areas_reversed = np.flip(cell_areas, axis=0)

# Multiply each cell in the heatmap by its corresponding area
total_power = np.sum(heatmap_W * cell_areas_reversed)

# Display the total power absorbed by the sphere
print("Total power absorbed by the sphere (calculated from the heatmap):", total_power, "W")

# Compare with the sum of ray_data['RayPower']
total_ray_power = ray_data['RayPower'].sum()
print("Total power absorbed by the sphere (sum of ray_data['RayPower']):", total_ray_power, "W")

# Check the difference
difference = abs(total_power - total_ray_power)
print("Difference between the two calculation methods:", difference, "W")

# Additional part to create the 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Create a mesh for theta and phi
theta_centers = (theta_edges[:-1] + theta_edges[1:]) / 2
phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2
theta, phi = np.meshgrid(theta_centers, phi_centers, indexing='ij')

# Convert spherical coordinates to Cartesian coordinates
X = sphere_radius * np.sin(theta) * np.cos(phi)
Y = sphere_radius * np.sin(theta) * np.sin(phi)
Z = sphere_radius * np.cos(theta) + center_z  # Adjust to center at z=center_z

# Normalize the values for colormap
heatmap_kW = np.flip(heatmap_reversed, axis=0)
heatmap_kW_normalized = (heatmap_kW - np.min(heatmap_kW)) / (np.max(heatmap_kW) - np.min(heatmap_kW))

# Create a colormap
colors = cm.inferno(heatmap_kW_normalized)

# Create the 3D surface
ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, shade=False)

# Add a color bar
# mappable = cm.ScalarMappable(cmap='inferno', norm=norm)
# mappable.set_array(heatmap_kW)
# cbar = plt.colorbar(mappable, shrink=0.5, aspect=5)
# cbar.set_label('Power per unit area (kW/m$^2$)', fontsize=16)

# Adjust the camera perspective for better visualization
ax.view_init(elev=-20, azim=-20)

# Adjust labels and title
ax.set_xlabel('X [mm]', fontsize=18)
ax.set_ylabel('Y [mm]', fontsize=18)
ax.set_zlabel('Z [mm]', fontsize=18)
# plt.title('Heat Map 3D Representation', fontsize=16)

plt.savefig('heat_map_3D_Case_1.pdf', format='pdf')
plt.show()



