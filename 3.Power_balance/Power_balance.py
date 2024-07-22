import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import pandas as pd

sns.set_style("ticks",{"xtick.direction": "in","ytick.direction": "in"})
#------------------------------------------------------------------------------
#                                  FILES    
#------------------------------------------------------------------------------

file_name1 = r'results_set_2'
	
load_here = os.getcwd()

def load_from_txt_or_csv(file):
    df = pd.read_table(
        file,
        sep=None,
        index_col=0,
        header=None,
        comment='#',
        engine='python'
    ).interpolate(method='index')
    df.insert(0, 0, df.index)
    return df.values

file_1 = os.path.join(load_here,file_name1 + '.txt')

lines_1 = load_from_txt_or_csv(file_1)



# ------------------------------ PLOT ANALISIS -------------------------------#
sns.set_style("ticks",{"xtick.direction": "in","ytick.direction": "in"})
sns.color_palette("tab10",5)

sigma_SB = 5.67E-8
G_DNI = 900
A_ap = 16 * 4000000 + 1760 * 2000
emissivity = 0.33
T_w = 600+ 273.15


font = {'family': 'serif',
    'weight': 'normal',
    'size': 12,
    }



fig, axes = plt.subplots(2, 3, figsize=(18, 10))

optical_ef = lines_1[:,13]
D_sphere = lines_1[:,3]
A_r = np.pi * (D_sphere / 1000) ** 2 / 4.0
power = G_DNI * (A_ap / 1000000) * optical_ef - A_r * emissivity * sigma_SB * T_w ** 4
power = power/ 1000

sns.scatterplot(ax=axes[0,0], x=lines_1[:,0], y=power, color = 'blue')
axes[0,0].set_ylabel("Power (kW)",fontdict=font)
axes[0,0].set_xlabel('Paraboloid focal distance (mm)',fontdict=font)
axes[0,0].set_ylim([44, 45.6])
axes[0,0].set_xlim(1400,4600)
axes[0,0].set_xticks(list(range(1500, 4501,500)))
axes[0,0].tick_params(labelsize=10)

sns.scatterplot(ax=axes[0,1], x=lines_1[:,1], y=power, color = 'blue')
axes[0,0].set_ylabel("Power (kW)",fontdict=font)
axes[0,1].set_xlabel('Hyperbola focus [mm]',fontdict=font)
axes[0,1].set_ylim([44, 45.6])
axes[0,2].set_xlim(-250,450)
axes[0,2].set_xticks(list(range(-200, 401, 100)))
axes[0,2].tick_params(labelsize=10)

sns.scatterplot(ax=axes[0,2], x=lines_1[:,2], y=power, color = 'blue')
axes[0,0].set_ylabel("Power (kW)",fontdict=font)
axes[0,2].set_xlabel('Receiver height distance [mm]',fontdict=font)
axes[0,2].set_ylim([44, 45.6])
axes[1,0].set_xlim(210,250)
axes[1,0].set_xticks(list(range(220, 241,20)))
axes[1,0].tick_params(labelsize=10)

sns.scatterplot(ax=axes[1,0], x=lines_1[:,3], y=power, color = 'blue')
axes[1,0].set_ylabel("Power (kW)",fontdict=font)
axes[1,0].set_xlabel('Receiver aperture diameter [mm]',fontdict=font)
axes[1,0].set_ylim([44, 45.6])
axes[1,0].tick_params(labelsize=10)

sns.scatterplot(ax=axes[1,1], x=lines_1[:,4], y=power, color = 'blue')
axes[0,0].set_ylabel("Power (kW)",fontdict=font)
axes[1,1].set_xlabel('Acceptance angle [deg.]',fontdict=font)
axes[1,1].set_ylim([44, 45.6])
axes[1,1].tick_params(labelsize=10)

sns.scatterplot(ax=axes[1,2], x=lines_1[:,5], y=power, color = 'blue')
axes[0,0].set_ylabel("Power (kW)",fontdict=font)
axes[1,2].set_xlabel('Truncation factor',fontdict=font)
axes[1,2].set_ylim([44, 45.6])
axes[1,2].set_xlim(0.25,0.75)
axes[1,2].set_xticks(list(np.arange(0.3, 0.71,0.1)))
axes[1,2].tick_params(labelsize=10)

plt.savefig('Power_balance'+'.png', format='png', dpi=500, bbox_inches='tight')

best_indexes = np.argpartition(power, -10)[-10:]
best_powers = power[best_indexes]
best_efficiency = optical_ef[best_indexes]
best_Pf = lines_1[best_indexes,0]
best_Hf = lines_1[best_indexes,1]
best_Rh = lines_1[best_indexes,2]
best_Dr = lines_1[best_indexes,3]
best_alpha = lines_1[best_indexes,4]
best_Tau = lines_1[best_indexes,5]
combined_array = np.column_stack((best_powers, best_efficiency, best_Pf, best_Hf, best_Rh, best_Dr, best_alpha, best_Tau))
sorted_indices = np.argsort(-best_powers)
sorted_combined_array = combined_array[sorted_indices]

# LaTeX code for the Table
table_code = r'''
\begin{table}
\centering
\label{table_results}
\begin{tabular}{cccccccc}
\toprule
$P_u$ [\si{kW}] & $\eta_{opt}$ & $F_p$ [\si{mm}] & $F_h$ [\si{mm}] & $H_r$ [\si{mm}] & $D_r$ [\si{mm}] & $\alpha$ [-] & $\tau$ [-] \\
\midrule
'''
for row in sorted_combined_array:
    table_code += f"{row[0]:.2f} & {row[1]*100:.2f} & {int(row[2])} & {int(row[3])} & {int(row[4])} & {int(row[5])} & {int(row[6])} & {row[7]:.1f} \\\\\n"
table_code += r'''
\bottomrule
\end{tabular}
\caption{Descripción de la tabla.}
\end{table}
'''

print(table_code)

#------------------------------------------------------------------------------     
