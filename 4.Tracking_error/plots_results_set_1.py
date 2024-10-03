import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import pandas as pd


#    sns.set_palette("BuGn_d",2)
sns.set_style("ticks",{"xtick.direction": "in","ytick.direction": "in"})
#    sns.set_palette("Set2")

#------------------------------------------------------------------------------
#                                  FILES    
#------------------------------------------------------------------------------

file_name1 = r'results_set_1'
	
load_here = os.getcwd()

def load_from_txt_or_csv(file):
    df = pd.read_table(
        file,
        delimiter=',',
        index_col=0,
        header=None,
        comment='#',
        engine='python'
    ).interpolate(method='index')
    df.insert(0, 0, df.index)
    return df.values

# Cargar todas las líneas del archivo
file_1 = os.path.join(load_here,file_name1 + '.txt')

lines_1 = load_from_txt_or_csv(file_1)



# ------------------------------ PLOT ANALISIS -------------------------------#
sns.set_style("ticks",{"xtick.direction": "in","ytick.direction": "in"})
sns.color_palette("tab10",5)



font = {'family': 'serif',
    'weight': 'normal',
    'size': 12,
    }

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

sns.scatterplot(ax=axes[0,0], x=lines_1[:,0], y=lines_1[:,13], color = 'blue')
axes[0,0].set_ylabel("Optical efficiency",fontdict=font)
axes[0,0].set_xlabel('Paraboloid focal distance [mm]',fontdict=font)
axes[0,0].set_ylim([0.66, 0.77])
axes[0,0].set_xlim(900,4600)
axes[0,0].set_xticks(list(range(1000, 4501,500)))
# axes[0,0].set_xlim(400,6100)
# axes[0,0].set_xticks(list(range(4500, 6001,500)))
axes[0,0].tick_params(labelsize=10)

sns.scatterplot(ax=axes[0,1], x=lines_1[:,1], y=lines_1[:,13], color = 'blue')
axes[0,0].set_ylabel("Optical efficiency",fontdict=font)
axes[0,1].set_xlabel('Hyperbola focus [mm]',fontdict=font)
axes[0,1].set_ylim([0.66, 0.77])
axes[0,1].set_xlim(50,950)
axes[0,1].set_xticks(list(range(100, 901,200)))
axes[0,1].tick_params(labelsize=10)

sns.scatterplot(ax=axes[0,2], x=lines_1[:,2], y=lines_1[:,13], color = 'blue')
axes[0,0].set_ylabel("Optical efficiency",fontdict=font)
axes[0,2].set_xlabel('Receiver height distance [mm]',fontdict=font)
axes[0,2].set_ylim([0.66, 0.77])
axes[0,2].set_xlim(-450,650)
axes[0,2].set_xticks(list(range(-400, 601,200)))
axes[0,2].tick_params(labelsize=10)

sns.scatterplot(ax=axes[1,0], x=lines_1[:,3], y=lines_1[:,13], color = 'blue')
axes[1,0].set_ylabel("Optical efficiency",fontdict=font)
axes[1,0].set_xlabel('Receiver aperture diameter [mm]',fontdict=font)
axes[1,0].set_xlim(190,250)
axes[1,0].set_xticks(list(range(200, 241,20)))
axes[1,0].set_ylim([0.66, 0.77])
axes[1,0].tick_params(labelsize=10)

sns.scatterplot(ax=axes[1,1], x=lines_1[:,4], y=lines_1[:,13], color = 'blue')
axes[1,0].set_ylabel("Optical efficiency",fontdict=font)
axes[1,1].set_xlabel('Acceptance angle [deg.]',fontdict=font)
axes[1,1].set_ylim([0.66, 0.77])
axes[1,1].tick_params(labelsize=10)

sns.scatterplot(ax=axes[1,2], x=lines_1[:,5], y=lines_1[:,13], color = 'blue')
axes[1,0].set_ylabel("Optical efficiency",fontdict=font)
axes[1,2].set_xlabel('Truncation factor',fontdict=font)
axes[1,2].set_ylim([0.66, 0.77])
axes[1,2].set_xlim(0.15,0.85)
axes[1,2].set_xticks(list(np.arange(0.2, 0.81,0.2)))
axes[1,2].tick_params(labelsize=10)




plt.savefig('Optical_Efficiency_set_1'+'.png', format='png', dpi=500, bbox_inches='tight')


#------------------------------------------------------------------------------     
