## This script is meant to be a simple example of how to process the simulations.

# Import the libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import _pickle as pkl
import seaborn as sns
import os 

# Use a custom matplotlib style for the plots
if os.path.exists("./paper.mplstyle"):
    plt.style.use("./paper.mplstyle")
    print("Using paper style")
else:
    print("Using default style")

###########################################
############# Data processing #############
###########################################

# Read data
values = [1, 5, 9, 13, 17] # Values of the scan
J23  = {} # Dictionary to store the results

for v in values:
    # Read the data
    df = pd.read_csv('./ScanJ23_' + str(v))
    # Group by time and cell type
    df = df.groupby(['time', 'cell_type']).size().reset_index(name='count')

    res = {}
    for c in [1, 2, 3]:
        time = df[df['cell_type'] == c]['time']
        count = df[df['cell_type'] == c]['count']
        res["Type"+str(c)] = {"time": time, "count": count}

    J23[v] = res

# Save the data
with open('../Results/scan_J23.pkl', 'wb') as f:
    pkl.dump(J23, f)


###########################################
################# Plotting ################
###########################################
# When it's needed a plot with three cell types
with open('../Results/scan_J23.pkl', 'rb') as f:
    data = pkl.load(f)

fig, ax = plt.subplots(1, 3, figsize=(17, 5))

data_legend = r"$J_{23}$"
fontsize=20
xlim = 3000
colors = sns.color_palette("Set2", len(data.keys())) # Mechanical competition
colors = sns.color_palette("husl", len(data.keys())) # Biochemical competition

for i, v in enumerate(data.keys()):
    c1 = data[v]["Type1"]["count"].iloc[0]
    c2 = data[v]["Type2"]["count"].iloc[0]
    c3 = data[v]["Type3"]["count"].iloc[0]
    v_lab = v

    time = data[v]["Type1"]["time"]
    count = data[v]["Type1"]["count"]
    count = count / count.iloc[0]
    ax[0].plot(time, count, label=v_lab, color=colors[i])
    ax[0].set_title("Type 1", fontsize=fontsize)
    ax[0].legend(title=data_legend)
    ax[0].set_xlabel("Time", fontsize=fontsize)
    ax[0].set_ylabel("Relative cell count", fontsize=fontsize)
    ax[0].set_xlim(0, xlim)
    ax[0].hlines(1, 0, xlim, "black", "--")
    ax[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
    ax[0].tick_params(axis='both', which='major', labelsize=fontsize)
    ax[0].tick_params(axis='both', which='minor', labelsize=fontsize)

    time = data[v]["Type2"]["time"]
    count = data[v]["Type2"]["count"]
    count = count / count.iloc[0]
    ax[1].plot(time, count, label=v, color=colors[i])
    ax[1].set_title("Type 2", fontsize=fontsize)
    ax[1].set_xlabel("Time", fontsize=fontsize)
    ax[1].set_ylabel("Relative cell count", fontsize=fontsize)
    ax[1].set_xlim(0, xlim)
    ax[1].hlines(1, 0, xlim, "black", "--")
    ax[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
    ax[1].tick_params(axis='both', which='major', labelsize=fontsize)
    ax[1].tick_params(axis='both', which='minor', labelsize=fontsize)

    time = data[v]["Type3"]["time"]
    count = data[v]["Type3"]["count"]
    count = count / count.iloc[0]
    ax[2].plot(time, count, label=v, color=colors[i])
    ax[2].set_title("Type 3", fontsize=fontsize)
    ax[2].set_xlabel("Time", fontsize=fontsize)
    ax[2].set_ylabel("Relative cell count", fontsize=fontsize)
    ax[2].set_xlim(0, xlim)
    ax[2].hlines(1, 0, xlim, "black", "--")
    ax[2].grid(True, linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
    ax[2].tick_params(axis='both', which='major', labelsize=fontsize)
    ax[2].tick_params(axis='both', which='minor', labelsize=fontsize)


ax[0].legend(title=data_legend, loc=[3.5, 0], ncols=1, 
             fontsize=fontsize, title_fontsize=fontsize)
plt.show()


###########################################
# When it's needed a plot with two cell types
with open('../Results/ref_Jhetero.pkl', 'rb') as f:
    data = pkl.load(f)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

data_legend = r"$J_{hetero}$"
xlim = 3000
fontsize=20
colors = sns.color_palette("Set2", len(data.keys())) # Mechanical competition
colors = sns.color_palette("husl", len(data.keys())) # Biochemical competition

for i, v in enumerate(data.keys()):
    time = data[v]["Type1"]["time"]
    count = data[v]["Type1"]["count"]
    count = count / count.iloc[0]
    ax[0].plot(time, count, label=v, color=colors[i])
    ax[0].set_title("Type 1", fontsize=fontsize)
    ax[0].legend(title=data_legend)
    ax[0].set_xlabel("Time", fontsize=fontsize)
    ax[0].set_ylabel("Relative cell count", fontsize=fontsize)
    ax[0].set_xlim(0, xlim)
    ax[0].hlines(1, 0, 3000, "black", "--")
    ax[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
    ax[0].tick_params(axis='both', which='major', labelsize=fontsize)
    ax[0].tick_params(axis='both', which='minor', labelsize=fontsize)

    time = data[v]["Type2"]["time"]
    count = data[v]["Type2"]["count"]
    count = count / count.iloc[0]
    ax[1].plot(time, count, label=v, color=colors[i])
    ax[1].set_title("Type 2", fontsize=fontsize)
    ax[1].set_xlabel("Time", fontsize=fontsize)
    ax[1].set_ylabel("Relative cell count", fontsize=fontsize)
    ax[1].set_xlim(0, xlim)
    ax[1].hlines(1, 0, 3000, "black", "--")
    ax[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.7, color='gray')
    ax[1].tick_params(axis='both', which='major', labelsize=fontsize)
    ax[1].tick_params(axis='both', which='minor', labelsize=fontsize)

ax[0].legend(title=data_legend, loc=[2.3, 0], ncols=1, 
             fontsize=fontsize, title_fontsize=fontsize)
plt.show()