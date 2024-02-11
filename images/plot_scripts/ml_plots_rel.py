# Runs relative plots, but requires changing evaluate_nn

import numpy as np
import matplotlib.pyplot as plt
import subprocess

# call script
Means = []
Maxs = []
sizes = []
for run in range(6, 16):
    sizes.append(2**run)
    run = 'run_' + str(2**run)
    out = subprocess.run(args=['python', 'scripts/09_nn_training/evaluate_nn.py', '-run', run], stdout = subprocess.PIPE)
    print(out)
    out = out.stdout.decode().split('\n')
    Means.append(out[0].split(': ')[1])
    Maxs.append(out[1].split(': ')[1])


def to_float(row):
    return [float(r) for r in row]

Maxs = np.array([to_float(row.split(" | ")) for row in Maxs])
Means = np.array([to_float(row.split(" | ")) for row in Means])

fs=12
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('text', usetex=True)

s = 30
fig1, ax1 = plt.subplots(2, 2, figsize=(10, 6),)
ax1[0][1].set_xscale("log", base=2)
markers = ['o','>','D']
colors = ['b','r', 'g']

# Pressures
# ax1[0][1].scatter(sizes, Maxs[:, 3], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[0][1].scatter(sizes, Maxs[:, 4], s = s, color = "r", marker='>', label = 'Mean')
ax1[0][1].scatter(sizes, Maxs[:, 5], s = s, color = "g", marker='D', label = 'Systolic')

ax1[0][1].set_title("Max Relative Errors",fontsize=fs)
ax1[0][1].set_ylabel("Pressure [mmHg]",fontsize=fs)
ax1[0][1].set_xlabel("Train Size",fontsize=fs)
ax1[0][1].set_xscale("log", base=2)
#ax1[0][1].set_yscale("log")
ax1[0][1].legend(fontsize=fs-2)
ax1[0][1].tick_params(axis='both', which='major', labelsize=fs)    

# ax1[0][0].scatter(sizes, Means[:, 3], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[0][0].scatter(sizes, Means[:, 4], s = s, color = "r", marker='>', label = 'Mean')
ax1[0][0].scatter(sizes, Means[:, 5], s = s, color = "g", marker='D', label = 'Systolic')

ax1[0][0].set_title("Mean Relative Errors",fontsize=fs)
ax1[0][0].set_ylabel("Pressure [mmHg]",fontsize=fs)
ax1[0][0].set_xlabel("Train Size",fontsize=fs)
ax1[0][0].set_xscale("log", base=2)
#ax1[0][0].set_yscale("log")
ax1[0][0].legend(fontsize=fs-2)
ax1[0][0].tick_params(axis='both', which='major', labelsize=fs)    


# Flows
ax1[1][1].set_title("Max Relative Errors",fontsize=fs)
# ax1[1][1].scatter(sizes, Maxs[:, 0], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[1][1].scatter(sizes, Maxs[:, 1], s = s, color = "r", marker='>', label = 'Mean')
ax1[1][1].scatter(sizes, Maxs[:, 2], s = s, color = "g", marker='D', label = 'Systolic')

ax1[1][1].set_ylabel(r"Flow [cc$^3$/sec]",fontsize=fs)
ax1[1][1].set_xlabel("Train Size",fontsize=fs)
ax1[1][1].set_xscale("log", base=2)
#ax1[1][1].set_yscale("log")
ax1[1][1].legend(fontsize=fs-2)
ax1[1][1].tick_params(axis='both', which='major', labelsize=fs)    

# ax1[1][0].scatter(sizes, Means[:, 0], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[1][0].scatter(sizes, Means[:, 1], s = s, color = "r", marker='>', label = 'Mean')
ax1[1][0].scatter(sizes, Means[:, 2], s = s, color = "g", marker='D', label = 'Systolic')

ax1[1][0].set_title("Mean Relative Errors",fontsize=fs)
ax1[1][0].set_ylabel(r"Flow [cc$^3$/sec]",fontsize=fs)
ax1[1][0].set_xlabel("Train Size",fontsize=fs)
ax1[1][0].set_xscale("log", base=2)
#ax1[1][0].set_yscale("log")
ax1[1][0].legend(fontsize=fs-2)
ax1[1][0].tick_params(axis='both', which='major', labelsize=fs)    

plt.tight_layout()
plt.savefig('images/paper/06_ann/11_ml_results_rel.pdf')