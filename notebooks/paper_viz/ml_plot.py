import numpy as np
import matplotlib.pyplot as plt



Means = ["0.0014 | 0.0295 | 0.0828 | 0.0058 | 0.0918 | 0.2339",
         "0.0004 | 0.008 | 0.0234 | 0.0017 | 0.0266 | 0.0688",
         "0.0002 | 0.003 | 0.0082 | 0.0007 | 0.0101 | 0.026",
         "0.0001 | 0.0022 | 0.0061 | 0.0005 | 0.0075 | 0.0187",
         "0.0001 | 0.0014 | 0.0038 | 0.0003 | 0.0046 | 0.0113",
         "0.0 | 0.0009 | 0.0024 | 0.0002 | 0.003 | 0.0075",
         "0.0 | 0.0006 | 0.0017 | 0.0002 | 0.0021 | 0.0052",
         "0.0 | 0.0004 | 0.001 | 0.0001 | 0.0013 | 0.0031",
         "0.0 | 0.0002 | 0.0006 | 0.0001 | 0.0008 | 0.002",
         "0.0 | 0.0002 | 0.0005 | 0.0001 | 0.0006 | 0.0016",
         
         
          ]
Maxs = ["13.1345 | 6.3374 | 17.8532 | 0.1035 | 1.9102 | 5.3067",
        "6.9166 | 2.5306 | 6.7035 | 0.046 | 0.8758 | 3.0448",
        "2.0724 | 1.2345 | 2.5609 | 0.0218 | 0.4211 | 1.0173",
        "1.9111 | 0.7476 | 2.3511 | 0.0229 | 0.3624 | 1.0413",
        "1.2331 | 0.6951 | 1.2411 | 0.0188 | 0.2844 | 0.561",
        "1.0272 | 0.1799 | 0.5851 | 0.0101 | 0.1146 | 0.3229",
        "0.5603 | 0.0983 | 0.2713 | 0.0074 | 0.0843 | 0.2042",
        "0.4759 | 0.0734 | 0.4816 | 0.0057 | 0.0615 | 0.2212",
        "0.1687 | 0.0798 | 0.1368 | 0.0043 | 0.0426 | 0.1324",
        "0.1501 | 0.0385 | 0.3103 | 0.004 | 0.0229 | 0.0679",
        ]
sizes = [64,
         128,
         256,
         512,
         1024,
         2048,
         4096,
         8192,
         16384,
         32768
         ]

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
fig1, ax1 = plt.subplots(2, 2, figsize=(10, 8),)
ax1[0][1].set_xscale("log", base=2)
markers = ['o','>','D']
colors = ['b','r', 'm']

# Pressures
ax1[0][1].scatter(sizes, Maxs[:, 3], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[0][1].scatter(sizes, Maxs[:, 4], s = s, color = "r", marker='>', label = 'Mean')
ax1[0][1].scatter(sizes, Maxs[:, 5], s = s, color = "m", marker='D', label = 'Systolic')

ax1[0][1].set_title("Max Absolute Errors",fontsize=fs)
ax1[0][1].set_ylabel("Pressure [mmHg]",fontsize=fs)
ax1[0][1].set_xlabel("Train Size",fontsize=fs)
ax1[0][1].set_xscale("log", base=2)
#ax1[0][1].set_yscale("log")
ax1[0][1].legend(fontsize=fs-2)
ax1[0][1].tick_params(axis='both', which='major', labelsize=fs)    

ax1[0][0].scatter(sizes, Means[:, 3], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[0][0].scatter(sizes, Means[:, 4], s = s, color = "r", marker='>', label = 'Mean')
ax1[0][0].scatter(sizes, Means[:, 5], s = s, color = "m", marker='D', label = 'Systolic')

ax1[0][0].set_title("Mean Absolute Errors",fontsize=fs)
ax1[0][0].set_ylabel("Pressure [mmHg]",fontsize=fs)
ax1[0][0].set_xlabel("Train Size",fontsize=fs)
ax1[0][0].set_xscale("log", base=2)
#ax1[0][0].set_yscale("log")
ax1[0][0].legend(fontsize=fs-2)
ax1[0][0].tick_params(axis='both', which='major', labelsize=fs)    


# Flows
ax1[1][1].scatter(sizes, Maxs[:, 0], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[1][1].scatter(sizes, Maxs[:, 1], s = s, color = "r", marker='>', label = 'Mean')
ax1[1][1].scatter(sizes, Maxs[:, 2], s = s, color = "m", marker='D', label = 'Systolic')

ax1[1][1].set_ylabel(r"Flow [cc$^3$/sec]",fontsize=fs)
ax1[1][1].set_xlabel("Train Size",fontsize=fs)
ax1[1][1].set_xscale("log", base=2)
#ax1[1][1].set_yscale("log")
ax1[1][1].legend(fontsize=fs-2)
ax1[1][1].tick_params(axis='both', which='major', labelsize=fs)    

ax1[1][0].scatter(sizes, Means[:, 0], s = s, color = "b", marker='o', label = 'Diastolic')
ax1[1][0].scatter(sizes, Means[:, 1], s = s, color = "r", marker='>', label = 'Mean')
ax1[1][0].scatter(sizes, Means[:, 2], s = s, color = "m", marker='D', label = 'Systolic')

ax1[1][0].set_ylabel(r"Flow [cc$^3$/sec]",fontsize=fs)
ax1[1][0].set_xlabel("Train Size",fontsize=fs)
ax1[1][0].set_xscale("log", base=2)
#ax1[1][0].set_yscale("log")
ax1[1][0].legend(fontsize=fs-2)
ax1[1][0].tick_params(axis='both', which='major', labelsize=fs)    

plt.tight_layout()
plt.savefig('plots/ml_results.pdf')
