import matplotlib.pyplot as plt
from numpy import *
# plot MaxCap, ER, PA in one figure based on their datafile

# get the sorted node rate from the datafile of each algorithm
MAX = loadtxt("MaxCap_plot.dat")
ER = loadtxt("ER_plot.dat")
PA = loadtxt("PA_plot.dat")

# sorted index from 1 to 10
sorted_index = [(i+1) for i in range(10)]

# plt.figure(1)
plt.xlabel("sorted node index")
plt.ylabel("Rate(kb/s)")
# ax1 = plt.plot(sorted_index, MAX, '-s',label = "MaxCap", sorted_index, ER, '-D', sorted_index, PA, '-o')
# plt.show()

# plot
fig, ax = plt.subplots()
plt.xlabel("sorted node index")
plt.ylabel("Rate(kb/s)")
ax.plot(sorted_index, MAX, '-s', label='MaxCap')
ax.plot(sorted_index, ER, '-D', label='SLP-ER')
ax.plot(sorted_index, PA, '-o', label='SLP-PA')
legend = ax.legend()

# legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
# Put a nicer background color on the legend.
# legend.get_frame().set_facecolor('C0')

plt.show()