
import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt("output_active.txt")
data2 = np.loadtxt("output_passive.txt")

plt.plot(data1[:,1])
