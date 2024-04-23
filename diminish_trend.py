import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
p = np.array([2,3,4,5,6])
R = np.array([0.66, 0.4, 0.285, 0.206, 0.140])

fig, ax = plt.subplots()
ax.plot(p, R)
plt.axis((1, 8, 0, 1))
ax.set(xlabel='No. of input modes', ylabel=' R',
       title='Diminishing Trend - Bunching of Photons')
ax.grid()

plt.show()