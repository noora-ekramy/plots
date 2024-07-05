import numpy as np
import matplotlib.pyplot as plt

# Load data from se3.dat
data = np.loadtxt('se3.dat')

# Extract columns 2 and 5 (Python is 0-indexed, so these are columns 1 and 4)
x = data[:, 1]
y = data[:, 4]

# Plot the data
plt.plot(x, y, marker='o')
plt.xlabel('Column 2')
plt.ylabel('Column 5')
plt.title('Cyclic Stress-Strain Curves')
plt.grid(True)
plt.show()
