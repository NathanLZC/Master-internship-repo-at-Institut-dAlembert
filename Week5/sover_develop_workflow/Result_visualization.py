import numpy as np
import matplotlib.pyplot as plt

L = 2  # Domain size

# Load the displacement data from C++ and Python solvers
displacement_cpp = np.loadtxt("displacement.dat")
displacement_python = np.loadtxt("displacement_python.dat")


# Plot the displacement field
# Plot the displacement field
plt.figure(figsize=(10, 8))
plt.imshow(displacement_cpp, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement_cpp (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement Field with C++ Solver')


plt.figure(figsize=(10, 8))
plt.imshow(displacement_python, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement_python (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement Field with Python Solver')

# Plot the difference between the displacement fields
plt.figure(figsize=(10, 8))
plt.imshow(displacement_cpp - displacement_python, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Difference (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Difference between solvers (C++ - Python)')

plt.show()