import numpy as np
import matplotlib.pyplot as plt
import time

def pressure_distribution(x, y, x0, y0, L, p0):
    return p0 * np.cos(2 * np.pi * (x - x0) / L) * np.cos(2 * np.pi * (y - y0) / L)

L = 1.0
p0 = 1.0
x0 = 0
y0 = 0
n = 100
m = 100

x, y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, m))
test_pressure = pressure_distribution(x, y, x0, y0, L, p0)
pressure_fourier = np.fft.fft2(test_pressure, norm='ortho')

# Visualization of the Fourier-transformed pressure
fig, ax = plt.subplots()
extent = [x.min(), x.max(), y.min(), y.max()]
im = ax.imshow(np.log(np.abs(pressure_fourier)), extent=extent, origin='lower')#not so clear with abs here, that means we are taking the absolute value of the pressure_fourier?
plt.colorbar(im, ax=ax)
plt.title('Log of Absolute Value of Fourier-transformed Pressure')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Material parameters
E = 2.1e2
nu = 0.3
E_star = E / (1 - nu**2)

# Kernel in Fourier space
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)
kernel_fourier = np.zeros_like(QX) #not so convincing with the zeros_like here, this works for a 2D problem?
non_zero_indices = (QX**2 + QY**2) != 0
kernel_fourier[non_zero_indices] = 2 / (E_star * (QX**2 + QY**2)[non_zero_indices])

# Displacement field in Fourier and real space
displacement_fourier = pressure_fourier * kernel_fourier
displacement_real = np.fft.ifft2(displacement_fourier, norm='ortho')

# Analytical solution
def analytical_solution(x, y, x0, y0, L, E, nu):
    E_star = E / (1 - nu**2)
    return np.cos(2*np.pi*(x-x0)/L)*np.cos(2*np.pi*(y-y0)/L) / (np.sqrt(2)*np.pi*E_star)

start_time_analytical = time.process_time()
displacement_analytical = analytical_solution(x, y, x0, y0, L, E, nu)
end_time_analytical = time.process_time()

cpu_time_analytical = end_time_analytical - start_time_analytical
print("CPU time for calling analytical_solution:", cpu_time_analytical, "seconds")

# Since displacement_analytical is a large array, we won't print it directly to avoid clutter.
# Instead, let's visualize the real part of the displacement in real space for comparison.
plt.figure()
plt.imshow(np.real(displacement_real), extent=extent, origin='lower')
plt.colorbar()
plt.title('Real Part of Displacement Field in Real Space')
plt.xlabel('x')
plt.ylabel('y')
plt.show()




from mpl_toolkits.mplot3d import Axes3D

# Calculate the error between the real part of the displacement obtained through FFT and the analytical solution
error = np.abs(np.real(displacement_real) - displacement_analytical)

# Create a 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Create a surface plot of the error
X, Y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, m))
surf = ax.plot_surface(X, Y, error, cmap='viridis', edgecolor='none')

# Add a color bar which maps values to colors
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

# Labels and title
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Error')
ax.set_title('Error between Analytical and Fourier Method Displacement')

plt.show()
