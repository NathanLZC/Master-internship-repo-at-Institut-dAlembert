#######################################################
#####Contact problem, complementary energy minimization
#######################################################
import numpy as np
import matplotlib.pyplot as plt

# Define the parameters
W = 1e2  # Total load
R = 1  # Radius of demi-sphere
L = 2  # Domain size
S = L**2  # Domain area

# Material parameters
E = 1e3  # Young's modulus
nu = 0.3  # Poisson's ratio
E_star = E / (1 - nu**2)  # Plane strain modulus

# Generate a 2D coordinate space
n = 100
m = 100



x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))



x0 = 1
y0 = 1


# Define the separation h
h_matrix = - ((x-x0)**2 + (y-y0)**2) / (2 * R)





# Initial pressure distribution
P = np.full((n, m), W / S)  # Initial guess for the pressure

# Frequency components for the Fourier transform
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

# Define the kernel in the Fourier domain
kernel_fourier = np.zeros_like(QX)
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Avoid division by zero at the zero frequency

# Initialize variables for the iteration
tol = 1e-6  # Tolerance for convergence
iter_max = 10000  # Maximum number of iterations
k = 0  # Iteration counter
error = np.inf  # Initialize error

