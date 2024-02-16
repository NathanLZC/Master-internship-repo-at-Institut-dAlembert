import numpy as np


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

# Define the separation h
h_matrix = - (x**2 + y**2) / (2 * R)

# Initial pressure distribution
P = np.full((n, m), W / S)  # Initial guess for the pressure

# Frequency components for the Fourier transform
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

# Define the kernel in the Fourier domain
kernel_fourier = np.zeros_like(QX)
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = np.inf  # Avoid division by zero at the zero frequency

# Initialize variables for the iteration
tol = 1e-6  # Tolerance for convergence
iter_max = 1000  # Maximum number of iterations
k = 0  # Iteration counter
error = np.inf  # Initialize error

while np.abs(error) > tol and k < iter_max:
    # Calculate the gradient G in the Fourier domain and transform it back to the spatial domain
    P_fourier = np.fft.fft2(P, norm='ortho')
    G_fourier = P_fourier * kernel_fourier
    G = np.fft.ifft2(G_fourier, norm='ortho').real - h_matrix
    
    # Update P by subtracting G
    P = P - alpha_0 * G
    
    # Ensure P is non-negative
    P = np.maximum(P, 0)
    
    # Adjust P to satisfy the total load constraint
    alpha_0 = (W - np.sum(P)) / np.size(P)
    P += alpha_0
    
    # Calculate the error for convergence checking
    error = np.linalg.norm(G) / np.linalg.norm(h_matrix)
    
    k += 1  # Increment the iteration counter

'''
# Ensure a positive gap by updating G
G = G - np.min(G)
'''

displacement_fourier = P_fourier * kernel_fourier
displacement = np.fft.ifft2(displacement_fourier, norm='ortho').real

# Plot the results
import matplotlib.pyplot as plt
plt.imshow(displacement, cmap='jet', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement Field')
plt.show()
