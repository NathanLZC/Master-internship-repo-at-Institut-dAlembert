import numpy as np
import matplotlib.pyplot as plt

# Define input parameters
F = np.linspace(0, 1e6, 1000)  # Total loading range
R = 1  # Radius of the elastic sphere
E1 = 1e5  # Elastic modulus of the sphere
E2 = 1e3  # Elastic modulus of the surface
v1 = 0.33  # Poisson ratio of the sphere
v2 = 0.33  # Poisson ratio of the surface
E_star = 1 / ((1 - v1**2) / E1 + (1 - v2**2) / E2)  # Reduced elastic modulus

L = 2  # Define the size of the half-space domain
n = 100  # Grid size in both dimensions

# Generate a 2D coordinate space
x, y = np.meshgrid(np.linspace(-L/2, L/2, n), np.linspace(-L/2, L/2, n))

# Distance from the center of the sphere
r = np.sqrt(x**2 + y**2)

# Corrected pressure distribution function for use with arrays
def pressure_distribution(x, y, a, p0):
    r2 = x**2 + y**2
    pressure = np.zeros_like(r2)
    within_contact_area = r2 <= a**2
    pressure[within_contact_area] = p0 * np.sqrt(1 - r2[within_contact_area] / a**2)
    return pressure

# Fourier transform method for calculating displacement field
def calculate_displacement_in_fourier_space(p0, a, E_star, L, n):
    test_pressure = pressure_distribution(x, y, a, p0)

    # Fourier transform of the pressure distribution
    pressure_fourier = np.fft.fft2(test_pressure)
    
    # Frequency coordinates
    q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
    q_y = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
    QX, QY = np.meshgrid(q_x, q_y)
    Q = np.sqrt(QX**2 + QY**2)
    
    # Define the Green's function in Fourier space
    G = np.zeros_like(Q)
    nonzero = Q != 0
    G[nonzero] = 1 / (E_star * Q[nonzero]**2)
    
    # Calculate the displacement field in Fourier space
    displacement_fourier = pressure_fourier * G
    
    # Inverse Fourier transform to get the displacement field in real space
    displacement_real = np.fft.ifft2(displacement_fourier).real
    
    return displacement_real

# Example calculation for a single F value
F_example = 1e5  # Example load value
a_example = (3*F_example*R/(4*E_star))**(1/3)  # Contact radius for the example
p0_example = (6*F_example*E_star**2/(np.pi**3*R**2))**(1/3)  # Reference pressure for the example

# Calculate the displacement field for the example F value
displacement_real_example = calculate_displacement_in_fourier_space(p0_example, a_example, E_star, L, n)

# Visualization of the real part of the displacement field
plt.figure(figsize=(6, 5))
plt.imshow(displacement_real_example, extent=(-L/2, L/2, -L/2, L/2), origin='lower', cmap='viridis')
plt.colorbar(label='Displacement')
plt.title(f'Real Part of Displacement Field for F={F_example:.1e} N')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
