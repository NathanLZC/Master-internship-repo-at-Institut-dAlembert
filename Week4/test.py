import numpy as np
import matplotlib.pyplot as plt

L = 2

# Constants for the piecewise function
C = 1 # represents phi_0
q_l = 2*np.pi/L
q_r = 2*np.pi/L
q_s = 2*np.pi*25/L
H = 0.8 # Hurst exponent, represents the roughness of the surface

# Defining the piecewise function
'''
def phi(q):
    if q_l <= q < q_r:
        return C
    elif q_r <= q < q_s:
        return C * q ** (-2 * (H + 1))
    else:
        return 0
'''
def phi(q):
    # Initialize the result with zeros
    result = np.zeros_like(q)
    
    # Apply conditions
    mask1 = (q_l <= q) & (q < q_r)  # Condition for C
    mask2 = (q_r <= q) & (q < q_s)  # Condition for the power-law decay
    
    result[mask1] = C
    result[mask2] = C * (q[mask2] / q_r) ** (-2 * (H + 1))
    
    return result

#We generate a 2D coordinate space
n = 300
m = 300

#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

q_values = np.sqrt(QX**2 + QY**2)



phi_values = phi(q_values)

# Generate random phase
gen = np.random.default_rng()
theta = gen.uniform(0, 2*np.pi, size=phi_values.shape)

# Generate white noise and apply PSD and phase
white_noise = gen.normal(size=phi_values.shape)
fft_noise = np.fft.fft2(white_noise)
filtered_noise = fft_noise * np.sqrt(phi_values) #* np.exp(1j * theta)
surface = np.fft.ifft2(filtered_noise).real

# Plotting
plt.figure(figsize=(10, 8))
plt.imshow(surface, cmap='viridis')
plt.colorbar()
plt.title("Synthesized Rough Surface")
plt.show()



