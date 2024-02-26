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
def phi(q):
    if q_l <= q < q_r:
        return C
    elif q_r <= q < q_s:
        return C * q ** (-2 * (H + 1))
    else:
        return 0

#We generate a 2D coordinate space
n = 30
m = 30

#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
#QX, QY = np.meshgrid(q_x, q_y)

q_values = np.logspace(np.log10(q_l/10), np.log10(q_s*10), 500)#just as an example, we will overload the range of q values later



phi_values = np.sqrt(phi(q_x)**2 + phi(q_y)**2)

gen = np.random.default_rng()
theta = gen.uniform(0, 2*np.pi, size=phi_values.shape)





# Creating the plot
plt.figure(figsize=(8, 6))
plt.loglog(q_values, phi_values, label=r'$\phi(q)$', color='blue')

# Adding titles and labels
plt.title('Log-Log plot of the function $\\phi(q)$')
plt.xlabel('$q$')
plt.ylabel('$\\phi(q)$')

# Adding a grid
plt.grid(which='both', linestyle='--', linewidth=0.5)

# Show the plot
plt.show()