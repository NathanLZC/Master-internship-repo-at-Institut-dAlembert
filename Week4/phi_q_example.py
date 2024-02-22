import numpy as np
import matplotlib.pyplot as plt

L = 2

# Constants for the piecewise function
C = 1 # Just as an example, C can be any constant
q_l = 2*np.pi/L
q_r = 2*np.pi/L
q_s = 2*np.pi*25/L
H = 0.8 # Hurst exponent, just as an example

# Defining the piecewise function
def phi(q):
    if q_l <= q < q_r:
        return C
    elif q_r <= q < q_s:
        return C * (q / q_r) ** (-2 * (H + 1))
    else:
        return 0

# Generating a range of q values for the plot
q_values = np.logspace(np.log10(q_l/10), np.log10(q_s*10), 500)
phi_values = np.array([phi(q) for q in q_values])

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
