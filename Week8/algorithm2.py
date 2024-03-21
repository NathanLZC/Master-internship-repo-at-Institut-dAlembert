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
n = 300
m = 300
x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))

x0 = 1
y0 = 1

##################################################################
#####First just apply for demi-sphere and compare with Hertz######
##################################################################

# We define the distance from the center of the sphere
r = np.sqrt((x-x0)**2 + (y-y0)**2)

#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

# Define the kernel in the Fourier domain
kernel_fourier = np.zeros_like(QX)
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Avoid division by zero at the zero frequency


h_profile = -(r**2)/(2*R)



# Initial pressure distribution
P = np.full((n, m), W / S)  # Initial guess for the pressure


#initialize the search direction
T = np.zeros((n, m))

#set the norm of surface(to normalze the error)
h_rms = np.std(h_profile)

#initialize G and G_old
G = np.zeros((n, m))
G_old = np.ones((n, m))

#initialize delta
delta = 0

# Initialize variables for the iteration
tol = 1e-6  # Tolerance for convergence
iter_max = 10000  # Maximum number of iterations
k = 0  # Iteration counter
error = np.inf  # Initialize error





while np.abs(error) > tol and k < iter_max:
    
    P_fourier = np.fft.fft2(P, norm='ortho')
    G_fourier = P_fourier * kernel_fourier
    G = np.fft.ifft2(G_fourier, norm='ortho').real - h_profile













    # Calculate the error for convergence checking
    error = np.vdot(P, (G - np.min(G))) / (h_profile.size*h_rms*W) #/ np.linalg.norm(h_matrix)
    print(error, k, np.mean(P))
    
    k += 1  # Increment the iteration counter


# Ensure a positive gap by updating G
G = G - np.min(G)




















#######################################
###Hertzian contact theory reference
#######################################

#Here we define p0 as the reference pressure
p0 = (6*W*E_star**2/(np.pi**3*R**2))**(1/3)
a = (3*W*R/(4*E_star))**(1/3)

u_z = -(r**2)/(2*R)

# Correctly applying the displacement outside the contact area
outside_contact = r > a
u_z_outside = -(r[outside_contact]**2)/(2*R) + a * np.sqrt(r[outside_contact]**2 - a**2)/(np.pi*R) + (r[outside_contact]**2-2*a**2)*np.arccos(a/r[outside_contact])/(np.pi*R)
u_z[outside_contact] = u_z_outside

#plot u_z
plt.imshow(u_z, cmap='jet', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement_Hertz Field')

#plot the difference between the two displacement fields
plt.figure(figsize=(10, 8))
plt.imshow(displacement - np.max(displacement) - u_z, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement Field Difference')


plt.show()