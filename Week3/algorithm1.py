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




def find_alpha_0(P, W, alpha_l, alpha_r, tol):
    # approximates a root, R, of f bounded 
    # by a and b to within tolerance 
    # | f(m) | < tol with m the midpoint 
    # between a and b Recursive implementation
    
    # check if a and b bound a root
    def f(alpha):
        return np.mean(P + alpha) - W

    while np.sign(f(alpha_l)) == np.sign(f(alpha_r)):
        alpha_r *= 2
        #raise Exception(f"The scalars alpha_l and alpha_r do not bound a root {np.sign(f(alpha_l))} {np.sign(f(alpha_r))}")
    # get midpoint
    alpha_c = (alpha_l + alpha_r)/2
    

    if np.abs(f(alpha_c)) < tol:
        # stopping condition, report alpha_c as root
        return alpha_c
    elif np.sign(f(alpha_l)) == np.sign(f(alpha_c)):
        # case where m is an improvement on a. 
        # Make recursive call with a = m
        return find_alpha_0(P, W, alpha_c, alpha_r, tol)
    elif np.sign(f(alpha_r)) == np.sign(f(alpha_c)):
        # case where m is an improvement on b. 
        # Make recursive call with b = m
        return find_alpha_0(P, W, alpha_l, alpha_c, tol)


while np.abs(error) > tol and k < iter_max:
    # Calculate the gradient G in the Fourier domain and transform it back to the spatial domain
    P_fourier = np.fft.fft2(P, norm='ortho')
    G_fourier = P_fourier * kernel_fourier
    G = np.fft.ifft2(G_fourier, norm='ortho').real - h_matrix
    
    # Update P by subtracting G
    P = P - G
    
    # Ensure P is non-negative
    P = np.maximum(P, 0)
    
    # Adjust P to satisfy the total load constraint
    alpha_0 = find_alpha_0(P, W/S, -np.max(P), W, tol)
    #alpha_0 = find_alpha_0(P, W, -1e2, 1e2, 1e-6)
    P += alpha_0
    P[P < 0] = 0
    
    # Calculate the error for convergence checking
    error = np.vdot(P, (G - np.min(G))) / (P.size*W) #/ np.linalg.norm(h_matrix)
    print(error, k)
    
    k += 1  # Increment the iteration counter


# Ensure a positive gap by updating G
G = G - np.min(G)


displacement_fourier = P_fourier * kernel_fourier
displacement = np.fft.ifft2(displacement_fourier, norm='ortho').real



'''
# Plot the results
import matplotlib.pyplot as plt
plt.imshow(displacement, cmap='jet', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Displacement (u_z)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement Field')
plt.show()
'''
plt.imshow(P, cmap='jet', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label='Pressure (P)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Pressure Field')
plt.show()



#######################################################
## The following is the code for week 2, hertz solution
#######################################################


import matplotlib.pyplot as plt

F_value = 1e2 

# We define the radius of the elastic sphere as R
R = 1

#define material parameters
E = 1e3  # Young's modulus
nu = 0.3  # Poisson's ratio
E_star = E / (1 - nu**2)  # plane strain modulus 

# We define the half-space domain is L^2
L = 2

#We generate a 2D coordinate space
#n = 100
#m = 100
'''
x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))#notice here that we use endpoint=False to avoid having the last point
'''
x0 = 1
y0 = 1

# We define the distance from the center of the sphere
r = np.sqrt((x-x0)**2 + (y-y0)**2)

#Here we define p0 as the reference pressure
p0 = (6*F_value*E_star**2/(np.pi**3*R**2))**(1/3)
a = (3*F_value*R/(4*E_star))**(1/3)

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
plt.title('Displacement Field')
plt.show()





#######################################################
## The following is comparasion of the two solutions
#######################################################



from mpl_toolkits.mplot3d import Axes3D

# Calculate the error between the real part of the displacement obtained through FFT and the analytical solution
error = np.abs(displacement - np.max(displacement) - u_z)

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
