import numpy as np
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
n = 100
m = 100

x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))#notice here that we use endpoint=False to avoid having the last point

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