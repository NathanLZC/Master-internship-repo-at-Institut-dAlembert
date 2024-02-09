import numpy as np
import matplotlib.pyplot as plt
import time


# We define the total loading P as F, we give F as input
F = np.linspace(0, 1e6, 1000)


# We define the radius of the elastic sphere as R
R = 1
# We define the elastic modulus of the sphere as E1, and the elastic modulus of the surface as E2
E1 = 1e5
E2 = 1e3
# We define the Poisson ratio of the sphere as v1, and the Poisson ratio of the surface as v2
v1 = 0.33
v2 = 0.33
# We define the reduced elastic modulus as E_star
E_star = 1 / ((1 - v1**2) / E1 + (1 - v2**2) / E2)
# We define the contact radius as a
a = 3*F*R/(4*E_star)

# We define the half-space domain is L^2
L = 2

#We generate a 2D coordinate space
n = 100
m = 100

x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))#notice here that we use endpoint=False to avoid having the last point

# We define the distance from the center of the sphere
r = np.sqrt(x**2 + y**2)

'''
#Here we define p0 as the reference pressure
p0 = (1/np.pi) * (6*F*(E_star**2)/(R**2))**(1/3) #this is the reference pressure

## We define the analytical displacement field

Analytical_displacement = np.pi*p0*(2*a**2-r**2)/(4*E_star*a) #this is the analytical displacement field


d_total_deformation = a**2 / R

#d_total_deformation

#d_total_deformation2 = (9 * F**2 / (16 * E_star**2 * R))**(1/3)

#d_total_deformation2
'''

# 修正：对每个 F 值计算并可视化位移场
for F_value in F[::100]:  # 示例：每100个F值取一个进行处理和可视化
    a = (3*F_value*R/(4*E_star))**(1/3)  # 接触半径
    p0 = (6*F_value*E_star**2/(np.pi**3*R**2))**(1/3)  # 参考压力
    u_z = np.pi*p0*(2*a**2-r**2)/(4*a*E_star)  # 位移场
    u_z[r > a] = 0  # 超出接触区域的位移设置为0

    # 可视化位移场
    plt.figure(figsize=(6, 5))
    plt.imshow(u_z, extent=(0, L, 0, L), origin='lower')
    plt.colorbar(label='Displacement')
    plt.title(f'Analytical Displacement Field for F={F_value:.1e} N')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()





#Here we define the pressure distribution
def pressure_distribution(x, y, x0, y0, a, p0):
    r2 = (x-x0)**2 + (y-y0)**2
    if r2 <= a**2:        
        return p0*np.sqrt(1-r2/a**2)
    else:
        return 0


x0 = 0
y0 = 0

#We calculate test pressure on the grid
test_pressure = pressure_distribution(x, x0, p0)
#We generate the pressure distribution in Fourier space
pressure_fourier = np.fft.fft2(test_pressure, norm='ortho')
# Visualization of the Fourier-transformed pressure
fig, ax = plt.subplots()
extent = [x.min(), x.max(), y.min(), y.max()]
im = ax.imshow(np.log(np.abs(pressure_fourier)), extent=extent, origin='lower')#not so clear with abs here # That means we take the norm of the complex number
plt.colorbar(im, ax=ax)
plt.title('Log of Absolute Value of Fourier-transformed Pressure')
plt.xlabel('x')
plt.ylabel('y')
plt.show()



#plot the pressure distribution
pressure_img = plt.imshow(P/p0, extent=[x_min + dx/2, x_max, y_min + dx/2, y_max], 
                            origin='lower', interpolation='bicubic', cmap=plt.cm.inferno, vmax = 1)
fig.colorbar(pressure_img, ax=ax[0], label='Pressure', orientation="horizontal")
plt.set_title('Pressure Distribution')
plt.set_xlabel('x')
plt.set_ylabel('y')



#we define the frequency with q_x and q_y
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

kernel_fourier = np.zeros_like(QX)  # Initialize the kernel array
#non_zero_indices = (QX**2 + QY**2) != 0  # Avoid division by zero
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Set the zero frequency component to zero



#Calculate the displacement field in fourrier space
displacement_fourrier = pressure_fourier * kernel_fourier

#We calculate the displacement field in real space
displacement_real = np.fft.ifft2(displacement_fourrier, norm='ortho')


# Since displacement_analytical is a large array, we won't print it directly to avoid clutter.
# Instead, let's visualize the real part of the displacement in real space for comparison.
plt.figure()
plt.imshow(np.real(displacement_real), extent=extent, origin='lower')
plt.colorbar()
plt.title('Real Part of Displacement Field in Real Space')
plt.xlabel('x')
plt.ylabel('y')
plt.show()