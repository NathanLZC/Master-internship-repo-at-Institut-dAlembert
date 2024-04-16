### This script is for the Maxwell one branch model.
### Deduce process is in generalized_Maxwell_backward_Euler.ipynb

import numpy as np
import matplotlib.pyplot as plt
#define input parameters
##time
t0 = 0
t1 = 1
dt = (t1 - t0)/50
##load(constant)
W = 1e0  # Total load

#domain size
R = 1  # Radius of demi-sphere
L = 2  # Domain size
Radius = 0.5
S = L**2  # Domain area

# Generate a 2D coordinate space
n = 300
m = 300

x, y = np.meshgrid(np.linspace(0, L, n, endpoint=False), np.linspace(0, L, m, endpoint=False))

x0 = 1
y0 = 1

E = 3  # Young's modulus
nu = 0.5
E_star = E / (1 - nu**2)  # Plane strain modulus

##################################################################
#####First just apply for demi-sphere and compare with Hertz######
##################################################################

# We define the distance from the center of the sphere
r = np.sqrt((x-x0)**2 + (y-y0)**2)

# Define the kernel in the Fourier domain
q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
QX, QY = np.meshgrid(q_x, q_y)

kernel_fourier = np.zeros_like(QX)
kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
kernel_fourier[0, 0] = 0  # Avoid division by zero at the zero frequency

h_profile = -(r**2)/(2*Radius)

def apply_integration_operator(Origin, kernel_fourier, h_profile):
    # Compute the Fourier transform of the input image
    Origin2fourier = np.fft.fft2(Origin, norm='ortho')

    Middle_fourier = Origin2fourier * kernel_fourier

    Middle = np.fft.ifft2(Middle_fourier, norm='ortho').real

    Gradient = Middle - h_profile

    return Gradient, Origin2fourier#true gradient

##define our elastic solver with constrained conjuagte gradient method

## E_star seems no use here
def contact_solver(n, m, W, S, h_profile, tol=1e-6, iter_max=200):
    

    # Initial pressure distribution
    P = np.full((n, m), W / S)  # Initial guess for the pressure

    #initialize the search direction
    T = np.zeros((n, m))

    #set the norm of surface(to normalze the error)
    h_rms = np.std(h_profile)

    #initialize G_norm and G_old
    G_norm = 0
    G_old = 1

    #initialize delta
    delta = 0

    # Initialize variables for the iteration
    k = 0  # Iteration counter
    error = np.inf  # Initialize error
    h_rms = np.std(h_profile)

    while np.abs(error) > tol and k < iter_max:
        S = P > 0

        G, P_fourier = apply_integration_operator(P, kernel_fourier, h_profile)

        G -= G[S].mean()

        G_norm = np.linalg.norm(G[S])**2

        # Calculate the search direction
        T[S] = G[S] + delta * G_norm / G_old * T[S]
        T[~S] = 0  ## out of contact area, dont need to update

        # Update G_old
        G_old = G_norm

        # Set R
        R, T_fourier  = apply_integration_operator(T, kernel_fourier, h_profile)
        R += h_profile
        R -= R[S].mean()

        # Calculate the step size tau
        tau = np.vdot(G[S], T[S]) / np.vdot(R[S], T[S])

        # Update P
        P -= tau * T        
        P *= P > 0

        # identify the inadmissible points
        R = (P == 0) & (G < 0)

        if R.sum() == 0:
            delta = 1
        else:
            delta = 0#change the contact point set and need to do conjugate gradient again

        # Enforce the applied force constraint
        P = W * P / np.mean(P) / L**2  

        # Calculate the error for convergence checking
        error = np.vdot(P, (G - np.min(G))) / (P.sum()*h_rms) 
        print(delta, error, k, np.mean(P), np.mean(P>0), tau)
        
        k += 1  # Increment the iteration counter

    # Ensure a positive gap by updating G
    G = G - np.min(G)

    displacement_fourier = P_fourier * kernel_fourier
    displacement = np.fft.ifft2(displacement_fourier, norm='ortho').real

    return displacement, P

##################################################################
#####shear modulus for one branch Maxwell model###################
##################################################################
G_inf = 2.75 #elastic branch
G_1 = 2.75 

print('G_inf:', G_inf, ' G_1:', G_1)

# Define the relaxation times
tau_1 = .1
eta_1 = G_1 * tau_1

##################################################################
#####define G_tilde for one-branch Maxwell model #################
##################################################################

G_tilde = tau_1 * G_1 / (tau_1 + dt)
#for multi-branch, we have $\tilde{G}=\sum_k \frac{\tau^k}{\tau^k+\Delta t} G^k$
#G_tilde = 0
#for k in range(n):
#    G_tilde += tau[k] / (tau[k] + dt) * G[k]


# Define parameters for updating the surface profile
alpha = G_inf + G_tilde
beta = G_tilde
gamma_k = tau_1 / (tau_1 + dt)#k=1 for one-branch model
#gamma_k = tau[k]/(tau[k] + dt)

Surface = h_profile

U = np.zeros((n, m))
M = np.zeros((n, m))

Ac=[]

for t in np.arange(t0, t1, dt):
    #Update the surface profile
    H_new = alpha*Surface - beta*U + gamma_k*M

    #main step1: Compute $P_{t+\Delta t}^{\prime}$
    #M_new, P = contact_solver(n, m, W, S, H_new, tol=1e-6, iter_max=200)
    







    #main step2: Update the l

    Ac.append(np.mean(P > 0)*S)

    U_new = (1/alpha2)*(M_new - gamma*M + beta*U)

    M = M_new
    #main step4: Update the total displacement field
    U = U_new




#######################################
###Hertzian contact theory reference
#######################################
##Hertz solution at t0    
E_effective_t0 = 2*G_inf*(1+nu)/(1-nu**2)

p0_t0 = (6*W*(E_effective_t0)**2/(np.pi**3*Radius**2))**(1/3)
a_t0 = (3*W*Radius/(4*(E_effective_t0)))**(1/3)
##Hertz solution at t_inf
E_effective_inf = 2*(G_inf+G_1)*(1+nu)/(1-nu**2)

p0_t_inf = (6*W*(E_effective_inf)**2/(np.pi**3*Radius**2))**(1/3)
a_t_inf = (3*W*Radius/(4*(E_effective_inf)))**(1/3)


plt.plot(x[n//2], P[n//2])
plt.plot(x[n//2], p0_t0*np.sqrt(1 - (x[n//2]-x0)**2 / a_t0**2))
plt.plot(x[n//2], p0_t_inf*np.sqrt(1 - (x[n//2]-x0)**2 / a_t_inf**2))
plt.legend(["Numerical", "Hertz at t=0", "Hertz at t=inf"])
plt.xlabel("x")
plt.ylabel("Pressure distribution")
plt.show()

Ac_hertz_t0 = np.pi*a_t0**2
Ac_hertz_t_inf = np.pi*a_t_inf**2

print("Analytical contact area radius at t0:", a_t0)
print("Analytical contact area radius at t_inf:", a_t_inf)
print("Analytical maximum pressure at t0:", p0_t0)
print("Analytical maximum pressure at t_inf:", p0_t_inf)
print("Numerical contact area at t0:", Ac[0])
print("Numerical contact area at t_inf",  Ac[-1])
print("Analyical contact area at t0:", Ac_hertz_t0)
print("Analyical contact area at t_inf:", Ac_hertz_t_inf)
plt.plot(np.arange(t0, t1, dt), Ac)
plt.axhline(Ac_hertz_t0, color='red', linestyle='dotted')
plt.axhline(Ac_hertz_t_inf, color='blue', linestyle='dotted')
plt.xlabel("Time(s)")
plt.ylabel("Contact area($m^2$)")
plt.legend(["Numerical", "Hertz at t=0", "Hertz at t=inf"])
#define a title that can read parameter tau_0
plt.title("Contact area vs time for a single branch viscoelastic model with tau_0 = " + str(tau_0) + "s")
#plt.axhline(Ac_hertz_t_inf, color='blue')
plt.show()