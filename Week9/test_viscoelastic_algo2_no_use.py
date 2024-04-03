#Script to implement the viscoelastic contact with algorithm 2, one branch model
#First compare with the hertz solution
import numpy as np
import matplotlib.pyplot as plt

# Define the material properties
#first we only consider one branch
G_0 = 2.75  # MPa
G_1 = 2.75  # MPa

G_inf = 1/(1/G_0 + 1/G_1)  # MPa

tau_0 = 0.5  # s
eta_1 = 0#G_1 * tau_0  # Characteristic time

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


#not to set E and nu for our elastic solver in Laplace domain
'''
E = 1.0  # Young's modulus
nu = 0.5
E_star = E / (1 - nu**2)  # Plane strain modulus
'''

nu = 0.3




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




h_profile = -(r**2)/(2*Radius)

def apply_integration_operator(Origin, kernel_fourier, h_profile):
    # Compute the Fourier transform of the input image
    Origin2fourier = np.fft.fft2(Origin, norm='ortho')

    Middle_fourier = Origin2fourier * kernel_fourier

    Middle = np.fft.ifft2(Middle_fourier, norm='ortho').real

    Gradient = Middle - h_profile

    return Gradient, Origin2fourier#true gradient


def contact_solver(n, m, W, S, E_star, h_profile, tol=1e-6, iter_max=200):
    

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

    # Define the kernel in the Fourier domain
    q_x = 2 * np.pi * np.fft.fftfreq(n, d=L/n)
    q_y = 2 * np.pi * np.fft.fftfreq(m, d=L/m)
    QX, QY = np.meshgrid(q_x, q_y)

    kernel_fourier = np.zeros_like(QX)
    kernel_fourier = 2 / (E_star * np.sqrt(QX**2 + QY**2))
    kernel_fourier[0, 0] = 0  # Avoid division by zero at the zero frequency


    while np.abs(error) > tol and k < iter_max:
    # try np.where(P > 0) to find the contact area
        S = P > 0

        #encapsulate into a function

        # Calculate the gap G(as Gradient, see Lucas(2020)) in the Fourier domain and transform it back to the spatial domain
        #P_fourier = np.fft.fft2(P, norm='ortho')
        #G_fourier = P_fourier * kernel_fourier
        #G = np.fft.ifft2(G_fourier, norm='ortho').real - h_profile
        ##function

        G, P_fourier = apply_integration_operator(P, kernel_fourier, h_profile)

        G -= G[S].mean()

        G_norm = np.linalg.norm(G[S])**2

        # Calculate the search direction
        T[S] = G[S] + delta * G_norm / G_old * T[S]
        T[~S] = 0  ## out of contact area, dont need to update
        ## size dont match

        # Update G_old
        G_old = G_norm

        # Set R
        R, T_fourier  = apply_integration_operator(T, kernel_fourier, h_profile)
        R += h_profile
        R -= R[S].mean()

        # Calculate the step size tau
        #######
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

        # Apply positive pressure on inadmissible points       
        #P[R] -= tau * G[R]


        # Enforce the applied force constraint
        P = W * P / np.mean(P) / L**2  ## be wise here#############
        ##############################



        # Calculate the error for convergence checking
        error = np.vdot(P, (G - np.min(G))) / (P.sum()*h_rms) 
        print(delta, error, k, np.mean(P), np.mean(P>0), tau)
        
        k += 1  # Increment the iteration counter


    # Ensure a positive gap by updating G
    G = G - np.min(G)

    displacement_fourier = P_fourier * kernel_fourier
    displacement = np.fft.ifft2(displacement_fourier, norm='ortho').real

    return displacement, P





#######################################
###if we let k=1, we can compare the real contact area with hertz solution at t=0 and t>>\tau_0
#######################################
k = 1

alpha = G_inf + (G_1 + eta_1/dt)/(1 + G_1/G_0 + eta_1/G_0/dt)
beta = (eta_1/dt)/(1+G_1/G_0+eta_1/G_0/dt)
gamma = (eta_1/G_0/dt)/(1+G_1/G_0+eta_1/G_0/dt)

Surface = h_profile

U = np.zeros((n, m))
M = np.zeros((n, m))

Ac=[]

for t in np.arange(t0, t1, dt):
    #main step0: Update the effective modulus


    #effictive modulus
    G_t = G_inf + (G_0 - G_inf) * np.exp(-t/tau_0)

    E_star = 2*G_t*(1+nu)
    


    #main step1: Update the surface profile
    H_new = alpha*Surface - beta*U + gamma*M

    #main step2: Update the displacement field
    U_new, P = contact_solver(n, m, W, S, E_star, H_new, tol=1e-6, iter_max=200)
    ###const no need to update the loading field W_new


    Ac.append(np.mean(P > 0)*S)



    #main step3: Update the partial displacement field
    M = (G_0*dt/((G_0+G_1)*dt+eta_1))*(eta_1*M /G_0/dt + (G_1+eta_1/dt)*U_new -eta_1*U/dt)

    #main step4: Update the total displacement field
    U = U_new
    #######################################
    #check and need to correct backward_euler.py
    #######################################

#######################################
###Hertzian contact theory reference
#######################################
    

#########################
#Question!!!!t=0, we should apply the hertzian contact theory 
#with G = G_inf + G_1 as effective modulus???
#Why not just G_0 as handwritten notes in paper one?
#########################
    
#Here we define p0 as the reference pressure
##Hertz solution at t0    
E_effective_t0 = 2*G_0*(1+nu)

p0_t0 = (6*W*(E_effective_t0)**2/(np.pi**3*Radius**2))**(1/3)
a_t0 = (3*W*Radius/(4*(E_effective_t0)))**(1/3)

##Hertz solution at t_inf
##Hertz solution at t_inf follows the Prony series in https://en.wikipedia.org/wiki/Viscoelasticity#cite_note-VanVliet-5
##We can use the Prony series to calculate the effective modulus at t_inf
##the following is from Marques, Severino P. C., and Guillermo J. Creus. Computational Viscoelasticity. 
##Then, reference is made to a rheological model (generalized Maxwell) and Prony series are introduced as its representation.
##To clarify the different notations, 
E_effective_inf = 2*G_inf*(1+nu)


p0_t_inf = (6*W*(E_effective_inf)**2/(np.pi**3*Radius**2))**(1/3)
a_t_inf = (3*W*Radius/(4*(E_effective_inf)))**(1/3)


plt.plot(x[n//2], P[n//2])
plt.plot(x[n//2], p0_t0*np.sqrt(1 - (x[n//2]-x0)**2 / a_t0**2))
plt.plot(x[n//2], p0_t_inf*np.sqrt(1 - (x[n//2]-x0)**2 / a_t_inf**2))
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


print(Ac)



plt.plot(np.arange(t0, t1, dt), Ac)
#add different color for the hertz solution at t0 and t_inf
plt.axhline(Ac_hertz_t0, color='red')
plt.axhline(Ac_hertz_t_inf, color='blue')
#plt.axhline(Ac_hertz)
plt.show()

