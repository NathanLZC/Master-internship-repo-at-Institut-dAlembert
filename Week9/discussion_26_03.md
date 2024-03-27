[1] **elastic_solver.py**, we didn't normalize *white_noise* when we generate *fft_noise* by *np.fft.fft2*, so we should normalize the surface, and we need to check the surface plot scale before our presentation.

```bash
white_noise = gen.normal(size=phi_values.shape)
fft_noise = np.fft.fft2(white_noise)
filtered_noise = fft_noise * np.sqrt(phi_values) #* np.exp(1j * theta)
surface = np.fft.ifft2(filtered_noise).real*np.sqrt(n*m)/(n*m)# add this np.sqrt(n*m) to implement the filtering algorithm
```
**For C++ solver**, when we call FFTW library, it is better to keep in mind to normalize *backward*:

```bash
Eigen::MatrixXd generateRandomSurface(const Eigen::MatrixXd& phiValues, int n, int m){
    //...
    // FFTW plan for inverse transform
    fftw_plan p_backward = fftw_plan_dft_c2r_2d(n, m, fftOutput, surface.data(), FFTW_ESTIMATE);
    //...
    return surface*std::sqrt(n*m)/(n*m);
}
```

Here, we just normalize *surface* like Python solver.

[2] **Week9/test_algo2.py**, we verify our solution with Hertz, and we need to encapsulate this script into a function for viscoelastic solver.

[3] We are supposed to check our C++ solver(algorithm 1) with Hertz solution for a demi-sphere. We found the problem in surface generation for normalizing the noise(or surface). But now we still face two problems to solve, the first one is that our C++ solution seems not to match Python solution, the second one is Python version seems faster than C++ version. In this week, we propose to develop Python based visoelastic solver with algorithm2.

[4] For debugging, first I am supposed to get more familier with these common tools like *gdb*. Then based on our contact problem solving algorithm, which mainly contains energy minimization process, we summarize a logic to find bugs:

1. We can plot the pressure distribution along the central line of our domain, where **p** represents our solver solution and **p0** and **a** represent analytical maximum pressure and analytical contact area radius(Hertz solution). And comparing solution with Hertz first is always a good and efficient way to verify our solver.
```bash
p0 = (6*W*E_star**2/(np.pi**3*Radius**2))**(1/3)
a = (3*W*Radius/(4*E_star))**(1/3)
plt.plot(x[n//2], P[n//2])
plt.plot(x[n//2], p0*np.sqrt(1 - (x[n//2]-x0)**2 / a**2))
```
2. Plot the proper parameters of our algorithm. For example, in *algorithm2*, *delta* controls whether we need to change our contact point set or not. *np.mean(P>0)* represents the percentage of contact area. *tau* shows the step size and *T* represents searching direction. And the relationship between gradient and integration operator $\mathcal{M}$ is given as:
$$
\nabla I_p(p)=\mathcal{M}\left[p \boldsymbol{e}_3\right] \cdot \boldsymbol{e}_3-h=u_3-h=g
$$
Such that:
```bash
    R, T_fourier  = apply_integration_operator(T, kernel_fourier, h_profile)
    R += h_profile
```
3. We are supposed to be wise about the possibility to have errors, and check our scipt:
```
Which iteration step does error jumps?
How is searching direction updated?
Step size increase or not?
...
```


[5] **Always go back to original paper.** During our develop process, we were stucked several times by the details of algorithms. It is better to follow the algorithm