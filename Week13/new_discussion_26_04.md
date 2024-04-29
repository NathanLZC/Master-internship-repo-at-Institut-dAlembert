26/04/2024 Jussieu 44-54 518

Not exponential(not 63% for t=$\tau$, it is around 78%)




In **slope_auxiliary_line.py**, we expect to observe that our reference auxiliary line could fit the tangent line at $t=0$ and $t=t_{\infty}$, we do a better discritization for our 2D surface and set simulation time $t=\tau$:

```bash
t1 = 1 #If we do a long-term, say 10 seconds, we expect to see the slope change from 1/tau1 to 1/tau2 then to zero
n=512
m=512
```

To reduce the effect by other branches(possible superposition effect from other Maxwell branches), we play with only one branch:

```bash
G = [2.75]
tau = [1]
```
For the reference auxiliary lines, with respect to the following relationship:
$$
\begin{aligned}
A_c(t) & =A_{\infty}+A_0 e^{-t / \tau} \\
& =A_{\infty}\left(1+\frac{A_0}{A_{\infty}} e^{-t / \tau}\right)
\end{aligned}
$$
we can define:

```bash
slope_t0_ref = 1 / np.min(tau)*Ac_hertz_t_inf
slope_t1_ref = 1 / np.max(tau)*Ac_hertz_t_inf
```

To have a better view of our results, we should use *plt.axline* to avoid extending our plotting zone:

```bash
plt.axline([t0,Ac[0]], slope=slope_t0_ref, color='r', linestyle='dotted', label='Tangent at Start (Ref)')
plt.axline([t1,Ac[-1]], slope=slope_t1_ref, color='b', linestyle='dotted', label='Tangent at End (Ref)')
```




