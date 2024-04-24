22/04/2024 14:00 Jessieu 44-54 518

In this discussion, we first review our development process, and decided to focus on finishing Generalized Maxwell model(one-branch, multi-branch) instead of generalized zener model. In our report, we will start with viscoelasticity knowledge background and explain our algorithm and results for Generalized Maxwell model(one-branch, multi-branch).

#### First, we fix the bug in our codes for one-branch Maxwell model:

##### Update surface profile $H_{t+\Delta t}$:

$$
\mathbb{H}^\prime = (\underbrace{G_{\infty}+\tilde{G}}_\alpha) \mathbb{H}-\underbrace{\tilde{G}}_\beta \mathbb{U_t}+\sum_k \underbrace{\frac{\tau^k}{\tau^k+\Delta t}}_{\gamma^k} \mathbb{M}_t^k
$$

So in our codes for one-branch model, we should have $\gamma^k \mathbb{M}_t^k (k=1)$.

```bash
H_new = alpha*Surface - beta*U + gamma_k*M_k_new
```

Another thing that needs to be considered is, in Python, variable assignment essentially binds the variable name to an object in memory. When a variable is reassigned, the original object is unbound, and the variable name is linked to a new object. If the original object is no longer referenced by any variable, its reference count drops to zero, making it eligible for garbage collection. This object will be reclaimed by Python's garbage collector at an appropriate time. This mechanism enables Python to manage memory efficiently, avoiding the complexities of manual memory management and preventing memory leaks from objects that are no longer needed.

Here in our case, we initialize **M_k_new** = $A_{z z} P_t^k$ to be zero, and we can update this parameter by:
```bash
M_k_new = gamma_k*(M_k_new + G_1*(U_new - U))
```
Then in next step, we can reassign the variable **H_new** without the risk of abnormal accumulation of variables.

#### Second, we move to multi-branch Generalized Maxwell model

We start from two Maxwell branches and set them the same shear modulus but different relaxation time, such that we can verify our model easily with one equivalent Maxwell branch.

To clarify **M**, **M_new**, **M_maxwell**, we should explain the so-called *partial displacement* and *global displacement*.

##### Explaination for partial displacement and global displacement 

In equations, we can define $\mathcal{M}[p] = A_{z z} \sum_k P_{t+\Delta t}^k+P_{t+\Delta t}^e$, where $\mathcal{M}$ is an integration operator in Lucas(2020)[1], $\mathcal{M}[p]$ is the so-called *partial displacement*. Here in our codes, **M** and **M_new** computes this variable. The reason that we call "partial" is that the partial convolution operation only consider the horizontal coordinates $\tilde{\boldsymbol{x}}=\left(x_1, x_2\right)$ and this integral doesn't contain modulus(Young's modulus, shear modulus...). In one word, the *partial displacement* represents the reaction of normal pressure and not strictly equivalent to displacement.

Global displacement is always updated by:

$$
A_{z z}\{\underbrace{\sum_k P_{t+\Delta t}^k+P_{t+\Delta t}^e}_{\text{Generated from elastic solver}}\}-\sum_k {\frac{\tau^k}{\tau^k+\Delta t}} A_{z z} P_t^k+\tilde{G} U_t=\left(G_{\infty}+\tilde{G}\right) U_{t+\Delta t} \quad \text{(4)}
$$

which represents the reaction of elastic branch.

Since we are dealing with multi-branch Maxwell model, **M** should be a 3D numpy array:

```bash
M = np.zeros((len(G), n, m))
```

And for each branch, **M_maxwell** has the same data structure as **U**:

```bash
M_maxwell = np.zeros_like(U)
```
Since we talked about the strategy to assign variables in Python, we should set **M_maxwell** to be zero in for loop to avoid erroneously accumulating unnecessary variables:

```bash
M_maxwell[:] = 0
```
The use of [:] is a common Python idiom for applying an operation to all elements of an array. 

Then for each banch, we can update the pressure:

```bash
for k in range(len(G)):
    M[k] = gamma[k]*(M[k] + G[k]*(U_new - U))
```

The so-called sanity check is to check whether **M_new** = **M_k_new** + **$M^e$** is correct or not, where **$M^e$** is $A_{z z} P_{t+\Delta t}^e$, which is equivalent to $G_{\infty} U_{t+\Delta t}$ .



### Task list

[1] Add data point at $t=t_0$, we can compute one solution for $\eta=1e6$ or we can regard the generalized Maxwell model only has elastic elements on each branch.

[2] Follow the handnotes to add two auxiliary lines with slopes $\frac{1}{\tau_1}$ and $\frac{1}{\tau_2}$, where $\tau_1 < \tau_2$. And we should also have a log-scale plot to have a better visualization of slope change.

[3] Do longer simulation for two branches with different relaxation time(0.1s and 1s), and check whether the transition at t = $\tau $ could reach 63% of the Hertz solution for $t=t_{\infty}$ or not. We can also try (0.01s and 1s). We can explain with Ref[4] and the figure in handnotes that our numerical results can never reach the limitation by Hertz solution at $t=t_{\infty}$.



For next step, we have two trends: contribute Tamaas[2] and modify our solver to simulate stress mutation in Ref[3].

We plan to talk about how to play with Gitlab to contribute Tamaas project on Frdiay(26/04/2024). 

#### Reference

[1] Frérot, Lucas Henri Galilée. ‘Bridging Scales in Wear Modeling with Volume Integral Methods for Elastic-Plastic Contact’, n.d.

[2] https://tamaas.readthedocs.io/en/latest/index.html

[3] Dillavou, Sam, and Shmuel M. Rubinstein. ‘Nonmonotonic Aging and Memory in a Frictional Interface’. Physical Review Letters 120, no. 22 (1 June 2018): 224101. https://doi.org/10.1103/PhysRevLett.120.224101.

[4] Marques, Severino P. C., and Guillermo J. Creus. Computational Viscoelasticity. SpringerBriefs in Applied Sciences and Technology. Berlin, Heidelberg: Springer Berlin Heidelberg, 2012. https://doi.org/10.1007/978-3-642-25311-9.


