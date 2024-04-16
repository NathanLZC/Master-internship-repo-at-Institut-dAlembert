[1] One detail about writing report(viscoelasticity) and notebook(one_branch_euler_backward.ipynb, generalized_Maxwell_backward_Euler.ipynb):

$$
\mathbb{H}_{t+\Delta t} = (\underbrace{G_{\infty}+\tilde{G}}_\alpha) \mathbb{H}_t-\underbrace{\tilde{G}}_\beta \mathbb{U_t}+\sum_k \underbrace{\frac{\tau^k}{\tau^k+\Delta t}}_{\gamma^k} \mathbb{M}_t^k
$$

So?

$$
\mathbb{H}_{t+\Delta t} or \mathbb{H}^\prime
\mathbb{H}_{t} or \mathbb{H}
$$

[2] From elastic solver, we can return **displacement** and **pressure**, shall we get **diaplacement** as partial displacement $\mathbb{M_t^k}$ for updating global displacement and use **pressure** for sanity check?

[3] explain what exactly "Partial displacement" mean

[4] explain why we can not update global displacement by updating partial displacement(I think that could be: Using global displacement will cause some parts to exceed the partial displacement, while the backward euler method is unconditionally stable and will cause errors to accumulate.) 

[5] does finer time discretization could fix the gap at t_0?