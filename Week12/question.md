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