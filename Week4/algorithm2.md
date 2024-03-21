Here's the algorithm described in the provided LaTeX style, formatted into Markdown with LaTeX math expressions:

### Algorithm 2: Polonsky and Keer (1999b) Constrained Conjugate Gradient Algorithm.

**Inputs:** normal total load $W$, surface profile $\mathbf{H}$, tolerance $\varepsilon$, maximum number of iterations $N_{\text{max}}$.

1. Initialize $\mathbf{P}$ with the average constant load initial guess:
   $$ \mathbf{P} \leftarrow \frac{W}{|\partial \beta_p|} $$
2. Set search direction $\mathbf{T}$ to zero:
   $$ \mathbf{T} \leftarrow 0 $$
3. Set $h_{\text{norm}}$ to the norm of $\mathbf{H}$:
   $$ h_{\text{norm}} \leftarrow \|\mathbf{H}\| $$
4. Initialize $G$ to zero and $G_{\text{old}}$ to one:
   $$ G \leftarrow 0, \quad G_{\text{old}} \leftarrow 1 $$
5. Initialize $\delta$ to zero:
   $$ \delta \leftarrow 0 $$
6. Set iteration counter $k$ to zero:
   $$ k \leftarrow 0 $$

**Repeat until** $e < \varepsilon$ **or** $k > N_{\text{max}}$:

7. Identify the set of current points in contact $S$ where $\mathbf{P} > 0$:
   $$ S \leftarrow \text{where}(\mathbf{P} > 0) $$
8. Calculate the gradient $\mathbf{G}$:
   $$ \mathbf{G} \leftarrow \nabla I $$
9. Center on contact zone by updating $\mathbf{G}$:
   $$ \mathbf{G} \leftarrow \mathbf{G} - \text{mean}(\mathbf{G}|_S) $$
10. Normalize $\mathbf{G}$:
    $$ G \leftarrow \left\| \mathbf{G}|_S \right\|^2 $$
11. Update only in the current contact zone $\mathbf{T}|_S$:
    $$ \mathbf{T}|_S \leftarrow \mathbf{G}|_S + \delta \frac{G}{G_{\text{old}}} \mathbf{T}|_S $$
12. Set $\mathbf{R}$ as the product of the mobility matrix $\mathcal{M}$ and $\mathbf{T}$:
    $$ \mathbf{R} \leftarrow \mathcal{M}[\mathbf{T}] $$
13. Center on contact zone by updating $\mathbf{R}$:
    $$ \mathbf{R} \leftarrow \mathbf{R} - \text{mean}(\mathbf{R}|_S) $$
14. Compute the critical step size $\tau$:
    $$ \tau \leftarrow \frac{\mathbf{G}|_S \cdot \mathbf{T}|_S}{\mathbf{R}|_S \cdot \mathbf{T}|_S} $$
15. Update $\mathbf{P}$:
    $$ \mathbf{P} \leftarrow (\mathbf{P} - \tau \mathbf{T})_+ $$
16. Identify the set of inadmissible points $R$ where $\mathbf{P} = 0$ and $\mathbf{G} < 0$:
    $$ R \leftarrow \text{where}(\mathbf{P}=0 \text{ and } \mathbf{G}<0) $$
17. If $R$ is empty, set $\delta$ to one, otherwise to zero:
    $$ \text{if } R = \emptyset \text{ then } \delta \leftarrow 1 \text{ else } \delta \leftarrow 0 $$
18. Apply positive pressure on inadmissible points $\mathbf{P}|_R$:
    $$ \mathbf{P}|_R \leftarrow \mathbf{P}|_R - \tau \mathbf{G}|_R $$
19. Enforce the applied force constraint on $\mathbf{P}$:
    $$ \mathbf{P} \leftarrow \frac{W}{\sum \mathbf{P}} \mathbf{P} $$
20. Compute the error $e$:
    $$ e \leftarrow \frac{\mathbf{P} \cdot (\mathbf{G} - \min(\mathbf{G}))}{W h_{\text{norm}}} $$
21. Increment the iteration counter $k$:
    $$ k = k + 1$$

**After exiting the loop**:

22. Ensure a positive gap by updating \( G \):
    $$ G \leftarrow G - \min(G) $$

   







