## MAC Schemes 

Current MAC schemes, along with suggestions by George, was tried to be summarized here. (Updated at 14:02, Apr 12 by Jinghu)



The original continuous equation is:
$$
\frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0,\;\;\text{or }\nabla\cdot\bold{\overline{u}}=0
$$
Let's assume right now we have discretization methods $\bold{C}$ for the continuous equation featuring unknown velocity properties $\overline{\bold{u}}$

The original momentum equation can be simplified to the following terms:
$$
\frac{\partial \bold{u}}{\partial t}+\bold{A}(\bold{u})=-\nabla \bold{\pi}+\nu\Delta\bold{u}+\bold{f}
$$
where $\pi=\bold{P}/\rho$

Let's say we discretize the continuous equation to have matrix operator $\bold{C}$ onto the $\bold{\overline{u}}^{n+1}$ (we have info on the $n$ timestep):
$$
\bold{C}(\bold{\overline{u}}^{n+1})=0
$$
where $n+1$ indicates the next time step.



Taking momentum equation into account, let's say we use $\bold{D}$ as discretization matrix for diffusion term, $\bold{A}$ as discertization matrix for advection term, $\bold{M}$ as discretization matrix for pressure term:
$$
\frac{\bold{\overline{u}}^{n+1}-\bold{\overline{u}}^n}{\Delta t}=-\bold{M}\left(\pi^{n+1},\pi^n\right)+\bold{D}\left(\bold{\overline{u}}^{n+1},\bold{\overline{u}}^{n}\right)-\bold{A}\left(\bold{\overline{u}}^{n}\right)+\bold{f}^n
$$
Both continuous and momentum equations can be converted to a matrix type equations:
$$
\begin{bmatrix}
\bold{C} & 0 \\
\bold{X} & \bold{Y} 
\end{bmatrix}\begin{bmatrix}
\bold{\overline{u}}^{n+1} \\
\pi^{n+1} 
\end{bmatrix}=\begin{bmatrix}
0 \\
\bold{Z}
\end{bmatrix}
$$
where $\bold{X}, \bold{Y}, \bold{Z}$can be derived from the discretized momentum equations, and they are the function of the values at $n$ time step.Solving this matrix equation could allow us to get the properties at $n+1$ time step.







The governing equations that will allow us to solve for the semi-steady state 2D velocity distribution in the cylindrical domain include the continuity equation and two momentum equations. Refering to the proposal, the non-dimensional continuity equation is

$$ \frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0 $$

and the two continuity equations for the $r$ and $z$ components of velocity, respectively, are

$$\frac{1}{Pe}\frac{\partial\overline{u}_r}{\partial Fo} + \overline{u}_r\frac{\partial\overline{u}_r}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_r}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_r}{\partial \overline{r}}\right) - \frac{\overline{u}_r}{\overline{r}^2} + \frac{\partial^2 \overline{u}_r}{\partial \overline{z}^2}\right] - \frac{\partial \overline{p}}{\partial \overline{r}} $$

and

$$\frac{1}{Pe}\frac{\partial\overline{u}_z}{\partial Fo} + \overline{u}_r\frac{\partial\overline{u}_z}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_z}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_z}{\partial \overline{r}}\right) + \frac{\partial^2 \overline{u}_z}{\partial \overline{z}^2}\right] + \frac{Ga}{Re^2} - \frac{\partial \overline{p}}{\partial \overline{z}},$$

where the non-dimensional terms are defined as

$$\overline{u}_r := \frac{u_r}{U_\infty}, \;\; \overline{u}_z := \frac{u_z}{U_\infty}, \;\; \overline{r} := \frac{r}{H}, \;\; \overline{z} := \frac{z}{H}, \;\; \overline{p} := \frac{p}{\rho U_\infty^2}, \;\; Re := \frac{HU_\infty}{\nu}, \;\; Pe := \frac{HU_\infty}{\alpha}, \;\; Ga := \frac{gH^3}{\nu^2}.$$

