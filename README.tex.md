
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
The governing equations that will allow us to solve for the semi-steady state 2D velocity distribution in the cylindrical domain include the continuity equation and two momentum equations. Refering to the proposal, the non-dimensional continuity equation is

$$ \frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0 $$

and the two continuity equations for the $r$ and $z$ components of velocity, respectively, are

$$\frac{1}{Pe}\frac{\partial\overline{u}_r}{\partial Fo} + \overline{u}_r\frac{\partial\overline{u}_r}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_r}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_r}{\partial \overline{r}}\right) - \frac{\overline{u}_r}{\overline{r}^2} + \frac{\partial^2 \overline{u}_r}{\partial \overline{z}^2}\right] - \frac{\partial \overline{p}}{\partial \overline{r}} $$

and

$$\frac{1}{Pe}\frac{\partial\overline{u}_z}{\partial Fo} + \overline{u}_r\frac{\partial\overline{u}_z}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_z}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_z}{\partial \overline{r}}\right) + \frac{\partial^2 \overline{u}_z}{\partial \overline{z}^2}\right] + \frac{Ga}{Re^2} - \frac{\partial \overline{p}}{\partial \overline{z}},$$

where the non-dimensional terms are defined as

$$\overline{u}_r := \frac{u_r}{U_\infty}, \;\; \overline{u}_z := \frac{u_z}{U_\infty}, \;\; \overline{r} := \frac{r}{H}, \;\; \overline{z} := \frac{z}{H}, \;\; \overline{p} := \frac{p}{\rho U_\infty^2}, \;\; Re := \frac{HU_\infty}{\nu}, \;\; Pe := \frac{HU_\infty}{\alpha}, \;\; Ga := \frac{gH^3}{\nu^2}.$$

---

Updated on Apr 11 13:57 by Jinghu, to compare with the original dimensionless equations.



Original continuous equation at the cylindrical coordinates (here $u_\theta$ and $\partial/\partial\theta$ is cancelled out):
$$
\frac{1}{r} \frac{\partial\left(r u_{r}\right)}{\partial r}+\frac{\partial u_{z}}{\partial z}=0
$$
Original momentum equations:
$$
\begin{eqnarray}
\begin{array}{2}
\rho\left(\frac{\partial u_{r}}{\partial t}+u_{r} \frac{\partial u_{r}}{\partial r}+u_{z} \frac{\partial u_{r}}{\partial z}\right) = -\frac{\partial P}{\partial r}+\mu\left[\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial u_{r}}{\partial r}\right)-\frac{u_{r}}{r^{2}}+\frac{\partial^{2} u_{r}}{\partial z^{2}}\right]\\
\rho\left(\frac{\partial u_{z}}{\partial t}+u_{r} \frac{\partial u_{z}}{\partial r}+u_{z} \frac{\partial u_{z}}{\partial z}\right) = -\frac{\partial P}{\partial z}+\rho g_{z}+\mu\left[\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial u_{z}}{\partial r}\right)+\frac{\partial^{2} u_{z}}{\partial z^{2}}\right]\\
\end{array}
\end{eqnarray}
$$


The non-dimensionalization process takes the following conversions:
$$
\overline{u}_r := \frac{u_r}{U_\infty}, \;\; \overline{u}_z := \frac{u_z}{U_\infty}, \;\; \overline{r} := \frac{r}{H}, \;\; \overline{z} := \frac{z}{H}, \;\; \overline{p} := \frac{p}{\rho U_\infty^2}, \;\; Re := \frac{HU_\infty}{\nu}, \;\; \overline{t}=\frac{U_\infty t}{H},\;\;Fr := \frac{U_\infty}{\sqrt{gH}}.
$$
The dimensionless continuous equation then is:
$$
\frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0
$$
The dimensionless radial direction momentum equation is:
$$
\frac{\partial\overline{u}_r}{\partial \overline{t}} + \overline{u}_r\frac{\partial\overline{u}_r}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_r}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_r}{\partial \overline{r}}\right) - \frac{\overline{u}_r}{\overline{r}^2} + \frac{\partial^2 \overline{u}_r}{\partial \overline{z}^2}\right] - \frac{\partial \overline{p}}{\partial \overline{r}}
$$


And dimensionless axial direction momentum equation is:
$$
\frac{\partial\overline{u}_z}{\partial \overline{t}} + \overline{u}_r\frac{\partial\overline{u}_z}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_z}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_z}{\partial \overline{r}}\right) + \frac{\partial^2 \overline{u}_z}{\partial \overline{z}^2}\right] - \frac{1}{Fr^2} - \frac{\partial \overline{p}}{\partial \overline{z}},
$$

---



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



## Projection methods

Another possible method is called projectio method. The basic algorithm is shown here.



First the original momentum equation is the one to use:
$$
\frac{\partial \bold{u}}{\partial t}+(\bold{u}\cdot\nabla)(\bold{u})=-\nabla \bold{\pi}+\nu\Delta\bold{u}+\bold{f}
$$
First the intermediate velocity is calculated:
$$
\frac{\mathbf{u}^{*}-\mathbf{u}^{n}}{\Delta t}=-\left(\mathbf{u}^{n} \cdot \nabla\right) \mathbf{u}^{n}+\nu \nabla^{2} \mathbf{u}^{n}
$$
Remember that $\Delta=\nabla^2$



The superscript indicates the time step, and we assume we know everything at $n$ timestep. In this step, we solve the intermediate velocity based on the real boundary conditions.



In the following step, we want to correct the obtained velocity to satisfy the zero divergence condition. We first need to solve for the pressure at $n+1$ time step by solving the Possion equation:
$$
\nabla^{2} \pi^{n+1}=\frac{1}{\Delta t} \nabla \cdot \mathbf{u}^{*}
$$


Remember we assume the density to be kept constant all the time.



After we obtain the pressure, the velocity at $n+1$ time step is calculated by the following equation:
$$
\mathbf{u}^{n+1}=\mathbf{u}^{*}-\frac{\Delta t}{\rho} \nabla p^{n+1}=\mathbf{u}^{*}-\Delta t
\cdot\nabla \pi^{n+1}
$$

## Steps 1 of the Projection Method
The starting equation for the $u$ velocity:
$$
{\bar{u^*}-\bar{u^n} \over \Delta t}+(\bar{u^n} \nabla) \bar{u^n} = \Delta \bar{u^*} 
$$
The equation rewritten in 2D cylindrical coordinates ($z$ and $r$ used):
$$
{\bar{u^*_r}-\bar{u^n_r} \over \Delta t} + \bar{u^n_r}{\partial\bar{u^n_r} \over \partial \bar{r}} +\bar{u^n_z}{\partial\bar{u^n_z} \over \partial \bar{z}} = {1 \over Re} \left[{1 \over \bar{r}}{\partial \over \partial \bar{r}} \left(\bar{r} {\partial\bar{u^*_r} \over \partial \bar{r}} \right)+{\partial^2\bar{u^*_r} \over \partial \bar{z^2}}-{\bar{u^*_r} \over \bar{r^2}}\right]
$$
Discretization of the Terms:









[1]:	https://www.overleaf.com/read/hzzczmvjnnht
