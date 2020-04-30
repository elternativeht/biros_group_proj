
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
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

## Projection methods



Another possible method is called projectio method. The basic algorithm is shown here.



First the original momentum equation is the one to use:
$$
\frac{\partial \bold{u}}{\partial t}+(\bold{u}\cdot\nabla)(\bold{u})=-\nabla \bold{\pi}+\nu\Delta\bold{u}+\bold{f}
$$
First the intermediate velocity is calculated:
$$
\left\{\begin{array}{c}
\frac{\boldsymbol{u}^{*}-\boldsymbol{u}^{n}}{\Delta t}+\left(\boldsymbol{u}^{n} \cdot \nabla\right) \boldsymbol{u}^{n}=\Delta \boldsymbol{u}^{*}+\bold{\it{f}} \\
\text {Various B.C.       on } \partial \Omega
\end{array}\right.
$$
Remember that $\Delta=\nabla^2$. Here we have $\pi=p/\rho$.



The superscript indicates the time step. We assume properties at timestep $n$ were known. Step one method allowed us to solve the intermediate velocity based on the real boundary conditions. It is an implicit equation, and after linearization it would become a linear matrix equation.



The intermediate velocity profile obtained in step one does not satisfy the continuous equation. In step two, we want to correct the obtained velocity to satisfy the zero divergence condition. We first need to solve for the pressure at $n+1$ time step by doing the correction:
$$
\boldsymbol{u}^{*}=\boldsymbol{u}^{n+1}+\Delta t \nabla \pi^{n+1}
$$


Taking the divergence of the equation above, the Possion equation can be obtained:
$$
\nabla\cdot\boldsymbol{u}^{*} = \Delta t \Delta \pi^{n+1}
$$


 solving the Possion equation:
$$
\nabla^{2} \pi^{n+1}=\frac{1}{\Delta t} \nabla \cdot \mathbf{u}^{*}
$$


Remember we assume the density to be kept constant all the time.



After we obtain the pressure, the velocity at $n+1$ time step is calculated by the following equation:
$$
\mathbf{u}^{n+1}=\mathbf{u}^{*}-\Delta t
\cdot\nabla \pi^{n+1}
$$

## Step 1 of the Projection Method
The starting equation for the $u$ velocity:
$$
{\overline{u^*}-\overline{u^n} \over \Delta t}+(\overline{u^n} \nabla) \overline{u^n} = \Delta \overline{u^*} 
$$

The equation rewritten in 2D cylindrical coordinates ($z$ and $r$ used):
$$
{\overline{u^*_r}-\overline{u^n_r} \over \Delta t} + \overline{u^n_r}{\partial\overline{u^n_r} \over \partial \overline{r}} +\overline{u^n_z}{\partial\overline{u^n_z} \over \partial \overline{z}} = {1 \over Re} \left[{1 \over \overline{r}}{\partial \over \partial \overline{r}} \left(\overline{r} {\partial\overline{u^*_r} \over \partial \overline{r}} \right)+{\partial^2\overline{u^*_r} \over \partial \overline{z^2}}-{\overline{u^*_r} \over \overline{r^2}}\right]
$$



Discretization of each term using a central differnce scheme:
$$
{\partial\overline{u^n_r} \over \partial \overline{r}} \approx {\overline{u_r^n}_{,i(j+1)}-\overline{u_r^n}_{,i(j-1)} \over 2 \overline{\Delta r}}
$$

$$
{\partial\overline{u^*_r} \over \partial \overline{r}} \approx {\overline{u_r^*}_{,i(j+1)}-\overline{u_r^*}_{,i(j-1)} \over 2 \overline{\Delta r}}
$$

$$
{\partial^2\overline{u^*_r} \over \partial \overline{r^2}} \approx {\overline{u_r^*}_{,i(j+1)}-2 \overline{u_r^*}_{,ij}+\overline{u_r^*}_{,i(j-1)} \over (\overline{\Delta r})^2} 
$$

$$
{\partial\overline{u^n_r} \over \partial \overline{z}} \approx {\overline{u_r^n}_{,(i+1)j}-\overline{u_r^n}_{,(i-1)j} \over 2 \overline{\Delta z}}
$$

$$
{\partial\overline{u^*_r} \over \partial \overline{z}} \approx {\overline{u_r^*}_{,(i+1)j}-\overline{u_r^*}_{,(i-1)j} \over 2 \overline{\Delta z}} 
$$

$$
{\partial^2\overline{u^*_r} \over \partial \overline{z^2}} \approx {\overline{u_r^*}_{,(i+1)j}-2 \overline{u_r^*}_{,ij}+\overline{u_r^*}_{,(i-1)j} \over (\overline{\Delta z})^2}
$$





## Grid Naming and Notation



The grid has following naming notations:

<img src="./References/Document Sources/Notation.png" alt="./References/Document Sources/Notation.png" style="zoom:80%;" />



where $M$ is the grid number in axial direction, $N$ the grid number in radial direction. The notation is using 1-index rule as Matlab is also 1-index (starting the index from 1). Each row is representing an axial coordinate and each column representing a radial coordinate. The data is using [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). Note that we are using staggered grid and the properties are having different dimensions.



## Step 1 in staggered grids

<img src="./References/Document Sources/Staggered grid.jpg" style="zoom:50%;" />



The radial momentum equation discretization is shown below:




$$
\frac{\boldsymbol{u}^{*}(M,N)-\boldsymbol{u}^{n}(M,N)}{\Delta t}+(\boldsymbol{u}^{n}(M,N)\cdot\nabla)\boldsymbol{u}^{n}(M,N)=\Delta \boldsymbol{u}^{*}(M,N)
$$


The advection term can be discretized as:
$$
(\boldsymbol{u}^{n}(M,N)\cdot\nabla)\boldsymbol{u}^{n}(M,N)=u^n_r(M,N)\frac{u^n_r(M,N+1)-u^n_r(M,N-1)}{2\Delta r}\\+\frac{u^n_z(M,N-1)+u^n_z(M,N)+u^n_z(M+1,N-1)+u^n_z(M+1,N)}{4}\left[\frac{u^n_r(M+1,N)-u^n_r(M-1,N)}{2\Delta z}\right]
$$


The diffusion term is:


$$
\frac{1}{\text{Re}}\frac{1}{r(M,N)}[(r(M,N)+\frac{\Delta r}{2})(\frac{u^*_r(M,N+1)-u^*_r(M,N)}{\Delta r})
$$

$$
-(r(M,N)-\frac{\Delta r}{2})(\frac{u^*_r(M,N)-u^*_r(M,N-1)}{\Delta r})]\frac{1}{\Delta r}
$$

$$
-\frac{1}{\text{Re}}\frac{1}{r^2(M,N)}u^*_r(M,N)+\frac{1}{\text{Re}}\frac{u^*_r(M+1,N)-u^*_r(M-1,n)}{(\Delta z)^2}
$$



The intermediate variable can be summarized by the following linearized equation:


$$
Au^*_r(M,N)+Bu^*_r(M,N+1)+Cu^*_r(M,N-1)+Du^*_r(M+1,N)+Eu^*_r(M-1,N)=-F
$$


where
$$
A = \frac{1}{\Delta t}+\frac{1}{\text{Re}}\frac{2}{(\Delta r)^2}+\frac{1 }{\text{Re}}\frac{2}{r^2}
$$

$$
B = -\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}\left(r+\frac{\Delta r}{2}\right)\frac{1}{r}
$$



$$
C = -\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}\left(r-\frac{\Delta r}{2}\right)\frac{1}{r}
$$



$$
D =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2},E =  \frac{1}{\text{Re}}\frac{1}{(\Delta z)^2}
$$



$$
F= \frac{u^n_r(M,N)}{\Delta t}-u_r^n\left(\frac{u^n_r(M,N+1)-u^n-r(M,N-1)}{2\Delta r}\right)
$$

$$
+ \frac{[u^n_z(M,N-1)+u^n_z(M,N)+u^n_z(M+1,N-1)+u^n_z(M+1,N)]}{4}\left[\frac{u^n_r(M+1,N)-u^n_r(M-1,N)}{2\Delta z}\right]
$$







where $\text{where }r=r(M,N)$, where $M$ is the z-direction coordinate and $N$ r-direction coordinate. The equation only applies to internal points, with z positive direction pointing downward and r positive direction pointint rightward (outward).







## Boundary conditions



<img src="./References/Document Sources/Boundary condition.jpg" style="zoom:24%;" />



Boundary (1): set velocity (not pressure)



Boundary (2): you need to set $v_r=0$ but also $\text{stress}_r =0$. The stress includes velocity and pressure terms. 



Boundary (3):  normal stress =0. The stress has a viscous component and the pressure. 



Boundaries 4 and 5:  no slip.










[1]:	https://www.overleaf.com/read/hzzczmvjnnht
