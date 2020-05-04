
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
\frac{\partial\overline{u}_r}{\partial \overline{t}} \left[{u^2_\infty \over H} \right]+ \overline{u}_r\frac{\partial\overline{u}_r}{\partial\overline{r}}\left[{u^2_\infty \over H} \right] + \overline{u}_z\frac{\partial\overline{u}_r}{\partial\overline{r}}\left[{u^2_\infty \over H} \right] = -\frac{\partial \overline{p}}{\partial \overline{r}}\left[{u^2_\infty \over H} \right] +\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_r}{\partial \overline{r}}\right)\left[{\nu u_\infty \over H^2} \right] - \frac{\overline{u}_r}{\overline{r}^2}\left[{\nu u_\infty \over H^2} \right] + \frac{\partial^2 \overline{u}_r}{\partial \overline{z}^2} \left[{\nu u_\infty \over H^2} \right] 
$$

$$
\frac{\partial\overline{u}_r}{\partial \overline{t}} + \overline{u}_r\frac{\partial\overline{u}_r}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_r}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_r}{\partial \overline{r}}\right) - \frac{\overline{u}_r}{\overline{r}^2} + \frac{\partial^2 \overline{u}_r}{\partial \overline{z}^2}\right] - \frac{\partial \overline{p}}{\partial \overline{r}}
$$


And dimensionless axial direction momentum equation is:
$$
\frac{\partial\overline{u}_z}{\partial \overline{t}}\left[{u^2_\infty \over H} \right] + \overline{u}_r\frac{\partial\overline{u}_z}{\partial\overline{r}}\left[{u^2_\infty \over H} \right] + \overline{u}_z\frac{\partial\overline{u}_z}{\partial\overline{z}}\left[{u^2_\infty \over H} \right] = - \frac{\partial \overline{p}}{\partial \overline{z}}\left[{u^2_\infty \over H} \right] +\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_z}{\partial \overline{r}}\right)\left[{\nu u_\infty \over H^2} \right]  + \frac{\partial^2 \overline{u}_z}{\partial \overline{z}^2}\left[{\nu u_\infty \over H^2} \right] -g, 
$$

$$
\frac{\partial\overline{u}_z}{\partial \overline{t}} + \overline{u}_r\frac{\partial\overline{u}_z}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_z}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_z}{\partial \overline{r}}\right) + \frac{\partial^2 \overline{u}_z}{\partial \overline{z}^2}\right] + \frac{1}{Fr^2} - \frac{\partial \overline{p}}{\partial \overline{z}},
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
{\overline{u^*_r}-\overline{u^n_r} \over \Delta t} + \overline{u^n_r}{\partial\overline{u^n_r} \over \partial \overline{r}} +\overline{u^n_r}{\partial\overline{u^n_z} \over \partial \overline{z}} = {1 \over Re} \left[{1 \over \overline{r}}{\partial \over \partial \overline{r}} \left(\overline{r} {\partial\overline{u^*_r} \over \partial \overline{r}} \right)+{\partial^2\overline{u^*_r} \over \partial \overline{z^2}}-{\overline{u^*_r} \over \overline{r^2}}\right]
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

Howver, due to a staggered grid being used the discretizations require an extra step of complexity, detailed below. 



## Grid Naming and Notation



The grid has following naming notations:

<img src="./References/Document Sources/Notation.png" alt="./References/Document Sources/Notation.png" style="zoom:80%;" />



where $M$ is the grid number in axial direction, $N$ the grid number in radial direction. The notation is using 1-index rule as Matlab is also 1-index (starting the index from 1). Each row is representing an axial coordinate and each column representing a radial coordinate. The data is using [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). Note that we are using staggered grid and the properties are having different dimensions.



## Projection method step 1 implementation in staggered grids

<img src="./References/Document Sources/Staggered grid.jpg" style="zoom:50%;" />

### Radial direction momentum equation

The radial momentum equation discretization is shown below:




$$
\frac{\boldsymbol{u}^{*}_{(i,j)}-\boldsymbol{u}^{n}_{(i,j)}}{\Delta t}+\left(\boldsymbol{u}^{n}_{(i,j)}\cdot\nabla\right)\boldsymbol{u}^{n}_{(i,j)}=\Delta \boldsymbol{u}^{*}_{(i,j)}
$$

$$
{\partial u_r \over \partial r} \approx {{u_{r,i(j+1)} + u_{r,ij} \over 2}-{u_{r,ij} +u_{r,i(j-1)} \over 2} \over \Delta r} \approx {u_{r,i(j+1)} - u_{r,i(j-1)} \over 2 \Delta r}
$$

$$
{\partial^2 u_r \over \partial r^2} \approx {{u_{r,i(j+1)} + u_{r,ij} \over \Delta r}-{u_{r,ij} +u_{r,i(j-1)} \over \Delta r} \over \Delta r} \approx {u_{r,i(j+1)} - 2u_{r,ij} + u_{r,i(j-1)} \over (\Delta r)^2}
$$

$$
{\partial u_r \over \partial z} \approx {{u_{r,(i+1)j} + u_{r,ij} \over 2}-{u_{r,ij} +u_{r,(i-1)j} \over 2} \over \Delta z} \approx {u_{r,(i+1)j} - u_{r,(i-1)j} \over 2 \Delta z}
$$

$$
{\partial^2 u_r \over \partial z^2} \approx {{u_{r,(i+1)j} + u_{r,ij} \over \Delta z}-{u_{r,ij} +u_{r,(i-1)j} \over \Delta z} \over \Delta z} \approx {u_{r,(i+1)j} - 2u_{r,ij} + u_{r,(i-1)j} \over (\Delta z)^2}
$$


$$
u_z \approx {u^k_{z,ij}+u^k_{z,i(j-1)}+u^k_{z,(i+1)(j-1)}+u^k_{z,(i+1)j} \over 4}
$$


The advection term can be discretized as:
$$
\left( \boldsymbol{u}^n_{r,ij} \cdot \nabla \right) \boldsymbol{u}^n_{r,ij}=u^n_{r,ij}\frac{u^n_{r,i(j+1)}-u^n_{r,i(j-1)}}{2\Delta r}+\frac{u^n_{z,i(j-1)}+u^n_{z,ij}+u^n_{z,(i+1)(j-1)}+u^n_{z,(i+1)j}}{4}\left[\frac{u^{n}_{r,(i+1)j}-u^{n}_{r,(i-1)j}}{2\Delta z}\right]
$$



The diffusion term is:

$$
\frac{1}{\text{Re}}\left[\frac{1}{r_{ij}}\frac{u^*_{r,i(j+1)}-u^*_{r,i(j-1)}}{2\Delta r}+\frac{u^*_{r,i(j+1)}-2u^*_{r,ij}+u^*_{r,i(j-1)}}{(\Delta r)^2}-\frac{u^*_{r,ij}}{r^2_{ij}}+\frac{u^*_{r,(i+1)j}-2u^*_{r,ij}+u^*_{r,(i-1)j}}{(\Delta z)^2}\right]
$$



The intermediate variable can be summarized by the following linearized equation:


$$
Au^*_{r,ij}+Bu^*_{r,i(j+1)}+Cu^*_{r,i(j-1)}+Du^*_{r,(i+1)j}+Eu^*_{r,(i-1)j}=F
$$


where
$$
A = \frac{1}{\Delta t}+\frac{1}{\text{Re}}\frac{2}{(\Delta r)^2}+\frac{1 }{\text{Re}}\frac{2}{r_{ij}^2}+\frac{1}{\text{Re}}\frac{2}{(\Delta z)^2}
$$




$$
B = -\frac{1}{\text{Re}}\frac{1}{r_{ij}}\frac{1}{2\Delta r}-\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}
$$




$$
C = +\frac{1}{\text{Re}}\frac{1}{r_{ij}}\frac{1}{2\Delta r}-\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}
$$




$$
D =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2},E =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2}
$$




$$
F= \frac{u^n_{r,ij}}{\Delta t}-\left(u^n_{r,ij}\frac{u^n_{r,i(j+1)}-u^n_{r,i(j-1)}}{2\Delta r}+\frac{u^n_{z,i(j-1)}+u^n_{z,ij}+u^n_{z,(i+1)(j-1)}+u^n_{z,(i+1)j}}{4}\left[\frac{u^{n}_{r,(i+1)j}-u^{n}_{r,(i-1)j}}{2\Delta z}\right]\right)
$$


$$
A_r^* u_r^* = F(u^n)
$$


$$
A_r^*= 
$$

  \begin{array}{ccccc}
 -	& -	&-	&- & \Omega^n_1 \\
 -	& -	&-	&\Omega^n_2 & - \\
 -	& -	&\Omega^n_3	&- & - \\
 -	& \Omega^n_4	&-	&- & - \\
 \Omega^n_5	& -	&-	&- & - \\
  \end{array}
  




where $\text{where }r_{ij}=(j-1)\cdot \Delta r$, where $i$ is the z-direction coordinate and $j$ r-direction coordinate. The equation only applies to internal points, with z positive direction pointing downward and r positive direction pointint rightward (outward).



### Axial direction momentum equation



Like radial, axial momentum equation discretization is shown below:




$$
\frac{\boldsymbol{u}^{*}_{(i,j)}-\boldsymbol{u}^{n}_{(i,j)}}{\Delta t}+\left(\boldsymbol{u}^{n}_{(i,j)}\cdot\nabla\right)\boldsymbol{u}^{n}_{(i,j)}=\Delta \boldsymbol{u}^{*}_{(i,j)}
$$


$$
{\partial u_z \over \partial r} \approx {{u_{z,i(j+1)} + u_{z,ij} \over 2}-{u_{z,i(j+1)} +u_{z,ij} \over 2} \over \Delta r} \approx {u_{z,i(j+1)} - u_{z,i(j-1)} \over 2 \Delta r}
$$

$$
{\partial^2 u_z \over \partial r^2} \approx {{{u_{z,i(j+1)} + u_{z,ij} \over 2} - u_{z,ij} \over {1 \over 2 }\Delta r}-{u_{r,ij} +u_{r,i(j-1)} \over \Delta r} \over \Delta r} \approx {u_{r,i(j+1)} - 2u_{r,ij} + u_{r,i(j-1)} \over (\Delta r)^2}
$$

$$
{\partial u_r \over \partial z} \approx {{u_{r,(i+1)j} + u_{r,ij} \over 2}-{u_{r,ij} +u_{r,(i-1)j} \over 2} \over \Delta z} \approx {u_{r,(i+1)j} - u_{r,(i-1)j} \over 2 \Delta z}
$$

$$
{\partial^2 u_r \over \partial z^2} \approx {{u_{r,(i+1)j} + u_{r,ij} \over \Delta z}-{u_{r,ij} +u_{r,(i-1)j} \over \Delta z} \over \Delta z} \approx {u_{r,(i+1)j} - 2u_{r,ij} + u_{r,(i-1)j} \over (\Delta z)^2}
$$


$$
u_z \approx {u^k_{z,ij}+u^k_{z,i(j-1)}+u^k_{z,(i+1)(j-1)}+u^k_{z,(i+1)j} \over 4}
$$





The advection term can be discretized as:
$$
\left( \boldsymbol{u}^n_{z,ij} \cdot \nabla \right) \boldsymbol{u}^n_{z,ij}=\frac{u^n_{r,i(j+1)}+u^n_{r,ij}+u^n_{r,(i-1)(j+1)}+u^n_{r,(i-1)j}}{4}\left[\frac{u^{n}_{z,i(j+1)}-u^{n}_{z,i(j-1)}}{2\Delta r}\right]+u^n_{z,ij}\frac{u^n_{z,(i+1)j}-u^n_{z,(i-1)j}}{2\Delta z}
$$



The diffusion term is:

$$
\frac{1}{\text{Re}}\left[\frac{1}{r_{ij}+\frac{\Delta r}{2}}\frac{u^*_{z,i(j+1)}-u^*_{z,i(j-1)}}{2\Delta r}+\frac{u^*_{z,i(j+1)}-2u^*_{z,ij}+u^*_{z,i(j-1)}}{(\Delta r)^2}+\frac{u^*_{z,(i+1)j}-2u^*_{z,ij}+u^*_{z,(i-1)j}}{(\Delta z)^2}\right]
$$



The intermediate variable can be summarized by the following linearized equation:


$$
Au^*_{z,ij}+Bu^*_{z,i(j+1)}+Cu^*_{z,i(j-1)}+Du^*_{z,(i+1)j}+Eu^*_{z,(i-1)j}=F+\frac{1}{\text{Fr}^2}
$$


where in state matrix form we have:
$$
A = \frac{1}{\Delta t}+\frac{1}{\text{Re}}\frac{2}{(\Delta r)^2}+\frac{1}{\text{Re}}\frac{2}{(\Delta z)^2}
$$




$$
B = -\frac{1}{\text{Re}}\frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}
$$




$$
C = +\frac{1}{\text{Re}}\frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}
$$




$$
D =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2},E =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2}
$$




$$
F= \frac{u^n_{z,ij}}{\Delta t}-\left(\frac{u^n_{r,i(j+1)}+u^n_{r,ij}+u^n_{r,(i-1)(j+1)}+u^n_{r,(i-1)j}}{4}\left[\frac{u^{n}_{z,i(j+1)}-u^{n}_{z,i(j-1)}}{2\Delta r}\right]+u^n_{z,ij}\frac{u^n_{z,(i+1)j}-u^n_{z,(i-1)j}}{2\Delta z}\right)
$$





where $\text{where }r_{ij}=(j-1)\cdot \Delta r$, where $i$ is the z-direction coordinate and $j$ r-direction coordinate. The equation only applies to internal points, with z positive direction pointing downward and r positive direction pointint rightward (outward).



### Step 1 boundary conditions



<img src="./References/Document Sources/Boundary condition.jpg" style="zoom:24%;" />



- Boundary (1): velocity boundary; assume that the inflow velocity is equal to $U_\inf$.


$$
\overline{u_{z,1j}}=1,j=1,2,3...N
$$


- Boundary (2):  $v_r=0$ and $\text{stress}_r =0$. The stress includes velocity and pressure terms. However, as pressure term does not appear in the step one, the boundary condition is simplified to the following:


$$
\left\{\begin{array}{c}
\overline{u_{r,i1}}=0,i=1,2,3...M \\
\frac{\partial \overline{u_z}}{\partial \overline{r}}=0 \rightarrow \overline{u_{z,i2}}-\overline{u_{z,i1}}=0,i=1,2,3...M
\end{array}\right.
$$




- Boundary (3):  normal stress =0. The stress has a viscous component and the pressure. 


$$
\left\{\begin{array}{c}
\frac{\partial \overline{u_r}}{\partial \overline{z}}=0 \rightarrow \overline{u_{r,Mj}}-\overline{u_{r,(M-1)j}}=0,j=1,2,3...N+1\\
\frac{\partial \overline{u_z}}{\partial \overline{z}}=0 \rightarrow \overline{u_{z,(M+1)j}}-\overline{u_{z,Mj}}=0,j=1,2,3...N
\end{array}\right.
$$




- Boundaries 4 and 5:  no slip.


$$
\left\{\begin{array}{c}
\overline{u_{r,i(N+1)}}=0,i=1,2,3...M\\
\overline{u_{z,(M+1)j}}=0,j=S,S+1,S+2...N, \text{ where }S=L/\Delta r+1\text{, and }L\text{ is the outlet length}
\end{array}\right.
$$



## Projection method step 2 implementation in staggered grids



The Possion equation needs to be solved:
$$
\nabla^{2} \pi^{n+1}=\frac{1}{\Delta t} \nabla \cdot \mathbf{u}^{*}
$$
After the first step.





A divergence operator $\nabla \cdot \mathbf{u}^{*}$ has the following expression in 2d cylindrical coordinates:


$$
\nabla \cdot \mathbf{u}^{*} = \frac{1}{r} \frac{\partial\left(r u_{r}\right)}{\partial r}+\frac{\partial u_{z}}{\partial z}=\frac{1}{r}u_r+\frac{\partial u_r}{\partial r}+\frac{\partial u_{z}}{\partial z}
$$


The discretization would need to be also at the node point of  pressure $\pi = p/\rho$:
$$
\nabla \cdot \mathbf{u}^{*} \approx \frac{1}{r_{ij}+\Delta r/2}\frac{u_{r,ij}+u_{r,i(j+1)}}{2}+\frac{u_{r,i(j+1)}-u_{r,ij}}{\Delta r}+\frac{u_{z,(i+1)j}-u_{z,ij}}{\Delta z}
$$


The Laplace operator has the following equation in the 2d cylindrical coordinates:

 
$$
\nabla^2\pi=\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial \pi}{\partial r}\right)+\frac{\partial^{2} \pi}{\partial z^{2}}=\frac{1}{r}\frac{\partial \pi}{\partial r}+\frac{\partial^2\pi}{\partial r^2}+\frac{\partial^2 \pi}{\partial z^2}
$$


The discretization method would be:
$$
\frac{1}{r_{ij}+\Delta r/2}\frac{\pi_{i(j+1)}-\pi_{i(j-1)}}{2\Delta r}+\frac{\pi_{i(j+1)}-2\pi_{ij}+\pi_{i(j-1)}}{(\Delta r)^2}+\frac{\pi_{(i+1)j}-2\pi_{ij}+\pi_{(i-1)j}}{(\Delta z)^2}
$$
where, as in step 1, $r_{ij}=(j-1)\cdot \Delta r$.



The discretization allows us to get, for all internal nodes:
$$
A\pi_{ij}+B\pi_{i(j+1)}+C\pi_{i(j-1)}+B\pi_{(i+1)j}+B\pi_{(i-1)j}=F
$$


where
$$
A = \frac{2}{(\Delta r)^2}+\frac{2}{(\Delta z)^2}
$$




$$
B = -\frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{(\Delta r)^2}
$$




$$
C = \frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{(\Delta r)^2}
$$




$$
D = - \frac{1}{(\Delta z)^2},E =  -\frac{1}{(\Delta z)^2}
$$




$$
F=\frac{1}{\Delta t} \left(\frac{1}{r_{ij}+\Delta r/2}\frac{u_{r,ij}+u_{r,i(j+1)}}{2}+\frac{u_{r,i(j+1)}-u_{r,ij}}{\Delta r}+\frac{u_{z,(i+1)j}-u_{z,ij}}{\Delta z}\right)
$$



Df 

## Validation of the model
The results of our numerical scheme are compard to the results of K. F. Zhang and  J. Y. OOi (1998). 
### Establishing a relationship between viscosity and properties of solid grain particles
K. F. Zhang and  J. Y. OOi (1998) investigated the kinematic constant, B, as a function of various particle parameters. Has been shown to most closely related to particle size. The kinematic eqn they use is identical to the 1D unsteady heat conduction equation. 

[1]:	https://www.overleaf.com/read/hzzczmvjnnht
