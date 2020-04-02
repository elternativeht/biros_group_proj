
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
The governing equations that will allow us to solve for the semi-steady state 2D velocity distribution in the cylindrical domain include the continuity equation and two momentum equations. Refering to the proposal, the non-dimensional continuity equation is

$$ \frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0 $$

and the two continuity equations for the $r$ and $z$ components of velocity, respectively, are

$$ \overline{u}_r\frac{\partial\overline{u}_r}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_r}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_r}{\partial \overline{r}}\right) - \frac{\overline{u}_r}{\overline{r}} + \frac{\partial^2 \overline{u}_r}{\partial \overline{z}^2}\right] $$

and

$$\overline{u}_r\frac{\partial\overline{u}_z}{\partial\overline{r}} + \overline{u}_z\frac{\partial\overline{u}_z}{\partial\overline{z}} = \frac{1}{Re}\left[\frac{1}{\overline{r}}\frac{\partial}{\partial \overline{r}}\left(\overline{r}\frac{\partial \overline{u}_z}{\partial \overline{r}}\right) + \frac{\partial^2 \overline{u}_z}{\partial \overline{z}^2}\right] + \frac{Ga}{Re^2} - \frac{1}{\rho}\frac{\partial p}{\partial \overline{z}}.$$

## Finite difference method equation discretization 

The original continuous equation is:

$$\frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0$$



Different discretization schemes can be used for both radial direction ($r$ direction) and axial direction ($z$ direction). For simplificity, the z direction is aligned with the gravity direction, i.e. the increase in z coordinate corresponds to the decrease in actual height.



The solid flows inward and downward, so a possible discretization upwind scheme would be:

$$\frac{1}{\overline{r}_N}\frac{(\overline{r}_{N+1}\:\overline{u}_{r,N+1}-\overline{r}_{N}\:\overline{u}_{r,N})}{\Delta \overline{r}} + \frac{\overline{u}_{z,N}-\overline{u}_{z,N-1}}{\Delta \overline{z}} = 0$$

The equation is based on the first order upwind discretization at the location $N$.





[1]:	https://www.overleaf.com/read/hzzczmvjnnht
