
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Finite difference method equation discretization 

The original continuous equation is:

$$\frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0$$



Different discretization schemes can be used for both radial direction ($r$ direction) and axial direction ($z$ direction). For simplificity, the z direction is aligned with the gravity direction, i.e. the increase in z coordinate corresponds to the decrease in actual height.



The solid flows inward and downward, so a possible discretization upwind scheme would be:

$$\frac{1}{\overline{r}_N}\frac{(\overline{r}_{N+1}\:\overline{u}_{r,N+1}-\overline{r}_{N}\:\overline{u}_{r,N})}{\Delta \overline{r}} + \frac{\overline{u}_{z,N}-\overline{u}_{z,N-1}}{\Delta \overline{z}} = 0$$

The equation is based on the first order upwind discretization at the location $N$.





[1]:	https://www.overleaf.com/read/hzzczmvjnnht