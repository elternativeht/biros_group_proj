
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
The governing equations that will allow us to solve for the semi-steady state 2D velocity distribution in the cylindrical domain include the continuity equation and two momentum equations. Refering to the proposal, the non-dimensional continuity equation is

<p align="center"><img src="/tex/b20441acd018750e8d678889ae462554.svg?invert_in_darkmode&sanitize=true" align=middle width=145.5763914pt height=34.7253258pt/></p>

and the two continuity equations for the <img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> and <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> components of velocity, respectively, are

<p align="center"><img src="/tex/f8daaed4e9eb2fcf3e5617c88d8bd851.svg?invert_in_darkmode&sanitize=true" align=middle width=377.71099574999994pt height=40.11819404999999pt/></p>

and

<p align="center"><img src="/tex/3647032c548ef0c54997ed14680b8d34.svg?invert_in_darkmode&sanitize=true" align=middle width=436.44889034999994pt height=40.11819404999999pt/></p>

where the non-dimensional terms are defined as
<p align="center"><img src="/tex/2311a25815cac661f12b38a3c700ee3d.svg?invert_in_darkmode&sanitize=true" align=middle width=522.77943795pt height=79.38279195pt/></p>\frac{1}{\overline{r}}\frac{\partial(\overline{r}\:\overline{u}_r)}{\partial \overline{r}} + \frac{\partial \overline{u}_z}{\partial \overline{z}} = 0<p align="center"><img src="/tex/cd7a4bcbe57a7012a6ab9afdece38133.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2746371499999pt height=74.70320054999999pt/></p>\frac{1}{\overline{r}_N}\frac{(\overline{r}_{N+1}\:\overline{u}_{r,N+1}-\overline{r}_{N}\:\overline{u}_{r,N})}{\Delta \overline{r}} + \frac{\overline{u}_{z,N}-\overline{u}_{z,N-1}}{\Delta \overline{z}} = 0$<img src="/tex/88d8d0b2ed7239e247ee5dfe0633ee79.svg?invert_in_darkmode&sanitize=true" align=middle width=544.3849834499999pt height=45.84475499999998pt/>N$.





[1]:	https://www.overleaf.com/read/hzzczmvjnnht
