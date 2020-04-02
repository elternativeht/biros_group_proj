
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Finite difference method equation discretization 

The original continuous equation is:

<p align="center"><img src="/tex/c83032e3696a0cdea9bff801e9295c41.svg?invert_in_darkmode&sanitize=true" align=middle width=145.5763914pt height=34.7253258pt/></p>



Different discretization schemes can be used for both radial direction (<img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> direction) and axial direction (<img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> direction). For simplificity, the z direction is aligned with the gravity direction, i.e. the increase in z coordinate corresponds to the decrease in actual height.



The solid flows inward and downward, so a possible discretization upwind scheme would be:

<p align="center"><img src="/tex/2fdee062a0651d5d58918344d18945db.svg?invert_in_darkmode&sanitize=true" align=middle width=359.63558894999994pt height=37.1910528pt/></p>

The equation is based on the first order upwind discretization at the location <img src="/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/>.





[1]:	https://www.overleaf.com/read/hzzczmvjnnht