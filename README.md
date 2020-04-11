
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
The governing equations that will allow us to solve for the semi-steady state 2D velocity distribution in the cylindrical domain include the continuity equation and two momentum equations. Refering to the proposal, the non-dimensional continuity equation is

<p align="center"><img src="/tex/b20441acd018750e8d678889ae462554.svg?invert_in_darkmode&sanitize=true" align=middle width=145.5763914pt height=34.7253258pt/></p>

and the two continuity equations for the <img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> and <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> components of velocity, respectively, are

<p align="center"><img src="/tex/09c0d1dd7ca84874dafaa93357c35213.svg?invert_in_darkmode&sanitize=true" align=middle width=494.6479824pt height=40.11819404999999pt/></p>

and

<p align="center"><img src="/tex/66b2ce8d0178306244e2c799d6b1d858.svg?invert_in_darkmode&sanitize=true" align=middle width=513.4111356pt height=40.11819404999999pt/></p>

where the non-dimensional terms are defined as

<p align="center"><img src="/tex/5a3f6e46b91c5d8fe520b10c3c78ca87.svg?invert_in_darkmode&sanitize=true" align=middle width=708.6747392999999pt height=39.84127125pt/></p>

---

Updated on Apr 11 13:57 by Jinghu, to compare with the original dimensionless equations.



Original continuous equation at the cylindrical coordinates (here <img src="/tex/0a5a0c3d35e8e061b4b7c57b0c7c6add.svg?invert_in_darkmode&sanitize=true" align=middle width=16.02556724999999pt height=14.15524440000002pt/> and <img src="/tex/5d3d08286c873e3706833a05670e23d7.svg?invert_in_darkmode&sanitize=true" align=middle width=35.67354779999999pt height=24.65753399999998pt/> is cancelled out):
<p align="center"><img src="/tex/3a67de4e3043ba54fef558cb3b461de0.svg?invert_in_darkmode&sanitize=true" align=middle width=144.66317085pt height=34.7253258pt/></p>
Original momentum equations:
<p align="center"><img src="/tex/42bd284bf063cb92a04bb1fdcf8f4837.svg?invert_in_darkmode&sanitize=true" align=middle width=502.0140807pt height=59.1786591pt/></p>


The non-dimensionalization process takes the following conversions:
<p align="center"><img src="/tex/81b32feaf0fb5dc627be6b9e56f5029c.svg?invert_in_darkmode&sanitize=true" align=middle width=685.9271034pt height=37.693258349999994pt/></p>
The dimensionless continuous equation then is:
<p align="center"><img src="/tex/b059c5503a448993f80a667638aaf64d.svg?invert_in_darkmode&sanitize=true" align=middle width=145.5763914pt height=34.7253258pt/></p>
The dimensionless radial direction momentum equation is:
<p align="center"><img src="/tex/789eb34231ed01609a5f77aa34263bc0.svg?invert_in_darkmode&sanitize=true" align=middle width=466.07951610000003pt height=40.11819404999999pt/></p>


And dimensionless axial direction momentum equation is:
<p align="center"><img src="/tex/8a5672e369e2d17bd0726ef06747f8f9.svg?invert_in_darkmode&sanitize=true" align=middle width=485.60185409999997pt height=40.11819404999999pt/></p>

---



## Finite difference method equation discretization 

The original continuous equation is:

<p align="center"><img src="/tex/c83032e3696a0cdea9bff801e9295c41.svg?invert_in_darkmode&sanitize=true" align=middle width=145.5763914pt height=34.7253258pt/></p>



Different discretization schemes can be used for both radial direction (<img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> direction) and axial direction (<img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> direction). For simplificity, the z direction is aligned with the gravity direction, i.e. the increase in z coordinate corresponds to the decrease in actual height.



The solid flows inward and downward, so a possible discretization upwind scheme would be:

<p align="center"><img src="/tex/2fdee062a0651d5d58918344d18945db.svg?invert_in_darkmode&sanitize=true" align=middle width=359.63558894999994pt height=37.1910528pt/></p>

The equation is based on the first order upwind discretization at the location <img src="/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/>.





[1]:	https://www.overleaf.com/read/hzzczmvjnnht
