
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
The governing equations that will allow us to solve for the semi-steady state 2D velocity distribution in the cylindrical domain include the continuity equation and two momentum equations. Refering to the proposal, the non-dimensional continuity equation is

<p align="center"><img src="/tex/b20441acd018750e8d678889ae462554.svg?invert_in_darkmode&sanitize=true" align=middle width=145.5763914pt height=34.7253258pt/></p>

and the two continuity equations for the <img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> and <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> components of velocity, respectively, are

<p align="center"><img src="/tex/d0eca5d23c58fb85589a9c040e39ee00.svg?invert_in_darkmode&sanitize=true" align=middle width=494.6479824pt height=40.11819404999999pt/></p>

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
<p align="center"><img src="/tex/a1da397c75f34ceb78f6f505752659ba.svg?invert_in_darkmode&sanitize=true" align=middle width=485.60185409999997pt height=40.11819404999999pt/></p>

---



## MAC Schemes 

Current MAC schemes, along with suggestions by George, was tried to be summarized here. (Updated at 14:02, Apr 12 by Jinghu)



The original continuous equation is:
<p align="center"><img src="/tex/7c80891cb1f06588c97bd489a481a22e.svg?invert_in_darkmode&sanitize=true" align=middle width=248.36120429999997pt height=34.7253258pt/></p>
Let's assume right now we have discretization methods <img src="/tex/8e423496dc713a5ecc5b76be73dead1d.svg?invert_in_darkmode&sanitize=true" align=middle width=13.652895299999988pt height=22.55708729999998pt/> for the continuous equation featuring unknown velocity properties <img src="/tex/27d932569045a2b5876d0a94af9d2b74.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=19.871860799999983pt/>

The original momentum equation can be simplified to the following terms:
<p align="center"><img src="/tex/e9dd68830097b7908ae8ff890a887045.svg?invert_in_darkmode&sanitize=true" align=middle width=219.26648205pt height=33.81208709999999pt/></p>
where <img src="/tex/43a06e2f32cf18ced46a8183ceacdb3e.svg?invert_in_darkmode&sanitize=true" align=middle width=61.51811984999999pt height=24.65753399999998pt/>

Let's say we discretize the continuous equation to have matrix operator <img src="/tex/8e423496dc713a5ecc5b76be73dead1d.svg?invert_in_darkmode&sanitize=true" align=middle width=13.652895299999988pt height=22.55708729999998pt/> onto the <img src="/tex/5f95d6e24f7d43cd41f1bfe4d4dcec62.svg?invert_in_darkmode&sanitize=true" align=middle width=35.27217044999999pt height=27.89013150000002pt/> (we have info on the <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> timestep):
<p align="center"><img src="/tex/918cd2429068351b33ca7f6df1d75697.svg?invert_in_darkmode&sanitize=true" align=middle width=92.66924535pt height=18.312383099999998pt/></p>
where <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> indicates the next time step.



Taking momentum equation into account, let's say we use <img src="/tex/b2f745d6c8328e3502fbe96121a457b3.svg?invert_in_darkmode&sanitize=true" align=middle width=14.49764249999999pt height=22.55708729999998pt/> as discretization matrix for diffusion term, <img src="/tex/e1616d38ca198927a525c18fad3716cd.svg?invert_in_darkmode&sanitize=true" align=middle width=14.29216634999999pt height=22.55708729999998pt/> as discertization matrix for advection term, <img src="/tex/69c52b5b30089ca77e145e6352879529.svg?invert_in_darkmode&sanitize=true" align=middle width=17.94511949999999pt height=22.55708729999998pt/> as discretization matrix for pressure term:
<p align="center"><img src="/tex/ffa280aa7a2ad69e296a5a687ad32e57.svg?invert_in_darkmode&sanitize=true" align=middle width=430.22172434999993pt height=36.34162455pt/></p>
Both continuous and momentum equations can be converted to a matrix type equations:
<p align="center"><img src="/tex/1b192869e13171685d78330bf885b063.svg?invert_in_darkmode&sanitize=true" align=middle width=169.85411355pt height=39.520717499999996pt/></p>
where <img src="/tex/1d86976b96ee0384870352131f454b7f.svg?invert_in_darkmode&sanitize=true" align=middle width=55.221135749999995pt height=22.55708729999998pt/>can be derived from the discretized momentum equations, and they are the function of the values at <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> time step.Solving this matrix equation could allow us to get the properties at <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> time step.



## Projection methods

Another possible method is called projectio method. The basic algorithm is shown here.



First the original momentum equation is the one to use:
<p align="center"><img src="/tex/e8d6b30ca704507088ae1701de6e61f1.svg?invert_in_darkmode&sanitize=true" align=middle width=253.83263069999995pt height=33.81208709999999pt/></p>
First the intermediate velocity is calculated:
<p align="center"><img src="/tex/4fe61313acb4b338bcca4c3093c653d2.svg?invert_in_darkmode&sanitize=true" align=middle width=245.97131459999997pt height=33.715788149999995pt/></p>
Remember that <img src="/tex/89b1bb306081ee6493af9931b0511665.svg?invert_in_darkmode&sanitize=true" align=middle width=55.86752159999999pt height=26.76175259999998pt/>



The superscript indicates the time step, and we assume we know everything at <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> timestep. In this step, we solve the intermediate velocity based on the real boundary conditions.



In the following step, we want to correct the obtained velocity to satisfy the zero divergence condition. We first need to solve for the pressure at <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> time step by solving the Possion equation:
<p align="center"><img src="/tex/7fd490426be54f30979bad7c126b3048.svg?invert_in_darkmode&sanitize=true" align=middle width=144.93071999999998pt height=32.990165999999995pt/></p>


Remember we assume the density to be kept constant all the time.



After we obtain the pressure, the velocity at <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> time step is calculated by the following equation:
<p align="center"><img src="/tex/5ac251767d819596a328ca01755db3a1.svg?invert_in_darkmode&sanitize=true" align=middle width=307.30683884999996pt height=36.82577085pt/></p>

## Steps 1 of the Projection Method
The starting equation for the <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> velocity:
<p align="center"><img src="/tex/68f971431ab134929b6e7ff03b54ec5b.svg?invert_in_darkmode&sanitize=true" align=middle width=193.26470625pt height=35.1292854pt/></p>

The equation rewritten in 2D cylindrical coordinates (<img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> and <img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> used):
<p align="center"><img src="/tex/d46068a6c2652de6a896644115e7254c.svg?invert_in_darkmode&sanitize=true" align=middle width=462.40435395000003pt height=40.11819404999999pt/></p>


Discretization of each term using a central differnce scheme:
<p align="center"><img src="/tex/c96e606a6fbaaf94a07d9c8d40d26c4c.svg?invert_in_darkmode&sanitize=true" align=middle width=421.3320111pt height=255.96313545pt/></p>







[1]:	https://www.overleaf.com/read/hzzczmvjnnht
