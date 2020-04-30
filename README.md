
# biros\_group\_proj
The group project for Biros

## The current project proposal

Current proposal is at the [link  (read only)][1]. The simulation process will be carried out in the rectangular region representing a 2-dimensional axis-symmetric region. Although finite volume method is the method in the proposal, possible finite difference method can be used if possible.

## Governing Equations
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

## Projection methods



Another possible method is called projectio method. The basic algorithm is shown here.



First the original momentum equation is the one to use:
<p align="center"><img src="/tex/e8d6b30ca704507088ae1701de6e61f1.svg?invert_in_darkmode&sanitize=true" align=middle width=253.83263069999995pt height=33.81208709999999pt/></p>
First the intermediate velocity is calculated:
<p align="center"><img src="/tex/632a6f2cd5073754ecb284d1d3802314.svg?invert_in_darkmode&sanitize=true" align=middle width=252.2778489pt height=40.1768301pt/></p>
Remember that <img src="/tex/89b1bb306081ee6493af9931b0511665.svg?invert_in_darkmode&sanitize=true" align=middle width=55.86752159999999pt height=26.76175259999998pt/>. Here we have <img src="/tex/63e77d3dfbf71950a9be1ae6519ad4dd.svg?invert_in_darkmode&sanitize=true" align=middle width=56.866378799999985pt height=24.65753399999998pt/>.



The superscript indicates the time step. We assume properties at timestep <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> were known. Step one method allowed us to solve the intermediate velocity based on the real boundary conditions. It is an implicit equation, and after linearization it would become a linear matrix equation.



The intermediate velocity profile obtained in step one does not satisfy the continuous equation. In step two, we want to correct the obtained velocity to satisfy the zero divergence condition. We first need to solve for the pressure at <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> time step by doing the correction:
<p align="center"><img src="/tex/63346aef0a1b6a686943dee91328556f.svg?invert_in_darkmode&sanitize=true" align=middle width=165.6108168pt height=15.572667pt/></p>


Taking the divergence of the equation above, the Possion equation can be obtained:
<p align="center"><img src="/tex/846d362d70a47606c766a457353badb0.svg?invert_in_darkmode&sanitize=true" align=middle width=134.3036376pt height=14.202794099999998pt/></p>


 solving the Possion equation:
<p align="center"><img src="/tex/7fd490426be54f30979bad7c126b3048.svg?invert_in_darkmode&sanitize=true" align=middle width=144.93071999999998pt height=32.990165999999995pt/></p>


Remember we assume the density to be kept constant all the time.



After we obtain the pressure, the velocity at <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> time step is calculated by the following equation:
<p align="center"><img src="/tex/8e36ce664b32c0754eeb58f9370c348c.svg?invert_in_darkmode&sanitize=true" align=middle width=176.09766735pt height=15.572667pt/></p>

## Step 1 of the Projection Method
The starting equation for the <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> velocity:
<p align="center"><img src="/tex/68f971431ab134929b6e7ff03b54ec5b.svg?invert_in_darkmode&sanitize=true" align=middle width=193.26470625pt height=35.1292854pt/></p>

The equation rewritten in 2D cylindrical coordinates (<img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> and <img src="/tex/89f2e0d2d24bcf44db73aab8fc03252c.svg?invert_in_darkmode&sanitize=true" align=middle width=7.87295519999999pt height=14.15524440000002pt/> used):
<p align="center"><img src="/tex/d46068a6c2652de6a896644115e7254c.svg?invert_in_darkmode&sanitize=true" align=middle width=462.40435395000003pt height=40.11819404999999pt/></p>



Discretization of each term using a central differnce scheme:
<p align="center"><img src="/tex/f77c67f9736a5e4ddd412484c37b0477.svg?invert_in_darkmode&sanitize=true" align=middle width=195.65030595pt height=39.2184573pt/></p>

<p align="center"><img src="/tex/3a16e88dc0a48e57cc27d709956be675.svg?invert_in_darkmode&sanitize=true" align=middle width=191.47775955pt height=39.61800315pt/></p>

<p align="center"><img src="/tex/71ef133b5aa04b7b5d1a16fbec1673a7.svg?invert_in_darkmode&sanitize=true" align=middle width=259.61135475pt height=43.72759214999999pt/></p>

<p align="center"><img src="/tex/e1ca1b5924d6177eff4293c8478975f5.svg?invert_in_darkmode&sanitize=true" align=middle width=195.65030595pt height=39.2184573pt/></p>

<p align="center"><img src="/tex/74ae5ee5383c7cdcb51cc74a5074fc4f.svg?invert_in_darkmode&sanitize=true" align=middle width=191.47775955pt height=39.61800315pt/></p>

<p align="center"><img src="/tex/7c5c5607a353042c445e4bfef8ab9faf.svg?invert_in_darkmode&sanitize=true" align=middle width=259.61135475pt height=43.72759214999999pt/></p>





## Grid Naming and Notation



The grid has following naming notations:

<img src="./References/Document Sources/Notation.png" alt="./References/Document Sources/Notation.png" style="zoom:80%;" />



where <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> is the grid number in axial direction, <img src="/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/> the grid number in radial direction. The notation is using 1-index rule as Matlab is also 1-index (starting the index from 1). Each row is representing an axial coordinate and each column representing a radial coordinate. The data is using [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). Note that we are using staggered grid and the properties are having different dimensions.



## Step 1 in staggered grids

<img src="./References/Document Sources/Staggered grid.jpg" style="zoom:50%;" />



The radial momentum equation discretization is shown below:




<p align="center"><img src="/tex/34b408443ae916fb3c7edb9d22b1cf21.svg?invert_in_darkmode&sanitize=true" align=middle width=473.64803309999996pt height=34.7253258pt/></p>


The advection term can be discretized as:
<p align="center"><img src="/tex/ce70d322c09720fd3035f3d9292b0cdb.svg?invert_in_darkmode&sanitize=true" align=middle width=1204.8766999499999pt height=39.452455349999994pt/></p>


The diffusion term is:


<p align="center"><img src="/tex/502b1633fd4c52b15683527ce143f807.svg?invert_in_darkmode&sanitize=true" align=middle width=410.76756765pt height=38.83491479999999pt/></p>

<p align="center"><img src="/tex/447af874d0ad3c97d2a8b9f10578e240.svg?invert_in_darkmode&sanitize=true" align=middle width=361.98238065pt height=34.7253258pt/></p>

<p align="center"><img src="/tex/535deca738544baf1ce4fc37f423ce51.svg?invert_in_darkmode&sanitize=true" align=middle width=432.7234593pt height=38.83491479999999pt/></p>



The intermediate variable can be summarized by the following linearized equation:


<p align="center"><img src="/tex/63abbe84a18e005227f3fccb42117868.svg?invert_in_darkmode&sanitize=true" align=middle width=640.3252482pt height=16.438356pt/></p>


where
<p align="center"><img src="/tex/e5c8f03e052873e81f6ba801fe328b42.svg?invert_in_darkmode&sanitize=true" align=middle width=207.60861824999998pt height=37.099754999999995pt/></p>

<p align="center"><img src="/tex/d8785d54a99155cda937e99964e1503e.svg?invert_in_darkmode&sanitize=true" align=middle width=210.37781489999998pt height=39.452455349999994pt/></p>



<p align="center"><img src="/tex/c3c9a3cb9ed4ec26f522352975dfc35c.svg?invert_in_darkmode&sanitize=true" align=middle width=210.00904814999998pt height=39.452455349999994pt/></p>



<p align="center"><img src="/tex/96f966d682548d731a7888aee4804105.svg?invert_in_darkmode&sanitize=true" align=middle width=228.14831654999998pt height=37.099754999999995pt/></p>



<p align="center"><img src="/tex/e5ec928a2d98e737cb178dd797073cac.svg?invert_in_darkmode&sanitize=true" align=middle width=423.5558959499999pt height=39.452455349999994pt/></p>

<p align="center"><img src="/tex/d0aa239cc79909a5b241fb2e7dc1f3e2.svg?invert_in_darkmode&sanitize=true" align=middle width=725.4352346999999pt height=39.452455349999994pt/></p>







where <img src="/tex/5bf3cf1ae2a463f03e1cffe5d95f8514.svg?invert_in_darkmode&sanitize=true" align=middle width=137.11557749999997pt height=24.65753399999998pt/>, where <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> is the z-direction coordinate and <img src="/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/> r-direction coordinate. The equation only applies to internal points, with z positive direction pointing downward and r positive direction pointint rightward (outward).







## Boundary conditions



<img src="./References/Document Sources/Boundary condition.jpg" style="zoom:24%;" />



Boundary (1): set velocity (not pressure)



Boundary (2): you need to set <img src="/tex/4aee50c7fe12fe3e346db76f22ffd69e.svg?invert_in_darkmode&sanitize=true" align=middle width=45.384226799999986pt height=21.18721440000001pt/> but also <img src="/tex/f6cdc5360d2d944d54ab5675731bb521.svg?invert_in_darkmode&sanitize=true" align=middle width=77.00521784999998pt height=21.18721440000001pt/>. The stress includes velocity and pressure terms. 



Boundary (3):  normal stress =0. The stress has a viscous component and the pressure. 



Boundaries 4 and 5:  no slip.










[1]:	https://www.overleaf.com/read/hzzczmvjnnht
