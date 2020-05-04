
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
<p align="center"><img src="/tex/972e7eb1f96482a6f3f60aaab979c2d4.svg?invert_in_darkmode&sanitize=true" align=middle width=753.6464232pt height=40.11819404999999pt/></p>

<p align="center"><img src="/tex/789eb34231ed01609a5f77aa34263bc0.svg?invert_in_darkmode&sanitize=true" align=middle width=466.07951610000003pt height=40.11819404999999pt/></p>


And dimensionless axial direction momentum equation is:
<p align="center"><img src="/tex/c74d1dbecc9bdbb3fd85f6914a272a43.svg?invert_in_darkmode&sanitize=true" align=middle width=698.30142525pt height=40.11819404999999pt/></p>

<p align="center"><img src="/tex/8a5672e369e2d17bd0726ef06747f8f9.svg?invert_in_darkmode&sanitize=true" align=middle width=485.60185409999997pt height=40.11819404999999pt/></p>

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
<p align="center"><img src="/tex/b40f94183b2cb779583fb91ad95e349a.svg?invert_in_darkmode&sanitize=true" align=middle width=462.40435395000003pt height=40.11819404999999pt/></p>



Discretization of each term using a central differnce scheme:
<p align="center"><img src="/tex/f77c67f9736a5e4ddd412484c37b0477.svg?invert_in_darkmode&sanitize=true" align=middle width=195.65030595pt height=39.2184573pt/></p>

<p align="center"><img src="/tex/3a16e88dc0a48e57cc27d709956be675.svg?invert_in_darkmode&sanitize=true" align=middle width=191.47775955pt height=39.61800315pt/></p>

<p align="center"><img src="/tex/71ef133b5aa04b7b5d1a16fbec1673a7.svg?invert_in_darkmode&sanitize=true" align=middle width=259.61135475pt height=43.72759214999999pt/></p>

<p align="center"><img src="/tex/e1ca1b5924d6177eff4293c8478975f5.svg?invert_in_darkmode&sanitize=true" align=middle width=195.65030595pt height=39.2184573pt/></p>

<p align="center"><img src="/tex/74ae5ee5383c7cdcb51cc74a5074fc4f.svg?invert_in_darkmode&sanitize=true" align=middle width=191.47775955pt height=39.61800315pt/></p>

<p align="center"><img src="/tex/7c5c5607a353042c445e4bfef8ab9faf.svg?invert_in_darkmode&sanitize=true" align=middle width=259.61135475pt height=43.72759214999999pt/></p>

Howver, due to a staggered grid being used the discretizations require an extra step of complexity, detailed below. 



## Grid Naming and Notation



The grid has following naming notations:

<img src="./References/Document Sources/Notation.png" alt="./References/Document Sources/Notation.png" style="zoom:80%;" />



where <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> is the grid number in axial direction, <img src="/tex/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/> the grid number in radial direction. The notation is using 1-index rule as Matlab is also 1-index (starting the index from 1). Each row is representing an axial coordinate and each column representing a radial coordinate. The data is using [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). Note that we are using staggered grid and the properties are having different dimensions.



## Projection method step 1 implementation in staggered grids

<img src="./References/Document Sources/Staggered grid.jpg" style="zoom:50%;" />

### Radial direction momentum equation

The radial momentum equation discretization is shown below:




<p align="center"><img src="/tex/367c48c6b0933d14354735ee9dc78818.svg?invert_in_darkmode&sanitize=true" align=middle width=309.64548614999995pt height=37.1774601pt/></p>

<p align="center"><img src="/tex/9b92aae872efd9c1b35726db243c8c54.svg?invert_in_darkmode&sanitize=true" align=middle width=417.15232184999996pt height=39.8302311pt/></p>

<p align="center"><img src="/tex/fef2e3fbea22a1a08cc1e119045fb33f.svg?invert_in_darkmode&sanitize=true" align=middle width=483.47846415pt height=43.9398201pt/></p>

<p align="center"><img src="/tex/948260ab6f92cb8397a91b0602178826.svg?invert_in_darkmode&sanitize=true" align=middle width=417.15232184999996pt height=39.8302311pt/></p>

<p align="center"><img src="/tex/09c041205d15e51546ad04e7fa294567.svg?invert_in_darkmode&sanitize=true" align=middle width=483.47846415pt height=43.9398201pt/></p>


<p align="center"><img src="/tex/a8c7eb7c20a9657f030bc08ce698ef59.svg?invert_in_darkmode&sanitize=true" align=middle width=335.76136275pt height=39.8144505pt/></p>


The advection term can be discretized as:
<p align="center"><img src="/tex/e0a03c688828cae924df852a766f79cd.svg?invert_in_darkmode&sanitize=true" align=middle width=771.1729624499999pt height=41.11867485pt/></p>



The diffusion term is:

<p align="center"><img src="/tex/ed5a150e229f2fff26e109e94bcf18be.svg?invert_in_darkmode&sanitize=true" align=middle width=697.01142885pt height=49.315569599999996pt/></p>



The intermediate variable can be summarized by the following linearized equation:


<p align="center"><img src="/tex/eaeb8dc11f27a98826a9bd7865654e77.svg?invert_in_darkmode&sanitize=true" align=middle width=442.07507024999995pt height=19.490725649999998pt/></p>


where
<p align="center"><img src="/tex/8f33e697eb8c579e1f44870ec1afc556.svg?invert_in_darkmode&sanitize=true" align=middle width=300.9689595pt height=39.8096622pt/></p>




<p align="center"><img src="/tex/84f7a99adc51ab127d55223dd7740e9a.svg?invert_in_darkmode&sanitize=true" align=middle width=215.16971025pt height=37.6933392pt/></p>




<p align="center"><img src="/tex/5b6692a9ca8be0f642e8bb1f246d400e.svg?invert_in_darkmode&sanitize=true" align=middle width=214.8009435pt height=37.6933392pt/></p>




<p align="center"><img src="/tex/9dc225eafcdbab1765cff21cd5b61922.svg?invert_in_darkmode&sanitize=true" align=middle width=240.93375074999997pt height=37.099754999999995pt/></p>




<p align="center"><img src="/tex/dce73c6fa4e5e38c1f6b8687c2738a44.svg?invert_in_darkmode&sanitize=true" align=middle width=754.6750529999999pt height=41.11867485pt/></p>


<p align="center"><img src="/tex/1f80ff94ba3536bd1f442bad550ebbcc.svg?invert_in_darkmode&sanitize=true" align=middle width=102.76844655pt height=16.438356pt/></p>


   $$A_r^*=$$
   
  <p align="center"><img src="/tex/ca2734ad8d4b0eb5efe5550fbb32edba.svg?invert_in_darkmode&sanitize=true" align=middle width=131.77355565pt height=188.4931719pt/></p>





where <img src="/tex/1c7887c98fd16c36851ae4c9f576fef7.svg?invert_in_darkmode&sanitize=true" align=middle width=170.69537265pt height=24.65753399999998pt/>, where <img src="/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/> is the z-direction coordinate and <img src="/tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode&sanitize=true" align=middle width=7.710416999999989pt height=21.68300969999999pt/> r-direction coordinate. The equation only applies to internal points, with z positive direction pointing downward and r positive direction pointint rightward (outward).



### Axial direction momentum equation



Like radial, axial momentum equation discretization is shown below:




<p align="center"><img src="/tex/367c48c6b0933d14354735ee9dc78818.svg?invert_in_darkmode&sanitize=true" align=middle width=309.64548614999995pt height=37.1774601pt/></p>





The advection term can be discretized as:
<p align="center"><img src="/tex/ff8896c15d9393b60ae3a2ddfbc8d7aa.svg?invert_in_darkmode&sanitize=true" align=middle width=772.2085040999999pt height=41.11867485pt/></p>



The diffusion term is:

<p align="center"><img src="/tex/360d75162335835d8a43e45062bf8b8e.svg?invert_in_darkmode&sanitize=true" align=middle width=691.6483876499999pt height=49.315569599999996pt/></p>



The intermediate variable can be summarized by the following linearized equation:


<p align="center"><img src="/tex/f4b4b7c970717d8a9f0f662c3d13265f.svg?invert_in_darkmode&sanitize=true" align=middle width=492.32601164999994pt height=34.4904021pt/></p>


where in state matrix form we have:
<p align="center"><img src="/tex/46527258472a972637eabc5732735851.svg?invert_in_darkmode&sanitize=true" align=middle width=234.58740195pt height=37.099754999999995pt/></p>




<p align="center"><img src="/tex/7b152f2d4a7fecdcff4fe0d7461fe9bc.svg?invert_in_darkmode&sanitize=true" align=middle width=273.27090615pt height=37.6933392pt/></p>




<p align="center"><img src="/tex/8ef2b205931e0405244b19a7dfc39ec1.svg?invert_in_darkmode&sanitize=true" align=middle width=272.9021394pt height=37.6933392pt/></p>




<p align="center"><img src="/tex/9dc225eafcdbab1765cff21cd5b61922.svg?invert_in_darkmode&sanitize=true" align=middle width=240.93375074999997pt height=37.099754999999995pt/></p>




<p align="center"><img src="/tex/a2cf7aec727c950abf692d4ccf5ba033.svg?invert_in_darkmode&sanitize=true" align=middle width=756.6804844499999pt height=41.11867485pt/></p>





where <img src="/tex/1c7887c98fd16c36851ae4c9f576fef7.svg?invert_in_darkmode&sanitize=true" align=middle width=170.69537265pt height=24.65753399999998pt/>, where <img src="/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/> is the z-direction coordinate and <img src="/tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode&sanitize=true" align=middle width=7.710416999999989pt height=21.68300969999999pt/> r-direction coordinate. The equation only applies to internal points, with z positive direction pointing downward and r positive direction pointint rightward (outward).



### Step 1 boundary conditions



<img src="./References/Document Sources/Boundary condition.jpg" style="zoom:24%;" />



- Boundary (1): velocity boundary; assume that the inflow velocity is equal to <img src="/tex/98f42838386f103b3bc6a224d72ccb8a.svg?invert_in_darkmode&sanitize=true" align=middle width=27.19988204999999pt height=22.465723500000017pt/>.


<p align="center"><img src="/tex/573d17445e72da82a0685fcf3516d699.svg?invert_in_darkmode&sanitize=true" align=middle width=168.5845029pt height=15.936036599999998pt/></p>


- Boundary (2):  <img src="/tex/4aee50c7fe12fe3e346db76f22ffd69e.svg?invert_in_darkmode&sanitize=true" align=middle width=45.384226799999986pt height=21.18721440000001pt/> and <img src="/tex/f6cdc5360d2d944d54ab5675731bb521.svg?invert_in_darkmode&sanitize=true" align=middle width=77.00521784999998pt height=21.18721440000001pt/>. The stress includes velocity and pressure terms. However, as pressure term does not appear in the step one, the boundary condition is simplified to the following:


<p align="center"><img src="/tex/ba67f53f7960015c02f3a8b987e04644.svg?invert_in_darkmode&sanitize=true" align=middle width=322.49373585pt height=39.66415035pt/></p>




- Boundary (3):  normal stress =0. The stress has a viscous component and the pressure. 


<p align="center"><img src="/tex/2021ae14f98ddd70580c474d316fb019.svg?invert_in_darkmode&sanitize=true" align=middle width=392.38602149999997pt height=41.2237353pt/></p>




- Boundaries 4 and 5:  no slip.


<p align="center"><img src="/tex/3840a2e4877ce728c63faa696a8f9001.svg?invert_in_darkmode&sanitize=true" align=middle width=644.7930819pt height=39.452455349999994pt/></p>



## Projection method step 2 implementation in staggered grids



The Possion equation needs to be solved:
<p align="center"><img src="/tex/7fd490426be54f30979bad7c126b3048.svg?invert_in_darkmode&sanitize=true" align=middle width=144.93071999999998pt height=32.990165999999995pt/></p>
After the first step.





A divergence operator <img src="/tex/3715c81142a6f4a78a29106b9a6799bd.svg?invert_in_darkmode&sanitize=true" align=middle width=42.80807519999998pt height=22.63846199999998pt/> has the following expression in 2d cylindrical coordinates:


<p align="center"><img src="/tex/9d5797db3163640ffc6af6338a0c09ff.svg?invert_in_darkmode&sanitize=true" align=middle width=331.87325655pt height=34.7253258pt/></p>


The discretization would need to be also at the node point of  pressure <img src="/tex/ddd51eabbdc98d8592c0e3bcd64764d8.svg?invert_in_darkmode&sanitize=true" align=middle width=56.866378799999985pt height=24.65753399999998pt/>:
<p align="center"><img src="/tex/6a146366189bc34cd9c240c1b33dd481.svg?invert_in_darkmode&sanitize=true" align=middle width=523.51350315pt height=37.8237354pt/></p>


The Laplace operator has the following equation in the 2d cylindrical coordinates:

 
<p align="center"><img src="/tex/da3ce264e2ecf34f57dc6d9ec8d1c18e.svg?invert_in_darkmode&sanitize=true" align=middle width=353.62222004999995pt height=40.11819404999999pt/></p>


The discretization method would be:
<p align="center"><img src="/tex/bc226e8a1d8a15153cbb058ab254f6a1.svg?invert_in_darkmode&sanitize=true" align=middle width=575.6432979pt height=38.8282917pt/></p>
where, as in step 1, <img src="/tex/c3d1dea018f98330b6e64554fe26423a.svg?invert_in_darkmode&sanitize=true" align=middle width=123.16109354999999pt height=24.65753399999998pt/>.



The discretization allows us to get, for all internal nodes:
<p align="center"><img src="/tex/acdbb5809b2ca484849c367b1bfbc457.svg?invert_in_darkmode&sanitize=true" align=middle width=393.04493414999996pt height=17.0776386pt/></p>


where
<p align="center"><img src="/tex/623c908c1d52288e9a729cf5534f6db3.svg?invert_in_darkmode&sanitize=true" align=middle width=144.21308385pt height=37.099754999999995pt/></p>




<p align="center"><img src="/tex/b23318b4a3ea8437ba4e5c37294dff94.svg?invert_in_darkmode&sanitize=true" align=middle width=226.56772379999998pt height=37.6933392pt/></p>




<p align="center"><img src="/tex/caccbd4ece2f4f220e54ede23188888b.svg?invert_in_darkmode&sanitize=true" align=middle width=213.4135245pt height=37.6933392pt/></p>




<p align="center"><img src="/tex/dbafe675788d8f5585d9cfa02daf71fc.svg?invert_in_darkmode&sanitize=true" align=middle width=194.23056839999998pt height=37.099754999999995pt/></p>




<p align="center"><img src="/tex/97ed8e70975e6798ee64e100a62efbe3.svg?invert_in_darkmode&sanitize=true" align=middle width=545.2306843499999pt height=39.814869599999994pt/></p>



Df 

## Validation of the model
The results of our numerical scheme are compard to the results of K. F. Zhang and  J. Y. OOi (1998). 
### Establishing a relationship between viscosity and properties of solid grain particles
K. F. Zhang and  J. Y. OOi (1998) investigated the kinematic constant, B, as a function of various particle parameters. Has been shown to most closely related to particle size. The kinematic eqn they use is identical to the 1D unsteady heat conduction equation. 

[1]:	https://www.overleaf.com/read/hzzczmvjnnht
