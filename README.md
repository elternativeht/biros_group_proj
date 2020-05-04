
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


   $$A_r^*=
  \left[ <p align="center"><img src="/tex/2b8dbe58dfa7f2bfb08f5c1a38e8b521.svg?invert_in_darkmode&sanitize=true" align=middle width=131.77355565pt height=188.4931719pt/></p> \right]<p align="center"><img src="/tex/374927edf3ea0e0fc8964a30b9d195e2.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2746965499999pt height=94.42922775pt/></p>
\frac{\boldsymbol{u}^{*}_{(i,j)}-\boldsymbol{u}^{n}_{(i,j)}}{\Delta t}+\left(\boldsymbol{u}^{n}_{(i,j)}\cdot\nabla\right)\boldsymbol{u}^{n}_{(i,j)}=\Delta \boldsymbol{u}^{*}_{(i,j)}
<p align="center"><img src="/tex/5174e3e2cc66af333d6ec31bf1b0be4e.svg?invert_in_darkmode&sanitize=true" align=middle width=300.2747649pt height=11.4155283pt/></p>
\left( \boldsymbol{u}^n_{z,ij} \cdot \nabla \right) \boldsymbol{u}^n_{z,ij}=\frac{u^n_{r,i(j+1)}+u^n_{r,ij}+u^n_{r,(i-1)(j+1)}+u^n_{r,(i-1)j}}{4}\left[\frac{u^{n}_{z,i(j+1)}-u^{n}_{z,i(j-1)}}{2\Delta r}\right]+u^n_{z,ij}\frac{u^n_{z,(i+1)j}-u^n_{z,(i-1)j}}{2\Delta z}
<p align="center"><img src="/tex/d352c2fe29961288bbd34dc377e06478.svg?invert_in_darkmode&sanitize=true" align=middle width=155.02324199999998pt height=11.4155283pt/></p>
\frac{1}{\text{Re}}\left[\frac{1}{r_{ij}+\frac{\Delta r}{2}}\frac{u^*_{z,i(j+1)}-u^*_{z,i(j-1)}}{2\Delta r}+\frac{u^*_{z,i(j+1)}-2u^*_{z,ij}+u^*_{z,i(j-1)}}{(\Delta r)^2}+\frac{u^*_{z,(i+1)j}-2u^*_{z,ij}+u^*_{z,(i-1)j}}{(\Delta z)^2}\right]
<p align="center"><img src="/tex/c3da89483fe22ee76470d058eabb93e8.svg?invert_in_darkmode&sanitize=true" align=middle width=595.7094pt height=14.611878599999999pt/></p>
Au^*_{z,ij}+Bu^*_{z,i(j+1)}+Cu^*_{z,i(j-1)}+Du^*_{z,(i+1)j}+Eu^*_{z,(i-1)j}=F+\frac{1}{\text{Fr}^2}
<p align="center"><img src="/tex/a38577c4da2a0986caf37e1891e2f11f.svg?invert_in_darkmode&sanitize=true" align=middle width=260.50297184999994pt height=11.4155283pt/></p>
A = \frac{1}{\Delta t}+\frac{1}{\text{Re}}\frac{2}{(\Delta r)^2}+\frac{1}{\text{Re}}\frac{2}{(\Delta z)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
B = -\frac{1}{\text{Re}}\frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
C = +\frac{1}{\text{Re}}\frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{\text{Re}}\frac{1}{(\Delta r)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
D =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2},E =  -\frac{1}{\text{Re}}\frac{1}{(\Delta z)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
F= \frac{u^n_{z,ij}}{\Delta t}-\left(\frac{u^n_{r,i(j+1)}+u^n_{r,ij}+u^n_{r,(i-1)(j+1)}+u^n_{r,(i-1)j}}{4}\left[\frac{u^{n}_{z,i(j+1)}-u^{n}_{z,i(j-1)}}{2\Delta r}\right]+u^n_{z,ij}\frac{u^n_{z,(i+1)j}-u^n_{z,(i-1)j}}{2\Delta z}\right)
<p align="center"><img src="/tex/51201ec0ed4a30a6f8b2c36fd65c65fd.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2746965499999pt height=115.06849364999998pt/></p>
\overline{u_{z,1j}}=1,j=1,2,3...N
<p align="center"><img src="/tex/f01ed562c9500c4932aae825e4a5b3ca.svg?invert_in_darkmode&sanitize=true" align=middle width=700.27463055pt height=54.9771717pt/></p>
\left\{<p align="center"><img src="/tex/41adc9f6b482e55c3b37fe026b4f79a6.svg?invert_in_darkmode&sanitize=true" align=middle width=299.973168pt height=37.5498585pt/></p>\right.
<p align="center"><img src="/tex/c67a4f3ce389602857ebf6cc46aeec09.svg?invert_in_darkmode&sanitize=true" align=middle width=639.4081039499999pt height=16.438356pt/></p>
\left\{<p align="center"><img src="/tex/36e376aa8ea769431ddd08c1027e87f6.svg?invert_in_darkmode&sanitize=true" align=middle width=369.86545364999995pt height=41.2237353pt/></p>\right.
<p align="center"><img src="/tex/8efb6f604f0265c0f55230cf84c600e2.svg?invert_in_darkmode&sanitize=true" align=middle width=214.6124112pt height=14.611878599999999pt/></p>
\left\{<p align="center"><img src="/tex/e064c8371fb15edf68c3ed36e28030fe.svg?invert_in_darkmode&sanitize=true" align=middle width=624.2451105pt height=36.8036658pt/></p>\right.
<p align="center"><img src="/tex/189b81a6bae2042f436d3faa8671a324.svg?invert_in_darkmode&sanitize=true" align=middle width=433.06055595pt height=34.3379058pt/></p>
\nabla^{2} \pi^{n+1}=\frac{1}{\Delta t} \nabla \cdot \mathbf{u}^{*}
<p align="center"><img src="/tex/a3872426bb1565e6ba7c604d9a961301.svg?invert_in_darkmode&sanitize=true" align=middle width=628.2891615pt height=44.2009194pt/></p>
\nabla \cdot \mathbf{u}^{*} = \frac{1}{r} \frac{\partial\left(r u_{r}\right)}{\partial r}+\frac{\partial u_{z}}{\partial z}=\frac{1}{r}u_r+\frac{\partial u_r}{\partial r}+\frac{\partial u_{z}}{\partial z}
<p align="center"><img src="/tex/6b141d720bb72fde33c0690a0a945747.svg?invert_in_darkmode&sanitize=true" align=middle width=572.4384765pt height=16.438356pt/></p>
\nabla \cdot \mathbf{u}^{*} \approx \frac{1}{r_{ij}+\Delta r/2}\frac{u_{r,ij}+u_{r,i(j+1)}}{2}+\frac{u_{r,i(j+1)}-u_{r,ij}}{\Delta r}+\frac{u_{z,(i+1)j}-u_{z,ij}}{\Delta z}
<p align="center"><img src="/tex/b6140b7e6872492331e2b0c77292a439.svg?invert_in_darkmode&sanitize=true" align=middle width=591.46284615pt height=14.611878599999999pt/></p>
\nabla^2\pi=\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial \pi}{\partial r}\right)+\frac{\partial^{2} \pi}{\partial z^{2}}=\frac{1}{r}\frac{\partial \pi}{\partial r}+\frac{\partial^2\pi}{\partial r^2}+\frac{\partial^2 \pi}{\partial z^2}
<p align="center"><img src="/tex/f4579c56b81f6fa0a57e39a3c87e35b3.svg?invert_in_darkmode&sanitize=true" align=middle width=264.52128120000003pt height=11.4155283pt/></p>
\frac{1}{r_{ij}+\Delta r/2}\frac{\pi_{i(j+1)}-\pi_{i(j-1)}}{2\Delta r}+\frac{\pi_{i(j+1)}-2\pi_{ij}+\pi_{i(j-1)}}{(\Delta r)^2}+\frac{\pi_{(i+1)j}-2\pi_{ij}+\pi_{(i-1)j}}{(\Delta z)^2}
<p align="center"><img src="/tex/6d2a00b89723c2a30d8abb2b668861bb.svg?invert_in_darkmode&sanitize=true" align=middle width=436.1198094pt height=63.92694825pt/></p>
A\pi_{ij}+B\pi_{i(j+1)}+C\pi_{i(j-1)}+B\pi_{(i+1)j}+B\pi_{(i-1)j}=F
<p align="center"><img src="/tex/f6851476009182f828d80792e4d78020.svg?invert_in_darkmode&sanitize=true" align=middle width=42.05492115pt height=11.4155283pt/></p>
A = \frac{2}{(\Delta r)^2}+\frac{2}{(\Delta z)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
B = -\frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{(\Delta r)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
C = \frac{1}{r_{ij}+\Delta r/2}\frac{1}{2\Delta r}-\frac{1}{(\Delta r)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
D = - \frac{1}{(\Delta z)^2},E =  -\frac{1}{(\Delta z)^2}
<p align="center"><img src="/tex/f8bfc3a483b0e785ac7c71bc777479cb.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>
F=\frac{1}{\Delta t} \left(\frac{1}{r_{ij}+\Delta r/2}\frac{u_{r,ij}+u_{r,i(j+1)}}{2}+\frac{u_{r,i(j+1)}-u_{r,ij}}{\Delta r}+\frac{u_{z,(i+1)j}-u_{z,ij}}{\Delta z}\right)
$$



Df 

## Validation of the model
The results of our numerical scheme are compard to the results of K. F. Zhang and  J. Y. OOi (1998). 
### Establishing a relationship between viscosity and properties of solid grain particles
K. F. Zhang and  J. Y. OOi (1998) investigated the kinematic constant, B, as a function of various particle parameters. Has been shown to most closely related to particle size. The kinematic eqn they use is identical to the 1D unsteady heat conduction equation. 

[1]:	https://www.overleaf.com/read/hzzczmvjnnht
