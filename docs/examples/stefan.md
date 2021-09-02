# Stefan problems

By setting `activeT = 0` in `bou.in`, we can decouple the temperature and velocity fields and model the class of problems involving only the diffusion of heat coupled to phase changes, known as Stefan problems.

**Note:** Due to the strong coupling of the temperature field to the phase-field variable $\phi$, the phase-field method is *not* well suited to modelling the growth of solids in a supercooled liquid.

## 1-D solidification from a cooled boundary
First, we consider the classical textbook example of a growing solidifaction front such as that presented in [Davis (2001)](https://doi.org/10.1017/CBO9780511546747).
As a standard Stefan problem, this example does not involve any fluid flow, so the dynamics are entirely determined by the 1-D diffusion equation for temperature

$$
\partial_t T =\frac{1}{Pe}\partial_{xx} T,
$$

and the boundary conditions.
The initial condition considers an isothermal fluid at the melting temperature $T_m$, which in our dimensionless framework we can set to $T=1$.
At $t=0$, the lower boundary at $x=0$ is set to a lower temperature $T=0$, and the fluid begins to solidify from the bottom.
At the moving solidification front $h(t)$, the boundary conditions are

$$
T = T_m = 1, \qquad \mathcal{S}\frac{dh}{dt} = \frac{1}{Pe}\frac{\partial T}{\partial x} ,
$$

where $\mathcal{S}$ is the Stefan number.
Using the similarity variable $\eta = \frac{x}{2\sqrt{t/Pe}}$, this problem can be solved analytically as

$$
T(x,t) = \frac{\mathrm{erf}(\eta)}{\mathrm{erf}(\Lambda)} \quad \textrm{for} \quad 0<x<h(t), \qquad h(t) = 2\Lambda \sqrt{t/Pe} ,
$$

where $\Lambda$ depends on the Stefan number through the condition

$$
\sqrt{\pi}\Lambda e^{\Lambda^2} \mathrm{erf}(\Lambda) = \mathcal{S}^{-1} .
$$

When validating our numerical solver against this solution, it makes sense to avoid the singularity at $t=0$, so we can instead start the simulation when the solidification front is already a finite distance from the lower wall.
We take the initial position of the interface to be $x=0.1$, which corresponds to a start time of $t=Pe\left(\frac{h}{2\Lambda}\right)^2 = 0.0025 Pe/\Lambda^2$, and an initial condition for the simulation of

$$
T(x,t) = \frac{1}{\mathrm{erf}(\Lambda)} \mathrm{erf}\left(\frac{\Lambda x}{0.1}\right) \quad \textrm{for} \quad 0<x<0.1, \qquad T = 1 \quad \textrm{for} \quad x>0.1.
$$

Input files to run this example are provided in [`examples/1DFreezing`](https://github.com/chowland/AFiD-MuRPhFi/examples/1DFreezing).
For this example, we use $\mathcal{S}=1$, which requires $\Lambda=0.62$.
The value of $\Lambda$ is currently hard-coded into `CreateInitialConditions`, so if you want to test other values of $\mathcal{S}$ then this will need modifying.
Since the domain is 1-D, this example must be run in serial (i.e. with `mpiexec -n 1`.

### 1-D melting from a hot boundary
For a pure substance, i.e. one where melting and freezing solely depends on the temperature, the inverse problem is mathematically identical.
Therefore we can also validate against the case of an initially isothermal solid at the melting temperature $T=T_m = 0$ melting due to a boundary whose temperature is higher at $T=1$.
The above initial condition can be used by simply mapping $T \rightarrow 1 - T$.

The input file for this example is provided in [`examples/1DMelting`](https://github.com/chowland/AFiD-MuRPhFi/examples/1DMelting).

## Axisymmetric melting of a solid disc in 2-D
This example is taken from appendix A.2 of [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773).
Here, we consider an isothermal solid disc of radius 0.1 that is immersed in an isothermal fluid of higher temperature.
The melting temperature is taken to be the half way between the temperature of the solid and that of the liquid.
This example aims to (a) test the solver in two dimensions, and (b) confirm that curvature has a minimal effect on the melting temperature.

The phase-field method cannot set the surface energy (which appears in our dimensionless framework as the parameter $A$) as small as in the physical problem.
The Gibbs-Thomson effect, by which surface curvature lowers the local melting point, is only significant at very small length scales or high curvatures, which are unlikely to be relevant for the problems we are considering.

Following [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773), we take the intial temperature field to be

$$
T(r,t) = \frac{1 + \tanh(100(r - 0.1))}{2} ,
$$

where $r = \sqrt{(x-0.5)^2 + (y - 0.5)^2}$, so the temperature in the solid disc is $T=0$ and the temperature in the liquid is $T=1$.
The melting temperature is set to $T_m=0.5$.
We provide a smooth initial condition for the phase-field variable as

$$
\phi = \frac{1 - \tanh( (r-0.1)/2\varepsilon)}{2} ,
$$

where $\varepsilon=1/N_x^r$ is the diffuse interface thickness, equal to one wall-normal grid space.

The input file for this example is provided in [`examples/2DAxisymMelting`](https://github.com/chowland/AFiD-MuRPhFi/examples/2DAxisymMelting).
Here we take $\mathcal{S}=1$ to match the example presented in [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773).