# Stefan problems

By setting `activeT = 0` in `bou.in`, we can decouple the temperature and velocity fields and model the class of problems involving only the diffusion of heat coupled to phase changes, known as Stefan problems.
Some of the results obtained from running these examples can be found on the [validation page](../phasefield_validation.md).

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

Input files to run this example are provided in [`examples/1DFreezing`](https://github.com/chowland/AFiD-MuRPhFi/tree/main/examples/1DFreezing).
For this example, we use $\mathcal{S}=1$, which requires $\Lambda=0.62$.
The value of $\Lambda$ is currently hard-coded into `CreateInitialConditions`, so if you want to test other values of $\mathcal{S}$ then this will need modifying.
Since the domain is 1-D, this example must be run in serial (i.e. with `mpiexec -n 1`.

### 1-D melting from a hot boundary
For a pure substance, i.e. one where melting and freezing solely depends on the temperature, the inverse problem is mathematically identical.
Therefore we can also validate against the case of an initially isothermal solid at the melting temperature $T=T_m = 0$ melting due to a boundary whose temperature is higher at $T=1$.
The above initial condition can be used by simply mapping $T \rightarrow 1 - T$.

The input file for this example is provided in [`examples/1DMelting`](https://github.com/chowland/AFiD-MuRPhFi/tree/main/examples/1DMelting).

## 1-D freezing into a supercooled melt
[Davis (2001)](https://doi.org/10.1017/CBO9780511546747) also presents an example where one-dimensional freezing occurs into a liquid that is initially at a temperature *below* the melting temperature.
In this example, the liquid is initially at $T=0$ and the melting temperature is higher at $T_m=1$.
As the system evolves, the solid occupies the domain $z<h(t)$ and remains isothermal at the melting temperature.
In the liquid, a similarity solution can be derived by imposing the boundary condition $T\rightarrow 0$ as $x\rightarrow \infty$, using the same similarity variable $\eta=x/2\sqrt{\kappa t}$ as above:

$$
T(\eta) = \frac{\mathrm{erfc}(\eta)}{\mathrm{erfc}(\Lambda)} , \qquad
\sqrt{\pi} \Lambda e^{\Lambda^2} \mathrm{erfc}(\Lambda) = \mathcal{S}^{-1} .
$$

Here, $\mathrm{erfc}(x) = 1 - \mathrm{erf}(x)$ is the complementary error function.
Note that the second equation above, which determines $\Lambda$, has **no** solutions for $\mathcal{S}\leq 1$.
For such small Stefan numbers, this highlights a breakdown of the thermodynamic equilibrium assumption at the interface.
We are not concerned with such problems, which could be addressed by including the effect of kinetic undercooling.
For the validation of this case, we shall simply set the Stefan number higher at $\mathcal{S}=10$, at which $\Lambda=0.060314$.
Setting the initial interface position to $x=0.1$ in this problem would result in a temperature significantly above zero at the upper boundary.
We therefore set the intial interface position to be $x=0.02$, for which the effect of enforcing $T=0$ at $x=1$ will be minimal.

The input file for this example is provided in [`examples/1DSupercooling`](https://github.com/chowland/AFiD-MuRPhFi/tree/main/examples/1DSupercooling).

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

The input file for this example is provided in [`examples/2DAxisymMelting`](https://github.com/chowland/AFiD-MuRPhFi/tree/main/examples/2DAxisymMelting).
Here we take $\mathcal{S}=1$ to match the example presented in [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773).

## 2-D axisymmetric growth of a disc from a supercooled melt
Finally, we consider the problem originally considered by [Frank (1950)](https://doi.org/10.1098/rspa.1950.0080) of the radially symmetric growth of a solid due to diffusion in a supercooled melt.
In 2-D the problem considers the growth of a disc at the origin in an unbounded domain.
The solid that grows is isothermal at the melting temperature $T_m=1$, and the liquid phase is supercooled such that it tends to a far field temperature of $T_\infty=0$ as $r\rightarrow\infty$.
The solution for the radial temperature profile is then given by

$$
T(r,t) = \frac{F(s)}{F(\Lambda)}, \quad \textrm{where } F(s) = E_1\left(\frac{s^2}{4}\right), \quad \textrm{and } E_1(z) = \int_z^\infty \frac{e^{-t}}{t} \mathrm{d}t .
$$

As in the above examples, $s=r\sqrt{Pe/t}$ is a similarity variable, and $\Lambda$ is the constant value of $s$ associated with the phase boundary.
$\Lambda$ satisfies the following solvability condition:

$$
\mathcal{S} \frac{\Lambda^2}{4}e^{\Lambda^2/4} F(\Lambda) = 1 .
$$

We fix the Stefan number to be 2.5, giving a value for $\Lambda$ of 1.2012.
The initial radius for the solid in the simulation is chosen to be $r=0.1$ such that the unbounded solution for the temperature field is close to zero at the boundaries of our finite domain.
This is associated with an initial time of $t_0=25.3$.

## 1-D multicomponent melting ([Martin and Kauffman, 1977](https://doi.org/10.1175/1520-0485(1977)007%3C0272:AEATSO%3E2.0.CO;2))

We consider the melting of a solid ice phase ($x>h(t)$) by a liquid phase ($x<h(t)$) that has a far-field temperature $T_0$ and a far-field concentration of salinity $C_0$.
The ice phase contains no salt, so has concentration $C=0$, and is uniform in temperature.
The ice temperature is not known *a priori* but can be derived from the interfacial boundary conditions.

We consider the problem as purely diffusive, so the governing equations in both phases are simply diffusion equations.
We nondimensionalise the temperature and concentration fields as

$$
T = T_e + \Delta T \ T^\ast, \qquad
C = C_0 C^\ast ,
$$

where $T_e=0 \degree\mathrm{C}$ is the equilibrium melting temperature and $\Delta T=T_0 - T_e$.

Dropping the asterisks, the dimensionless governing equations then read

$$
\partial_t T = \kappa_T \partial_{xx} T, \qquad
\partial_t C = \kappa_S \partial_{xx} C,
$$

using $\kappa_T$ as shorthand for $Pe^{-1}$, with the far-field boundary conditions

$$
T\rightarrow 1 \quad C\rightarrow 1 \quad \textrm{as} \quad x\rightarrow -\infty ,
$$

and the interfacial boundary conditions

$$
T+\Lambda C = 0, \qquad \mathcal{S}h'=-\kappa_T \partial_x T, \qquad
Ch' = -\kappa_S \partial_x C,
$$

at $x=h(t)$.
Here $\Lambda=\lambda \Delta T/C_0$ is the dimensionless liquidus slope where $\lambda=5.71\times 10^{-2} \ \degree\mathrm{C} (\mathrm{psu})^{-1}$, and $\mathcal{S}=L/c_p\Delta T$ is the Stefan number.

Note that although $C=0$ in the solid phase, $C$ can be non-zero in the liquid phase at the interface, meaning that the concentration field can jump at $x=h(t)$.

The solution reads

$$
T(x,t) = \begin{cases}
1 - A\mathrm{erfc}(-\eta) & x<h(t) \\
1 - A\mathrm{erfc}(-\alpha) &  x\geq h(t)
\end{cases}
$$

$$
C(x,t) = \begin{cases}
1 - B\mathrm{erfc}(-\eta/\sqrt{\tau}) & x<h(t) \\
0 &  x\geq h(t)
\end{cases}
$$

where

$$
\eta = \frac{x}{2\sqrt{\kappa_T t}}
$$

is the similarity variable and $\tau=\kappa_S/\kappa_T$ is the diffusivity ratio.
The interface moves as $h(t)=2\alpha \sqrt{\kappa_T t}$ for the constant value $\alpha$.
The coefficients $\alpha$, $A$, and $B$ can be determined through the boundary conditions.

Specifically, they satisfy the following relations:

$$
\alpha e^{\alpha^2} \mathrm{erfc}(-\alpha) = \frac{1+\Lambda}{\mathcal{S}\sqrt{\pi}}\left[ 1 - \frac{\Lambda}{1+\Lambda} f(\alpha/\sqrt{\tau})\right] , \quad \textrm{where} \ f(\gamma) = \frac{\gamma \sqrt{\pi} \mathrm{erfc}(-\gamma)}{e^{-\gamma^2} + \gamma \sqrt{\pi}\mathrm{erfc}(-\gamma)} ,
$$

$$
\alpha e^{\alpha^2} = \frac{A}{\mathcal{S}\sqrt{\pi}} , \qquad A = \frac{1 + \Lambda - \Lambda B \mathrm{erfc}(-\alpha/\sqrt{\tau})}{\mathrm{erfc}(-\alpha)} .
$$

The coefficients are fully determined once values of $\Lambda$, $\tau$ and $\mathcal{S}$ are specified.
For our validation case, we will take $\mathcal{S}=2.5$, $\tau=0.1$ and $\Lambda =0.4$, giving the following values for the solution coefficients:

$$
A=0.90954, \quad B=0.44748, \quad \alpha = 0.19742 .
$$

An animation of the solution is shown below, where temperature is plotted as a red line, concentration is in blue, and the interface position is shown by a vertical black line.
The origin for the similarity variable has been shifted such that the initial position of the interface is at $x=0.8$.
Note the constant values of $T$ and $C$ and the interface, as well as the significant jump in $C$ at the interface.

<video width="100%" controls>
  <source src="../../assets/1D_DD_analytic.mp4" type="video/mp4">
</video>