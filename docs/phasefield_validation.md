# Phase-field model validation

## 1-D Freezing from a cold boundary

First, we shall validate the phase-field implementation against the classical textbook example of a growing solidifaction front such as that presented in [Davis (2001)](https://doi.org/10.1017/CBO9780511546747).
This example does not involve any fluid flow, so the dynamics are entirely determined by the 1-D diffusion equation for temperature

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

### 1-D melting from a hot boundary
For a pure substance, i.e. one where melting and freezing solely depends on the temperature, the inverse problem is mathematically identical.
Therefore we can also validate against the case of an initially isothermal solid at the melting temperature $T=T_m = 0$ melting due to a boundary whose temperature is higher at $T=1$.
The above initial condition can be used by simply mapping $T \rightarrow 1 - T$.