# SOLVED: Ice melting in salt water

<!-- Currently, the phase-field implementation is unable to replicate the analytic diffusive solution for the melting of ice in salt water.
On this page, I outline the multicomponent Stefan problem that we wish to replicate, the phase-field model being implemented, and the erroneous results that we obtain. -->
The phase-field implementation now successfully reproduces the analytic diffusive solution for the melting of ice in salt water.
The problem and results are outlined below.

## 1-D diffusive melting ([Martin and Kauffman, 1977](https://doi.org/10.1175/1520-0485(1977)007%3C0272:AEATSO%3E2.0.CO;2))

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
\alpha e^{\alpha^2} \mathrm{erfc}(\alpha) = \frac{1+\Lambda}{\mathcal{S}\sqrt{\pi}}\left[ 1 - \frac{\Lambda}{1+\Lambda} f(\alpha/\sqrt{\tau})\right] , \quad \textrm{where} \ f(\gamma) = \frac{\gamma \sqrt{\pi} \mathrm{erfc}(-\gamma)}{e^{-\gamma^2} + \gamma \sqrt{\pi}\mathrm{erfc}(-\gamma)} ,
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
  <source src="../../assets/diffusive_S2_5.mp4" type="video/mp4">
</video>

## Phase-field model

We have implemented the phase-field model of [Hester et al. (2020)](https://doi.org/10.1098/rspa.2020.0508) to simulate the melting of ice in salty water.
For the validation of this Stefan problem, we decouple the velocity field from the scalars, leaving only the following three equations:

$$
\partial_t T = \kappa_T \nabla^2 T + \mathcal{S} \partial_t \phi ,
$$

$$
\partial_t C = \kappa_S \left(\nabla^2 C - \frac{\boldsymbol{\nabla} \phi \cdot \boldsymbol{\nabla} C}{1 - \phi + \delta}\right) + \frac{C\partial_t \phi}{1 - \phi + \delta},
$$

$$
\partial_t \phi = \mathcal{A}  \nabla^2 \phi - \frac{\mathcal{A}}{\varepsilon^2} \phi (1 - \phi) [1 - 2\phi + \Gamma^{-1}(T +\Lambda C)] .
$$

Here, $\mathcal{A}=6\Gamma\kappa_T/5\mathcal{S}$ is the phase-field diffusion coefficient, $\varepsilon$ is the diffuse interface thickness set to the width of the grid spacing, and $\Gamma$ is a numerical parameter that is linked physically to the concept of interfacial surface energy.
The small factor $\delta$ is used to prevent terms in the concentration equation from blowing up in the solid phase where $\phi=1$.

<!-- Note that these equations differ slightly from those used in the later paper of [Hester et al. (2021)](https://doi.org/10.1103/PhysRevFluids.6.023802) to model the melting of ice blocks.
In that paper,  -->

## Model run results

Below, we present the results obtained from a run of the phase-field model with $Pe_T=10^3$, $Pe_S=10^4$, $\mathcal{S}=2.5$, $\Gamma=1$, and $\Lambda=0.4$.
No-flux boundary conditions are applied to $T$ and $C$ at $x=0$ and $x=1$, and a resolution of $512$ points is used for $T$, whereas $C$ and $\phi$ are resolved on a grid of $1024$ points.

The simulation is initialised with the analytic solution given above at time $t=1$, and an initial phase-field profile of

$$
\phi = \frac{1}{2}\left[ 1 + \tanh \left(\frac{x - h_0}{2\varepsilon} \right) \right]
$$

where $h_0 = 2\alpha /\sqrt{Pe_T}$ is the corresponding interface position.

First, we inspect the temperature profile over time.
<!-- Although the temperature in the solid phase remains lower than zero, we see that it increases over time.
This is in contrast to the analytic solution, which prescribes that the temperature in the solid (and at the interface) should be constant. -->
The profiles match well, and the temperature in the solid remains at the constant value predicted by the analytic solution.

{% include "temp_validation.html" %}

<!-- Inspection of the salinity field reveals the cause for this temperature difference, and the greatest issue with the current phase-field implementation.
The persistent jump in concentration in the analytic solution is smeared out almost straight away.
The condition of zero salinity is not imposed in the solid phase, and this in turn affects the salinity profile in the liquid. -->
Inspection of the salinity field reveals good agreement once again in the liquid, although somewhat spurious behaviour in the solid.
Salinity must be continuous at the phase boundary, which does not match the physical solution, but we already know that the salinity must be zero in the solid phase.
The phase-field model is set up such that the correct boundary conditions are applied at the interface, making the solution in the liquid reliable.

{% include "sal_validation.html" %}

As further evidence of the model's good agreement, we plot the salinity profiles such that the concentration is set to zero in the solid phase ($\phi>0.5$):

{% include "sal_validation_trunc.html" %}

The evolution of the phase-field interface over time also recovers the analytic $t^{1/2}$ expression near-perfectly:

{% include "phi_validation.html" %}

## Notes

- One of the boundary condition expressions in Martin & Kauffman (1977) contains an extra factor of $1/2$ that I am pretty sure should not be there. This only changes the prefactors in the analytic solution. Out of caution, I have tried initial conditions with and without the factor of 2 - there is no significant change to the above results.
- Various values of $\delta$ have been tried, from $10^{-14}$ to $10^{-2}$. None have produced the desired solution.