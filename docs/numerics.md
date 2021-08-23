# Numerical schemes

## Time stepping: 3rd order Runge-Kutta
For the time-stepping, the evolution of each state variable ($u$, $v$, $w$, $T$, $S$, $\phi$) is essentially treated in the same way.
The wall-normal diffusion (associated with a $\partial_xx$ term) is treated semi-implicitly, and all other terms are treated explicitly.
A third-order Runge-Kutta scheme is used:

$$
\frac{f^{l+1} - f^l}{\Delta t} = \gamma_l H^l + \rho_l H^{l-1} + \alpha_l \mathcal{L} \left(\frac{f^{l+1} + f^l}{2}\right),
$$

where $H^l$ denotes the explicit terms of the equation at the beginning of Runge-Kutta substep $l$.
The coefficients of the time stepper are those given by [Rai and Moin (1991)](https://doi.org/10.1016/0021-9991(91)90264-L):

$$
\gamma_1 = 8/15, \quad \gamma_2 = 5/12, \quad \gamma_3 = 3/4, \\
\rho_1 = 0, \quad \rho_2 = -17/60, \quad \rho_3 = -5/12 ,
$$

and $\alpha_l = \gamma_l + \rho_l$.
The pressure term in the momentum equation is treated using a split time step to ensure the velocity field remains divergence free.
This is detailed below.

### Crank-Nicolson (semi-implicit) diffusion
The only term in each equation that is not treated explicitly is the wall-normal component of the diffusion.
This term is treated semi-implicitly, which we illustrate below with a simple 1D diffusion scheme to discretise $\partial f/\partial t = \kappa \partial^2 f/\partial x^2$:

$$
\frac{f_i^{n+1} - f_i^n}{\Delta t} = \kappa \mathcal{L}\left(\frac{f_i^{n+1} + f_i^n}{2}\right) = \kappa \mathcal{L}\left(\frac{f_i^{n+1} - f_i^n}{2} + f_i^n\right) .
$$

Here $\mathcal{L}$ is the discrete form of the second spatial derivative.
Writing the scheme as above allows us to rearrange and solve for the increment $\Delta f_i = f_i^{n+1} - f_i$:

$$
\left(1 - \frac{\kappa \Delta t}{2}\mathcal{L}\right) \Delta f_i = \kappa \Delta t \mathcal{L} f_i^n
$$

This amounts to solving a tridiagonal matrix equation, for which AFiD uses a fast LAPACK routine `dgttrs`.
This step is performed in the subroutines named `SolveImpEqnUpdate_XX`.
Note that for a Runge-Kutta substep, $\Delta t$ in the above equation should be replaced with $\alpha \Delta t$.

### Pressure solver - split time step
For the momentum equations, each Runge-Kutta substep is split into two steps - one to evolve the velocity due to advection, diffusion, buoyancy etc., and one to solve for the pressure correction needed to ensure that $\boldsymbol{\nabla} \cdot \boldsymbol{u}^{l+1} = 0$.
An intermediate velocity $\hat{u}_i$ is obtained by the RK step

$$
\frac{\hat{u}_i - u_i^l}{\Delta t} = \gamma_l H^l_i + \rho_l H^{l-1}_i - \alpha_l \mathcal{G}_ip^l + \alpha_l \mathcal{L} \left(\frac{u_i^{l+1} + u_i^l}{2}\right),
$$

where $\mathcal{G}_i$ is the discrete form of the gradient operator and $p^l$ is the pressure at the start of the substep.
The continuity equation $\boldsymbol{\nabla}\cdot\boldsymbol{u}=0$ is then enforced by

$$
u_i^{l+1} - \hat{u}_i = -\alpha_l \Delta t \mathcal{G}_i \Phi^{l+1} ,
$$

where $\Phi$ is the pressure correction.

## Staggered grid layout

