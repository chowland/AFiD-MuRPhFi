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
\frac{f_k^{n+1} - f_k^n}{\Delta t} = \kappa \mathcal{L}\left(\frac{f_k^{n+1} + f_k^n}{2}\right) = \kappa \mathcal{L}\left(\frac{f_k^{n+1} - f_k^n}{2} + f_k^n\right) .
$$

Here $\mathcal{L}$ is the discrete form of the second spatial derivative.
Writing the scheme as above allows us to rearrange and solve for the increment $\Delta f_k = f_k^{n+1} - f_k$:

$$
\left(1 - \frac{\kappa \Delta t}{2}\mathcal{L}\right) \Delta f_k = \kappa \Delta t \mathcal{L} f_k^n
$$

This amounts to solving a tridiagonal matrix equation, for which AFiD uses a fast LAPACK routine `dgttrs`.
This step is performed in the subroutines named `SolveImpEqnUpdate_XX`.
Note that for a Runge-Kutta substep, $\Delta t$ in the above equation should be replaced with $\alpha \Delta t$.

### Pressure solver - split time step
For the momentum equations, each Runge-Kutta substep is split into two steps - one to evolve the velocity due to advection, diffusion, buoyancy etc., and one to solve for the pressure correction needed to ensure that $\boldsymbol{\nabla} \cdot \boldsymbol{u}^{l+1} = 0$.
An intermediate velocity $\hat{u}_k$ is obtained by the RK step

$$
\frac{\hat{u}_k - u_k^l}{\Delta t} = \gamma_l H^l_k + \rho_l H^{l-1}_k - \alpha_l \mathcal{G}_kp^l + \alpha_l \mathcal{L} \left(\frac{u_k^{l+1} + u_k^l}{2}\right),
$$

where $\mathcal{G}_k$ is the discrete form of the gradient operator and $p^l$ is the pressure at the start of the substep.
The continuity equation $\boldsymbol{\nabla}\cdot\boldsymbol{u}=0$ is then enforced by

$$
u_k^{l+1} - \hat{u}_k = -\alpha_l \Delta t \mathcal{G}_k \Phi^{l+1} ,
$$

where $\Phi$ is the pressure correction.

The equation for the pressure correction arises by taking the divergence of the above equation, such that

$$
\boldsymbol{\nabla}\cdot\boldsymbol{u}^{l+1} = 0 \approx \boldsymbol{\nabla}\cdot\boldsymbol{\hat{u}} - \alpha_l \Delta t \nabla^2 \Phi^{l+1} .
$$

The derivatives in the periodic directions are computed using Fourier transforms, whereas the wall-normal derivatives use a finite difference.
Denoting $\boldsymbol{\tilde{u}}$ and $\tilde{\Phi}$ as the ($yz$-)Fourier transforms of $\boldsymbol{\hat{u}}$ and $\Phi^{l+1}$ respectively, we can write

$$
\mathcal{D}\tilde{u} - ik_y \tilde{v} - ik_z\tilde{w} = \alpha_l \Delta t\left( \mathcal{DG}\tilde{\Phi} - (k_y^2 + k_z^2)\tilde{\Phi}\right) ,
$$

where $\mathcal{D}$ is the wall-normal discrete divergence operator.
Note that since we are aiming to set the divergence of $\boldsymbol{u}^{l+1}$ to zero, we use the composite discrete operator $\mathcal{DG}$ for the second derivative of $\Phi$ rather than the discrete Laplacian $\mathcal{L}$ used above in the equations of motion.
These are *not* the same operator.
The above equation can be expressed as a tridiagonal matrix problem for each Fourier mode and solved using similar LAPACK routines as described above for the semi-implicit step.
This is performed in the `SolvePressureCorrection` subroutine.

For completeness, the composite operator $\mathcal{DG}$ takes the form

$$
\mathcal{DG}\Phi_k = a_+ (\Phi_{k+1} - \Phi_k) - a_- (\Phi_k - \Phi_{k-1}) ,
$$

$$
a_+ = [(x^m_{k+1} - x^m_k)(x^c_{k+1} - x^c_k)]^{-1}, \qquad
a_+ = [(x^m_k - x^m_{k-1})(x^c_{k+1} - x^c_k)]^{-1},
$$

where $x^c$ is the wall-normal velocity grid, and $x^m$ is the cell-centre grid.

## Staggered grid layout

