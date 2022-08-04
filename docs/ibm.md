# Immersed boundary method

## Modification to the implicit solver

When using the immersed boundary method, we use the direct forcing method as described by [Fadlun et al. (2000)](https://doi.org/10.1006/jcph.2000.6484).
We show below the modifications to the implicit solve for the wall-normal velocity component $u$ (or `vx`) but the other components are analogous.
The new equation to solve is

$$
\left[ 1 - \frac{\alpha_l \Delta t}{2} \nu F(\boldsymbol{x}) \mathcal{L}\right] \Delta u = RHS^* .
$$

Recall that $\Delta u = u^{l+1} - u^l$ is the velocity increment over the RK substep.
For all points away from the solid immersed boundaries, $RHS^*$ is simply the sum of the explicit terms and the explicit portion of the Crank-Nicolson term, and $F(\boldsymbol{x})=1$ so that the implicit step is unchanged in the liquid portion of the domain.

### Solid points

Inside the solid, $F(\boldsymbol{x})=0$ and $RHS^*=-u^l$ are imposed.
This enforces zero velocity in the solid phase as we see by writing out

$$
[1 - 0]\Delta u = u^{l+1} - u^l = RHS^* = -u^l \ \Rightarrow u^{l+1}=0
$$

### Boundary points

The boundary points are where we must take the most care.
These points are defined as being those grid points adjacent (in the wall-normal direction `x`) to the immersed solids.
Again, we set $F(\boldsymbol{x})=0$ at these points, so do not solve the Navier--Stokes equations.
Rather, we aim to interpolate the velocity linearly between the boundary and the second point into the liquid:

$$
u_k^{l+1} \approx \frac{\delta_2}{\delta_1 + \delta_2} u_{k+1}^{l+1} ,
$$

where $\delta_1$ is the $x$-distance between the first two grid points in the fluid, and $\delta_2$ is the $x$-distance between the first fluid grid point and the solid boundary.

<!-- We also have to approximate the temporal evolution of $V_{k+1}$ to calculate this.
To do this, we record the previous time step and assume the acceleration remains approximately constant close to the boundary

$$
\frac{u_{k+1}^{l+1} - u_{k+1}^l}{\alpha_l \Delta t} \approx \frac{u_{k+1}^l - u_{k+1}^{l-1}}{\alpha_{l-1} \Delta t} \ \Rightarrow u_{k+1}^{l+1} \approx u_{k+1}^l + \frac{\alpha_l}{\alpha_{l-1}} (u_{k+1}^l - u_{k+1}^{l-1})
$$

Putting all this together, the expression used for $RHS^*$ is

$$
\Delta u = RHS^* = -u_k^l + \frac{\delta_2}{\delta_1 + \delta_2} \left(u_{k+1}^l + \frac{\alpha_l}{\alpha_{l-1}} (u_{k+1}^l - u_{k+1}^{l-1}) \right)
$$ -->
Since the interpolation involves the adjacent velocity at the next time step, this requires a further modification to the implicit solver.
Rewriting in terms of the velocity increments, the condition we want to impose is

$$
\Delta u_k - \frac{\delta_2}{\delta_1 + \delta_2} \Delta u_{k+1} = \frac{\delta_2}{\delta_1 + \delta_2} u_{k+1}^l - u_k^l
$$

Note also that $k+1$ will be replaced by $k-1$ for boundary points where the liquid phase is below the solid boundary.
We modify the matrix solve to include this by adding, for example

$$
- (1 - F(x_k)) F(x_{k+1}) \frac{\delta_2}{\delta_1 + \delta_2}
$$

to the tridiagonal term.
Recall that $F(x)$ is zero in the solid *and* at the boundary points, so this term is only nonzero when $x_k$ is a boundary point and $x_{k+1}$ is in the fluid.
An analogous term is added to the lower diagonal.
This step is performed in the `SolveImpEqnUpdate*_ibm` routines, where the distance ratio in the above expressions is stored in the explicit term storage array (e.g. `qcap`, `dph`, `dq`) which is no longer needed during the implicit solve.

### Temperature at the boundary

The above treatment of the boundary points is a special case enforcing a fixed velocity of zero.
For the temperature field, sometimes we want to enforce a specific non-zero temperature in the solid.
This can be achieved simply by prescribing

$$
T_k = T_b + \frac{T_{k+1} - T_b}{\delta_1 + \delta_2} \delta_2 = \frac{\delta_2}{\delta_1 + \delta_2} T_{k+1} + \left( 1 - \frac{\delta_2}{\delta_1 + \delta_2} \right) T_b
$$

where $T_b$ is the fixed temperature enforced in the solid boundaries.

### No-flux conditions

Imposing zero-flux conditions at the immersed boundary is a bit more troublesome than enforcing a fixed value.
Due to the complex nature of the three-dimensional semi-implicit time stepping on a staggered grid, we enforce zero gradients in *every* direction for points across the boundaries rather than only the normal gradient.
In the wall-normal direction, we simply apply a similar approach to that above.
Considering the concentration field $C$, we assume that $\partial_x C=0$ across the boundary by setting $C_{k}$ equal to $C_{k+1}$:

$$
C_k^{l+1} = C_{k+1}^{l+1} \ \Rightarrow \Delta C_k - \Delta C_{k+1} = C_{k+1}^l - C_k^l
$$

This is implemented in the same way as above into the implicit solver, just without any length ratios since they are not needed here.

In the horizontal directions, there is no implicit step, so we cannot apply direct forcing for the immersed boundary method.
To apply the zero-gradient condition, we locate the points adjacent to the boundary.
At these points, we can write out the diffusive term as

$$
\kappa \partial_{yy} C = \frac{\kappa (C_{j+1} - C_j) - \kappa (C_j - C_{j-1})}{(\Delta y)^2}
$$

and set $C_{j-1}=C_j$ to eliminate half of the expression (where we consider point $y_{j-1}$ to be in the solid phase) and impose $\partial_y C=0$ across the boundary.
For diffusion in the $x$-direction, we do not need to hard-code the modified diffusion since the direct immersed boundary forcing overrides the diffusion at the relevant points.