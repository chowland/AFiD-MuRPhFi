# Immersed boundary method

When using the immersed boundary method, we use the direct forcing method as described by [Fadlun et al. (2000)](https://doi.org/10.1006/jcph.2000.6484).
This must be incorporated into the solution of the implicit step of the time stepper, and is taken care of by the subroutines `SolveImpEqnUpdateXX_ibm` stored in the `ibm` subdirectory.
We distinguish points as either solid/fluid/boundary using the arrays `ibmaskXX`.

## Modification to the implicit solver
We show below the modifications to the implicit solve for the wall-normal velocity component $u$ (or `vx`) but the other components are analogous.
In solving the implicit step, we solve the matrix system

$$
a_{ij} (\Delta u)_j = RHS^*_{i} ,
$$

where in the standard fluid domain $a_{ij}$ defines the wall-normal diffusion operator (see [the time stepper documentation](numerics.md#crank-nicolson-semi-implicit-diffusion) for more details).
Recall that $\Delta u = u^{l+1} - u^l$ is the velocity increment over the RK substep.
For all points away from the solid immersed boundaries, $RHS^*$ is simply the sum of the explicit terms and the explicit portion of the Crank-Nicolson term.

### Solid points (`ibmask==0`)

Inside the solid, we enforce $u^{l+1}=0$.
If we consider the velocity at position $x_k$, this is achieved by setting $a_{kk}=1$, and $a_{kj}=0$ for $j\neq k$, with $RHS^*_k=-u^l_k$.

$$
(1)\Delta u = u^{l+1} - u^l = RHS^* = -u^l \quad \Rightarrow u^{l+1}=0 .
$$

Ideally, in the solid $u^l$ would be zero, but the pressure correction can introduce small velocities in the solid, so it is good practice to ensure the right hand side exactly cancels this previous velocity.

### Boundary points (`ibmask==±1`)

The boundary points are where we must take the most care.
These points are defined as being those grid points adjacent (in the wall-normal direction `x`) to the immersed solids.
For the points just above the interface, we set `ibmask=1` and for points just below a solid we set `ibmask=-1`, and we do not solve the Navier--Stokes equations.
Rather, we interpolate the velocity linearly between the boundary and the second point into the liquid:

$$
u_k^{l+1} = \frac{\delta_2}{\delta_1 + \delta_2} u_{k+1}^{l+1} ,
$$

where $\delta_1$ is the $x$-distance between the first two grid points in the fluid, and $\delta_2$ is the $x$-distance between the first fluid grid point and the solid boundary.

Since the interpolation involves the adjacent velocity at the next time step, this requires a further modification to the implicit solver.
Rewriting in terms of the velocity increments, the condition we want to impose is

$$
\Delta u_k - \frac{\delta_2}{\delta_1 + \delta_2} \Delta u_{k+1} = \frac{\delta_2}{\delta_1 + \delta_2} u_{k+1}^l - u_k^l
$$

Note also that $k+1$ will be replaced by $k-1$ for boundary points where the liquid phase is below the solid boundary.
The matrix for the implicit solver therefore reads

$$
a_{kj} = \begin{cases} 1 & j=k \\ -\delta_2/(\delta_1+\delta_2) & j= k+1 \\ 0 & \textrm{otherwise} \end{cases}
$$

with $RHS^*$ calculated as the right hand side of the equation above.
This step is performed in the `SolveImpEqnUpdate*_ibm` routines, where the distance ratio in the above expressions is stored in the arrays `distx`, `disty`,... which are defined in `topogr_ibm`.

### Temperature at the boundary

The above treatment of the boundary points is a special case enforcing a fixed velocity of zero.
For the temperature field, sometimes we want to enforce a specific non-zero temperature in the solid.
This can be achieved simply by prescribing

$$
T_k = T_b + \frac{T_{k+1} - T_b}{\delta_1 + \delta_2} \delta_2 = \frac{\delta_2}{\delta_1 + \delta_2} T_{k+1} + \left( 1 - \frac{\delta_2}{\delta_1 + \delta_2} \right) T_b
$$

where $T_b$ is the fixed temperature enforced in the solid boundaries.
Compared to the above case of zero velocity, all that is required is and addition of $T_b$ to the $RHS^*$ in the solid points, and an addition of $(1 - \delta_2/(\delta_1 + \delta_2)) T_b$ to the $RHS^*$ for the boundary points.

### Concentration field: No-flux conditions and multi-resolution considerations

Imposing zero-flux conditions at the immersed boundary is a bit more troublesome than enforcing a fixed value.
Due to the complex nature of the three-dimensional semi-implicit time stepping on a staggered grid, we enforce zero gradients in *every* direction for points across the boundaries rather than only the normal gradient.
In the wall-normal direction, we apply a similar approach to that above.

A key difference here, is that the boundary points denoted `ibmask=±1` are now the closest points to the boundary *inside* the solid.
Consider the point $x_k$ being one of these such boundary points.
For the concentration field $C$, we then assume that $\partial_x C=0$ across the boundary by setting $C_{k}$ equal to $C_{k+1}$:

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

Since the concentration field inside the solid has no physical meaning, we can do whatever is most convenient numerically for the points inside the solid (away from the boundary).
This has no effect on the surrounding fluid region, since we have prescribed all the boundary points.
In the current implementation, we choose to allow diffusion of $C$ inside the solid, such that the global solution is relatively smooth, which leads to a more reliable interpolation of $C$ onto the coarse velocity grid.