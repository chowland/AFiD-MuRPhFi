# Hermite interpolation between grids

To ensure efficient and accurate interpolation between the multiple grids used in our simulations, we implement Hermite interpolation.
As handily described on [this Wikipedia page](https://en.wikipedia.org/wiki/Cubic_Hermite_spline), Hermite interpolation uses the value and gradient of a function at two points to calculate the value of the function at any intermediate point.

Say we wish to interpolate a function $f(x)$ from the coarse grid to the refined grid.
At each point on the refined grid $x^r$, we find the index $k$ such that $x_k < x^r < x_{k+1}$.
The Hermite interpolation is then given by

$$
f(x^r) = h_{00}(t)f(x_k) + \Delta_k h_{10}(t)f'(x_k) + h_{01}(t)f(x_{k+1}) + \Delta_k h_{11}(t)f'(x_{k+1}),
$$

where $t = (x^r-x_k)/(x_{k+1}-x_k)$, $\Delta_k = x_{k+1} - x_k$ and $h_{ij}$ are the polynomial basis functions

$$
h_{00}(t) = (1+2t)(1-t)^2, \quad h_{10}(t) = t(1-t)^2,
$$

$$
h_{01}(t) = t^2(3-2t), \quad h_{11}(t) = t^2(t-1) .
$$

We compute the derivatives using a second order finite difference approximation, writing $f_k = f(x_k)$ and so on:

$$
f'(x_k) = \frac{
    \left[\delta_-^2 \left(f_{k+1} - f_k\right) + \delta_+^2 \left(f_k - f_{k-1}\right)\right]
}{\delta_+ \delta_-(\delta_+ + \delta_-)}
$$

where $\delta_+ = x_{k+1} - x_k$ and $\delta_- = x_k - x_{k-1}$ are the grid spacings either side of $x_k$.
For each point on the refined grid, we therefore need a four-point stencil containing the values $\{f_{k-1}, f_k, f_{k+1}, f_{k+2}\}$.

Exactly the same technique can be applied when interpolating from the refined grid to the coarse grid, such as for adding salinity to the buoyancy term.
The module `HermiteInterpolations` contains the relevant subroutines for constructing the stencils.
These subroutines are also used for the interpolation of the initial condition, and work regardless of whether the resolution is being upscaled or downscaled.

## Boundary points
When the point we are interpolating to is between one of the solid boundaries and the adjacent point, we do not have sufficient information to perform such a detailed interpolation.
In this case, we use a simple linear interpolation between the boundary and the first grid point.
This will not significantly impact the accuracy of the simulation if the fields are adequately resolved.

## Minimising the velocity divergence

In order to preserve the divergence-free nature of the velocity field when interpolating it to the refined grid, we actually interpolate the velocity gradients to the refined grid, and integrate over space to reconstruct the velocity components on the refined grid.

For example, for the wall normal velocity `vx`, we compute

$$
v_x(x,y,z) = \int_0^x \frac{\partial v_x}{\partial x}(x',y,z) \ \mathrm{d}x' ,
$$

on the refined grid, having interpolated $\partial_x v_x$ to higher resolution.

For the components of velocity in the other (homogeneous) directions, we have to first interpolate a single plane of the velocity field and then integrate over the desired direction to recover, for example

$$
v_y(x,y,z) = v_y(x,0,z) + \int_0^y \frac{\partial v_y}{\partial y} (x,y',z) \ \mathrm{d}y' .
$$

We choose to use the interpolated plane from $y=0$ here since there is no interpolation in $y$ necessary: the first index of both grids (`yc` and `ycr`) for $v_y$ corresponds to this plane.

The velocity interpolation is performed in the subroutine [`InterpVelMgrd`](https://github.com/chowland/AFiD-MuRPhFi/blob/main/src/multires/InterpVelMgrd.F90).