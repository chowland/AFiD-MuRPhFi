# Input file overview

At present only one input file is needed to run a simulation using AFiD-MuRPhFi.
For historical reasons unknown to me, this file is named `bou.in`.
The only other (optional) input file at the moment is `spectra.in` which relates to the calculation of power spectra.
If `spectra.in` is absent from the simulation directory then no power spectra are computed.

### bou.in

`bou.in` provides a simple way to define various parameters for running the simulation.
These are detailed below.

#### Grid size, time step algorithm, and multi-resolution flags
- `NXM` and `NXMR` are the number of computational cells in the wall-normal direction for the coarse velocity grid and the refined scalar grid respectively.
- `NYM(R)` and `NZM(R)` are the number of computational cells in each of the periodic directions for the refined grid.
    - *Since the code is pencil-parallelized in the periodic directions, `NYM(R)` and `NZM(R)` should be divisible by the number of decompositions in that direction.
    These variables should also be products of powers of low-value primes (e.g. $N=2^n 3^m$) to speed up the FFT used to calculate the pressure field.*
    - *The refined variables (e.g. `NXMR`) should satisfy this condition even when `MULTIRES = 0` since the initialization routine for the domain decomposition is currently run for both grids regardless of `MULTIRES`.*
    - *The variables `NX`, `NY`, and `NZ` are also used throughout the source code, defined as `NX = NXM + 1` etc.*
- `NSST` prescribes the time-stepping algorithm:
    - `NSST = 3` uses a third-order Runge-Kutta integrator (recommended)
    - `NSST = 1` uses an Adams-Bashforth integrator
- `MULTIRES` is a logical (set to `0` or `1`) that determines whether or not the second, refined grid is to be used.  
**Note**: *`MULTIRES = 1` is required if one wants to simulate the second scalar (salinity) or use the phase-field method.*
- `FLAGSAL` is a logical that determines whether or not to evolve a second scalar variable, namely salinity.
- `FLAGPF` is a logical that determines whether or not to evolve the phase-field variable used to model melting solid objects.

#### Simulation initialization flags
- `NREAD` flags whether we are continuing an old simulation
    - `NREAD = 1` starts the simulation by reading the flow field from existing `outputdir/continua*` files
    - `NREAD = 0` starts the simulation from the initial conditions prescribed in `src/CreateInitialConditions.F90`
- `IRESET` flags whether to reset time in the simulation to zero or not
    - `IRESET = 0` leaves the simulation unchanged if reading from `continua` files. 
    - The value of `IRESET` has no effect when `NREAD = 0`.

#### Simulation time and I/O parameters
- `NTST` is the maximum number of time steps for the simulation and must be an integer
- `WALLTIMEMAX` is the maximum wall time for the simulation in seconds and must be an integer
- `TMAX` is the time limit for the simulation in flow units (typically the number of free-fall time scale)
- `TOUT` prescribes the simulation time interval between successive statistics calculations saved to `means.h5`, and if in use `spectra.h5`.
- `TFRAME` sets the simulation time interval between saving 2D slices for generating movies in `outputdir/flowmov`.
- `SAVE_3D` sets the equivalent time interval between successive saving of the full 3D flow fields to the directory `outputdir/fields`.

#### Domain size and grid stretching parameters
- `ALX3D` is the dimensionless length of the wall-bounded direction. It is strongly recommended to always set this to `1.0`.
- `YLEN` and `ZLEN` are the dimensionless lengths of each periodic direction.
- `ISTR3` sets the type of grid stretching applied in the wall-bounded direction. `ISTR3R` does the same for the refined scalar grid.
    - `ISTR3 = 0` produces uniform spacing (no grid stretching)
    - `ISTR3 = 4` uses hyperbolic tangent-type clustering
    - `ISTR3 = 6` uses clipped Chebychev-type clustering
- `STR3` sets the clipping parameter used when `ISTR3 > 0`. This is applied to both the base grid and the refined grid.

#### Physical parameters
- `RAYT` and `RAYS` are the Rayleigh numbers (e.g. $Ra_T = g\alpha \Delta T H^3 /\nu \kappa_T$) for temperature and salinity respectively.
    - Note: positive values for these parameters impose mean gradients $T_x < 0, \ S_x > 0$ through the boundary conditions. The sign of these gradients can be reversed by setting negative values here.
- `PRAT` and `PRAS` prescribe the Prandtl (or Schmidt) numbers for the two scalars (e.g. $Pr_T = \nu/\kappa_T$).
- `FFscaleS` determines which time scale to non-dimensionalize the equations by:
    - `FFscaleS = 0` scales the governing equations using the free-fall velocity for the temperature $U_T = \sqrt{g\alpha \Delta T H}$, and the corresponding inertial time scale $[t] = H/U$
    - `FFscaleS = 1` scales the equations using the free-fall velocity associated with the salinity field $U_S = \sqrt{g\beta \Delta S H}$ and its corresponding inertial time scale.  
    **Note:** *if `FLAGSAL = 0` then this option is ignored, and the temperature free-fall scale is used.*

#### Time stepping parameters
- `IDTV` flags whether or not to use variable time stepping
- `DT` sets the size of the first time step used in the simulation. If variable timestepping is not used (i.e. if `IDTV = 0`), then `DT` remains the fixed time step for the whole simulation
- `RESID` prescribes the maximum allowed residual velocity divergence for mass conservation.
- `CFLMAX` sets the maximum CFL number for the simulation.  
**Note:** *if `MULTIRES = 1` then the CFL condition uses the grid spacing from the refined grid.*
- `DTMIN` sets the minimum time step for variable time stepping
- `DTMAX` sets the maximum time step for variable time stepping

#### Boundary conditions and gravity
- `inslwS` and `inslwN` prescribe the velocity boundary conditions on the lower (`S`) and upper (`N`) walls respectively:
    - `inslw(s/n) = 1` imposes no-slip $v = 0, w = 0$ at the wall.
    - `inslw(s/n) = 0` imposes free-slip $\partial_x v = 0, \partial_x w = 0$ at the wall.
- `Tfix(S/N)` and `Sfix(S/N)` prescribe the boundary conditions for temperature and salinity at each wall:
    - `(T/S)fix(S/N) = 1` imposes fixed values, e.g. $T = \pm 0.5$. These values are set in the subroutines `SetTempBCs` and `SetSalBCs` respectively.
    - `(T/S)fix(S/N) = 0` imposes no flux, e.g. $\partial_x T = 0$.
- `active_T(S)` is a flag determining whether each scalar is active or passive. That is, if `active_T = 1` then the buoyancy effect of temperature is added to the momentum equation.
- `gAxis` sets the direction of the vertical axis:
    - e.g. if `gAxis = 1` then gravity acts in the negative $x$ direction.

#### Shear, pressure gradient, and melting condition
- `xplusU` and `xminusU` specify the velocities (in the $y$-direction) of the upper and lower walls respectively.
- `dPdy` and `dPdz` set mean pressure gradients in the $y$ and $z$ directions respectively.
    - `dPdy` is interpreted as a target shear Reynolds number $Re_\tau = v_\tau H/2\nu$ where $v_\tau$ is the friction velocity in the $y$-direction. The pressure gradient is kept to a constant value such that in a statistically steady state $Re_\tau$ measured at the wall will match the target, i.e. that the dimensionless pressure gradient in the equations satsifies $G=8 Re_\tau^2/Gr$. This constant pressure gradient is simply added to the $y$-momentum equation in `ExplicitTermsVY`.
    - `dPdz` is interpreted as a target bulk Reynolds number $Re_b=\langle w \rangle H/\nu$ where $\langle w\rangle$ is the volume-averaged $z$-component of velocity. The pressure gradient is adjusted at each time step to maintain this constant volume flux. We follow [Quadrio et al. (2016)](https://doi.org/10.1016/j.euromechflu.2015.09.005) in setting the wall-normal profile of the pressure gradient such that it remains consistent with the boundary conditions. This adjustment is performed in the subroutine `CorrectVelocity`.
- `MELT` is a logical that implements a dynamic boundary condition for temperature and salinity according to the three-equation model for a melting boundary. With `MELT = 1`, the boundary at $x=0$ is treated as a stationary, planar wall of ice. **This feature is highly experimental and is not recommended for current use.**

#### Phase-field parameters
As described in the documentation page on the [governing equations](equations.md#double-diffusive-melting), the evolution for the phase-field equation is

$$
\partial_t \phi = A\nabla^2 \phi - (A/\varepsilon^2) \phi (1-\phi)(1-2\phi+C(T-T_m + \Lambda S))
$$

The diffuse interface thickness $\varepsilon$ is always set equal to the mean wall-normal grid spacing for the refined grid $\varepsilon = 1/n_x^r$.

- `pf_A`$=A Pe$ is value of the diffusivity coefficient in the phase-field equation as a multiple of the Peclet number $Pe=(RaPr)^{-1/2}$;
- `pf_C`$=C$ is the value of the thermal coupling parameter;
- `pf_S`$=\mathcal{S}=L/c_p\Delta T$ is the Stefan number;
- `pf_Tm`$=T_m$ is the dimensionless melting temperature;
- `IBM` is a logical that determines which method is used to constrain the velocity field in the solid:
    - `IBM = 0` uses the volume penalty method, adding a forcing term of $-{Pe}^{-1} \phi \boldsymbol{u}/(1.5\varepsilon)^2$ to the momentum equation;
    - `IBM = 1` imposes an immersed boundary in the solid phase by simply setting $\boldsymbol{u} = 0$ wherever $\phi>0.5$ at each time step;
- `pf_IC` determines the initial condition when the phase-field method is used:
    - `pf_IC = 1` provides the initial condition for the 1-D melting/freezing example described [here](examples/stefan.md#1-d-solidification-from-a-cooled-boundary);
    - `pf_IC = 2` produces an initial condition with a solid disc at the centre of the domain, following the example described [here](examples/stefan.md#axisymmetric-melting-of-a-solid-disc-in-2-d);
    - `pf_IC = 3` provides an initial condition to match the validation example of Rayleigh-Bénard convection from Favier et al. (2019) as described [here](examples/coupled_flows.md#2-d-rayleigh-benard-with-a-melting-boundary). The velocity is zero everywhere, the phase boundary is at $x=0.5$, and the temperature profile is linear plus a wave-like perturbation in the liquid phase.