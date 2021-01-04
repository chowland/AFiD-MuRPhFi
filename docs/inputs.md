# Input file overview

### bou.in

`bou.in` provides a simple way to define various parameters for running the simulation.
These are detailed below.

#### Grid size, time step algorithm, and initial condition
- `NX` and `NXR` are the number of computational cells in the wall-normal direction for the coarse velocity grid and the refined scalar grid respectively.
- `NY(R)` and `NZ(R)` are the number of computational cells in each of the periodic directions.
    - _Since the code is pencil-parallelized in the periodic directions, `NYM(R)` and `NZM(R)` should be divisible by the number of decompositions in that direction.
    These variables should also be products of powers of low-value primes (e.g. $`N=2^n 3^m`$) to speed up the FFT used to calculate the pressure field.
    The variables `NX`, `NY`, and `NZ` are also used throughout the source code, defined as `NX=NXM+1` etc._
- `NSST` prescribes the time-stepping algorithm:
    - `NSST = 3` uses a third-order Runge-Kutta integrator (recommended)
    - `NSST = 1` uses an Adams-Bashforth integrator
- `NREAD` flags whether we are continuing an old simulation
    - `NREAD = 1` starts the simulation by reading the flow field from existing `outputdir/continua*` files
    - `NREAD = 0` starts the simulation from the initial conditions prescribed in `src/CreateInitialConditions.F90`

#### Time stepping details
- `NTST` is the maximum number of time steps for the simulation
- `WALLTIMEMAX` is the maximum wall time for the simulation
- `TOUT` prescribes the simulation time interval between successive statistics calculations
- `IRESET` flags whether to reset time in the simulation to zero or not
    - `IRESET = 0` leaves the simulation unchanged if reading from `continua` files
- `tscaleT` determines which time scale to non-dimensionalize the equations by:
    - `tscaleT = 1` scales the governing equations using the free-fall velocity for the temperature $`U_T = \sqrt{g\alpha \Delta T H}`$, and the corresponding inertial time scale $`[t] = H/U`$
    - `tscaleT = 0` scales the equations using the free-fall velocity associated with the salinity field $`U_S = \sqrt{g\beta \Delta S H}`$ and its corresponding inertial time scale.

#### Domain size and grid stretching parameters
- `ALX3D` is the dimensionless length of the wall-bounded direction. It is strongly recommended to always set this to `1.0`.
- `ISTR3` sets the type of grid stretching applied in the wall-bounded direction. `ISTR3R` does the same for the refined scalar grid.
    - `ISTR3 = 0` produces uniform spacing (no grid stretching)
    - `ISTR3 = 4` uses hyperbolic tangent-type clustering
    - `ISTR3 = 6` uses clipped Chebychev-type clustering
- `STR3` sets the clipping parameter used when `ISTR3 > 0`. This is applied to both the base grid and the refined grid.
- `YLEN` and `ZLEN` are the dimensionless lengths of each periodic direction.

#### Physical and stability parameters
- `RAYT` and `RAYS` are the Rayleigh numbers (e.g. $`Ra_T = g\alpha \Delta T H^3 /\nu \kappa_T`$) for temperature and salinity respectively.
    - Note: positive values for these parameters impose mean gradients $`T_x < 0, \ S_x > 0`$ through the boundary conditions. The sign of these gradients can be reversed by setting negative values here.
- `PRAT` and `PRAS` prescribe the Prandtl (or Schmidt) numbers for the two scalars (e.g. $`Pr_T = \nu/\kappa_T`$).
- `DT` sets the size of the first time step used in the simulation. If variable timestepping is not used (i.e. if `IDTV = 0`), then `DT` remains the fixed time step for the whole simulation
- `RESID` prescribes the maximum allowed residual velocity divergence for mass conservation.
- `CFLMAX` sets the maximum CFL number for the simulation.

#### Statistics flags
- `STATON` flags whether to calculate statistics during the simulation
- `BALANCEON` flags whether to calculate the Nusselt number through global balance equations relating dissipation and heat transport
- `TSTA` sets the time to start the statistics
- `STAREAD` flags whether to read in statistics from `continua` files

#### Boundary conditions and gravity
- `inslws` and `inslwn` prescribe the velocity boundary conditions on the lower (`s`) and upper (`n`) walls respectively:
    - `inslw(s/n) = 1` imposes no-slip $`v = 0, w = 0`$ at the wall.
    - `inslw(s/n) = 0` imposes free-slip $`\partial_x v = 0, \partial_x w = 0`$ at the wall.
- `Tfix(S/N)` and `Sfix(S/N)` prescribe the boundary conditions for temperature and salinity at each wall:
    - `(T/S)fix(S/N) = 1` imposes fixed values, e.g. $`T = \pm 0.5`$
    - `(T/S)fix(S/N) = 0` imposes no flux, e.g. $`\pd_x T = 0`$.
- `gAxis` sets the direction of the vertical axis:
    - e.g. if `gAxis = 1` then gravity acts in the negative $`x`$ direction.

#### Time stepping details
- `IDTV` flags whether or not to use variable time stepping
- `DTMIN` sets the minimum time step for variable time stepping
- `DTMAX` sets the maximum time step for variable time stepping
- `VLIM` sets the maximum allowable local velocity in the simulation

#### Slab output (outdated)
- `SLABDUMP(STST3)` flags whether to output wall-normal slabs at locations specified in stst3.in. If this option is active, you need to create a directory named `stst3`.

#### Shear, pressure gradient, and movie writing
- `xplusU` and `xminusU` specify the velocities (in the $`y`$-direction) of the upper and lower walls respectively.
- `dPdz` sets the mean pressure gradient in the $`z`$ direction.
- `TFRAME` sets the simulation time interval between saving 2D slices for generating movies.