# Moist convection module - Rainy-Bénard

The `afid_moisture` module provides additional features needed to simulate a simple model of moist convection.
This follows the setup of Vallis et al (2019) which is outlined below.

## Quick start: `humid.in` input file

If you already know the setup of the Rainy-Bénard model, and just want to get started running simulations, you need to use/modify the `humid.in` file in `examples/RainyBenard`, which defines the following parameters (given by their internal names in the code)

- `alpha_q`: $\alpha = L\Delta T/R_v T_0^2$ the dimensionless parameter determining the temperature dependence of the saturation humidity
- `beta_q`: $\beta = gH/c_p \Delta T$ the dimensionless lapse rate parameter
- `gamma_q`: $\gamma = L q_0/c_p \Delta T$ the dimensionless parameter determining the strength of heating due to condensation
    - Setting `gamma_q < 0` overrides the input value by setting $\gamma = \beta \Delta q$, where $\Delta q$ is the moisture difference between the top and bottom boundaries. This ensures that 
- `tau_q`: $\tau$ the condensation timescale used in the fast relaxation term
    - `ifdiff==.true.` implies that $\tau$ is defined in diffusive time units $H^2/\kappa$ (e.g. to match Vallis et al)
    - `ifdiff==.false.` implies that $\tau$ is defined in buoyancy-scaled time units $H/\sqrt{g\Delta T H/T_0}$
    - **N.B.** for stability, we must enforce that the time step $\Delta t < 0.1 \tau$. If $\tau$ is sufficiently small, then `dtmax` will be overwritten to enforce this constraint.
- `qfixN` sets the top boundary condition:
    - `qfixN==0` sets a no-flux condition $\partial_z q=0$
    - `qfixN==1` sets a fixed, saturated value at the top $q=q_s=1-e^{-\alpha}$
    - `qfixN==2` sets a fixed, dry value at the top $q=0$
- `qfixS` sets the bottom boundary condition in the same way, and can be set as saturated (`1`) or as no flux (`0`).
- `Sm`: $S_m = \kappa_q/\kappa$ the diffusivity ratio of moisture and heat.

The model for moist convection will only be run if both `bou.in` and `humid.in` are in the current directory.
The logical variable `moist` in the code determines which additional subroutines to call when running a simulation of moist convection.

## Governing equations and variables

### Buoyancy and temperature

In the Rainy-Bénard model, we solve buoyancy $b$ in place of temperature, where buoyancy is defined as

$$
b =-\frac{g(\rho-\rho_0)}{\rho_0} =\frac{g}{T_0} \left(T - T_0  + \frac{g}{c_p} z \right).
$$

The final term in the above expression represents the dry lapse rate, i.e. the rate at which air of the same density decreases in temperature with height due to the changes in pressure.
The Rainy-Bénard setup, like classical Rayleigh-Bénard, enforces fixed values of temperature at the top and bottom of the domain, differing by $\Delta T$.
Buoyancy can therefore be non-dimensionalised using $g\Delta T/T_0$, and satisfies the equation

$$
D_t b^* = (RaPr)^{1/2} \nabla^2 b^* + Q_c,
$$

where $Q_c$ represents the heating due to condensation of water vapour.
In the code, we re-use the variable `temp` for $b^*$ and add the condensation term (defined below) in the subroutine `AddCondensation`.

In dimensionless form, the boundary conditions on temperature are

$$
\theta = 0, \ z=0, \qquad \theta = -1, \ z = 1,
$$

leading to the following boundary conditions for the dimensionless buoyancy

$$
b = 0, \ z = 0, \qquad b = \beta - 1, \ z=1 .
$$

Here, the dimensionless parameter $\beta$, which can be set by the input file `humid.in`, is defined as

$$
\beta = \frac{gH}{c_p \Delta T} .
$$

### Specific humidity
The additional scalar field introduced by the Rainy-Bénard model is the specific humidity $q$.
This field also satisfies an advection-diffusion equation with an additional term to represent the effect of condensation.
The condensation term takes the form of a fast relaxation term that damps humidity above saturation back to the saturated value.

The saturation humidity is defined as

$$
q_s = q_0\exp [\alpha \theta] = q_0 \exp [\alpha (b^* - \beta z)], \qquad \alpha = \frac{L\Delta T}{R_v T_0^2}.
$$

Here, an additional dimensionless control parameter $\alpha$ has been introduced, that can be set by the input file `humid.in`.

In the solver, specific humidity is scaled by $q_0$, such that in dimensionless form a saturated bottom boundary takes the value $q=1$.
The governing equation for humidity is as follows

$$
D_t q = \frac{S_m}{\sqrt{Ra Pr}} \nabla^2 q + \frac{q - q_s}{\tau} \mathcal{H}(q - q_s) ,
$$

where $S_m=\kappa_q/\kappa$ accounts for the fact that moisture and heat may not diffuse at the same rates.
Both $S_m$ and the condensation time-scale $\tau$ can be set in the `humid.in` input file.
This equation is solved in an identical manner to the temperature equation.

Finally, the heating due to condensation can be expressed by invoking the first law of thermodynamics, that $L \Delta q = c_p \Delta T$.
In dimensionless form, this leads to the expression

$$
Q_c = \gamma \frac{q - q_s}{\tau} \mathcal{H}(q - q_s), \qquad \gamma = \frac{L q_0}{c_p \Delta T} .
$$

As with the other dimensionless control parameters, $\gamma$ can be set in `humid.in`.

## Additional notes

### The drizzle solution

The Rainy-Bénard system has a static solution known as the 'drizzle' solution where the moist static energy $m=b+\gamma q$ is linear.
We can optionally start a simulation from this state (with small-amplitude random noise) by providing an input `drizzle.h5` in the run directory.
This contains the individual profiles for $b(z)$ and $q(z)$, and can be created using the python script `tools/drizzle.py`.

### Statistics written out

The following profiles are added to `outputdir/means.h5` if the moisture model is running:

- `qbar`: $\overline{q}(z,t)$ horizontally-averaged specific humidity
- `qrms`: $\sqrt{\overline{q^2}}(z,t)$ root-mean-squared humidity
- `qsat`: $\overline{q_s}(z,t)$ saturation humidity
- `qrel`: $\overline{(q/q_s)}(z,t)$ *relative* humidity
- `vxq`: $\overline{wq}(z,t)$ vertical advective moisture flux
- `vyq`: $\overline{uq}(z,t)$ horizontal advective moisture flux (in $x$)
- `vzq`: $\overline{vq}(z,t)$ horizontal advective moisture flux (in $y$)

Cuts of $q$ are also written out to the `movie_Xcut.h5` files under the name `qhum`.