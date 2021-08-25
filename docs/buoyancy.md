# Buoyancy - equation of state

**Note:** Only the linear equation of state is currently implemented in AFiD.

## Boussinesq approximation

## Linear equation of state

The buoyancy term is commonly related linearly to temperature and salinity with constant expansion/contraction coefficients $\alpha$ and $\beta$:

$$
b = -\frac{g\rho'}{\rho_0} = g \left( \alpha T - \beta S\right) .
$$

This definition of buoyancy has units of acceleration, and its dimensionless form in AFiD is scaled by $U^2/H$, where $U$ is the free-fall velocity scale for either temperature or salinity.

The density ratio $R_\rho$ is determined by the input parameters `RAYT`, `RAYS`, `PRAT`, and `PRAS` through

$$
R_\rho = \frac{Ra_T Pr_S}{Ra_S Pr_T} .
$$

In reality, the thermal expansion coefficient $\alpha=\rho_0^{-1} \partial \rho/\partial T$ depends quite significantly on temperature.
For example in fresh water (with no salinity component), there is a density maximum at $4°\mathrm{C}$, below which $\alpha$ takes negative values.
In salt water with a typical oceanic salinity of 35psu, the value of $\alpha$ varies from approximately $10^{-5}\ \mathrm{K}^{-1}$ at the liquidus temperature of $-2°\mathrm{C}$ to $\alpha\approx 10^{-4}\ \mathrm{K}^{-1}$ at $20°\mathrm{C}$.
The linear assumption is therefore only reasonable when considering small temperature differences or temperatures far from the density maximum.
We can include the effects of the density maximum, and its dependence on salinity, by considering a nonlinear equation of state as described below.

## "Realistic" seawater equation of state

We follow the work of [Roquet et al. (2015)](https://doi.org/10.1175/JPO-D-15-0080.1) in defining a simple, polynomial expression to approximate the [TEOS-10](http://www.teos-10.org) International Thermodynamic Equation of Seawater.
Equation (18) of that paper defines the following form for density perturbations due to temperature and salinity:

$$
\rho' = \frac{-C_b}{2} \left( \Theta - \Theta_0 - \varepsilon S_A \right)^2 - T_h Z \Theta + b_0 S_A ,
$$

where $\Theta$ is conservative temperature, $S_A$ is absolute salinity, $Z$ is geopotential depth, and the coefficients have values

$$
C_b = 0.011 \mathrm{kg\ m^{-3}\ K^{-2}}, \quad
T_h = 2.5\times 10^{-5} \mathrm{kg\ m^{-4}\ K^{-1}}, \quad
b_0 = 0.77 \mathrm{kg\ m^{-3}\ (g\ kg^{-1})^{-1}} \\
\Theta_0 = 4°\mathrm{C} \quad
\varepsilon = -0.25 \mathrm{K\ (g\ kg^{-1})^{-1}} .
$$

When performing DNS, we can neglect the second term in the equation of state that quantifies variations in the temperature dependence with depth (due to pressure changes).
Non-dimensionalizing the temperature and salinity by $\Delta T$ and $\Delta S$ respectively then gives the dimensionless buoyancy as

$$
R_\rho^\ast \left( T - T_0 + \varepsilon^\ast S \right)^2 - S
$$

if velocities are scaled by the solutal free-fall scale $U_S = \sqrt{g\beta \Delta S H}$, where the modified density ratio $R_\rho^\ast$ and the density maximum coefficient $\varepsilon^\ast$ are defined as

$$
R_\rho^\ast = \frac{C_b (\Delta T)^2}{2b_0 \Delta S}, \quad \varepsilon^\ast = \frac{\varepsilon \Delta S}{\Delta T} .
$$

This definition of the density ratio also motivates an alternative Rayleigh number, based on the scale of the thermal component of the buoyancy term:

$$
\mathcal{R}_T = \frac{g C_b (\Delta T)^2 H^3}{\rho_0 \nu \kappa_T} .
$$

This can also be used to define a new free-fall velocity based on temperature as

$$
U_T = \sqrt{\frac{gC_b H}{2\rho_0}} \Delta T .
$$

In the case where this scale is used to non-dimensionalize the velocity, the buoyancy term will be given by

$$
\left( T - T_0 + \varepsilon^\ast S \right)^2 - {R_\rho^\ast}^{-1} S.
$$