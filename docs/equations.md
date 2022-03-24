# Governing equations

## Dimensional system
The velocity field $\boldsymbol{u}$ is determined by the Navier-Stokes equation subject to the Oberbeck-Boussinesq (OB) approximation such that density changes are related linearly to temperature changes and are negligible everywhere except the buoyancy term:

$$
\partial _t \boldsymbol{u} + \boldsymbol{\nabla} \cdot (\boldsymbol{u} \boldsymbol{u}) = -\rho_0^{-1} \boldsymbol{\nabla} p + g\alpha T \mathbf{\hat{e}_g} + \nu \nabla^2 \boldsymbol{u} .
$$

The pressure $p$ is determined to make the velocity field divergence-free ($\boldsymbol{\nabla}\cdot\boldsymbol{u}=0$), consistent with the OB approximation.
The nonlinear term $\boldsymbol{\nabla} \cdot (\boldsymbol{u} \boldsymbol{u}) = \partial_j(u_j u_i) \mathbf{e_i}$ is computed in conservative form to improve energy conservation properties.
The gravitational axis $\mathbf{\hat{e}_g}$ can be specified in the input file `bou.in` using the variable `gAxis`.
The temperature field $T$ satisfies an advection-diffusion equation:

$$
\partial_t T + \boldsymbol{\nabla} \cdot (\boldsymbol{u} T) = \kappa_T \nabla^2 T .
$$

Parameters featured in the above equations are defined in the table below.

Parameter    | Definition
:----------: | :-----------------------------------------
$\rho_0$   | (Constant) Mean fluid density
$g$        | Gravitational acceleration
$\alpha$   | Thermal expansion coefficient
$\nu$      | Kinematic viscosity / momentum diffusivity
$\kappa_T$ | Thermal diffusivity


## Dimensionless system
AFiD solves these equations in dimensionless form.
This dimensionless system is obtained by scaling the variables by the following:

$$
[\boldsymbol{x}] = H, \quad [T] = \Delta T, \quad [\boldsymbol{u}] = U_T = \sqrt{g\alpha H\Delta T}, \quad [t] = H/U_T, \quad [p] = \rho_0 U_T^2 .
$$

Here, $H$ is the distance between the two walls of the domain, $\Delta T$ is the temperature difference between the walls, and $U_T$ is known as the free-fall velocity scale, obtained by assuming a leading-order balance between the buoyancy and inertial terms.

After rescaling, the equations read

$$
\partial_t \boldsymbol{u} + \boldsymbol{\nabla} \cdot (\boldsymbol{u} \boldsymbol{u}) = - \boldsymbol{\nabla} p + T \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + \boldsymbol{\nabla} \cdot (\boldsymbol{u} T) = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T .
$$

The dimensionless velocity and temperature fields are therefore solely determined by the Rayleigh number $Ra_T$ and the Prandtl number $Pr$, which act as the control parameters of the system and are defined as

$$
Ra_T = \frac{g\alpha H^3 \Delta T}{\nu \kappa_T}, \quad Pr = \frac{\nu}{\kappa_T} .
$$

These parameters are prescribed by the user in the input file `bou.in` as the variables `RAYT` and `PRAT`.
AFiD can also treat $T$ as a passive scalar if the variable `active_T` is set to `0` in the input file.
In this case, the equations read the same as above, except that the buoyancy term $T\mathbf{\hat{e}_g}$ is set to zero.

## Salinity as a second scalar

AFiD-MuRPhFi also allows for the evolution of a second scalar field on a grid of finer resolution.
This option can be enabled by setting `MULTIRES` and `FLAGSAL` equal to 1 in the input file.
Using the refined grid is useful for scalars with low diffusivity, such as salinity, so we proceed by labelling this scalar $C$ for concentration.
The second scalar satisfies an advection-diffusion equation, but with a different diffusivity $\kappa_C$:

$$
\partial_t C + \boldsymbol{\nabla} \cdot (\boldsymbol{u} C) = \kappa_C \nabla^2 C .
$$

If $C$ is not a passive scalar, then it will also contribute towards the buoyancy of the system.
For temperature and salinity, we can assume a linear equation of state such that the (dimensional) buoyancy term becomes $g(\alpha T - \beta C)\mathbf{\hat{e}_g}$, where $\beta$ is the haline contraction coefficient.

We now have two choices for the velocity scale used to non-dimensionalise the system, depending on whether the dominant contribution to the buoyancy is from temperature or salinity.
The relevant scale can be prescribed in the input file by the variable `FFscaleS`.

If salinity is the dominant scalar in the buoyancy, we set `FFscaleS = 1` and take the free-fall scale based on $C$, that is $[\boldsymbol{u}] = U_C = \sqrt{g\beta H \Delta C}$.
Analogous to the case above, $\Delta C$ is the difference in salinity between the walls.
The dimensionless equations now read

$$
\partial_t \boldsymbol{u} + \boldsymbol{\nabla} \cdot (\boldsymbol{u} \boldsymbol{u}) = - \boldsymbol{\nabla} p + (R_\rho T - C) \mathbf{\hat{e}_g} + \sqrt{\frac{Sc}{Ra_C}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + \boldsymbol{\nabla} \cdot (\boldsymbol{u} T) = \frac{1}{\tau\sqrt{Ra_C Sc}} \nabla^2 T ,
$$

$$
\partial_t C + \boldsymbol{\nabla} \cdot (\boldsymbol{u} C) = \frac{1}{\sqrt{Ra_C Sc}} \nabla^2 C .
$$

The solutal Rayleigh number $Ra_C$ and Schmidt number $Sc$ can be defined in the input file through the variables `RAYS` and `PRAS` and are defined similarly to the thermal Rayleigh number and Prandtl number as

$$
Ra_C = \frac{g\beta H^3 \Delta C}{\nu \kappa_C}, \quad Sc = \frac{\nu}{\kappa_C} .
$$

Two new dimensionless parameters also appear in the above equations: the density ratio $R_\rho$ and the diffusivity ratio $\tau$.
These parameters are defined as follows and are uniquely determined by prescribing ($Ra_T$, $Pr$, $Ra_C$, $Sc$).

$$
R_\rho = \frac{\alpha \Delta T}{\beta \Delta C} = \frac{Ra_T Sc}{Ra_C Pr}, \quad \tau = \frac{\kappa_C}{\kappa_T} = \frac{Pr}{Sc} .
$$

If temperature is the dominant scalar, we set `FFscaleS = 0` and take the free-fall scale as $U_T$.
In this case, the dimensionless equations take the form

$$
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + (T - {R_\rho}^{-1} C) \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T ,
$$

$$
\partial_t C + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) C = \frac{\tau}{\sqrt{Ra_T Pr}} \nabla^2 C .
$$

As with the temperature field, we can treat $C$ as a passive scalar by setting `active_S = 0` in the input file.

## Phase-field method
The latest addition to AFiD-MuRPhFi is the implementation of a phase-field method to simulate melting or dissolving solid objects immersed in a fluid flow.
This method evolves an indicator function $\phi$ which takes the values 0/1 in the fluid/solid, approximating the fluid-solid interface by a finite, diffuse region.
We follow the formulation of [Hester et al (2020)](https://doi.org/10.1098/rspa.2020.0508), appropriately scaled to match the dimensionless system detailed above.

To start with, we consider the case with no salinity effects for simplicity.
The (dimensionless) equations for the evolution of the velocity and temperature gain additional term due to volume penalisation in the solid and latent heat respectively, and now take the form

$$
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + T \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\left(\nabla^2 \boldsymbol{u} - \frac{\phi \boldsymbol{u}}{\eta}\right) ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T + \mathcal{S}\partial_t \phi ,
$$

where $\eta = (\beta \varepsilon)^2$ is the volume penalty coefficient for a dimensionless diffusive interface thickness $\varepsilon$ and $\mathcal{S}=L/c_p\Delta T$ is the Stefan number.
$L$ is the latent heat of solidification, and $c_p$ is the specific heat capacity which, like the diffusivity of heat $\kappa_T$, is assumed to be equal in the solid and liquid phases.
Following [Hester et al (2020)](https://doi.org/10.1098/rspa.2020.0508), we take $\beta=1.51044385$ as the optimal coefficient for the smooth volume penalty term.

As an alternative to the volume penalty method, an extremely simple immersed boundary method is also implemented for imposing zero velocity in the solid.
For this implementation, the the velocity is simply reset to zero where $\phi<0.5$.
This restriction is imposed after the implicit step of the time-stepper but before the pressure correction to ensure that the velocity field remains divergence free at all times.
The immersed boundary technique is activated when `IBM` is set to 1 in the input file, otherwise the volume penalty method is used.

In dimensional form, the evolution equation for the phase variable $\phi$ is

$$
\hat{\varepsilon} \frac{5}{6} \frac{L}{c_p\kappa_T} \partial_t \phi - \Gamma \nabla^2 \phi = -\frac{1}{\hat{\varepsilon}^2} \phi (1-\phi)(\Gamma (1-2\phi) + \hat{\varepsilon} (T-T_m)) .
$$

Here $\hat{\varepsilon}=H\varepsilon$ is the dimensional interface thickness, $T_m$ is the equilibrium melting temperature, and $\Gamma$ is the coefficient arising from the Gibbs-Thomson relationship linking the interface temperature to the curvature $K$ of the interface:

$$
T - T_m = \Gamma K = \frac{\gamma T_m}{\rho_s L} K.
$$

$\gamma$ is the surface energy and $\rho_s$ is the density of the solid.

We now non-dimensionalise the phase-field equation using the same scales as before.
This results in the evolution equation

$$
\partial_t \phi = D \nabla^2 \phi - \frac{D}{\varepsilon^2} \phi (1 - \phi) (1 - 2\phi + A(T - T_m)) ,
$$

where the dimensionless coefficients are given by

$$
D = \frac{6}{5\varepsilon} \frac{\kappa}{U_T H} \frac{c_p \Gamma}{LH}, \quad A = \frac{\varepsilon H \Delta T}{\Gamma} .
$$

Note that we can typically set $\varepsilon=1/n_x^r$ to be the mean refined grid spacing.
Furthermore, $D$ can be written as $D=1.2 (Pe_T \mathcal{S} A)^{-1}$, so choosing $A$ also determines $D$.
In this case, we can write down the dimensionless phase-field equation in terms of the physical parameters and a single model parameter $A$:

$$
\frac{\partial\phi}{\partial t} = \frac{6}{5 Pe_T \mathcal{S} A} \left\{ \nabla^2 \phi - \frac{1}{\varepsilon^2} \phi (1 - \phi) \left[1 - 2\phi + A(T - T_m)\right]\right\}.
$$

In the absence of temperature variations, $\phi$ has a steady state solution of

$$
\phi = \frac{1}{2}\left(1 + \tanh \left(\frac{x-x_i}{2\varepsilon}\right)\right) ,
$$

where $x_i$ is the position of the fluid-solid interface.

When running phase-field simulations, the Stefan number $\mathcal{S}$ and the parameter $A$ must be specified in the input file as the variables `pf_S` and `pf_A` respectively.
The dimensionless melting temperature must also be set as the variable `pf_Tm`.

<!-- ####  -->
<!-- 
### Parameter choices
A wide range of parameters have appeared in the phase-field implementation, so it is important to identify what values are relevant in the current setup.
Two examples of parameters used in the literature are detailed below.

Both Hester et al (2021) and Couston et al (2021) use A=1 without further justification.-->

### Double-diffusive melting
Following [Hester et al (2020)](https://doi.org/10.1098/rspa.2020.0508), we can also simulate the ablation of ice in salt water.
In this case, the velocity of the ice-water boundary depends on the gradients of both temperature and salinity at the interface.
In our non-dimensionalization (using the free-fall velocity scale for temperature), the full collection of equations reads as follows:

$$
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + (T - {R_\rho}^{-1} C) \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T + \mathcal{S} \partial_t \phi ,
$$

$$
\partial_t C + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) C = \frac{\tau}{\sqrt{Ra_T Pr}} \left(\nabla^2 C - \frac{\boldsymbol{\nabla} \phi \cdot \boldsymbol{\nabla} C}{1 - \phi + \delta}\right) + \frac{C\partial_t \phi}{1 - \phi + \delta},
$$

$$
\partial_t \phi = D \nabla^2 \phi - \frac{D}{\varepsilon^2} \phi (1 - \phi) (1 - 2\phi + A(T - T_m +\Lambda C)) .
$$

Two new dimensionless parameters have appeared in the equations.
$\delta \ll 1$ is a small parameter that is solely present to stabilise the terms on the right hand side of the salinity equation.
In the phase-field equation, we obtain a new physical parameter defined as

$$
\Lambda = \frac{\lambda \Delta C}{\Delta T} ,
$$

where $\Delta C$ and $\Delta T$ are the previously used scales for salinity and temperature, and $\lambda$ is the liquidus slope.
This refects how the presence of salinity at the ice-water interface lowers the local melting temperature.

<!-- ## Boundary conditions -->