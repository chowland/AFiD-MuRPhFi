# Governing equations

## Dimensional system
The velocity field $\boldsymbol{u}$ is determined by the Navier-Stokes equation subject to the Oberbeck-Boussinesq (OB) approximation such that density changes are related linearly to temperature changes and are negligible everywhere except the buoyancy term:

$$
\partial _t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = -\rho_0^{-1} \boldsymbol{\nabla} p + g\alpha T \mathbf{\hat{e}_g} + \nu \nabla^2 \boldsymbol{u} .
$$

The pressure $p$ is determined to make the velocity field divergence-free ($\boldsymbol{\nabla}\cdot\boldsymbol{u}=0$), consistent with the OB approximation.
The gravitational axis $\mathbf{\hat{e}_g}$ can be specified in the input file `bou.in` using the variable `gAxis`.
The temperature field $T$ satisfies an advection-diffusion equation:

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \kappa_T \nabla^2 T .
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
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + T \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T .
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
Using the refined grid is useful for scalars with low diffusivity, such as salinity, so we proceed by labelling this scalar $S$.
The second scalar satisfies an advection-diffusion equation, but with a different diffusivity $\kappa_S$:

$$
\partial_t S + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) S = \kappa_S \nabla^2 S .
$$

If $S$ is not a passive scalar, then it will also contribute towards the buoyancy of the system.
For temperature and salinity, we can assume a linear equation of state such that the (dimensional) buoyancy term becomes $g(\alpha T - \beta S)\mathbf{\hat{e}_g}$, where $\beta$ is the haline contraction coefficient.

We now have two choices for the velocity scale used to non-dimensionalise the system, depending on whether the dominant contribution to the buoyancy is from temperature or salinity.
The relevant scale can be prescribed in the input file by the variable `FFscaleS`.

If salinity is the dominant scalar in the buoyancy, we set `FFscaleS = 1` and take the free-fall scale based on $S$, that is $[\boldsymbol{u}] = U_S = \sqrt{g\beta H \Delta S}$.
Analogous to the case above, $\Delta S$ is the difference in salinity between the walls.
The dimensionless equations now read

$$
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + (R_\rho T - S) \mathbf{\hat{e}_g} + \sqrt{\frac{Sc}{Ra_S}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\tau\sqrt{Ra_S Sc}} \nabla^2 T ,
$$

$$
\partial_t S + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) S = \frac{1}{\sqrt{Ra_S Sc}} \nabla^2 S .
$$

The solutal Rayleigh number $Ra_S$ and Schmidt number $Sc$ can be defined in the input file through the variables `RAYS` and `PRAS` and are defined similarly to the thermal Rayleigh number and Prandtl number as

$$
Ra_S = \frac{g\beta H^3 \Delta S}{\nu \kappa_S}, \quad Sc = \frac{\nu}{\kappa_T} .
$$

Two new dimensionless parameters also appear in the above equations: the density ratio $R_\rho$ and the diffusivity ratio $\tau$.
These parameters are defined as follows and are uniquely determined by prescribing ($Ra_T$, $Pr$, $Ra_S$, $Sc$).

$$
R_\rho = \frac{\alpha \Delta T}{\beta \Delta S} = \frac{Ra_T Sc}{Ra_S Pr}, \quad \tau = \frac{\kappa_S}{\kappa_T} = \frac{Pr}{Sc} .
$$

If temperature is the dominant scalar, we set `FFscaleS = 0` and take the free-fall scale as $U_T$.
In this case, the dimensionless equations take the form

$$
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + (T - {R_\rho}^{-1} S) \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T ,
$$

$$
\partial_t S + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) S = \frac{\tau}{\sqrt{Ra_T Pr}} \nabla^2 S .
$$

As with the temperature field, we can treat $S$ as a passive scalar by setting `active_S = 0` in the input file.

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
\partial_t \phi = A \nabla^2 \phi - B \phi (1 - \phi) (1 - 2\phi + C(T - T_m)) ,
$$

where the dimensionless coefficients are given by

$$
A = \frac{6}{5\varepsilon} \frac{\kappa}{U_T H} \frac{c_p \Gamma}{LH}, \quad B = \frac{A}{\varepsilon^2}, \quad C = \frac{\varepsilon H \Delta T}{\Gamma} .
$$

Note that we can typically set $\varepsilon=1/n_x^r$ to be the grid spacing, so that choosing $A$ determines $B$.
Furthermore, $A$ can be written as $A=1.2 (Pe_T \mathcal{S} C)^{-1}$, so choosing $C$ also determines $A$.
In this case, we can write down the dimensionless phase-field equation in terms of the physical parameters and a single model parameter $C$:

$$
\frac{\partial\phi}{\partial t} = \frac{6}{5 Pe_T \mathcal{S} C} \left\{ \nabla^2 \phi - \frac{1}{\varepsilon^2} \phi (1 - \phi) \left[1 - 2\phi + C(T - T_m)\right]\right\}.
$$

In the absence of temperature variations, $\phi$ has a steady state solution of

$$
\phi = \frac{1}{2}\left(1 + \tanh \left(\frac{x-x_i}{2\varepsilon}\right)\right) ,
$$

where $x_i$ is the position of the fluid-solid interface.

### Parameter choices
A wide range of parameters have appeared in the phase-field implementation, so it is important to identify what values are relevant in the current setup.
Two examples of parameters used in the literature are detailed below.

#### Hester et al. (2021)
This study defines all the parameters in dimensional units.
The physical parameters used are

$$
\Delta T = 20°\mathrm{C}, \quad L = 3.34\times 10^2 \mathrm{Jg^{-1}}, \quad c_p = 4.2 \mathrm{Jg^{-1}°C^{-1}}, \\
H=15 \mathrm{cm}, \quad \nu = 1.3\times 10^{-2}\mathrm{cm^2 s^{-1}}, \quad \kappa_T = 1.86\times 10^{-3} \mathrm{cm^2s^{-1}} .
$$

Rather than using a linear equation of state, Hester et al. use a more realistic equation appropriate for seawater.
Given the temperature range of $0-20°\mathrm{C}$, we can approximate the thermal expansion coefficient as $\alpha \approx 10^{-4}$.
This gives dimensionless physical parameters as

$$
Ra_T \approx 2.7\times 10^{8}, \quad Pr = 7, \quad \mathcal{S}\approx 4.
$$

The numerical parameters for the phase field equation are also given in dimensional units:

$$
\hat{\varepsilon} = 0.01 \mathrm{cm}, \quad \Gamma = 0.2 \mathrm{cm °C}, \quad \hat{\eta} = 0.1(\beta \hat{\varepsilon})^2 = 2.275\times 10^{-5} \mathrm{cm^2} .
$$

**N.B.** The value of $\hat{\eta}$ quoted in the paper is actually the value for the damping time scale $\hat{\eta}/\nu = 1.75\times 10^{-3} \mathrm{s}$.
The damping term used is not "optimal" in the sense that an additional prefactor of 0.1 is used.
The reason that this was used is unclear.

Scaling by the domain height gives $\varepsilon \approx 6.6\times 10^{-4}$ and $\eta \approx 10^{-7}$.
By comparison, the dimensionless grid spacing in the vertical is $\Delta z = 1/n_z = 6.5\times 10^{-4}$ and the dimensional time step is $\Delta t = 1.6\times 10^{-4} \mathrm{s}$.
Finally, we can recover the coefficients of the phase-field equation:

$$
A = \frac{1.2}{\mathcal{S} Pe_T C} \approx \frac{0.3}{Pe_T} , \quad B = \frac{A}{\varepsilon^2} \approx 11.8, \quad C = 1.
$$

#### Couston et al. (2021)
Opposite to the system presented above, in this study the phase-field variable $\phi$ takes the value 1 in the fluid and 0 in the solid.
This results in a slight modification of the equations, but the relevant coefficients are still comparable.
Couston et al. present their results in dimensionless form, but where the time scale used is the diffusive time scale $H^2/\kappa_T$ rather than the free-fall time scale $H/U_T$.
We therefore need to account for a factor of $Pe_T=HU_T/\kappa_T=\sqrt{RaPr}$ when scaling some of their coefficients to our dimensionless system.

The convective simulation presented by Couston et al uses the following dimensionless physical parameters:

$$
Ra = 4.5\times 10^5, \quad Pr = 1, \quad \mathcal{S} = 1.
$$

The "diffuse interface thickness" discussed by Couston et al is related to our parameter $\varepsilon$ by $\delta = 4\varepsilon$.
The interface thickness $\varepsilon$ is therefore equal to *half* the grid spacing stated, giving $\varepsilon \approx 8\times 10^{-4}$.
The remaining (appropriately scaled) numerical parameters for the phase field they use are

$$
A = \frac{6}{5\mathcal{S}\sqrt{RaPr}} \approx 1.7 \times 10^{-3}, \quad B = \frac{A}{\varepsilon^2} = \frac{A}{(8\times 10^{-4})^2} \approx 2656, \quad C = 1,
$$

and the volume penalty prefactor used is

$$
\eta = (\beta \varepsilon)^2 = 1.46 \times 10^{-6} .
$$

Couston et al remark that $\eta/2 = 7\times 10^{-7}$ is used as an upper bound for the time step, although due to their use of the diffusive time scale, this is equivalent to $\Delta t \leq 4.6\times 10^{-4}$.

Instead of scaling by the free-fall time scale, it may be more sensible to rescale the parameters by the Peclet number associated with the bulk velocity of the flow $Pe=RePr$, since the flow is primarily driven by a large-scale pressure gradient.
Only $A$ and $B$ are affected by this change, and taking $Re\approx 2000$ results in

$$
A = \frac{6}{5\mathcal{S}RePr} \approx 6\times 10^{-4}, \quad B = \frac{A}{\varepsilon^2} \approx 937 .
$$

Such a rescaling would also modify Couston et al's time step restriction to $\Delta t \leq 1.4\times 10^{-3}$.

### Double-diffusive melting
Following Hester et al. (2020), we can also simulate the ablation of ice in salt water.
In this case, the velocity of the ice-water boundary depends on the gradients of both temperature and salinity at the interface.
In our non-dimensionalization (using the free-fall velocity scale for temperature), the full collection of equations reads as follows:

$$
\partial_t \boldsymbol{u} + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) \boldsymbol{u} = - \boldsymbol{\nabla} p + (T - {R_\rho}^{-1} S) \mathbf{\hat{e}_g} + \sqrt{\frac{Pr}{Ra_T}}\nabla^2 \boldsymbol{u} ,
$$

$$
\partial_t T + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) T = \frac{1}{\sqrt{Ra_T Pr}} \nabla^2 T + \mathcal{S} \partial_t \phi ,
$$

$$
\partial_t S + (\boldsymbol{u} \cdot \boldsymbol{\nabla}) S = \frac{\tau}{\sqrt{Ra_T Pr}} \left(\nabla^2 S - \frac{\boldsymbol{\nabla} \phi \cdot \boldsymbol{\nabla} S}{1 - \phi + \delta}\right) + \frac{S\partial_t \phi}{1 - \phi + \delta},
$$

$$
\partial_t \phi = A \nabla^2 \phi - B \phi (1 - \phi) (1 - 2\phi + C(T - T_m +\Lambda S)) .
$$

Two new dimensionless parameters have appeared in the equations.
$\delta \ll 1$ is a small parameter that is solely present to stabilise the terms on the right hand side of the salinity equation.
In the phase-field equation, we obtain a new physical parameter defined as

$$
\Lambda = \frac{\lambda \Delta S}{\Delta T} ,
$$

where $\Delta S$ and $\Delta T$ are the previously used scales for salinity and temperature, and $\lambda$ is the liquidus slope.
This refects how the presence of salinity at the ice-water interface lowers the local melting temperature.

## Boundary conditions