# AFiD_MultiRes_staggered
[Documentation](https://physicsoffluids.github.io/AFiD_MultiRes_staggered/)

AFiD is a highly parallel application for simulating canonical flows in a channel domain.
More technical details can be found in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007).
The code is developed jointly by the University of Twente, SURFsara, and the University of Rome "Tor Vergata".

In addition to the method described in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007), this *multi-resolution* version of AFiD evolves a second scalar field.
This scalar field is simulated on a refined grid to allow for low diffusivity values or high Schmidt numbers.
The multi-resolution method applied is detailed in [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031) and was implemented into AFiD by Chong Shen Ng.
Interpolation between the two grids is performed using a two-point Hermite interpolation scheme.

One key difference from the original version of AFiD is that the scalar fields are evaluated at the mid-points of the computational cells.
This staggered grid was implemented by Chris Howland and allows one to change the direction of gravity.

## Contributing
If you would like to contribute to bug fixing/feature development/documentation, please create a new branch, commit your changes and then submit a pull request.
