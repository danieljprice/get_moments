# get_momemts
Fortran 90 module for maximum entropy reconstruction of a function from its moments, essentially a pure Fortran version of PyMaxEnt

# how to use
compile the code (you will need the gfortran compiler):
```
make
```
run the code:
```
./recon
```
this should produce two ascii files, exact.out and recon.out which can be used to plot the exact function (exact.out) and the reconstructed function (recon.out)

# how to use in another code
just use the module reconstruct_from_moments.f90 in your Fortran code. You will need the helper module integrate.f90 which implements a simple integrator using the trapezoidal rule. Feel free to swap this out for something more sophisticated. Also fsolve.f90 which is originally from MINPACK.

# credits
Saad & Ruai 2019 (https://doi.org/10.1016/j.softx.2019.100353)

fsolve.f90: Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN90 version by John Burkardt. This module distributed under the GNU LGPL license


