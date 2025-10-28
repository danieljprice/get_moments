FC=gfortran
FFLAGS=-Wall -fopenmp -fcheck=all -fdefault-real-8 -fdefault-double-8
SRC=fsolve.f90 integrate.f90 reconstruct_moments_gamma.f90 reconstruct_moments_maxent.f90 test_moments.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

recon: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm -f *.o *.mod *.f90~
