EXEC = calc_sigma_pol
FC = ifort
FFLAGS = -g -C

MODULES=module_data.f90

SRC=\
main.f90 \
input.f90 \
read_wfc_PEtot.f90 \
read_wfc_CPMD.f90 \
read_wfc_VASP.f90 \
cartesian_to_spherical.f90 \
integral.f90

OBJMOD=$(MODULES:.f90=.o)
OBJ=$(SRC:.f90=.o)

cart2intnl: $(OBJMOD) $(OBJ)
	$(FC) -o $(EXEC) $(OBJMOD) $(OBJ)

$(OBJMOD): %.o: %.f90
	$(FC) -c $(FFLAGS) $<

$(OBJ): %.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f *.o *.mod *~
