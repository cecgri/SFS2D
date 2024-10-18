.SUFFIXES: .F90 .o .mod .F90.o

FC = /usr/local/bin/gfortran

##########################
EXECT = SFS
default: $(EXECT)
##########################

LINK = $(FC)  -v

# MPIFLAGS = -L/usr/local/mpi/lib
FLAGS = -I mod/

OPTIONS = -extend_source -O4 -Wall

MOD_FILES = inc/general_mod.F90 inc/mylib_mod.F90 inc/grids_mod.F90 inc/io_mod.F90 inc/relax_mod.F90
F90_FILES = src/main.F90 src/adi.F90 src/boundaries.F90 src/init.F90 src/geometry.F90 src/mpdata.F90 src/outputs.F90 src/solvers.F90 src/time.F90 src/viscosity.F90

MODS = $(patsubst inc/%.F90, mod/%.o, $(MOD_FILES))
OBJS = $(patsubst src/%.F90, obj/%.o, $(F90_FILES))

#BINARY : link compiled files
$(EXECT) : $(MODS) $(OBJS)
	$(LINK) -o $@ $(MODS) $(OBJS)

#MODULES
mod/general_mod.o: inc/general_mod.F90
	$(FC) $(FLAGS) $(OPTIONS) -c $< -o $@
	mv -f general.mod mod/

mod/mylib_mod.o: inc/mylib_mod.F90 mod/general_mod.o
	$(FC) $(FLAGS) $(OPTIONS) -c $< -o $@
	mv -f mylib.mod mod/

mod/grids_mod.o: inc/grids_mod.F90 mod/general_mod.o mod/mylib_mod.o
	$(FC) $(FLAGS) $(OPTIONS) -c $< -o $@
	mv -f grids.mod mod/

mod/io_mod.o: inc/io_mod.F90 mod/general_mod.o mod/grids_mod.o
	$(FC) $(FLAGS) $(OPTIONS) -c $< -o $@
	mv -f io.mod mod/

mod/relax_mod.o: inc/relax_mod.F90 mod/general_mod.o mod/grids_mod.o mod/io_mod.o
	$(FC) $(FLAGS) $(OPTIONS) -c $< -o $@
	mv -f relax.mod mod/

#OTHER .F90 - compile without linking
obj/%.o: src/%.F90 $(MODS)
	$(FC) $(FLAGS) $(OPTIONS) -c $< -o $@

all: $(EXECT)

clean :
	rm -f obj/*.o mod/*.o mod/*.mod *.mod

allclean :
	rm -f obj/*.o mod/*.o mod/*.mod *.mod $(EXECT)
