# makefile for the finite elements source code

# activate lots of flags for debugging
# -fcheck=all and -fbounds-check: checks for arrays and array calling
# -Wall: warn about all, such at built-in names, wrong intent variables, etc.
# -Wextra: even more warnings, such as unused variables
# -Wconversion: checks for possible conversion errors
# -pedantic: checks backwards compatibility with F95
# -fbacktrace: produce backtrace when crashing
FFLAGS   := -fcheck=all -fbounds-check -Wall -Wextra -Wconversion \
					 -pedantic -fbacktrace
FC  		 := gfortran

SRCDIR   := src
BINDIR   := bin
MODDIR   := mod

OUTPUT   := program.exe

objects  := $(addprefix $(BINDIR)/, 			  		  \
								mod_physical_constants.o					\
								mod_global_variables.o 						\
								mod_gravity.o											\
								mod_cooling_curves.o							\
								mod_radiative_cooling.o						\
								mod_thermal_conduction.o					\
								mod_spline_functions.o						\
								mod_grid.o												\
								mod_equilibrium.o									\
								mod_equilibrium_derivatives.o			\
								mod_make_subblock.o								\
								mod_setup_matrix_b.o							\
								mod_setup_matrix_a.o							\
								main.o														\
							)

all: compile_modules		\
		 compile_main

compile_modules: $(objects)

compile_main:
	$(FC) $(FFLAGS) $(objects) -o $(OUTPUT)

# check directories
$(BINDIR):
	mkdir -p $(BINDIR)
	mkdir -p $(MODDIR)

# root modules
$(BINDIR)/%.o: $(SRCDIR)/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ -o $@ -J $(MODDIR)

# physics modules
$(BINDIR)/%.o: $(SRCDIR)/physics/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ -o $@ -J $(MODDIR)

clean:
	rm -rf $(BINDIR)
	rm -rf $(MODDIR)
	rm -f $(OUTPUT)
	make

purge:
	rm -rf $(BINDIR)
	rm -rf $(MODDIR)
	rm -f $(OUTPUT)
