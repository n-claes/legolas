# makefile for the finite elements code LEGOLAS

# activate lots of flags for debugging
# -fcheck=all and -fbounds-check: checks for arrays and array calling
# -Wall: warn about all, such at built-in names, wrong intent variables, etc.
# -Wextra: even more warnings, such as unused variables
# -Wconversion: checks for possible conversion errors
# -pedantic: checks backwards compatibility with F95
# -fbacktrace: produce backtrace when crashing
FFLAGS   := -fcheck=all -fbounds-check -Wall -Wextra -Wconversion \
            -pedantic -fbacktrace
FC       := gfortran

LIBS    := -lblas -llapack

SRCDIR   := src
BINDIR   := bin
MODDIR   := mod

OUTDIR   := output
OUTPUT   := legolas

objects  := $(addprefix $(BINDIR)/,         \
              mod_physical_constants.o      \
              mod_global_variables.o        \
              mod_types.o                   \
              mod_output.o                  \
              mod_check_values.o            \
              mod_spline_functions.o        \
              mod_input.o                   \
              mod_resistivity.o             \
              mod_cooling_curves.o          \
              mod_radiative_cooling.o       \
              mod_thermal_conduction.o      \
              mod_grid.o                    \
              mod_equilibrium_derivatives.o \
              mod_equilibrium.o             \
              mod_make_subblock.o           \
              mod_boundary_conditions.o     \
              mod_matrix_creation.o         \
              mod_solvers.o                 \
              mod_eigenfunctions.o          \
             )

main_objects := $(addprefix $(BINDIR)/, \
                  main.o          \
                 )


all: compile_modules      \
     compile_main_program \
     main_exec

compile_modules: $(objects)

compile_main_program: $(main_objects)

main_exec:
	$(FC) $(FFLAGS) $(objects) $(main_objects) $(LIBS) -o $(OUTPUT)

# check directories
$(BINDIR):
	mkdir -p $(BINDIR)
	mkdir -p $(MODDIR)
	mkdir -p $(OUTDIR)
	mkdir -p $(OUTDIR)/eigenfunctions
	mkdir -p $(OUTDIR)/figures

# root modules
$(BINDIR)/%.o: $(SRCDIR)/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@  -J $(MODDIR)

# physics modules
$(BINDIR)/%.o: $(SRCDIR)/physics/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@ -J $(MODDIR)

# dataIO modules
$(BINDIR)/%.o: $(SRCDIR)/dataIO/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@ -J $(MODDIR)

clean:
	rm -rf $(BINDIR)
	rm -rf $(MODDIR)
	rm -f $(OUTPUT)
	rm -rf $(OUTDIR)

output_clean:
	rm -rf $(OUTDIR)
	mkdir -p $(OUTDIR)/eigenfunctions
	mkdir -p $(OUTDIR)/figures

run:
	make
	./legolas -i legolas_config.par
