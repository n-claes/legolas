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

objects  := $(addprefix $(BINDIR)/, \
              mod_global_variables.o \
							mod_equilibrium_params.o \
							mod_logging.o \
              mod_physical_constants.o \
							mod_units.o	\
							mod_grid.o \
              mod_types.o \
              mod_check_values.o \
              mod_spline_functions.o \
              mod_input.o \
              mod_resistivity.o \
              mod_cooling_curves.o \
              mod_radiative_cooling.o \
              mod_thermal_conduction.o \
							mod_inspections.o \
              mod_equilibrium.o     \
							smod_equil_adiabatic_homo.o \
							smod_equil_constant_current.o \
							smod_equil_discrete_alfven.o \
							smod_equil_flow_driven_instabilities.o \
							smod_equil_gravito_acoustic.o \
							smod_equil_gravito_mhd.o \
							smod_equil_ideal_quasimodes.o \
							smod_equil_interchange_modes.o \
							smod_equil_interface_modes.o \
							smod_equil_internal_kink_instability.o \
							smod_equil_kelvin_helmholtz_cd.o \
							smod_equil_kelvin_helmholtz.o \
							smod_equil_magneto_rotational_instability.o \
							smod_equil_nonuniform_conduction.o \
							smod_equil_resistive_homo.o \
							smod_equil_resistive_tearing_flow.o \
							smod_equil_resistive_tearing.o \
							smod_equil_rotating_plasma_cylinder.o \
							smod_equil_rotating_theta_pinch.o \
							smod_equil_suydam_cluster.o \
							smod_equil_uniform_conduction.o \
							smod_equil_resonant_absorption.o \
							smod_equil_magnetothermal_instabilities.o \
							smod_equil_photospheric_flux_tube.o \
							smod_equil_coronal_flux_tube.o \
							smod_equil_gold_hoyle.o \
							smod_test_beta0.o \
							smod_test_hydro.o \
							smod_user_defined.o \
              mod_make_subblock.o \
              mod_boundary_conditions.o \
              mod_matrix_creation.o \
              mod_solvers.o \
              mod_eigenfunctions.o \
              mod_output.o \
             )

main_objects := $(addprefix $(BINDIR)/, \
                  main.o \
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

# root modules
$(BINDIR)/%.o: $(SRCDIR)/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@  -J $(MODDIR)

# physics modules
$(BINDIR)/%.o: $(SRCDIR)/physics/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@ -J $(MODDIR)

# dataIO modules
$(BINDIR)/%.o: $(SRCDIR)/dataIO/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@ -J $(MODDIR)

# equilibrium submodules
$(BINDIR)/%.o: $(SRCDIR)/equilibria/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(LIBS) -o $@ -J $(MODDIR)

clean:
	rm -rf $(BINDIR)
	rm -rf $(MODDIR)
	rm -f $(OUTPUT)
	rm -rf $(OUTDIR)

output_clean:
	rm -rf $(OUTDIR)

run:
	make
	./legolas -i legolas_config.par
