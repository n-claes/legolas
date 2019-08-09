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
TSTFLAGS := -Wno-integer-division
FC  		 := gfortran

MKLCOMP := -m64 -I${MKLROOT}/include
MKLLINK := ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group \
					 ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
					 ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \
					 ${MKLROOT}/lib/intel64/libmkl_core.a \
					 ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
					 						-Wl,--end-group -lgomp -lpthread -lm -ldl

LIBS 		 := -lblas -llapack

SRCDIR   := src
BINDIR   := bin
MODDIR   := mod

TSTDIR   := tests
OUTDIR   := output

OUTPUT   := esonas.exe
TEST_OUT := run_tests.exe

objects  := $(addprefix $(BINDIR)/, 			  		  \
								mod_physical_constants.o					\
								mod_global_variables.o 						\
								mod_timer.o												\
								mod_types.o												\
								mod_check_values.o								\
								mod_spline_functions.o						\
								mod_info.o												\
								mod_input.o												\
								mod_io.o													\
								mod_resistivity.o									\
								mod_gravity.o											\
								mod_cooling_curves.o							\
								mod_radiative_cooling.o						\
								mod_thermal_conduction.o					\
								mod_grid.o												\
								mod_equilibrium_derivatives.o			\
								mod_equilibrium.o									\
								mod_make_subblock.o								\
								mod_boundary_conditions.o					\
								mod_setup_matrix_b.o							\
								mod_setup_matrix_a.o							\
								mod_solvers.o											\
								mod_eigenfunctions.o							\
							)

main_objects := $(addprefix $(BINDIR)/,						\
								main.o														\
								)

test_objects := $(addprefix $(BINDIR)/,		     		\
								testmod_assert.o									\
								testmod_core_tests.o							\
								testmod_homo_adiabatic.o					\
								testmod_homo_gravity.o						\
								tests_main.o											\
								)

# If compiling for single core
# MKLLINK := $(LIBS)
# MKLCOMP := $(LIBS)

all: compile_modules		  \
		 compile_main_program \
		 main_exec

test: compile_modules			 \
		 	compile_main_program \
		 	compile_test_program \
			test_exec

compile_modules: $(objects)

compile_main_program: $(main_objects)

compile_test_program: $(test_objects)

main_exec:
	$(FC) $(FFLAGS) $(objects) $(main_objects) $(MKLLINK) -o $(OUTPUT)

test_exec:
	$(FC) $(FFLAGS) $(TSTFLAGS) $(objects) $(test_objects) $(MKLLINK) -o $(TEST_OUT)

# check directories
$(BINDIR):
	mkdir -p $(BINDIR)
	mkdir -p $(MODDIR)
	mkdir -p $(OUTDIR)
	mkdir -p $(OUTDIR)/eigenfunctions
	mkdir -p $(OUTDIR)/figures
	mkdir -p $(OUTDIR)/tests
	mkdir -p $(OUTDIR)/tests/figures

# root modules
$(BINDIR)/%.o: $(SRCDIR)/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(MKLCOMP) -o $@  -J $(MODDIR)

# physics modules
$(BINDIR)/%.o: $(SRCDIR)/physics/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(MKLCOMP) -o $@ -J $(MODDIR)

# dataIO modules
$(BINDIR)/%.o: $(SRCDIR)/dataIO/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(MKLCOMP) -o $@ -J $(MODDIR)

# test modules
$(BINDIR)/%.o: $(TSTDIR)/%.f08 | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ $(MKLCOMP) -o $@ -J $(MODDIR)

clean:
	rm -rf $(BINDIR)
	rm -rf $(MODDIR)
	rm -f $(OUTPUT)
	rm -f $(TEST_OUT)
	rm -rf $(OUTDIR)
	make

output_clean:
	rm -rf $(OUTDIR)
	mkdir -p $(OUTDIR)/eigenfunctions
	mkdir -p $(OUTDIR)/figures
	mkdir -p $(OUTDIR)/tests
	mkdir -p $(OUTDIR)/tests/figures

run:
	make
	./esonas.exe -i config_esonas.par
