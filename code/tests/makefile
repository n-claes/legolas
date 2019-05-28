# makefile for the finite elements source code
# order:
# 1. mod_global_variables
# 2. mod_setup_grid
# 3. mod_setup_equilibrium
# 4. mod_setup_matrix_b
# 5. main

COMPILER  := gfortran

SRC := /media/niels/DATA/School/PhD/finite_elements_code/code/src/

MODULE_FLAGS	:= -c -fcheck=all -fbounds-check
# activate lots of flags for debugging
# -fcheck=all and -fbounds-check: checks for arrays and array calling
# -Wall: warn about all, such at built-in names, wrong intent variables, etc.
# -Wextra: even more warnings, such as unused variables
# -Wconversion: checks for possible conversion errors
# -pedantic: checks backwards compatibility with F95
# -fbacktrace: produce backtrace when crashing
MAIN_FLAGS		:= -fcheck=all -fbounds-check -Wall -Wextra -Wconversion -pedantic -fbacktrace
MODULES			  :=	$(SRC)mod_global_variables.f08 		\
								  $(SRC)mod_setup_grid.f08 				\
								  $(SRC)mod_setup_equilibrium.f08
TEST_MODULES  := test_grid.f08 \
								 mod_assert.f08

MAIN	 		:= run_tests.f08

OUTPUT		:= execute_tests.exe


all:	compile_source_modules	\
			compile_test_modules 		\
			compile_main	\
			clean

compile_source_modules:
	$(COMPILER) $(MODULE_FLAGS) $(MODULES)
compile_test_modules:
	$(COMPILER) $(MODULE_FLAGS) $(TEST_MODULES)
compile_main:
	$(COMPILER) $(MAIN_FLAGS) $(MODULES) $(TEST_MODULES) $(MAIN) -o $(OUTPUT)
clean:
	if [ mod_global_variables.o ]; then  \
		rm mod_global_variables.o;		 \
	fi
	if [ mod_setup_grid.o ]; then		 \
		rm mod_setup_grid.o;			 \
	fi
	if [ mod_setup_equilibrium.o ]; then \
		rm mod_setup_equilibrium.o;		 \
	fi
	if [ test_grid.o ]; then \
		rm test_grid.o;				\
	fi
	if [ mod_assert.o ]; then \
		rm mod_assert.o;			\
	fi
