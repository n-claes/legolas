# makefile for the finite elements source code
# order:
# 1. mod_global_variables
# 2. mod_setup_grid
# 3. mod_setup_equilibrium
# 4. mod_setup_matrix_b
# 5. main

COMPILER  := gfortran

MODULE_FLAGS	:= -c
MODULES			:=	mod_global_variables.f08 		\
					mod_setup_grid.f08 				\
					mod_setup_equilibrium.f08		\
					mod_setup_matrix_b.f08
	
MAIN	 		:= main.f08

OUTPUT			:= program.exe


all:	compile_modules	\
		compile_main	\
		clean		

compile_modules:
	$(COMPILER) $(MODULE_FLAGS) $(MODULES)
compile_main:
	$(COMPILER) $(MODULES) $(MAIN) -o $(OUTPUT)
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
	if [ mod_setup_matrix_b.o ]; then	 \
		rm mod_setup_matrix_b.o;		 \
	fi

