SHELL = /bin/sh
#
# To run under unix systems including cygwin or mingw on windows make the software with the following:
# 'make'
#
# To make a native windows executable under mingw use
# 'make OS=WINDOWS'
#
# To clean all object files, module files and executables use
# 'make clean' 
#
# Note that you may need to change the name of the fortran compiler (FC below)
# and FLAGS depending on your system. The setup below is suitable for gnufortran
#
FC=gfortran
#
ifeq ($(OS),WINDOWS)
	MODE=WINDOWS
else
	MODE=UNIX
endif

# Compilation flags, Option 1. Here we compile with bounds check - slow but good for development testing...

FLAGS=  -D$(MODE) -cpp -std=f2008 -fbounds-check -finit-real=nan -finit-integer=nan -J$(OBJ_MOD_DIR) -I$(OBJ_MOD_DIR) -freal-4-real-8

# Compilation flags, Option 1b. as for option 1 but we compile with -static (for Windows)

#FLAGS=  -D$(MODE) -static -cpp -std=f2008 -fbounds-check -finit-real=nan -finit-integer=nan -J$(OBJ_MOD_DIR) -I$(OBJ_MOD_DIR)

# Compilation flags, Option 2.  Compilation flags, Here we compile with optimisation - fastest

#FLAGS= -D$(MODE) -cpp -J$(OBJ_MOD_DIR) -I$(OBJ_MOD_DIR) -O4

# Compilation flags, Option 2b. as for option 2 but we compile with -static (for Windows)

#FLAGS= -D$(MODE) -static -cpp -J$(OBJ_MOD_DIR) -I$(OBJ_MOD_DIR) -O4

SRC_DIRECTORY:= $(shell pwd)

TOP_LEVEL_DIRECTORY:=$(shell dirname ${SRC_DIRECTORY})

EXECUTABLE_DIR=$(TOP_LEVEL_DIRECTORY)/bin

OBJ_MOD_DIR=$(TOP_LEVEL_DIRECTORY)/obj_mod

# export the directories, compiler name and flags

export OBJ_MOD_DIR
export FLAGS
export FC

#  Specify the object files to compile in the sub-directories and the
#  modules which are required to specify the dependencies correctly for a F90 project.
#  Then we can export them so that they can be used in dependency lists.
#  Note that the order of compilation has been arranged carefully
#  such that compilation of files is only dependent on what has happened upstream in the make process
#  If a fortran source file uses a module then the dependency must be on both the .mod and the .o file
#  as the .mod file timestamp does not necessarily change when the module is recompiled if the contents of the mod file do not change
#  see for example: https://gcc.gnu.org/ml/fortran/2007-11/msg00013.html
# Note also that gfortran creates module files with lower case names whatever the case of the module name in the .F90 file...

COMPILATION_DATE:= $(shell date)

include ../version_information.inc

TYPE_SPEC_OBJS= $(OBJ_MOD_DIR)/type_specifications.o 
export TYPE_SPEC_OBJS
TYPE_SPEC_MODULE=$(OBJ_MOD_DIR)/type_specifications.mod  $(TYPE_SPEC_OBJS)
export TYPE_SPEC_MODULE

GENERAL_MODULES_OBJS= $(OBJ_MOD_DIR)/constants.o \
                      $(OBJ_MOD_DIR)/general_module.o \
                      $(OBJ_MOD_DIR)/frequency_spec.o 
export GENERAL_MODULES_OBJS
GENERAL_MODULE=  $(OBJ_MOD_DIR)/general_module.o 
export GENERAL_MODULE
CONSTANTS_MODULE=$(OBJ_MOD_DIR)/constants.mod $(OBJ_MOD_DIR)/constants.o
export CONSTANTS_MODULE
FREQ_MODULE=$(OBJ_MOD_DIR)/frequency_spec.mod $(OBJ_MOD_DIR)/frequency_spec.o
export FREQ_MODULE

EISPACK_OBJS=$(OBJ_MOD_DIR)/eispack.o 
export EISPACK_OBJS
EISPACK_MODULE=$(OBJ_MOD_DIR)/eispack.o $(OBJ_MOD_DIR)/eispack.mod
export EISPACK_MODULE

MATHS_OBJS=$(OBJ_MOD_DIR)/maths.o
export MATHS_OBJS
MATHS_MODULE=$(OBJ_MOD_DIR)/maths.o $(OBJ_MOD_DIR)/maths.mod
export MATHS_MODULE

POLYNOMIAL_AND_FILTER_MODULES_OBJS=$(OBJ_MOD_DIR)/filter_module.o
export POLYNOMIAL_AND_FILTER_MODULES_OBJS
POLYNOMIAL_AND_FILTER_MODULE=$(OBJ_MOD_DIR)/filter_module.o $(OBJ_MOD_DIR)/filter_module.mod
export POLYNOMIAL_AND_FILTER_MODULES


OBJS= \
$(TYPE_SPEC_OBJS) \
$(GENERAL_MODULES_OBJS) \
$(EISPACK_OBJS) \
$(MATHS_OBJS) \
$(POLYNOMIAL_AND_FILTER_MODULES_OBJS) 

MODS= \
$(TYPE_SPEC_MODULE) \
$(CONSTANTS_MODULE) \
$(GENERAL_MODULE) \
$(FREQ_MODULE) \
$(EISPACK_MODULE) \
$(MATHS_MODULE) \
$(POLYNOMIAL_AND_FILTER_MODULES)

default: \
MAKE_COMPILATION_DATE \
MAKE_TYPE_SPECIFICATION_MODULE \
MAKE_GENERAL_MODULES \
MAKE_EISPACK \
MAKE_MATHS_MODULES \
MAKE_POLYNOMIAL_AND_FILTER_MODULES \
network_synthesis \
compare_results
	
MAKE_COMPILATION_DATE:
	echo "NETWORK_SYNTHESIS_compilation_date='$(COMPILATION_DATE)' " > compilation_date.inc 

MAKE_TYPE_SPECIFICATION_MODULE:
	( cd TYPE_SPECIFICATION_MODULE ; $(MAKE) )

MAKE_GENERAL_MODULES:
	( cd GENERAL_MODULES ; $(MAKE) )
	
MAKE_EISPACK:
	( cd EISPACK ; $(MAKE) )

MAKE_MATHS_MODULES:
	( cd MATHS_MODULES ; $(MAKE) )
	
MAKE_POLYNOMIAL_AND_FILTER_MODULES:
	( cd POLYNOMIAL_AND_FILTER_MODULES ; $(MAKE) )
	        
include_f90_files=\
include_local_filter_subroutines.F90 \
include_test_positive_real.F90 \
include_minimum_resistance_function.F90 \
include_write_spice_model.F90 \
include_read_and_plot_function.F90 \
include_write_frequency_response.F90 \
include_pole_residue_test_functions.F90 \
include_RLC_test.F90 \
include_RC_test.F90 \
include_RL_test.F90 \
include_LC_test.F90 \
include_C_test.F90 \
include_L_test.F90 \
include_R_test.F90 \
include_R2_test.F90 \
include_BRUNE_test.F90

network_synthesis: network_synthesis.F90 $(include_f90_files) $(MODS)
	$(FC) $(FLAGS) -o $(EXECUTABLE_DIR)/network_synthesis network_synthesis.F90 $(OBJS) $(LIBS)

compare_results: compare_results.F90 $(MODS)
	$(FC) $(FLAGS) -o $(EXECUTABLE_DIR)/compare_results compare_results.F90 $(OBJS) $(LIBS)

clean:
	( rm -f $(EXECUTABLE_DIR)/network_synthesis* )
	( rm -f $(EXECUTABLE_DIR)/compare_results* )
	( rm -f *.o )
	( rm -f $(OBJ_MOD_DIR)/*.mod )
	( rm -f $(OBJ_MOD_DIR)/*.o )
# note that on windows machines the executables will have a .exe extension
