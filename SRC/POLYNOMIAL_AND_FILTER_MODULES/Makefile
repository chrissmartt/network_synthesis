default: $(POLYNOMIAL_AND_FILTER_MODULES_OBJS)
# note we add the include files to the dependencies
$(OBJ_MOD_DIR)/%.o: %.F90 $(TYPE_SPEC_MODULE) $(CONSTANTS_MODULE) \
                    polynomial_types.F90 polynomial_functions.F90 polynomial_subroutines.F90 polynomial_operators.F90 \
                    filter_types.F90 filter_functions.F90 filter_operators.F90 filter_subroutines.F90 findroots.F90
	$(FC) $(FLAGS) -c -o $@ $<
