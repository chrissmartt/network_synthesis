default: $(GENERAL_MODULES_OBJS)
#
# note we have added the include files to the dependencies
#
$(OBJ_MOD_DIR)/constants.o: constants.F90 $(TYPE_SPEC_MODULE)
	$(FC) $(FLAGS) -c -o $@ $<
#
$(OBJ_MOD_DIR)/general_module.o: general_module.F90 $(TYPE_SPEC_MODULE) program_status.F90 \
                                 read_write_subroutines.F90 read_version.F90 
	$(FC) $(FLAGS) -c -o $@ $<
#
$(OBJ_MOD_DIR)/frequency_spec.o: frequency_spec.F90 $(TYPE_SPEC_MODULE)
	$(FC) $(FLAGS) -c -o $@ $<

