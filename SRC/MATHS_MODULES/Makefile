default: $(MATHS_OBJS)
#
$(OBJ_MOD_DIR)/%.o: %.F90 $(TYPE_SPEC_MODULE)  $(CONSTANTS_MODULE) dmatrix.F90 cmatrix.F90 FFT.F90
	$(FC) $(FLAGS) -c -o $@ $<
