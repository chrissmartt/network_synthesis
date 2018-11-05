default: $(SFILTER_FIT_OBJS)
#
$(OBJ_MOD_DIR)/%.o: %.F90 $(TYPE_SPEC_MODULE)  $(CONSTANTS_MODULE) Weiner_Hopf.F90 Calculate_Sfilter.F90 \
                    filter_MSE.F90
	$(FC) $(FLAGS) -c -o $@ $<
