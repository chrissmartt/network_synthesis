default:  $(EISPACK_OBJS)

# eispack has no dependencies beyond the type specifications

$(OBJ_MOD_DIR)/%.o: %.F90 $(TYPE_SPEC_MODULE)
	$(FC) $(FLAGS) -c -o $@ $<
