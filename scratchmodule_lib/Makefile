AR     = /usr/bin/ar ruv
RANLIB = /usr/bin/ranlib
FFLAGS =  -g -traceback -xHOST -qopt-zmm-usage=high -unroll -funroll-loops -O2 -IPO
CPPFLAGS =

FC = ifort

SCRATCHMODULE_LIB = libscratch_module.a
SRC_SCRATCHMODULE = scratch_data.f90 scratch_utils.f90 scratch_managment.f90 scratch_segment_managment.f90 scratch_interface.f90
.SUFFIXES: .f90

%.mod:

%.o:    %.f90 
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(@F:.o=.f90)


OBJ_SCRATCHMODULE     = $(SRC_SCRATCHMODULE:%.f90=%.o)
$(SCRATCHMODULE_LIB): $(OBJ_SCRATCHMODULE)
	$(AR) $(SCRATCHMODULE_LIB) $(OBJ_SCRATCHMODULE)
	$(RANLIB) $(SCRATCHMODULE_LIB)


scratch_kinds.mod: scratch_kinds.o
	@true
scratch_data.mod: scratch_data.o
	@true
scratch_utils.mod: scratch_utils.o
	@true
scratch_managment.o: scratch_data.mod scratch_managment.f90
scratch_segment_mangament.o: scratch_data.mod scratch_utils.mod scratch_managment.mod scratch_segment_managment.f90
scratch_interface.o: scratch_segment_managment.o scratch_interface.f90


clean: 
	rm -f *.o *.mod *.a
