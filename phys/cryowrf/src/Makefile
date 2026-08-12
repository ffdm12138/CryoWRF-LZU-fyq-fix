FC = mpif90 # mpif90 (for MPI compatibility)
CXX = mpicxx # mpicxx (for MPI compatibility)

COMPNAME := $(shell hostname -s | sed 's/[0-9]//g' )
FCFLAGS = -O2 -ftree-vectorize -funroll-loops -w -ffree-form -ffree-line-length-none -fconvert=big-endian -frecord-marker=4 -fopenmp
CCFLAGS = -O2 -w 
ifeq ($(COMPNAME),cryospc)
   $(info UNAME is $(COMPNAME))
   FCFLAGS = -O3 -g -fbacktrace -ffree-line-length-512
   CCFLAGS = -O3 -g 
endif

ifeq ($(COMPNAME),daint)
   FC = ftn
   CXX = CC
   FCFLAGS = -O3 -g -traceback
   CCFLAGS = -O3 -g -traceback
endif

# Set include paths based on context
ifdef WRF_LIBS_DIR
# Building as part of WRF - centralized location
INCLUDES := -I$(WRF_LIBS_DIR)/include
else
# Standalone build - local installation
INCLUDES := -I../../../external/cryowrf/include
endif

all: libcoupler.a

libcoupler.so : coupler_mod.o global_f_contracts.o funcs.o coupler_capi.o Coupler.o
	gcc -shared $^ -o $@

libcoupler.a : coupler_mod.o global_f_contracts.o funcs.o coupler_capi.o Coupler.o
	ar rcs $@ $^

global_f_contracts.o: global_f_contracts.f90
	${FC} ${FCFLAGS} -c $^ -o $@ $(INCLUDES)

funcs.o : coupler_mod.o global_f_contracts.o

%.o : %.f90
	${FC} ${FCFLAGS} -c $^ -o $@ $(INCLUDES)

%.o : %.cpp
	${CXX} ${CCFLAGS} -c $< -o $@ $(INCLUDES)

.PHONY : clean

clean :
	${RM} -rf *.o *.mod test.x *.a *.so
