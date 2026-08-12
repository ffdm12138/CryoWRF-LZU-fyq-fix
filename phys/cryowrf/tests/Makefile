.PHONY: all run clean

FC = mpif90

# Use WRF centralized libs if building as part of WRF, otherwise use local libs
ifdef WRF_LIBS_DIR
INCLUDES = -I../src -I$(WRF_LIBS_DIR)/include
LIBDIRS = -L../src -L$(WRF_LIBS_DIR)/lib
RPATH = -Wl,-rpath,$(WRF_LIBS_DIR)/lib
else
INCLUDES = -I../src -I../external/meteoio/lib/include -I../external/snowpack/lib/include
LIBDIRS = -L../src -L../external/snowpack/lib -L../external/meteoio/lib/lib
RPATH = -Wl,-rpath,../external/snowpack/lib -Wl,-rpath,../external/meteoio/lib/lib
endif

FCFLAGS = -O2 -g $(INCLUDES) -ffpe-summary=none
LIBS = -lcoupler -lsnowpack -lmeteoio -lstdc++
LDFLAGS = $(LIBDIRS) $(LIBS) $(RPATH)

all: test_coupler

test_coupler: test_coupler.f90 ../src/libcoupler.a
	$(FC) $(FCFLAGS) $< -o $@ $(LDFLAGS)

run: test_coupler
	./test_coupler

clean:
	rm -f test_coupler test_funcs *.o *.mod
