.PHONY: all deps meteoio snowpack test clean clean-all distclean

NPROC := $(shell nproc 2>/dev/null || echo 4)
CMAKE_MIN_VERSION := 3.5

all: deps
	$(MAKE) -C src libcoupler.a

deps: meteoio snowpack

meteoio:
	@cmake_version=$$(cmake --version | head -1 | cut -d' ' -f3 2>/dev/null || echo "0.0"); \
	cmake_major=$$(echo $$cmake_version | cut -d'.' -f1); \
	if [ $$cmake_major -gt 3 ]; then \
		echo "WARNING: CMake $$cmake_version detected. Using compatibility flag."; \
	fi
	git submodule update --init --recursive external/meteoio
	mkdir -p external/meteoio/build
	cd external/meteoio/build && \
		cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=$(CMAKE_MIN_VERSION) \
		-DCMAKE_INSTALL_PREFIX=../lib -DBUILD_SHARED_LIBS=ON \
		-DBUILD_PLUGIN=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TOOLS=OFF -DBUILD_TESTS=OFF \
		-DENABLE_DBO=OFF -DENABLE_IMIS=OFF -DENABLE_MYSQL=OFF -DENABLE_ODBC=OFF \
		-DENABLE_ORACLE=OFF -DENABLE_POSTGRESQL=OFF -DENABLE_SQLITE=OFF \
		-DENABLE_SMET=ON -DENABLE_GSN=ON -DENABLE_BORMA=ON -DENABLE_GRASS=ON \
		-DENABLE_GEOTOP=ON -DENABLE_PROJ=ON -DENABLE_CURL=ON && \
		$(MAKE) -j $(NPROC) && $(MAKE) install

snowpack: meteoio
	git submodule update --init --recursive external/snowpack
	mkdir -p external/snowpack/build
	cd external/snowpack/build && \
		cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=$(CMAKE_MIN_VERSION) \
		-DCMAKE_INSTALL_PREFIX=../lib -DBUILD_SHARED_LIBS=ON \
		-DENABLE_CRYOWRF=ON -DBUILD_PLUGIN=OFF -DBUILD_EXAMPLES=OFF \
		-DBUILD_TOOLS=OFF -DBUILD_TESTS=OFF -DOPERATIONAL=ON \
		-DMETEOIO_BASE=../../meteoio/lib && \
		$(MAKE) -j $(NPROC) && $(MAKE) install

test: all
	$(MAKE) -C tests
	cd tests && \
		export LD_LIBRARY_PATH=$(PWD)/external/snowpack/lib:$(PWD)/external/meteoio/lib/lib:$$LD_LIBRARY_PATH && \
		$(MAKE) run

clean:
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	@if [ -d external/meteoio/build ]; then $(MAKE) -C external/meteoio/build clean; fi
	@if [ -d external/snowpack/build ]; then $(MAKE) -C external/snowpack/build clean; fi

clean-all: clean
	rm -rf external/meteoio/build external/meteoio/lib
	rm -rf external/snowpack/build external/snowpack/lib

distclean: clean-all
	@echo "WARNING: This will remove ALL untracked files in submodules!"
	@echo "Press Ctrl-C to cancel, or Enter to continue..."
	@read dummy
	git submodule foreach --recursive git clean -fdx
	git submodule foreach --recursive git reset --hard
