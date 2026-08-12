# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.4] - 2026-06-27

### Fixed

- Coupler.cpp: Removed thread-unsafe std::cout debug print statements from OpenMP parallel regions (I,J grid loops). Concurrent writes to stdout from multiple threads were causing heap corruption and simulation crashes in mountain domains. Affected prints: per-layer SNOWPACK initialisation, height_meteo diagnostic, TSK/forcing variable dumps, and compute/write counter.
- WRF/phys/physics_mmm/sf_sfclayrev.F90: Added lower bound clamp (nzol = max(0, nzol)) in psim_stable, psih_stable, psim_unstable, and psih_unstable look-up table functions. In steep mountain terrain, near-neutral z/L values can yield slightly negative floating-point inputs that truncate to index -1, causing out-of-bounds array access and numerical instability.

## [1.2.3] - 2026-04-27

### Changed

- WRF/phys/module_sf_snowpacklsm.f90: added a dt to e_budg to take into account dt<1s.
- Snowpack.cc: changes to not start at printing at nE-1 when throwing compTemperatureProfile.
- Coupler.cpp: changes to recalculate hfx every wrf timestep and not only every snowpack timestep in line 1124: HFX = -1.0 * SnLaws::compSensibleHeatCoefficient(Mdata, vecXdata, meteo_height) * Mdata.density_air * Constants::specific_heat_air * (Mdata.ta - TSK); with meteo_height is first WRF level height.
- module_sf_snowpacklsm.F: change Z0WRF from 0.002 to 0.01 in ALPS which is hardcoded.
- module_sf_snowpacklsm.F: changed SWDOWN to SWNORM, which is the sw rad incoming which includes slope projection and shading.
- SWDOWN and SWNORM should change depending if slope projection and shading is switched on.
- snowpack and meteoio plugins CMakeLists.txt: from CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/meteoio/IOHandler.cmake.cc ${PROJECT_SOURCE_DIR}/meteoio/IOHandler.cc ) to CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/meteoio/IOHandler.cmake.cc ${CMAKE_CURRENT_BINARY_DIR}/IOHandler.cc).

## [1.2.0] - 2023-11-15

### Added

- Repository restructured to work as external module inside WRF.
- Dynamic roughness length (z0) from Snowpack to WRF layers.
- Test data and simple tests for build verification.
- Makefile for build and test automation.
- CRYOWRF-WRF build trigger in CI.
- Air density as 3D timeseries output variable.
- Blowing snow model enhancements with modified exponential profile of saltation mass flux.
- Stability-corrected log profiles for temperature and humidity in fine mesh.
- Daniela's parameterization for particle speed.
- Fine mesh vertical levels stretched.
- Output for wind speed of fine mesh and reciprocal Obukhov length.

### Changed

- Updated to WRF v4.7.1.
- Replaced external libraries with submodules (WRF, WPS, MeteoIO, Snowpack).
- Snow density in Antarctic mode now defined in io.ini file.
- Modified discretization in blowing snow prognostic equations.
- Set h_salt to constant value of 0.15 m.
- Set average diameter at h_salt to 140 microns.

### Fixed

- Bug fix for Snowpack layers filled with -999 at restart.
- Correct number of output variables for TS timeseries output.
- Use actual air density instead of 1 kg m^-3 in calculations.
- Sea issue resolved to enter all gridpoints except sea in land mode.

## [1.1.0] - 2022-06-01

### Added

- Validation datasets published (Gerber et al., 2022).
- Test file for validation purposes.

### Changed

- Updated to match version used in Sharma et al., 2023 GMD publication.

## [1.0.0] - 2021-12-15

### Added

- Initial release of CRYOWRF coupled atmosphere-snow cover model.
- WRF v4.2.1 as atmospheric core.
- SNOWPACK and MeteoIO as land-surface snow cover model.
- Coupler library with Fortran and C++ interfaces.
- Blowing snow physics implementation.
- Time series output for multiple variables.
- Support for air density, sensible heat flux, and surface temperature outputs.
- Documentation and installation instructions.
- License under GNU General Public License v3.
- References to Sharma et al., 2023 GMD paper.
