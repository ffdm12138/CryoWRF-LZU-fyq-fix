# CRYOWRF MODEL version v1.2.2

_Note: This repository contains only the **coupler** code. The full CRYOWRF implementation also includes
[CRYOWRF-WRF](https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF), and [CRYOWRF-WPS](https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WPS)
._

# CRYOWRF Coupler

CRYOWRF is a coupled atmosphere-snow cover model with
[WRF](https://github.com/wrf-model/WRF) as the atmospheric core and
[SNOWPACK](https://code.wsl.ch/snow-models/snowpack) /
[MeteoIO](https://code.wsl.ch/snow-models/meteoio) as an added land-surface snow
cover model.

Since v1.2.0, this repository is designed to be used as an external module
inside WRF, following the pattern of [NOAH-MP](https://github.com/NCAR/noahmp), the exisiting land surface model in WRF. SNOWPACK/MeteoIO can be found in the
/phys directory of the
[CRYOWRF-WRF](https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF)
repository where they sit with associated WRF modifications, but following the
upstream WRF build system.

## Table of Contents

- [Versioning](#versioning)
- [Setup](#setup)
  - [CRYOWRF with WRF](#cryowrf-with-wrf)
    - [Environment Variables](#environment-variables)
    - [Compiling WRF and WPS](#compiling-wrf-and-wps)
  - [Coupler development](#coupler-development)
    - [Prerequisites](#prerequisites)
    - [Build Instructions](#build-instructions)
- [Repository Structure](#repository-structure)
- [Contributing and Development](#contributing-and-development)
  - [Understanding Submodules](#understanding-submodules)
  - [Development Workflow](#development-workflow)
  - [Commit Order](#commit-order)
- [CI/CD Pipeline](#cicd-pipeline)
- [Attribution](#attribution)
  - [Related Repositories](#related-repositories)
  - [References](#references)
    - [Developed by](#developed-by)
    - [Original authors](#original-authors)
    - [Current maintainers](#current-maintainers)
  - [Funding](#funding)
  - [License](#license)

## Versioning

CRYOWRF releases are versioned to match specific versions of WRF, WPS, MeteoIO, and SNOWPACK:

| CRYOWRF Version                            | WRF    | WPS  | MeteoIO              | SNOWPACK             |
| ------------------------------------------ | ------ | ---- | -------------------- | -------------------- |
| v1.0.0 (Published in Sharma et al., 2023)  | v4.2.1 | v4.2 | ~v2.9.0 @ `00b76b37` | ~v3.7.0 @ `e68810d6` |
| v1.1.0 (Validation by Gerber et al., 2022) | v4.2.1 | v4.2 | ~v2.9.0 @ `00b76b37` | ~v3.7.0 @ `3d1fa02f` |
| v1.2.2 (Current)                           | v4.7.1 | v4.6 | ~v2.9.0 @ `00b76b37` | ~v3.7.0 @ `3d1fa02f` |

## Setup

### CRYOWRF with WRF

To use CRYOWRF,

clone [CRYOWRF-WPS](https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WPS), the WRF Preprocessing System

```bash
git clone https://git.wsl.ch/atmospheric-models/CRYOWRF-WPS.git
```

Then clone the
[CRYOWRF-WRF](https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF)
repository, which includes this coupler as a git submodule.

```bash
git clone https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF.git
cd CRYOWRF-WRF
git submodule update --init --recursive
```

Then build WRF as per normal (some general steps are detailed below, but may be
specific to the hardware you are using). The repository includes a `Dockerfile`
that can help to understand a build on ubuntu:24.04. To build this, whilst in
the `CRYOWRF-WRF` folder, run:

```bash
docker build -t cryowrf:latest .
```

#### Environment Variables

Before compiling WRF, you need to set the following environment variables. It is
recommended to add these to your `.bashrc` or shell configuration file:

```bash
# NetCDF and parallel I/O libraries
export NETCDF=/path/to/netcdf        # Path to NetCDF installation
export PNETCDF=/path/to/pnetcdf      # Path to Parallel NetCDF installation
export HDF5=/path/to/hdf5            # Path to HDF5 installation

# WRF core configuration
export WRF_EM_CORE=1                 # Enable Eulerian Mass core
export WRF_NMM_CORE=0                # Disable NMM core
export WRF_DA_CORE=0                 # Disable Data Assimilation core

# Chemistry options (disabled for CRYOWRF)
export WRF_CHEM=0
export WRF_KPP=0

# NetCDF-4 support
export NETCDF4=1
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export WRFIO_NCD_NO_LARGE_FILE_SUPPORT=0
```

#### Compiling WRF and WPS

After setting up environment variables and initializing submodules:

1. **Compile WRF:**

   ```bash
   cd CRYOWRF-WRF
   ./configure
   # Option 35: GNU (gfortran/gcc) dmpar + smpar
   # Choose option 35, then option 1 for basic nesting
   ./compile em_real -j 20 &> compile.log
   cd ..
   ```

2. **Compile WPS:**
   ```bash
   cd CRYOWRF-WPS
   ./configure
   # Option 3: GNU compiler (gfortran/gcc), serial
   ./compile &> compile.log
   cd ..
   ```

**Note:** The above options (35 for WRF, 3 for WPS) are confirmed working on
Ubuntu 24.04. If
attempting to use a `-j` value higher than 20, WRF defaults it back to 2.

Refer to the `Dockerfile` in the CRYOWRF-WRF repository for the most up-to-date
compilation steps and dependencies.

Additional useful information for the installation on different clusters can be found in our [WIKI](https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF/-/wikis/home).

### Coupler development

#### Prerequisites

- MPI Fortran/C++ compilers (`mpif90`, `mpicxx`)
- CMake ≥ 3.5
- Build essentials: `gfortran`, `build-essential`
- Libraries: `libproj-dev`, `libopenmpi-dev`

#### Build Instructions

The `Makefile` can do most of this for you. From the root of the repository,
run:

```bash
make              # Build external libraries and coupler
make test         # Build and run unit tests
make clean        # Clean build artifacts
make distclean    # Reset and clean submodules (be careful!)
```

## Repository Structure

```
CRYOWRF/
├── src/                     # Coupler library source
│   ├── coupler_mod.f90      # Fortran module interface
│   ├── funcs.f90            # WRF global array interface
│   ├── Coupler.cpp/.hpp     # C++ implementation
│   ├── coupler_capi.cpp     # C API bridge
│   └── Makefile
├── external/                # Git submodules
│   ├── meteoio/
│   └── snowpack/
├── tests/                   # Unit tests
│   ├── test_coupler.f90     # Direct coupler tests
│   └── Makefile
└── Makefile                 # Root build orchestration
```

## Contributing and Development

This section explains how to make changes to CRYOWRF and its components. Because
CRYOWRF uses nested git submodules, there is a specific workflow to follow when
contributing code.

### Understanding Submodules

CRYOWRF uses git submodules to manage its dependencies. Think of submodules as
"pointers" to specific commits in other repositories. This design allows us to:

- **Track upstream changes**: WRF, MeteoIO, and SNOWPACK each have their own
  development cycles. Submodules let us integrate new versions when ready.
- **Pin exact versions**: Each CRYOWRF release points to specific, tested
  commits of all dependencies (see the [Versioning](#versioning) table).
- **Keep histories separate**: Changes to SNOWPACK stay in the SNOWPACK
  repository, making it easier to merge upstream updates.

The repository structure looks like this:

```
CRYOWRF-WRF/                          (main repository)
└── phys/cryowrf/                     (submodule -> CRYOWRF coupler)
    └── external/
        ├── meteoio/                  (submodule -> MeteoIO)
        └── snowpack/                 (submodule -> SNOWPACK)
```

Currently, the `cryowrf` branch in both MeteoIO and SNOWPACK holds the
CRYOWRF-specific modifications. These branches have diverged significantly from
their upstream master branches, so for CRYOWRF development purposes, treat them
as the working branches.

**A note on submodule isolation**: Because submodules point to specific commits,
changes you make in one clone won't automatically appear in another. For example,
if you clone SNOWPACK independently and make changes there, your CRYOWRF-WRF
clone won't see those changes until you:

1. Push your SNOWPACK changes to the remote
2. In your CRYOWRF-WRF clone, navigate to the SNOWPACK submodule and pull the
   new commits (`git pull origin <branch>`)
3. Update the submodule pointer in the parent repositories (see
   [Commit Order](#commit-order))

This is also true for the CRYOWRF coupler itself — if you work on it
independently, those changes need to be pushed and then pulled into your
CRYOWRF-WRF working copy.

_Note_: For pushing changes to SNOWPACK and MeteoIO, make sure you have the
correct privileges on those repositories. This may not always be the case —
contact the maintainers if you need access.

### Development Workflow

**Important**: For development work, always clone CRYOWRF-WRF (not this coupler
repository directly). This ensures you have the full context and can test your
changes with WRF.

#### 1. Clone and initialize

```bash
git clone https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF.git
cd CRYOWRF-WRF
git submodule update --init --recursive
```

Your code is now available at:

- Coupler: `CRYOWRF-WRF/phys/cryowrf/`
- MeteoIO: `CRYOWRF-WRF/phys/cryowrf/external/meteoio/`
- SNOWPACK: `CRYOWRF-WRF/phys/cryowrf/external/snowpack/`

#### 2. Navigate to the component you want to modify

For example, to work on SNOWPACK:

```bash
cd phys/cryowrf/external/snowpack
```

#### 3. Check your current state

Before making changes, verify where you are:

```bash
# See which commit you're on
git log --oneline -3

# Check which remote repository you're pointing to
git remote -v
```

#### 4. Create a feature branch and make changes

```bash
# Create a new branch for your feature
git checkout -b feature/my-new-feature

# Make your changes to the code...
# Then stage and commit
git add <files-you-changed>
git commit -m "Add description of your changes"
```

#### 5. Push your changes

```bash
git push -u origin feature/my-new-feature
```

### Commit Order

When you modify a submodule, the parent repository needs to be updated to point
to your new commit. This must be done **from the inside out** — commit and push
the deepest submodule first, then work your way up.

**Example**: You made changes to SNOWPACK

```
CRYOWRF-WRF                    (3. Update last)
└── phys/cryowrf               (2. Update second)
    └── external/snowpack      (1. Commit and push first)
```

Here is the step-by-step process:

#### Step 1: Commit in SNOWPACK (deepest level)

```bash
cd CRYOWRF-WRF/phys/cryowrf/external/snowpack
git checkout -b feature/my-snowpack-change
# ... make changes ...
git add <changed-files>
git commit -m "Describe your SNOWPACK changes"
git push -u origin feature/my-snowpack-change
```

#### Step 2: Update the CRYOWRF coupler (middle level)

```bash
cd ../..  # From external/snowpack -> phys/cryowrf (the coupler directory)
git status
# You'll see something like:
#   modified:   external/snowpack (new commits)

# This is a NEW branch in CRYOWRF (the coupler), not related to the branch
# you just created in SNOWPACK. The name can correspond to other changes
# you're making in the coupler, or simply describe what submodule update
# it includes.
git checkout -b feature/update-snowpack-pointer
git add external/snowpack
git commit -m "Update SNOWPACK pointer: my-snowpack-change"
git push -u origin feature/update-snowpack-pointer
```

#### Step 3: Update CRYOWRF-WRF (top level)

```bash
cd ../..  # From phys/cryowrf -> CRYOWRF-WRF root
git status
# You'll see:
#   modified:   phys/cryowrf (new commits)

# Again, this is a NEW branch in CRYOWRF-WRF, separate from the branches
# in CRYOWRF and SNOWPACK.
git checkout -b feature/integrate-snowpack-change
git add phys/cryowrf
git commit -m "Update CRYOWRF pointer: integrate snowpack changes

Includes [brief description]"
git push -u origin feature/integrate-snowpack-change
```

After pushing, create merge/pull requests for each level, merging from the
inside out (SNOWPACK first, then CRYOWRF, then CRYOWRF-WRF).

_Tip_: The `git status` command in a parent repository will show you when a
submodule has new commits that haven't been recorded yet. This is your reminder
to update the pointer.

## CI/CD Pipeline

CRYOWRF components use the CI pipeline in Gitlab to ensure all components build
correctly together:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         CRYOWRF Repository                              │
├─────────────────────────────────────────────────────────────────────────┤
│  1. Build dependencies stage:                                           │
│     └─ Build MeteoIO and SNOWPACK from submodules                       │
│  2. Build stage:                                                        │
│     └─ Compile CRYOWRF coupler library (libcoupler.a)                   │
│  3. Test stage:                                                         │
│     └─ Run unit tests                                                   │
│  4. Integration stage:                                                  │
│     └─ Trigger CRYOWRF-WRF pipeline                                     │
│        (passes: CRYOWRF_COMMIT_SHA, CRYOWRF_BRANCH, CRYOWRF_PIPELINE_ID)│
└────────────────────────────────┬────────────────────────────────────────┘
                                 │
                                 │ trigger (strategy: depend)
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       CRYOWRF-WRF Repository                            │
├─────────────────────────────────────────────────────────────────────────┤
│  1. Dependencies stage:                                                 │
│     ├─ Build Jasper library for GRIB2 support                           │
│     └─ If triggered: checkout CRYOWRF_COMMIT_SHA + submodules           │
│  2. WRF stage:                                                          │
│     ├─ Compile WRF with CRYOWRF (em_real configuration)                 │
│     ├─ Generate wrf.exe and real.exe                                    │
│     └─ Create minimal artifacts (configure.wrf, I/O libs, frame objects)│
│  3. Integration stage:                                                  │
│     └─ Trigger CRYOWRF-WPS pipeline                                     │
│        (passes: WRF_PIPELINE_ID, WRF_COMMIT_SHA, WRF_BRANCH +           │
│                 CRYOWRF_COMMIT_SHA, CRYOWRF_BRANCH, CRYOWRF_PIPELINE_ID)│
└────────────────────────────────┬────────────────────────────────────────┘
                                 │
                                 │ trigger (strategy: depend)
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       CRYOWRF-WPS Repository                            │
├─────────────────────────────────────────────────────────────────────────┤
│  1. Prepare stage:                                                      │
│     ├─ If triggered: use WRF_PIPELINE_ID artifacts                      │
│     └─ If standalone: check for CRYOWRF-WRF artifacts from master build │
│  2. Build stage:                                                        │
│     ├─ Download WRF artifacts from CRYOWRF-WRF pipeline                 │
│     ├─ Configure WPS against WRF build                                  │
│     ├─ Compile WPS tools                                                │
│     └─ Verify binaries: geogrid.exe, metgrid.exe, ungrib.exe            │
└─────────────────────────────────────────────────────────────────────────┘
```

**CRYOWRF Repository Pipeline:**

- On every push or merge request:
  - Builds MeteoIO and SNOWPACK dependencies from submodules
  - Compiles the CRYOWRF coupler library
  - Runs some basic unit tests to verify coupler functionality
  - Triggers CRYOWRF-WRF pipeline for integration testing

**CRYOWRF-WRF Repository Pipeline:**

- On every push or merge request:
  - When triggered from CRYOWRF, checks out the specific CRYOWRF commit being
    tested along with its submodules
  - Compiles WRF with CRYOWRF integration using the em_real configuration
  - Generates `wrf.exe` and `real.exe` executables
  - Produces minimal WRF artifacts (`configure.wrf`, I/O libraries, frame
    objects) for downstream WPS build
  - Triggers CRYOWRF-WPS pipeline for full integration verification

**CRYOWRF-WPS Repository Pipeline:**

- On every push or merge request:
  - Downloads WRF artifacts from the upstream CRYOWRF-WRF build
  - For standalone builds, it verifies WRF artifacts exist from the current
    master build of CRYOWRF-WRF
  - For pipeline triggers, such as CRYOWRF -> CRYOWRF-WRF -> CRYOWRF-WPS, it
    uses a provided commit to checkout the correct WRF build
  - Configures and compiles WPS (`geogrid.exe`, `metgrid.exe`, `ungrib.exe`)
    against the WRF build and is marked as success if the binaries exist.

The pipeline uses GitLab's trigger mechanism with `strategy: depend` to wait
for each stage to complete, ensuring the full integration chain passes before
marking the original CRYOWRF commit as successful. Variables like
`CRYOWRF_COMMIT_SHA`, `CRYOWRF_BRANCH`, and `CRYOWRF_PIPELINE_ID` are passed
through the chain to track the integration test context.

The pipelines ensure that any changes to CRYOWRF, WRF, WPS, MeteoIO, or
SNOWPACK changes that affect the library are tested together, maintaining
compatibility across all components.

_Note_: Despite the build passing, there are no tests on processing or running
the models, therefore further validation is required before running the model.

## Attribution

### Related Repositories

- **CRYOWRF-WRF**: https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF
- **CRYOWRF-WPS**: https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WPS

### References

**Scientific Background:**

- Sharma, V., Gerber, F., and Lehning, M. (2023). _Introducing CRYOWRF v1.0:
  multiscale atmospheric flow simulations with advanced snow cover modelling_.
  **Geoscientific Model Development**, 16(2), 719–749.
  - [https://gmd.copernicus.org/articles/16/719/2023/](https://gmd.copernicus.org/articles/16/719/2023/)
  - DOI: [10.5194/gmd-16-719-2023](https://doi.org/10.5194/gmd-16-719-2023).

**Datasets:**

- Sharma, V., Gerber, F., and Lehning, M. (2021). _Reproducibility Dataset for
  CRYOWRF v1.0_. **EnviDat**.

  - https://www.envidat.ch/dataset/envidat_232
  - DOI: https://doi.org/10.16904/envidat.232

- Gerber, F., Sharma, V., and Lehning, M. (2022). _Reproducibility dataset for
  CRYOWRF validation_. **EnviDat**.
  - https://www.envidat.ch/dataset/reproducibility-dataset-for-cryowrf-validation
  - DOI: https://doi.org/10.16904/envidat.347

### Developed by

- [Laboratory of Cryospheric Sciences
  (CRYOS)](https://www.epfl.ch/labs/cryos/)

  EPFL Valais Campus, Sion,
  Switzerland

- [Snow Processes Group, Snow and Atmosphere Research
  Unit](https://www.slf.ch/en/about-the-slf/organisation/research-units/snow-and-atmosphere/snow-processes/)

  WSL Institute for Snow and Avalanche Research SLF, Davos, Switzerland

#### Original authors

- Varun Sharma - WSL SLF
- Franziska Gerber - WSL SLF

#### Current maintainers

- Evan Thomas - EPFL Valais
- Sergi Gonzàlez Herrero - WSL SLF
- Samuele Viaro - CRYOS

### Funding

- Swiss National Science Foundation: "From cloud to ground: Snow deposition in
  extreme environments" Grant ID: 200020-179130
- EPFL ENAC Big Call: Project LOSUMEA
- Swiss National Supercomputer Center (CSCS): Projects s873, s938, s1308

### License

CRYOWRF is released under GNU General Public License v3. See `LICENSE.txt`.
Each component (WRF, WPS, MeteoIO, SNOWPACK) retains its own license. See
respective directories for details.
