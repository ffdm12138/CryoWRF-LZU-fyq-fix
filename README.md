# CryoWRF-LZU-fyq-fix

**Unofficial personal development repository** based on CryoWRF (WRF V4.7.1 + CryoWRF V1.2.4).
This is a personal mirror with local bug fixes — **not** the official WSL/SLF repository.

- Upstream projects:
  - WRF superproject (CRYOWRF-WRF): <https://gitlabext.wsl.ch/atmospheric-models/CRYOWRF-WRF>
  - CryoWRF coupler: <https://git.wsl.ch/atmospheric-models/CRYOWRF>
- Upstream base: WSL `origin/master` (`7b097e5a`); all local fixes sit on top of it.
- `phys/cryowrf` was a private-GitLab submodule; in this mirror it is **vendored into the tree**
  (monorepo) so the full source, including local fixes, is self-contained. It corresponds to the
  working tree of the validated build (`v1.2.2-22-gca8ac9a` + final coupler tweaks); the
  submodule's original Git history remains only on the source machine (`.git/modules/`).

## Local fixes (12 commits over upstream master)

| Commit | Description |
|--------|-------------|
| `dac6a94e` | fix: guard QZ0_SEA blending by producer capability |
| `799766a6` | fix: harden MYNN psi lookup inputs and bounds |
| `e8753661` | fix: guard snowpack output interval step count |
| `9f72401f` | fix: initialize snowpack diagnostics before sea-ice early exit |
| `3d46ac71` | fix: return MYNN sea-side diagnostics to surface driver |
| `55212317` | fix: pass ordered (I,J) layout to initialize_f |
| `7174cf42` | fix: call mark_layout_initialized after init_snowpack loop |
| `35eae0ea`–`181072a2` | cryowrf submodule lifecycle fixes incl. FINALIZED terminal-state fix (`ca8ac9a`), STATE_UNALLOCATED transition, NDEBUG build, gfortran-9 compat, ordered layout init |

## Licenses and copyright

All original copyright and license notices are preserved unmodified:

- **WRF** — `LICENSE.txt`: NCAR/UCAR public-domain notice (attribution requested).
- **CryoWRF** — `phys/cryowrf/LICENSE.txt`: GNU GPL v3.
- **MeteoIO** — `phys/cryowrf/external/meteoio/`: `License.txt`, `COPYING`, `COPYING.LESSER`.
- **SNOWPACK** — `phys/cryowrf/external/snowpack/`: `License.txt`, `COPYING`, `COPYING.LESSER`.
- Other external components keep their notices in-tree:
  `external/esmf_time_f90/LICENSE`, `external/fftpack/fftpack5/LICENSE`,
  `tools/manage_externals/LICENSE.txt`, `var/da/makedepf90-2.8.8/COPYING`.

Redistribution must comply with the GPL-family and per-component licenses above.

## Layout notes

- Standard WRF V4.7.1 source layout with CryoWRF extensions.
- Fully vendored monorepo: all four former submodules (`phys/cryowrf`,
  `phys/noahmp`, `phys/MYNN-EDMF`, `.ci/hpc-workflows`) are inlined into the
  tree — no `git submodule` setup needed.
- Build products, run outputs and personal data are excluded from this repository.

## Contact

Maintained by ffdm12138 for research use.

---

### WRF-ARW Modeling System  ###

We request that all new users of WRF please register. This allows us to better determine how to support and develop the model. Please register using this form:[https://www2.mmm.ucar.edu/wrf/users/download/wrf-regist.php](https://www2.mmm.ucar.edu/wrf/users/download/wrf-regist.php).

For an overview of the WRF modeling system, along with information regarding downloads, user support, documentation, publications, and additional resources, please see the WRF Model Users' Web Site: [https://www2.mmm.ucar.edu/wrf/users/](https://www2.mmm.ucar.edu/wrf/users/).
 
Information regarding WRF Model citations (including a DOI) can be found here: [https://www2.mmm.ucar.edu/wrf/users/citing_wrf.html](https://www2.mmm.ucar.edu/wrf/users/citing_wrf.html).

The WRF Model is open-source code in the public domain, and its use is unrestricted. The name "WRF", however, is a registered trademark of the University Corporation for Atmospheric Research. The WRF public domain notice and related information may be found here: [https://www2.mmm.ucar.edu/wrf/users/public.html](https://www2.mmm.ucar.edu/wrf/users/public.html).
