# FastSpec

[FastSpec](https://arxiv.org/abs/2610.31763) is a fast and versatile stellar spectrum synthesis code, written in modern Fortran and designed for both raw computational speed and broad astrophysical applicability.
Working in local thermodynamic equilibrium (LTE), it computes emergent intensity and flux spectra from 1D stellar model atmospheres, and can solve the full polarised radiative transfer problem, yielding all four Stokes parameters (I, Q, U, V) in the presence of a magnetic field.

The code is applicable across essentially the entire Hertzsprung-Russell diagram, from cool M dwarfs to hot B-type stars.
It treats both plane-parallel and spherical model atmospheres, and synthesises atomic as well as molecular lines -- handling molecular line lists that reach hundreds of millions of transitions, as needed to model the spectra of cool stars.
Chemical and ionisation equilibrium is obtained from the corrected and extended SME equation-of-state package; continuum opacities follow ATLAS12; and hydrogen bound-bound and bound-free absorption use the occupation-probability formalism.
The radiative transfer equation is integrated with DELO-Bezier formal solvers in both scalar and Stokes-vector form.

Speed comes from aggressive OpenMP parallelisation, highly optimised adaptive per-line wavelength windowing, new efficient implementations of the Voigt and Faraday-Voigt profile algorithms, and a memory-capped block-processing scheme that lets arbitrarily large line lists be handled within a user-defined RAM budget.
On a modern multi-core laptop, FastSpec can synthesise cool-dwarf spectra from hundreds of millions of molecular lines over the full optical-infrared range in few minutes. Calculating the full 0.3 to 2.5 micron spectrum of a hotter solar-type star at 10^6 resolving power requires on the order of 10 seconds.

## Requirements

- A Fortran compiler with HDF5 support. The Makefile uses `h5fc` (the HDF5
  gfortran wrapper) for the main code and `gfortran` for the equation-of-state
  routines; OpenMP is used for parallel calculations.
- For preparing line lists and viewing results: Python3 with `numpy` and `matplotlib`.

## Compiling

Precompiled executables are included in `src/`: `fastspec_Darwin_arm64` for macOS on Apple Silicon and `fastspec_Linux_x86_64` for Linux on 64-bit Intel/AMD processors -- if one of them matches your platform you may be able to run it directly and skip this section.
Otherwise, you can rebuild from source by running `make` in `src/`.
This requires you to have the gfortran compiler with HDF5 and OpenMP support.

## Running

FastSpec has two modes of operation, selected by the number of command-line arguments.

**Full synthesis** (one argument) is driven by a plain-text control file:

    fastspec input.dat

The control file sets the wavelength range and sampling, the microturbulent velocity, the line lists, the model atmosphere, chemical abundances, and the number of limb angles. It also selects which outputs to save.
All keywords of the control file are described in [README_INPUT.md](README_INPUT.md).
See the `examples/` directory for ready-to-run setups.

Depending on the `save=` flags, FastSpec writes binary output files such as:

- `fs_spec.bin`  -- intensity spectra at the requested limb angles
- `fs_flux.bin`  -- disk-integrated flux spectrum
- `fs_linst.bin` -- per-line strengths

**Disk integration and convolution** (four arguments) reuses a pre-computed intensity file, performing only the disk-integration and convolution steps:

    fastspec fs_spec.bin <Vsini> <Vmacro> <Vgau|R>

Here `fs_spec.bin` is a multi-angle intensity spectrum saved by a previous full-synthesis run, and the three numbers are the projected rotational velocity `Vsini` (km/s), the radial-tangential macroturbulence `Vmacro` (km/s), and the Gaussian instrumental broadening -- given either as a FWHM in km/s or as a resolving power.
This mode writes only the broadened flux spectrum (`fs_flux.bin`).

## Input data

FastSpec needs two kinds of input, both referenced from the control file:

**Model atmospheres** are read in the [.krz format](https://marcs.oreme.org/documents/krzformat), whose header gives the effective temperature, surface gravity, chemical composition and depth structure.
The format is compatible with ATLAS9/12, MARCS and LLmodels atmospheres.
MARCS model atmospheres in this format can be downloaded from [MARCS homepage](https://marcs.oreme.org).
Both plane-parallel and spherical models are supported; a spherical model is recognised from its header and additionally carries the stellar radius and the geometric height at each depth.

**Line lists** are read in HDF5 format, for both atomic and molecular species, plus a dedicated magnetic format (`.z` / `.zh5`) that carries the Zeeman-splitting data required for polarised synthesis.
The Python tools in `python/` prepare these files: `vald2hdf.py` converts VALD atomic extractions to HDF5 (with `--zeeman` it writes the magnetic `.zh5` variant); `vald2z.py` writes the alternative magnetic `.z` format; and `mol2hdf.py` converts molecular lists in Turbospectrum format to single-species HDF5.

Line lists in VALD `extract stellar` short-format are also directly supported but are less efficient in terms of disk storage and computational speed and therefore are not recommended.

## Viewing results

The `python/` directory contains a reader and a plotting demo:

- `readfs.py`       -- functions `read_spec`, `read_flux`, `read_linst`
- `example_read.py` -- example that reads and plots the output files

## Author, citation and acknowledgements

FastSpec is developed by [Oleg Kochukhov](https://www.astro.uu.se/~oleg) (Department of Physics and Astronomy, Uppsala University, Sweden).

If you use FastSpec in your research, please cite [Kochukhov (2026)](https://arxiv.org/abs/2610.31763).

FastSpec incorporates, in modified form, several external packages, whose authors are gratefully acknowledged:

- **Equation of state**: the EOS routines from the SME package by
  N. Piskunov, J. Valenti, and P. Barklem -- <https://github.com/SpectroscopyMadeEasy/SMElib>
  (bug fixes, modified molecular data)
- **Hydrogen line and bound-free opacity**: the HLINOP / HBOP package by
  P. Barklem -- <https://github.com/barklem/hlinop>
  (re-written in Fortran90, optimised for speed, with corrected line profiles)
- **Continuum opacity**: routines from the ATLAS12 code by 
  R. Kurucz -- <http://kurucz.harvard.edu/programs/atlas12>
  (re-written in Fortran90 and optimised for speed)
