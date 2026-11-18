# Changelog

## changes since release 1.2

### Changed

- Fixed the `h5open_f failed` error when reading HDF5 line lists with HDF5 library versions 1.10.10, 1.10.11, 1.12.3 and 1.14.0.

## release 1.2 - code version 3.4 - 2026-11-16

### Added

- NH photodissociation continuum opacity (LTE cross-sections of P. C. Stancil, as distributed with Turbospectrum, smoothed with a 10 A Gaussian).
- Depth-dependent abundances in the model atmosphere layers (see `MODEL` in README_INPUT.md).
- Optional explicit isotopic factor in place of the species code for single-species molecular line lists (see `LINES` in README_INPUT.md).
- gf factor for atomic line lists, e.g. for isotopic scaling (see `LINES` in README_INPUT.md).
- Automatic isotopic scaling of gf values in atomic HDF5 line lists with an `ISOT` column (see `LINES` in README_INPUT.md).
- User-controlled `ISOTOPE` keyword changing the isotopic fractions used for the line lists with isotopic information (see README_INPUT.md).
- `PREFIX` keyword replacing the default `fs` prefix of the output file names (see `SAVE` in README_INPUT.md).

### Changed

- Faster scalar radiative transfer.
- An unsupported ASCII line list in VALD long format now stops the run with an error instead of being misread.
- Default isotopic fractions are adopted from [Lodders et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025SSRv..221...23L), corresponding to proto-solar values.
- Vacuum wavelengths of line lists are converted to air with the formula from [Birch & Downs (1994)](https://ui.adsabs.harvard.edu/abs/1994Metro..31..315B).
- Continuous opacity is computed at the vacuum wavelengths of the air wavelength grid.
- Faster hydrogen line opacity calculation and more accurate interpolation of the hydrogen line profiles for `HPROF` 1 to 3.
- Hydrogen bound-free cross-sections from [Karzas & Latter (1961)](https://ui.adsabs.harvard.edu/abs/1961ApJS....6..167K) Gaunt factors are used for the continuum opacity and for the occupation probability treatment of the hydrogen line opacities. The `BFKURUCZ` keyword, which selected the less accurate ATLAS9 cross-sections for the latter, is removed.
- Dissolved hydrogen bound-free opacity is cut redward of the alpha line of each series instead of tapered.
- Fixed overestimated hydrogen level populations in the occupation probability treatment. This increases the strength of H lines for hot stars, bringing them in agreement with other codes.
- Shorter output with default `VERBOSE=1`.

## release 1.1 - code version 3.3 - 2026-11-03

### Added

- `DISKQUAD` keyword selecting the disk-integration quadrature: equal-area annuli (reproduces v3.2), Gauss-Legendre (default) or Gauss-Radau.
- Disk integration of saved intensity files recognises Gauss-Legendre and Gauss-Radau &mu; sets.
- `HYDROSTAT` keyword for re-solving the hydrostatic structure at fixed T(&tau;) for abundances that differ from the model atmosphere.
- `VWMOL` keyword setting the van der Waals constant of all molecular lines, independently of `MOLBRD`.
- H<sub>2</sub>-H collision-induced absorption ([Gustafsson & Frommhold 2003](https://ui.adsabs.harvard.edu/abs/2003A%26A...400.1161G)).
- H<sub>2</sub><sup>-</sup> free-free absorption ([Bell 1980](https://ui.adsabs.harvard.edu/abs/1980JPhB...13.1859B)).

### Changed

- `H2BRD` is a scaling factor (0 to 10) of the H<sub>2</sub> van der Waals perturber contribution instead of an on/off flag.
- H<sub>2</sub>-H<sub>2</sub> and H<sub>2</sub>-He collision-induced absorption updated to [Borysow et al. (2001)](https://ui.adsabs.harvard.edu/abs/2001JQSRT..68..235B) and [J&oslash;rgensen et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000A%26A...361..283J); stimulated emission is no longer applied twice.
- Rotational broadening kernel of each disk-integration annulus integrated over the velocity pixels instead of sampled at their centres.

## release 1.0 - code version 3.2 - 2026-10-26

### Added

- `HPROF` keyword selecting the Stark profiles of hydrogen lines.
- Tabulated Stark profiles of [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L), [Stehl&eacute; & Hutcheon (1999)](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S) and [Stehl&eacute; & Fouquet (2010)](https://ui.adsabs.harvard.edu/abs/2010IJSp.2010E...1S), compiled into the code.
- Corrected, area-normalised HLINOP profiles for lines and plasma conditions not covered by the tables.

### Changed

- Hydrogen lines use the [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L) tables by default (`HPROF=2`); `HPROF=0` corresponding to the original uncorrected HLINOP reproduces v3.1.

## code version 3.1

Development version described in [Kochukhov (2026)](https://doi.org/10.1051/0004-6361/202662504).
