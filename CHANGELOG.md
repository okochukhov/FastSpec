# Changelog

## release 1.1 - code version 3.3 - 2026-11-03

### Added

- `DISKQUAD` keyword selecting the disk-integration quadrature: equal-area annuli (reproduces v3.2), Gauss-Legendre (default) or Gauss-Radau.
- Disk integration of saved intensity files recognises Gauss-Legendre and Gauss-Radau &mu; sets.
- `HYDROSTAT` keyword for re-solving the hydrostatic structure at fixed T(&tau;) for abundances that differ from the model atmosphere.
- `VWMOL` keyword setting the van der Waals constant of all molecular lines, independently of `MOLBRD`.
- H<sub>2</sub>-H collision-induced absorption ([Gustafsson & Frommhold 2003](https://ui.adsabs.harvard.edu/abs/2003A%26A...400.1161G/abstract)).
- H<sub>2</sub><sup>-</sup> free-free absorption ([Bell 1980](https://ui.adsabs.harvard.edu/abs/1980JPhB...13.1859B/abstract)).

### Changed

- `H2BRD` is a scaling factor (0 to 10) of the H<sub>2</sub> van der Waals perturber contribution instead of an on/off flag.
- H<sub>2</sub>-H<sub>2</sub> and H<sub>2</sub>-He collision-induced absorption updated to [Borysow et al. (2001)](https://ui.adsabs.harvard.edu/abs/2001JQSRT..68..235B/abstract) and [J&oslash;rgensen et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000A%26A...361..283J/abstract); stimulated emission is no longer applied twice.
- Rotational broadening kernel of each disk-integration annulus integrated over the velocity pixels instead of sampled at their centres.

## release 1.0 - code version 3.2 - 2026-10-26

### Added

- `HPROF` keyword selecting the Stark profiles of hydrogen lines.
- Tabulated Stark profiles of [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L/abstract), [Stehl&eacute; & Hutcheon (1999)](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S/abstract) and [Stehl&eacute; & Fouquet (2010)](https://ui.adsabs.harvard.edu/abs/2010IJSp.2010E...1S/abstract), compiled into the code.
- Corrected, area-normalised HLINOP profiles for lines and plasma conditions not covered by the tables.

### Changed

- Hydrogen lines use the [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L/abstract) tables by default (`HPROF=2`); `HPROF=0` corresponding to the original uncorrected HLINOP reproduces v3.1.

## code version 3.1

Development version described in [Kochukhov (2026)](https://ui.adsabs.harvard.edu/abs/2026arXiv261031763K/abstract).
