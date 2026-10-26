# Changelog

## release 1.0 - code version 3.2 - 2026-10-26

### Added

- `HPROF` keyword selecting the Stark profiles of hydrogen lines; see [README_INPUT.md](README_INPUT.md).
- Tabulated Stark profiles of Lemke (1997), Stehl&eacute; & Hutcheon (1999) and Stehl&eacute; & Fouquet (2010), compiled into the code.
- Corrected, area-normalised HLINOP profiles for lines and plasma conditions not covered by the tables.

### Changed

- Hydrogen lines use the Lemke (1997) tables by default (`HPROF=2`); `HPROF=0` corresponding to the original uncorrected HLINOP reproduces v3.1.

## code version 3.1

Development version described in [Kochukhov (2026)](https://arxiv.org/abs/2610.31763).
