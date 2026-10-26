# FastSpec input parameter file

A full synthesis run (`fastspec input.dat`) is controlled by a plain-text parameter file.
This page lists all keywords this file can contain.
It extends Table B.1 of [Kochukhov (2026)](https://arxiv.org/abs/2610.31763) and follows the current version of the code; where the two differ, this page is the reference.
Ready-to-run parameter files are in the `examples/` directory.

## Syntax

- Each setting is a `keyword=value` line. Keywords are case-insensitive.
- Lines starting with `#` or `!` are comments, and text after `#` or `!` on a keyword line is ignored.
- Several numbers in one value are separated by commas and/or spaces.
- `LINES`, `ABUNDANCE`, `BFIELD` and an explicit `MUS` list are blocks that may continue on the following lines. A block ends at an empty line or at the next keyword.
- Keywords are processed in the order they appear. `ABUNDANCE` and `[M/H]` must come after `MODEL`.
- An unknown keyword stops the run with an error message.

## Keywords

| Keyword | Meaning | Default |
|---|---|---|
| **Mandatory** | | |
| `WAVE` | Start and end wavelength, step, and mode (0: constant in wavelength, 1: constant in velocity) | -- |
| `SAVE` | Five flags selecting the outputs: line-strength list, line opacity, continuum opacity, intensity, flux | -- |
| `LINES` | Block of line-list entries, one per line: species code and file name | -- |
| `MODEL` | Model-atmosphere file name (`.krz`) | -- |
| **Optional** | | |
| `MUS` | Number of &mu; angles, or an explicit list of &mu; values | 7 |
| `VMIC` | Microturbulent velocity (km/s) | 0.0 |
| `CONVOL` | *v* sin *i*, *v*<sub>mac</sub> (km/s) and Gaussian width (km/s) or resolving power | none |
| `BFIELD` | Magnetic field vector (G): three numbers for `BFRAME` 0 or 1, or 2*N*<sub>&mu;</sub>+1 numbers for `BFRAME` 2 | none |
| `BFRAME` | Field vector frame: stellar (0), observer (1), or a list of inclination/azimuth pairs (2) | 0 |
| `[M/H]` | Metallicity scaling of the default solar abundance pattern | -- |
| `ABUNDANCE` | Individual element abundances overriding those in the model header/default solar pattern | -- |
| `OPACITY` | Precomputed line-opacity file to add | none |
| `LIMITRAM` | RAM limit (GB) for the large arrays | 16 |
| `EPSOPAC` | Line-to-continuum threshold &epsilon; for line inclusion | 10<sup>-4</sup> |
| `VWFUDGE` | Scaling factor of the Uns&ouml;ld van der Waals width for atomic lines without broadening parameter | 1.0 |
| `H2BRD` | Include H<sub>2</sub> in the van der Waals broadening (0/1) | 0 |
| `MOLBRD` | Apply Stark and van der Waals broadening to molecular lines (0/1) | 0 |
| `BFKURUCZ` | Choice of hydrogen bound-free cross-sections: ATLAS (1) or MARCS (0) | 1 |
| `HPROF` | Hydrogen line profiles: original HLINOP (0), corrected HLINOP (1), Lemke tables (2), Stehl&eacute; & Hutcheon tables (3) | 2 |

## Details

### `WAVE`

    wave = Wstart, Wend, step, mode

Wavelengths are in &Aring; (air).
With `mode=0` the step is a constant wavelength increment in &Aring;.
With `mode=1` it is a constant velocity increment in km/s, which gives a grid uniform in ln &lambda;.
The grid starts at `Wstart` and ends at the last grid point that does not exceed `Wend`; it must contain at least two points.

### `SAVE`

    save = f1, f2, f3, f4, f5

Exactly five flags are required:

| Flag | Output file | Values |
|---|---|---|
| f1 | `fs_linst.bin`, central line depths of the retained lines | 0: off; 1: output for atomic and magnetic line lists; 2: in addition, output for single-species molecular lists |
| f2 | `fs_lineop.bin`, line opacity | 0/1 |
| f3 | `fs_contop.bin`, continuum opacity | 0/1 |
| f4 | `fs_spec.bin`, intensity spectra at the &mu; angles (Stokes IQUV in the magnetic case) | 0/1 |
| f5 | `fs_flux.bin`, disk-integrated (and convolved) flux | 0/1 |

Line strengths and intensities need at least one &mu; angle, and the flux needs at least two.

### `LINES`

    lines =
    -1            path/to/VALD_H1.h5
    -1            path/to/VALD_atoms.h5
    0126.001056   path/to/56Fe1H.h5

Each entry is a species code followed by the file name. File paths are relative to the directory in which FastSpec is run.

- A code starting with `-` (e.g. `-1`) marks a multi-species atomic list: HDF5 from `vald2hdf.py`, or its magnetic `.zh5` variant (`vald2hdf.py --zeeman`). A file that is not HDF5 is read as an ASCII list in the `.z` format produced by `vald2z.py`, which may also contain molecular lines.
- A list containing only hydrogen lines is required by the hydrogen opacity routines. It must be the first entry.
- Any other code marks a single-species molecular HDF5 list produced by `mol2hdf.py`. The code is in the Turbospectrum convention (e.g. `0126.001056` for <sup>56</sup>FeH) and sets the isotopic fraction applied to the gf values. A code of one or two characters (e.g. `0`) uses the species code stored in the file instead.

### `MODEL`

The model atmosphere in `.krz` format. Plane-parallel and spherical models are supported; a spherical model is recognised from its header.

### `MUS`

- An integer *N* (no decimal point) gives *N* angles &mu;<sub>i</sub> = [1 - (i - 0.5)/*N*]<sup>1/2</sup>, i = 1...*N*; `mus=1` is the disk centre.
- A value with a decimal point starts an explicit list of &mu; values (0 < &mu; &le; 1), e.g. `mus=1.0, 0.5, 0.2`. The list may continue on the following lines. A single off-centre angle is written as e.g. `mus=0.5`.
- `mus=0` skips the radiative transfer; only line and continuum opacities (`SAVE` flags 2 and 3) can then be saved.

### `CONVOL`

    convol = Vsini, Vmacro, Vgau

The projected rotational velocity and the radial-tangential macroturbulence in km/s, and the Gaussian instrumental FWHM in km/s. A `Vgau` value larger than 1000 is taken as the resolving power *R*.
When the flux is saved, all three broadenings are applied to the disk-integrated flux.
When only intensity spectra are saved (`SAVE` flag 4 without flag 5), they are convolved with the Gaussian instrumental profile alone.

### `BFIELD` and `BFRAME`

Polarised radiative transfer is performed when `BFIELD` is present (even with all components zero) and at least one line list is magnetic (`.zh5` or `.z` line list).

| `BFRAME` | `BFIELD` values |
|---|---|
| 0 | Field components (*B*<sub>r</sub>, *B*<sub>m</sub>, *B*<sub>p</sub>) in G, in the stellar spherical coordinate system with *B*<sub>r</sub> aligned with the line of sight |
| 1 | \|*B*\| (G), inclination to the line of sight and azimuth (degrees), the same at every disk position |
| 2 | \|*B*\| (G) followed by one (inclination, azimuth) pair per &mu; angle, in the order of the `MUS` list: 2*N*<sub>&mu;</sub>+1 values, which may continue on the following lines |

### `ABUNDANCE` and `[M/H]`

    abundance = Cr:-6.41, Fe:-4.57

Abundances are given as log(*N*<sub>el</sub>/*N*<sub>tot</sub>), the convention of the `.krz` header. The quotes `'Fe: -4.57'` are optional, and the block may continue on the following lines.

Without either keyword, the abundances of the model header are used.
`[M/H]` replaces them by the solar abundances of [Asplund et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...653A.141A/abstract), with all elements from Li onward scaled by [M/H].
Since keywords are processed in order, `[M/H]` should precede `ABUNDANCE`; an `[M/H]` placed after `ABUNDANCE` overwrites it.

### `OPACITY`

A line-opacity file (`fs_lineop.bin`) from an earlier run on the same wavelength grid and model depth points. Its opacity is added to the line opacity computed in this run.

### `LIMITRAM`, `EPSOPAC`, `VWFUDGE`

- `LIMITRAM`: memory budget (GB, at least 2) for the large arrays; longer line lists are processed in blocks to stay within it.
- `EPSOPAC`: lines whose opacity stays below &epsilon; times the continuum opacity are skipped. Allowed range 10<sup>-6</sup> to 0.1.
- `VWFUDGE`: multiplies the van der Waals width of lines computed with the Uns&ouml;ld approximation. Allowed range 0.1 to 100.

### `H2BRD`, `MOLBRD`, `BFKURUCZ`

- `H2BRD=1` adds molecular hydrogen to the perturbers of the van der Waals broadening.
- `MOLBRD=1` applies Stark and van der Waals broadening to molecular lines, when these constants are given in the line list.
- `BFKURUCZ` selects the hydrogen bound-free cross-sections used with the occupation-probability treatment of hydrogen (ATLAS: 1, MARCS: 0). It matters only when a hydrogen line list is given.

### `HPROF`

This keyword selects the profiles of the hydrogen lines and matters only when a hydrogen line list is given.
The level populations, occupation probabilities and bound-free opacity do not depend on it.

| `HPROF` | Stark profiles of the hydrogen lines |
|---|---|
| 0 | HLINOP profile routine by P. Barklem, <https://github.com/barklem/hlinop> |
| 1 | Corrected HLINOP for all lines: the quasistatic profile is convolved with the electron-impact Lorentzian and normalised to unit area |
| 2 (default) | [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L/abstract) tables (Lyman to Brackett series, upper level up to 22, except Br11 and Br12) where they cover the plasma conditions; corrected HLINOP otherwise |
| 3 | [Stehl&eacute; & Hutcheon (1999)](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S/abstract) tables (Lyman, Balmer and Paschen series, upper level up to 30) and [Stehl&eacute; & Fouquet (2010)](https://ui.adsabs.harvard.edu/abs/2010IJSp.2010E...1S/abstract) tables (Br&alpha;, Br&beta; and Br&gamma;) where they cover the plasma conditions, then [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L/abstract), then corrected HLINOP |

For `HPROF` 1 to 3, the profile of each line at each depth is the convolution of the Stark profile with the Doppler profile (thermal and microturbulent, including fine structure) and the Lorentz profile (natural, self- and helium broadening). Lyman &alpha; keeps the original HLINOP treatment in all modes.

The tables are compiled into the code from the CDS catalogues [J/A+AS/122/285](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+AS/122/285) ([Lemke 1997](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L/abstract)) and [VI/98A](https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/98A) ([Stehl&eacute; & Hutcheon 1999](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S/abstract); with the 2007 correction of the Paschen data), and from the [STARK-H](https://stark-h.obspm.fr/) database ([Stehl&eacute; & Fouquet 2010](https://ui.adsabs.harvard.edu/abs/2010IJSp.2010E...1S/abstract)).
