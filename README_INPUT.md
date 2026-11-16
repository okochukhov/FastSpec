# FastSpec input parameter file

A full synthesis run (`fastspec input.dat`) is controlled by a plain-text parameter file.
This page lists all keywords this file can contain.
It extends Table B.1 of [Kochukhov (2026)](https://doi.org/10.1051/0004-6361/202662504) and follows the current version of the code; where the two differ, this page is the reference.
Ready-to-run parameter files are in the `examples/` directory.

## Syntax

- Each setting is a `keyword=value` line. Keywords are case-insensitive.
- Lines starting with `#` or `!` are comments, and text after `#` or `!` on a keyword line is ignored.
- Several numbers in one value are separated by commas and/or spaces.
- `LINES`, `ABUNDANCE`, `ISOTOPE`, `BFIELD` and an explicit `MUS` list are blocks that may continue on the following lines. A block ends at an empty line or at the next keyword.
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
| `VMIC` | Microturbulent velocity (km/s) | 0.0 |
| `[M/H]` | Metallicity scaling of the default solar abundance pattern | -- |
| `ABUNDANCE` | Individual element abundances overriding those in the model header/default solar pattern | -- |
| `ISOTOPE` | Isotopic fractions overriding the default values | -- |
| `CONVOL` | *v* sin *i*, *v*<sub>mac</sub> (km/s) and Gaussian width (km/s) or resolving power | none |
| `BFIELD` | Magnetic field vector (G): three numbers for `BFRAME` 0 or 1, or 2*N*<sub>&mu;</sub>+1 numbers for `BFRAME` 2 | none |
| `BFRAME` | Field vector frame: stellar (0), observer (1), or a list of inclination/azimuth pairs (2) | 0 |
| `MUS` | Number of &mu; angles, or an explicit list of &mu; values | 7 |
| `OPACITY` | Precomputed line-opacity file to add | none |
| `PREFIX` | Prefix of the output file names | `fs` |
| `DISKQUAD` | Disk-integration quadrature for a number of &mu; angles: equal-area annuli (0), Gauss-Legendre (1), Gauss-Radau with a node at &mu; = 1 (2) | 1 |
| `LIMITRAM` | RAM limit (GB) for the large arrays | 16 |
| `EPSOPAC` | Line-to-continuum threshold &epsilon; for line inclusion | 10<sup>-4</sup> |
| `VWFUDGE` | Scaling factor of the Uns&ouml;ld van der Waals width for atomic lines without broadening parameter | 1.0 |
| `H2BRD` | Scaling factor of the H<sub>2</sub> contribution to the van der Waals broadening | 0.0 |
| `MOLBRD` | Apply Stark and van der Waals broadening to molecular lines (0/1) | 0 |
| `VWMOL` | Van der Waals broadening constant for all molecular lines | 0.0 |
| `HPROF` | Hydrogen line profiles: original HLINOP (0), corrected HLINOP (1), Lemke tables (2), Stehl&eacute; & Hutcheon tables (3) | 2 |
| `HYDROSTAT` | Re-solve the hydrostatic structure when the abundances differ from those of the model (0/1) | 0 |

## Details

### `WAVE`

    wave = wstart, wend, step, mode

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
| f2 | `fs_lineop.bin`, line opacity (not saved with polarised radiative transfer) | 0/1 |
| f3 | `fs_contop.bin`, continuum opacity | 0/1 |
| f4 | `fs_spec.bin`, intensity spectra at the &mu; angles (Stokes IQUV in the magnetic case) | 0/1 |
| f5 | `fs_flux.bin`, disk-integrated (and convolved) flux | 0/1 |

Line strengths and intensities need at least one &mu; angle, and the flux needs at least two.

The default prefix `fs` can be modified with the keyword `PREFIX` in the input file, e.g. `prefix=star1` writes `star1_flux.bin`. The prefix may include a directory, which must exist. The disk integration of a saved intensity file (see `DISKQUAD`) operates without an input parameter file and always writes `fs_flux.bin`.

### `LINES`

    lines =
    -1            path/to/VALD_H1.h5
    -1            path/to/VALD_atoms.h5
    -0.5          path/to/151Eu.lin
    0126.001056   path/to/56Fe1H.h5
    0.0909        path/to/13C16O.h5

Each entry is a species code followed by the file name. File paths are relative to the directory in which FastSpec is run.

- A code starting with `-` (e.g. `-1`) marks a multi-species atomic list: HDF5 from `vald2hdf.py`, or its magnetic `.zh5` variant (written with `vald2hdf.py --zeeman`). A file that is not HDF5 is read as an ASCII list, either in the `.z` format produced by `vald2z.py` or in the short format of a VALD "Extract Stellar" request. Both ASCII formats may also contain molecular lines.
- An atomic HDF5 list may carry an `ISOT` column with the isotope mass number of each line, written by `vald2hdf.py` for VALD lists whose gf values are not scaled by the solar isotopic fractions. The gf value of every line with a given isotope is then multiplied by the proto-solar fraction of that isotope given by [Lodders et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025SSRv..221...23L).
- A number after the `-` is a factor 0 &le; *s* &le; 1 by which all gf values of the list are multiplied (e.g. `-0.5`). A code `-` or `-1` leaves them unchanged. This factor can be also used for isotopic scaling, provided that the lines of each isotope are in a separate list and their gf values do not already include the isotopic fraction.
- A list containing only hydrogen lines is required by the hydrogen opacity routines. It must be the first entry and does not accept a gf factor other than 1.
- Any other code marks a single-species molecular HDF5 list produced by `mol2hdf.py`. The code sets the isotopic factor by which all gf values of the list are multiplied:

| Code | Isotopic factor |
|---|---|
| One or two characters without a decimal point (e.g. `0`) | Product of the default solar isotopic fractions for the species code stored in the file |
| Species code in the Turbospectrum convention, with two digits per atom before the point and three after (e.g. `0126.001056` for <sup>56</sup>FeH) | Product of the default solar isotopic fractions of the isotopes in the code. An isotope `000` contributes 1, so `0126.000000` gives *s* = 1 resulting in no gf scaling |
| Decimal number with one digit before the point (e.g. `0.0909`) | This number, 0 &le; *s* &le; 1. A zero factor is written `0.0`, since `0` selects the species code stored in the file |

The default isotopic fractions are those of [Lodders et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025SSRv..221...23L); `ISOTOPE` changes them. A species code naming an isotope absent from this table stops the run.

### `MODEL`

The model atmosphere in `.krz` format. Plane-parallel and spherical models are supported; a spherical model is recognised from its header.

A layer may end with depth-dependent abundances in the format of `ABUNDANCE`, each enclosed in quotes, which replace the header abundances of these elements in that layer:

    4.69674E-08,  4222.5, 7.60800E+06, 1.26921E+09, 2.85327E-15,'Ca:-6.50','Fe:-7.0'

An element given in some layers only keeps its header abundance in the others.

### `MUS`

- An integer *N* (no decimal point) gives *N* angles placed according to `DISKQUAD`; `mus=1` is the disk centre.
- A value with a decimal point starts an explicit list of &mu; values (0 < &mu; &le; 1), e.g. `mus=1.0, 0.5, 0.2`. The list may continue on the following lines. A single off-centre angle is written as e.g. `mus=0.5`. The flux is integrated over annuli whose boundaries lie halfway in projected area between adjacent angles, each weighted by its area. An explicit list of more than one value requires `DISKQUAD=0`.
- `mus=0` skips the radiative transfer; only line and continuum opacities (`SAVE` flags 2 and 3) can then be saved.

### `DISKQUAD`

This keyword selects the angles and weights used to integrate the flux over the stellar disk when `MUS` gives a number of angles *N*.
The angles are set after the whole file has been read, so `DISKQUAD` may appear before or after `MUS`.

| `DISKQUAD` | &mu; angles and flux weights |
|---|---|
| 0 | Equal-area annuli: &mu;<sub>i</sub> = [1 - (i - 0.5)/*N*]<sup>1/2</sup>, i = 1...*N*, each with weight 1/*N* |
| 1 (default) | Gauss-Legendre nodes &mu;<sub>i</sub> on [0, 1] with weights *w*<sub>i</sub>; the flux weight is 2*w*<sub>i</sub>&mu;<sub>i</sub> |
| 2 | Gauss-Radau nodes on [0, 1] including &mu; = 1, with flux weights 2*w*<sub>i</sub>&mu;<sub>i</sub> |

For a single angle (`mus=1`) all three give the disk centre.
For rotational broadening, each angle represents an annulus enclosing the projected area of its flux weight.

When a saved intensity file is disk-integrated afterwards (`fastspec fs_spec.bin Vsini Vmacro Vgau`), the quadrature is recognised from the stored &mu; values: Gauss-Legendre and Gauss-Radau sets receive their own weights, and any other set is integrated over annuli as for `DISKQUAD=0`.

### `CONVOL`

    convol = Vsini, Vmacro, Vgau

The projected rotational velocity and the radial-tangential macroturbulence in km/s, and the Gaussian instrumental FWHM in km/s. A `Vgau` value larger than 1000 is taken as the resolving power *R*.
When the flux is saved, all three broadenings are applied to the disk-integrated flux.
When only intensity spectra are saved (`SAVE` flag 4 without flag 5), they are convolved with the Gaussian instrumental profile alone.

### `BFIELD` and `BFRAME`

Polarised radiative transfer is performed when `BFIELD` is present (even with all components zero) and at least one line list is magnetic. Lines from `.zh5` and `.z` lists have anomalous Zeeman splitting patterns, whereas lines from a VALD short-format list are split as normal triplets with the effective Land&eacute; factor.

| `BFRAME` | `BFIELD` values |
|---|---|
| 0 | Field components (*B*<sub>r</sub>, *B*<sub>m</sub>, *B*<sub>p</sub>) in G, in the stellar spherical coordinate system with *B*<sub>r</sub> aligned with the line of sight |
| 1 | \|*B*\| (G), inclination to the line of sight and azimuth (degrees), the same at every disk position |
| 2 | \|*B*\| (G) followed by one (inclination, azimuth) pair per &mu; angle, in the order of the `MUS` list: 2*N*<sub>&mu;</sub>+1 values, which may continue on the following lines |

### `ABUNDANCE` and `[M/H]`

    abundance = Cr:-6.41, Fe:-4.57

Abundances are given as log(*N*<sub>el</sub>/*N*<sub>tot</sub>), the convention of the `.krz` header. The quotes `'Fe: -4.57'` are optional, and the block may continue on the following lines.

Without either keyword, the abundances of the model header are used.
`[M/H]` replaces them by the solar abundances of [Asplund et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...653A.141A), with all elements from Li onward scaled by [M/H].
Since keywords are processed in order, `[M/H]` should precede `ABUNDANCE`; an `[M/H]` placed after `ABUNDANCE` overwrites it.
Both keywords set the same abundance in all layers, replacing any depth-dependent abundances of the model.

### `ISOTOPE`

    isotope = 13C:0.1, 18O:0.004

Each entry gives the mass number and element symbol of an isotope and its atom fraction 0 &le; *s* &le; 1. The quotes `'13C: 0.1'` are optional, and the block may continue on the following lines.

- The isotopes of an element that are not given share the rest of the element in proportion to their current fractions; e.g. `13C:0.1` sets <sup>12</sup>C to 0.9.
- The fractions given for an element may not sum to more than 1. If all isotopes of an element are given, their sum must be 1 within 10<sup>-3</sup>.
- An isotope absent from the table of [Lodders et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025SSRv..221...23L) stops the run.
- `ISOTOPE` may appear anywhere in the file; a second `ISOTOPE` applies on top of the first.
- The new fractions apply to molecular lists with a species code (including `0`) and to atomic HDF5 lists with an `ISOT` column. They do not change an explicit molecular factor (a warning is printed if its isotopes were changed), an isotope `000`, or atomic lists without `ISOT`.
- The changed isotopic fractions are printed at the start of the run.

### `OPACITY`

A line-opacity file (`fs_lineop.bin`) from an earlier run on the same wavelength grid and model depth points. Its opacity is added to the line opacity computed in this run.

### `LIMITRAM`, `EPSOPAC`, `VWFUDGE`

- `LIMITRAM`: memory budget (GB, at least 2) for the large arrays; longer line lists are processed in blocks to stay within it.
- `EPSOPAC`: lines whose opacity stays below &epsilon; times the continuum opacity are skipped. Allowed range 10<sup>-6</sup> to 0.1.
- `VWFUDGE`: multiplies the van der Waals width of lines computed with the Uns&ouml;ld approximation. Allowed range 0.1 to 100.

### `H2BRD`, `MOLBRD`, `VWMOL`

- `H2BRD`: multiplies the contribution of molecular hydrogen to the perturbers of the van der Waals broadening; 0 omits H<sub>2</sub>. Allowed range 0 to 10.
- `MOLBRD=1` applies Stark and van der Waals broadening to molecular lines, when these constants are given in the line list.
- `VWMOL`: van der Waals constant applied to all molecular lines, in the line-list convention: log<sub>10</sub>&gamma;<sub>6</sub> if negative, &gamma;<sub>6</sub> if positive. A non-zero value replaces the line-list constants and applies with or without `MOLBRD=1`. Values above 10 are not accepted.

### `HPROF`

This keyword selects the profiles of the hydrogen lines and matters only when a hydrogen line list is given.
The level populations, occupation probabilities and bound-free opacity do not depend on it.

| `HPROF` | Stark profiles of the hydrogen lines |
|---|---|
| 0 | HLINOP profile routine by P. Barklem, <https://github.com/barklem/hlinop> |
| 1 | Corrected HLINOP for all lines: the quasistatic profile is convolved with the electron-impact Lorentzian and normalised to unit area |
| 2 (default) | [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L) tables (Lyman to Brackett series, upper level up to 22, except Br11 and Br12) where they cover the plasma conditions; corrected HLINOP otherwise |
| 3 | [Stehl&eacute; & Hutcheon (1999)](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S) tables (Lyman, Balmer and Paschen series, upper level up to 30) and [Stehl&eacute; & Fouquet (2010)](https://ui.adsabs.harvard.edu/abs/2010IJSp.2010E...1S) tables (Br&alpha;, Br&beta; and Br&gamma;) where they cover the plasma conditions, then [Lemke (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L), then corrected HLINOP |

For `HPROF` 1 to 3, the profile of each line at each depth is the convolution of the Stark profile with the Doppler profile (thermal and microturbulent, including fine structure) and the Lorentz profile (natural, self- and helium broadening). Lyman &alpha; keeps the original HLINOP treatment in all modes.

The tables are compiled into the code from the CDS catalogues [J/A+AS/122/285](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+AS/122/285) ([Lemke 1997](https://ui.adsabs.harvard.edu/abs/1997A%26AS..122..285L)) and [VI/98A](https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/98A) ([Stehl&eacute; & Hutcheon 1999](https://ui.adsabs.harvard.edu/abs/1999A%26AS..140...93S); with the 2007 correction of the Paschen data), and from the [STARK-H](https://stark-h.obspm.fr/) database ([Stehl&eacute; & Fouquet 2010](https://ui.adsabs.harvard.edu/abs/2010IJSp.2010E...1S)).

### `HYDROSTAT`

With `HYDROSTAT=1` and abundances that differ from those of the model (set by `ABUNDANCE` or `[M/H]`), the model structure is adjusted to the new abundances.
Each layer keeps its temperature and its continuum optical depth at the standard wavelength of the model (`WLSTD` in the model header), and the gas pressure keeps its dependence on column mass.
The column mass of each layer is found such that the optical depth computed with the new abundances reproduces that of the model with its own abundances, and the equation of state is solved at the resulting gas pressures. Depth-dependent abundances stay with their layers.
The adjusted structure is used throughout the synthesis and is the one written to the opacity output files.
Without an abundance change the keyword has no effect. Only plane-parallel models on a column-mass scale (`MODEL TYPE= 0`) are supported.
