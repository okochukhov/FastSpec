#!/usr/bin/env python3.12
"""
mol2hdf.py -- convert a molecular line list in Turbospectrum ASCII format
into the HDF5 file read by FastSpec.

    group /Lines : one 1-D dataset per column {W:f8, E:f4, LGGF:f4, GU:f4,
                   GRAD:f4}, all sorted in wavelength.  Payload columns are
                   chunked + shuffled + gzip-compressed; W is left contiguous
                   (uncompressed) so the reader's full-column bound scan is a
                   fast straight read.  Each column dataset carries a UNIT
                   attribute.  Group attributes: NL (int32), ID, NAME, FILE
                   (null-terminated ASCII), ION, WL_VAC (int16), DATE
                   (creation date YYYY-MM-DD), PROVENANCE (input file name, or
                   --provenance text), FORMAT_VERSION, and WMIN/WMAX (f8, the
                   file's wavelength range, for fast file-level culling).

The input list must contain a single species.  It may be gzip-, bzip2- or
xz-compressed; compression is detected from the file content, not from the
file name.

Format notes
    Header line 1: "'species_id'  ion  nlines"; line 2 is a free-form
    comment.  Each data line carries at least 6 numbers: wavelength (A),
    E_low, log gf, fdamp (unused), g_up, gamma_rad.  A zero radiative
    damping is replaced by the classical estimate
        gamma_rad = 6.669e15 * 10**loggf / (g_up * wl**2).
    The molecule name stored in the NAME attribute must be given with
    molname, in the form understood by the FastSpec EOS: the names found
    in the original line lists follow different conventions and there is no
    systematic rule connecting the two.  SH is renamed to HS for 
    compatibility with the EOS. Wavelengths in files whose name contains 
    '.air' are converted to vacuum (SME air2vac formula); the WL_VAC 
    attribute is always 1.
"""

import argparse
import datetime
import os
import sys

import numpy as np
import h5py

FORMAT_VERSION = 'm-1.0'     # FastSpec HDF5 format version

# Columns of the line list, read in file order.  The output is a columnar
# (struct-of-arrays) HDF5 layout: one compressed dataset per column under a
# /Lines group.  The in-memory dtype below is only used to parse the input.
LINE_DTYPE = np.dtype([('W', '<f8'), ('E', '<f4'), ('LGGF', '<f4'),
                       ('GU', '<f4'), ('GRAD', '<f4')])

# Chunking + compression for the per-column datasets.  W is monotonic and the
# f4 columns are smooth, so shuffle+gzip compresses them well; chunking also
# lets the reader touch only the wavelength-window chunks.
CHUNK_SIZE = 1 << 18        # 262144 rows/chunk
GZIP_LEVEL = 4


class TSFormatError(ValueError):
    """The input file is not a readable Turbospectrum line list."""


# =============================================================================
#  Input handling (same transparent decompression as vald2z.py)
# =============================================================================
def open_ts(filename):
    """Open a TS text file, transparently decompressing gzip/bzip2/xz.

    Compression is detected from the file's magic bytes rather than its name,
    so a mis-named compressed file still opens.  Text is decoded as Latin-1
    (a superset of the ASCII these lists contain) so that stray bytes cannot
    raise a decode error; genuinely binary input is caught by the format
    validation instead.
    """
    with open(filename, 'rb') as f:
        magic = f.read(6)
    if magic.startswith(b'\x1f\x8b'):
        import gzip
        return gzip.open(filename, 'rt', encoding='latin-1')
    if magic.startswith(b'BZh') and magic[3:4].isdigit():
        import bz2
        return bz2.open(filename, 'rt', encoding='latin-1')
    if magic.startswith(b'\xfd7zXZ\x00'):
        import lzma
        return lzma.open(filename, 'rt', encoding='latin-1')
    return open(filename, 'r', encoding='latin-1')


def strip_extensions(filename):
    """Base name without a compression suffix and the format extension:
    ``path/list.vac.gz`` -> ``list``  (as in ts2hdf.pro).
    """
    base = os.path.basename(filename)
    root, ext = os.path.splitext(base)
    if ext.lower() in ('.gz', '.bz2', '.xz'):
        base = root
    root, ext = os.path.splitext(base)
    return root if root else base


# =============================================================================
#  Air -> vacuum wavelength conversion (port of SME air2vac.pro)
# =============================================================================
def air2vac(wl_air):
    """Convert an air wavelength (A) to vacuum; identity below 1999.352 A."""
    if wl_air <= 1999.352:
        return wl_air
    sigma2 = (1e4 / wl_air) ** 2          # wavenumber squared (um^-2)
    fact = (1.0 + 8.336624212083e-5
            + 2.408926869968e-2 / (1.301065924522e2 - sigma2)
            + 1.599740894897e-4 / (3.892568793293e1 - sigma2))
    return wl_air * fact


# =============================================================================
#  Turbospectrum reader (port of ReadTS)
# =============================================================================
def read_ts(infile):
    """Read a Turbospectrum line list; return (species_id, ion, lines)
    with ``lines`` a LINE_DTYPE record array in file order.
    """
    air_wave = '.air' in os.path.basename(infile)

    try:
        with open_ts(infile) as f:
            header = f.readline()
            if not header:
                raise TSFormatError('{}: input file is empty.'.format(infile))

            # "'0126.001056   '   1   111404" -> id, ion, number of lines
            tokens = header.replace("'", '').split()
            try:
                species_id = tokens[0]
                ion = int(tokens[1])
                nlin = int(tokens[2])
            except (IndexError, ValueError):
                raise TSFormatError(
                    "{}: not a Turbospectrum line list: the first line is "
                    "not a \"'species_id' ion nlines\" header: {!r}"
                    .format(infile, header.strip()[:70]))
            if nlin <= 0:
                raise TSFormatError(
                    '{}: the header reports no lines.'.format(infile))

            f.readline()               # free-form comment line

            lines = np.empty(nlin, dtype=LINE_DTYPE)
            for i in range(nlin):
                ln = f.readline()
                if not ln:
                    raise TSFormatError(
                        '{}: file truncated? The header declares {} lines '
                        'but the data end after {}.'.format(infile, nlin, i))
                t = ln.split()
                try:
                    w = float(t[0])
                    e = float(t[1])
                    lggf = float(t[2])
                    gu = float(t[4])
                    grad = float(t[5])
                except (IndexError, ValueError):
                    raise TSFormatError(
                        '{}: line {}: expected at least 6 numeric fields '
                        '(wl, E, loggf, fdamp, gu, grad): {!r}'
                        .format(infile, i + 3, ln.strip()[:70]))
                if air_wave:
                    w = air2vac(w)
                if grad == 0.0:
                    # Classical damping from the single-precision loggf and
                    # gu (as in IDL, where the struct fields are floats).
                    grad = 6.669e15 * 10.0 ** float(np.float32(lggf)) \
                        / (float(np.float32(gu)) * w * w)
                lines[i] = (w, e, lggf, gu, grad)
    except EOFError:
        raise TSFormatError(
            '{}: compressed file is truncated or corrupt.'.format(infile))

    return species_id, ion, lines


# =============================================================================
#  HDF5 output (matching IDL H5_CREATE datatypes exactly)
# =============================================================================
def _write_scalar_attr(dset, attname, value):
    """Write a scalar attribute the way IDL's H5_CREATE does.

    Strings become fixed-length null-terminated ASCII of size len+1 (h5py's
    high-level interface would write null-padded instead); integers keep the
    width of the numpy type passed in.  The dataspace is scalar.
    """
    if isinstance(value, str):
        raw = value.encode('latin-1')
        tid = h5py.h5t.C_S1.copy()
        tid.set_size(len(raw) + 1)
        tid.set_strpad(h5py.h5t.STR_NULLTERM)
        tid.set_cset(h5py.h5t.CSET_ASCII)
        space = h5py.h5s.create(h5py.h5s.SCALAR)
        aid = h5py.h5a.create(dset.id, attname.encode('ascii'), tid, space)
        try:
            aid.write(np.array(raw, dtype='S{}'.format(len(raw) + 1)),
                      mtype=tid)
        finally:
            aid.close()
    else:
        dset.attrs.create(attname, value)


# Physical units / conventions for each /Lines column, written as per-column
# <FIELD>_UNIT attributes so the file is self-describing.  Note GRAD is stored
# linearly here (unlike the atomic lists, where the damping columns are log10).
COLUMN_UNITS = {
    'W':    'Angstrom (air or vacuum per the WL_VAC attribute)',
    'E':    'eV (lower-level excitation energy)',
    'LGGF': 'log10(gf) (dimensionless)',
    'GU':   'dimensionless (upper-level statistical weight)',
    'GRAD': '1/s (radiative damping constant, linear)',
}


def _write_column(group, name, data, unit, compress=True):
    """Write one column dataset with a UNIT attribute.  Payload columns are
    chunked + shuffled + gzip-compressed; the sorted wavelength column W is
    written contiguous/uncompressed so the reader's full-column bound scan is a
    fast straight read (compress=False)."""
    n = len(data)
    if n > 0 and compress:
        d = group.create_dataset(name, data=data, chunks=(min(n, CHUNK_SIZE),),
                                 shuffle=True, compression='gzip',
                                 compression_opts=GZIP_LEVEL)
    else:
        d = group.create_dataset(name, data=data)
    _write_scalar_attr(d, 'UNIT', unit)


def write_hdf(outfile, lines, species_id, name, ion, base, provenance):
    """Write the FastSpec molecular-list HDF5 file (columnar layout: one
    compressed dataset per column under the /Lines group)."""
    n = len(lines)
    with h5py.File(outfile, 'w') as f:
        g = f.create_group('Lines')
        # Per-column datasets, each with its UNIT attribute (W uncompressed
        # for a fast bound scan; payload columns compressed).
        for field, unit in COLUMN_UNITS.items():
            _write_column(g, field, lines[field], unit, compress=(field != 'W'))
        # File-level metadata on the /Lines group.
        _write_scalar_attr(g, 'NL', np.int32(n))
        _write_scalar_attr(g, 'ID', species_id)
        _write_scalar_attr(g, 'NAME', name)
        _write_scalar_attr(g, 'ION', np.int16(ion))
        _write_scalar_attr(g, 'WL_VAC', np.int16(1))
        _write_scalar_attr(g, 'FILE', base)
        _write_scalar_attr(g, 'DATE', datetime.date.today().isoformat())
        _write_scalar_attr(g, 'PROVENANCE', provenance)
        _write_scalar_attr(g, 'FORMAT_VERSION', FORMAT_VERSION)
        # Wavelength range (sorted by W) for cheap file-level culling.
        _write_scalar_attr(g, 'WMIN', np.float64(lines['W'][0] if n > 0 else 0.0))
        _write_scalar_attr(g, 'WMAX', np.float64(lines['W'][-1] if n > 0 else 0.0))


# =============================================================================
#  Conversion driver
# =============================================================================
def convert(infile, molname, outfile=None, wmin=None, wmax=None,
            provenance=None):
    """Convert a (possibly compressed) TS list to the FastSpec HDF5 format.

    ``molname`` is the molecule name as understood by the FastSpec EOS
    (see molname_translate.txt); it cannot be inferred from the line list.
    Returns (outfile, n_written).  ``wmin``/``wmax`` restrict the output to
    lines with wavelength >= wmin and/or <= wmax (Angstroms, inclusive).
    ``provenance`` is stored in the PROVENANCE attribute; it defaults to the
    input file name.
    """
    base = strip_extensions(infile)
    if provenance is None:
        provenance = os.path.basename(infile)
    if outfile is None:
        outfile = os.path.join(os.path.dirname(infile), base + '.h5')
    print('{} -> {}'.format(infile, outfile), file=sys.stderr)

    species_id, ion, lines = read_ts(infile)

    # Sort in wavelength (stable: lines of equal wavelength keep file order).
    lines = lines[np.argsort(lines['W'], kind='stable')]

    # Limit lines to a specific wavelength interval.
    if wmin is not None or wmax is not None:
        keep = np.ones(len(lines), dtype=bool)
        if wmin is not None:
            keep &= lines['W'] >= wmin
        if wmax is not None:
            keep &= lines['W'] <= wmax
        nkeep = int(keep.sum())
        if nkeep > 0:
            print('Initial number of lines {}'.format(len(lines)),
                  file=sys.stderr)
            print('Final number of lines {}'.format(nkeep), file=sys.stderr)
            lines = lines[keep]

    # SH must be renamed HS for compatibility with the EOS.
    name = molname
    if name.strip().upper() == 'SH' and ion == 1:
        name = 'HS'

    write_hdf(outfile, lines, species_id, name, ion, base, provenance)
    return outfile, len(lines)


# =============================================================================
#  Command-line interface
# =============================================================================
def main(argv=None):
    p = argparse.ArgumentParser(
        description='Convert a molecular line list in Turbospectrum ASCII '
                    'format into the HDF5 file read by FastSpec')
    p.add_argument('infile', help='input TS line list, single species '
                                  '(plain or gzip/bzip2/xz compressed)')
    p.add_argument('molname',
                   help='molecule name stored in the NAME attribute, in the '
                        'form understood by the FastSpec EOS')
    p.add_argument('outfile', nargs='?',
                   help='output file (default: <infile stem>.h5)')
    p.add_argument('-wmin', dest='wmin', type=float, default=None,
                   metavar='WL',
                   help='output only lines with wavelength >= WL (Angstroms)')
    p.add_argument('-wmax', dest='wmax', type=float, default=None,
                   metavar='WL',
                   help='output only lines with wavelength <= WL (Angstroms)')
    p.add_argument('--provenance', dest='provenance', default=None,
                   metavar='TEXT',
                   help='text stored in the PROVENANCE attribute '
                        '(default: the input file name)')
    args = p.parse_args(argv)

    if args.outfile is not None and args.outfile == args.infile:
        p.error('Same input and output file names. Aborting.')

    try:
        outfile, nlin = convert(args.infile, args.molname, args.outfile,
                                wmin=args.wmin, wmax=args.wmax,
                                provenance=args.provenance)
    except (TSFormatError, OSError) as e:
        print('Error: {}'.format(e), file=sys.stderr)
        return 1

    print('{} lines written'.format(nlin), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
