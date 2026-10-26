#!/usr/bin/env python3.12
"""
vald2hdf.py -- convert a VALD3 extraction into the FastSpec atomic-line HDF5 file

  * All four VALD output formats are understood and auto-detected:
    extract stellar / extract all, each in short or long form.  The input
    may be gzip-, bzip2- or xz-compressed (detected from the file content).
  * Molecular lines are rejected and counted .  This converter writes a 
    multi-species atomic list (species carried per line in the EL column); 
    single-species molecular lists are produced by mol2hdf.py.  FastSpec 
    reads both with one routine (ReadHDF), which tells them apart by the 
    presence of the EL column.
  * Hydrogen lines are skipped and counted.  With --hlist (a dedicated 
    VALD hydrogen extraction) H lines are kept and the HFLAG attribute is 
    set to 1.
  * The air/vacuum wavelength scale (WL_VAC attribute) and the excitation
    energy unit (eV vs cm^-1, converted to eV) are detected from the VALD
    column headers; IDL assumed air and eV.
  * Ions above ION = 5 are excluded (EOS limit) and counted.
  * With --zeeman the output (written with a .zh5 extension) carries, in
    addition, everything needed to build Zeeman splitting patterns
    downstream: the lower/upper level J values (JLO, JUP) and Lande factors
    (GLO, GUP).  Missing Lande factors are estimated with LS/JJ/JK method; 
    a factor that cannot be computed keeps 99.0, and electric-dipole-forbidden 
    lines get 99.0 for both.  This mode requires a VALD long-format extraction 
    and sets FORMAT_VERSION to 'z-1.0'.

Output layout (columnar / struct-of-arrays):
    group /Lines : one 1-D dataset per column {EL:str, ION:i16, W:f8, E:f4,
                   LGGF:f4, GRAD:f4, GQST:f4, GVDW:f4}, all sorted in
                   wavelength.  With --zeeman four more columns follow
                   {JLO:f4, JUP:f4, GLO:f4, GUP:f4}.  Payload columns are
                   chunked + shuffled + gzip-compressed; W is left contiguous
                   (uncompressed) so the reader's full-column bound scan is a
                   fast straight read.  Each column dataset carries a UNIT
                   attribute.  Group attributes: NL (int32), NAME ('atom'),
                   WL_VAC, HFLAG (int16), FILE (null-terminated ASCII), DATE
                   (YYYY-MM-DD), PROVENANCE (input file name, or --provenance
                   text), FORMAT_VERSION ('a-1.0', or 'z-1.0' with --zeeman),
                   and WMIN/WMAX (f8, the file's wavelength range, for fast
                   file-level culling).
"""

import argparse
import datetime
import os
import re
import sys

import numpy as np
import h5py


CM_PER_EV = 8065.543937          # eV -> cm^-1 factor
MAXION = 5                       # highest ionization stage supported by the EOS
ATOMIC_FORMAT_VERSION = 'a-1.0'  # FastSpec HDF5 format version (atomic list)
ZEEMAN_FORMAT_VERSION = 'z-1.0'  # ... plus per-level J and Lande factors

# Chunking + compression for the per-column datasets (columnar output layout).
CHUNK_SIZE = 1 << 18         # 262144 rows/chunk
GZIP_LEVEL = 4

# All chemical element symbols (case-sensitive, as VALD writes them).
ELEMENTS = frozenset([
    'H',  'He', 'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',  'K', 'Ca',
    'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',  'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te',  'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
])


class VALDFormatError(ValueError):
    """The input file is not a readable VALD extraction."""


# =============================================================================
#  Lande-factor estimation (LS / JJ / JK coupling)
#
#  A Lande factor that VALD flags as missing (99.000) is estimated from the 
#  level's coupling scheme and its term/configuration string; parsing failures 
#  return the 99.0 sentinel.
# =============================================================================
# Letters used by VALD for orbital angular-momentum quantum numbers.
LVAL = ('S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
        'N', 'O', 'Q', 'R', 'T', 'U', 'V')
LVAL_LOWER = tuple(s.lower() for s in LVAL)


def clean_config(raw):
    """Strip surrounding quotes/whitespace from a VALD config/term line."""
    return raw.strip().strip("'").strip()


def _parse_fraction(text):
    """Convert a text fraction such as '3/2' (or an integer '5') to float.

    Returns None if the string cannot be parsed.
    """
    text = text.strip()
    try:
        if '/' in text:
            num, den = text.split('/')
            return float(num) / float(den)
        return float(text)
    except (ValueError, ZeroDivisionError):
        return None


def get_ls_numbers(config):
    """Return (L, S) for a level in LS coupling, or (None, None) on failure.

    ``config`` is the full VALD configuration+term string, e.g.
    ``"LS  3d6.(5D).4s.4p.(3P*) z7F*"``.  The term symbol is its last token,
    of the form (2S+1)L possibly decorated with a leading label and a trailing
    parity '*'.
    """
    tokens = config.split()
    if not tokens:
        return None, None
    term = tokens[-1]

    L = None
    S = None
    for i, ch in enumerate(term):
        if ch in LVAL:
            L = LVAL.index(ch)
            # The multiplicity 2S+1 precedes the letter; allow two digits.
            ii = i - 1
            if term[0:i].isnumeric():
                ii = 0
            frac = _parse_fraction(term[ii:i])
            if frac is None:
                return None, None
            S = (frac - 1.0) / 2.0
    if L is None or S is None or S < 0.0:
        return None, None
    return L, S


def _find_subgroup_ls(econfig, start):
    """Scan ``econfig[start::-1]`` for the innermost L (capital) or l (lower).

    Returns (index, L, S) where S = 1/2 for a single (lower-case) electron,
    or is derived from the preceding multiplicity for a subgroup (capital
    letter).  Returns (None, None, None) if nothing is found, or
    (index, None, None) if a multiplicity could not be parsed.
    """
    for i in range(start, 0, -1):
        ch = econfig[i]
        if ch in LVAL:
            L = LVAL.index(ch)
            ii = i - 1
            if econfig[i - 2:i].isnumeric():
                ii = i - 2
            frac = _parse_fraction(econfig[ii:i])
            if frac is None:
                return i, None, None
            return i, L, (frac - 1.0) / 2.0
        if ch in LVAL_LOWER:
            return i, LVAL_LOWER.index(ch), 0.5
    return None, None, None


def get_jj_numbers(config):
    """Return (L1, S1, J1, L2, S2, J2) for a level in JJ coupling.

    Any component being None signals a parse failure. Covers the J1-J2, J1-j and j-j sub-cases.
    """
    tokens = config.split()
    if len(tokens) < 2:
        return (None,) * 6
    term = tokens[-1]
    econfig = tokens[-2]

    # J1, J2 appear in parentheses in the term symbol, e.g. "(5,1/2)".
    if ',' not in term:
        return (None,) * 6
    try:
        sj1, sj2 = term.strip('*').strip('()').split(',')
    except ValueError:
        return (None,) * 6
    J1 = _parse_fraction(sj1)
    J2 = _parse_fraction(sj2)
    if J1 is None or J2 is None:
        return (None,) * 6

    # L2 (and possibly S2) from the innermost / last (sub)group.
    i2, L2, S2 = _find_subgroup_ls(econfig, len(econfig) - 1)
    if i2 is None or L2 is None or S2 is None:
        return (None,) * 6
    # L1 (and possibly S1) from the (sub)group preceding L2.
    _, L1, S1 = _find_subgroup_ls(econfig, i2 - 1)
    if L1 is None or S1 is None:
        return (None,) * 6
    return L1, S1, J1, L2, S2, J2


def get_jk_numbers(config):
    """Return (L1, S1, J1, L2, S2, K) for a level in JK coupling.

    Any component being None signals a parse failure.
    """
    tokens = config.split()
    if len(tokens) < 2:
        return (None,) * 6
    term = tokens[-1]
    econfig = tokens[-2]

    # K is given in square brackets in the term symbol, e.g. "(5,1/2)"->"[..]".
    if '[' not in term or ']' not in term:
        return (None,) * 6
    sK = term.split('[')[1].split(']')[0]
    K = _parse_fraction(sK)
    if K is None:
        return (None,) * 6

    # L2 (and possibly S2) from the last (sub)group / single electron.
    i2, L2, S2 = _find_subgroup_ls(econfig, len(econfig) - 1)
    if i2 is None or L2 is None or S2 is None:
        return (None,) * 6

    # Try to find L1 as a capital-letter subgroup preceding L2, along with J1
    # given in <> right after it.  Otherwise fall back to a lower-case l1.
    L1 = None
    S1 = None
    J1 = None
    for i1 in range(i2 - 1, 0, -1):
        if econfig[i1] in LVAL:
            L1 = LVAL.index(econfig[i1])
            ii = i1 - 1
            if econfig[i1 - 2:i1].isnumeric():
                ii = i1 - 2
            frac = _parse_fraction(econfig[ii:i1])
            if frac is None:
                return (None,) * 6
            S1 = (frac - 1.0) / 2.0
            if '<' in econfig[i1:i2] and '>' in econfig[i1:i2]:
                sJ1 = econfig[econfig.find('<') + 1:econfig.find('>')]
                J1 = _parse_fraction(sJ1)
                if J1 is None:
                    return (None,) * 6
            break
    if L1 is None:
        # j-l case: a single lower-case electron carries J1, given before [K].
        for i1 in range(i2 - 1, 0, -1):
            if econfig[i1] in LVAL_LOWER:
                L1 = LVAL_LOWER.index(econfig[i1])
                S1 = 0.5
                J1 = _parse_fraction(term.split('[')[0])
                if J1 is None:
                    return (None,) * 6
                break

    if L1 is None or L2 is None or S1 is None or S2 is None or J1 is None:
        return (None,) * 6
    return L1, S1, J1, L2, S2, K


def _triangular(a, b, c):
    """True if the angular momenta (a, b, c) satisfy the vector-coupling
    triangle condition |a - b| <= c <= a + b with an integer perimeter."""
    if c < abs(a - b) - 1e-6 or c > a + b + 1e-6:
        return False
    return abs((a + b + c) - round(a + b + c)) < 1e-6


def estimate_lande(J, config, schemes=('LS', 'JJ', 'JK')):
    """Estimate a level's Lande factor from its J and VALD term/config string.

    Supports LS, JJ and JK coupling (restricted to ``schemes``).  Returns
    99.0 when the coupling scheme is unsupported/disabled or the string
    cannot be parsed -- exactly the sentinel VALD uses for "no Lande factor".
    """
    if len(config.split()) < 2:
        return 99.0
    coupling = config.split()[0]

    if coupling not in schemes:
        return 99.0

    if coupling == 'LS':
        L, S = get_ls_numbers(config)
        if L is None or S is None:
            return 99.0
        if J == 0.0:
            return 0.0
        return 1.0 + (J * (J + 1.0) - L * (L + 1.0) + S * (S + 1.0)) \
            / (2.0 * J * (J + 1.0))

    if coupling == 'JJ':
        L1, S1, J1, L2, S2, J2 = get_jj_numbers(config)
        if None in (L1, S1, J1, L2, S2, J2):
            return 99.0
        if J == 0.0:
            return 0.0
        gam1 = (J * (J + 1.) + J1 * (J1 + 1.) - J2 * (J2 + 1.)) / (2. * J * (J + 1.))
        gam3 = (J * (J + 1.) + J2 * (J2 + 1.) - J1 * (J1 + 1.)) / (2. * J * (J + 1.))
        gam2 = 0.0 if J1 == 0 else \
            (J1 * (J1 + 1.) + S1 * (S1 + 1.) - L1 * (L1 + 1.)) / (2. * J1 * (J1 + 1.))
        gam4 = 0.0 if J2 == 0 else \
            (J2 * (J2 + 1.) + S2 * (S2 + 1.) - L2 * (L2 + 1.)) / (2. * J2 * (J2 + 1.))
        return 1.0 + gam1 * gam2 + gam3 * gam4

    if coupling == 'JK':
        L1, S1, J1, L2, S2, K = get_jk_numbers(config)
        if None in (L1, S1, J1, L2, S2, K):
            return 99.0
        # The final coupling is J = K + S2, so J must satisfy the triangle
        # condition with K and S2.  A VALD term whose (K, J) are inconsistent
        # (a strongly-mixed level given a nominal label) is not computable.
        if not _triangular(K, S2, J):
            return 99.0
        if J == 0.0:
            return 0.0
        # Exact two-step projection (g_s = 2): K = J1 + L2 with the orbital L2
        # carrying g = 1, then J = K + S2 with the spin S2 carrying g = 2.  The
        # compact (2J+1)(2K+1) form is exact only for S2 = 1/2 and is wrong for
        # multi-electron outer groups (S2 != 1/2), so use the projection.
        gam1 = 0.0 if J1 == 0.0 else \
            (J1 * (J1 + 1.) - L1 * (L1 + 1.) + S1 * (S1 + 1.)) / (2. * J1 * (J1 + 1.))
        cK = 0.0 if K == 0.0 else \
            (K * (K + 1.) + J1 * (J1 + 1.) - L2 * (L2 + 1.)) / (2. * K * (K + 1.))
        aJ = (J * (J + 1.) + K * (K + 1.) - S2 * (S2 + 1.)) / (2. * J * (J + 1.))
        gK = 1.0 + gam1 * cK
        return gK * aJ + 2.0 * (1.0 - aJ)

    return 99.0


# =============================================================================
#  Input handling (same transparent decompression as vald2z.py / mol2hdf.py)
# =============================================================================
def open_vald(filename):
    """Open a VALD text file, transparently decompressing gzip/bzip2/xz.

    Compression is detected from the file's magic bytes rather than its name,
    so a mis-named compressed file still opens.  Text is decoded as Latin-1
    (a superset of the ASCII VALD emits) so that stray bytes cannot raise a
    decode error; genuinely binary input is caught by the format validation
    instead.
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


def read_vald_lines(infile):
    """Read all lines of a (possibly compressed) VALD extraction."""
    try:
        with open_vald(infile) as f:
            return [ln.rstrip('\r\n') for ln in f]
    except EOFError:
        raise VALDFormatError(
            '{}: compressed file is truncated or corrupt.'.format(infile))


def strip_extensions(filename):
    """Base name without a compression suffix and the format extension:
    ``path/list.dat.gz`` -> ``list``.
    """
    base = os.path.basename(filename)
    root, ext = os.path.splitext(base)
    if ext.lower() in ('.gz', '.bz2', '.xz'):
        base = root
    root, ext = os.path.splitext(base)
    return root if root else base


# =============================================================================
#  Format detection
# =============================================================================
# Column positions of the used quantities in a comma-split data line.
# (Damping constants sit before the quoted reference string, so commas
# inside references cannot shift them.)
_COLUMNS = {                 # w, e, lggf, grad, gqst, gvdw
    ('stellar', 'long'):  (1, 3, 2, 10, 11, 12),
    ('stellar', 'short'): (1, 2, 4, 5, 6, 7),      # field 3 is Vmic
    ('all', 'long'):      (1, 3, 2, 10, 11, 12),
    ('all', 'short'):     (1, 2, 3, 4, 5, 6),
}
_MIN_FIELDS = {
    ('stellar', 'long'):  14,
    ('stellar', 'short'): 10,
    ('all', 'long'):      13,
    ('all', 'short'):     7,
}

# Extra long-format columns read for the --zeeman output: J and VALD Lande
# factor of the lower/upper level.  (Same field positions in both the stellar
# and all long layouts.)
_ZEEMAN_COLUMNS = (4, 6, 7, 8)   # jlo, jup, lande_lo, lande_up

_ENERGY_UNIT_RE = re.compile(r'(?:E[_ ]?low|Excit)\s*\(([^)]*)\)',
                             re.IGNORECASE)


def detect_format(infile, lines):
    """Identify the VALD output flavour.

    Returns (kind, layout, header_offset, nlin) where kind is 'stellar' or
    'all', layout is 'short' or 'long', and nlin is the number of
    transitions (declared in the extraction summary for extract stellar,
    obtained by scanning for extract all).
    """
    if not lines:
        raise VALDFormatError('{}: input file is empty.'.format(infile))

    # Extract stellar starts with the extraction summary; extract all with
    # the column headers (which mention the Lande/damping columns).
    fields = lines[0].split(',')
    try:
        float(fields[0])
        float(fields[1])
        nlin = int(fields[2])
        kind = 'stellar'
    except (ValueError, IndexError):
        if 'lande' in lines[0].lower() or 'damping' in lines[0].lower():
            kind = 'all'
            nlin = None
        else:
            raise VALDFormatError(
                '{}: not a VALD extraction: the first line is neither an '
                'extraction summary (extract stellar) nor a column header '
                '(extract all): {!r}'
                .format(infile, lines[0].strip()[:70]))

    header_offset = 3 if kind == 'stellar' else 2
    colhdr = lines[header_offset - 1] if len(lines) > header_offset - 1 \
        else ''
    layout = 'long' if 'j lo' in colhdr.lower() else 'short'

    if kind == 'stellar':
        if nlin <= 0:
            raise VALDFormatError(
                '{}: the extraction summary reports no selected lines.'
                .format(infile))
        group = 4 if layout == 'long' else 1
        need = header_offset + group * nlin
        if len(lines) < need:
            raise VALDFormatError(
                '{}: file truncated? The extraction summary declares {} '
                'transitions ({} lines including the header), but the file '
                'has only {} lines.'.format(infile, nlin, need, len(lines)))
    else:
        # Scan for the number of transitions: data (group-boundary) lines
        # start with a quoted species name; the trailing reference section
        # does not (and "'_" opens an in-group reference line).
        group = 4 if layout == 'long' else 1
        nlin = 0
        i = header_offset
        while i + group - 1 < len(lines):
            ln = lines[i].lstrip()
            if not ln.startswith("'") or ln.startswith("'_"):
                break
            nlin += 1
            i += group
        if nlin == 0:
            raise VALDFormatError(
                '{}: no transition data found after the extract-all column '
                'headers.'.format(infile))

    # Consistency check on the first data line.
    ndat = len(lines[header_offset].split(','))
    if ndat < _MIN_FIELDS[(kind, layout)]:
        raise VALDFormatError(
            '{}: line {}: an extract-{} {}-format transition line has at '
            'least {} comma-separated fields, found {}: {!r}'
            .format(infile, header_offset + 1, kind, layout,
                    _MIN_FIELDS[(kind, layout)], ndat,
                    lines[header_offset].strip()[:70]))

    return kind, layout, header_offset, nlin


def detect_energy_in_cm(header_lines):
    """True if the excitation energies are in cm^-1 (else eV, VALD default)."""
    for line in header_lines:
        m = _ENERGY_UNIT_RE.search(line)
        if m:
            unit = m.group(1).lower()
            if 'cm' in unit:
                return True
            if 'ev' in unit:
                return False
    return False


def detect_vacuum(header_lines):
    """True if the wavelength column header says vacuum (WL_vac)."""
    for line in header_lines:
        low = line.lower()
        if 'wl_vac' in low:
            return True
        if 'wl_air' in low:
            return False
    print('Warning: no WL_air/WL_vac column header found; assuming air',
          file=sys.stderr)
    return False


# =============================================================================
#  VALD parsing
# =============================================================================
class Stats:
    def __init__(self):
        self.processed = 0      # transitions in the input
        self.n_mol = 0          # molecular lines rejected
        self.n_h = 0            # hydrogen lines skipped
        self.n_ion = 0          # ion > MAXION lines excluded
        self.written = 0        # lines written to the output
        self.vald_lande = 0     # magnetic lines using VALD Lande factors
        self.computed_lande = 0  # magnetic lines with >=1 computed factor
        self.computed_ls = 0    # of those, lines that used LS coupling
        self.computed_jj = 0    # ... JJ coupling
        self.computed_jk = 0    # ... JK coupling
        self.n_forbidden = 0    # E1-forbidden lines (glo = gup = 99)
        self.n_failed = 0       # atomic lines with no computable Lande factor


def classify_lande(jlo, jup, lande_lo, lande_up, cfg_lo, cfg_up):
    """Lower/upper-level Lande factors and the transition's Lande category.

    Returns (glo, gup, category, schemes_used), mirroring vald2z.py with its
    default settings (LS/JJ/JK coupling, LS applied to every element).  A
    VALD-provided factor is kept; a factor VALD flags as missing (99.0) is
    estimated from the level's coupling scheme and left at 99.0 when it cannot
    be computed.  ``category`` is one of:

      'forbidden' : E1-forbidden line (Jlo = Jup = 0 or |Jlo - Jup| > 1);
                    glo = gup = 99.0.
      'vald'      : magnetic line, both factors taken from VALD.
      'computed'  : magnetic line, >=1 factor estimated (``schemes_used`` names
                    the coupling scheme(s) used).
      'failed'    : a missing factor could not be estimated (kept 99.0), or the
                    estimate gives a null (both-zero) split -- as vald2z.py, a
                    non-magnetic outcome.
    """
    if (jlo == 0.0 and jup == 0.0) or abs(jlo - jup) > 1:
        return 99.0, 99.0, 'forbidden', frozenset()

    if lande_up == 99.0 or lande_lo == 99.0:
        lo_missing = lande_lo == 99.0
        up_missing = lande_up == 99.0
        glo = estimate_lande(jlo, cfg_lo) if lo_missing else lande_lo
        gup = estimate_lande(jup, cfg_up) if up_missing else lande_up
        if glo == 99.0 or gup == 99.0:
            return glo, gup, 'failed', frozenset()
        schemes = set()
        if lo_missing:
            schemes.add(cfg_lo.split()[0])
        if up_missing:
            schemes.add(cfg_up.split()[0])
        # A pattern with both g = 0 does not split -> non-magnetic.
        if glo == 0.0 and gup == 0.0:
            return glo, gup, 'failed', frozenset()
        return glo, gup, 'computed', schemes

    return lande_lo, lande_up, 'vald', frozenset()


def read_vald(infile, hlist=False, zeeman=False):
    """Read any VALD extraction; return (records, stats, wl_vac).

    ``records`` is a list of (el, ion, w, e, lggf, grad, gqst, gvdw) tuples
    in file order, restricted to atomic species with ion <= MAXION;
    hydrogen is kept only when ``hlist`` is set.  With ``zeeman`` each tuple
    is extended by (jlo, jup, glo, gup) -- the lower/upper level J values and
    Lande factors (see classify_lande) -- and the Lande statistics in
    ``stats`` are filled in; this requires the long format.
    """
    lines = read_vald_lines(infile)
    kind, layout, header_offset, nlin = detect_format(infile, lines)
    print('Format: extract {}, {} format'.format(kind, layout),
          file=sys.stderr)

    if zeeman and layout != 'long':
        raise VALDFormatError(
            '{}: --zeeman requires a VALD long-format extraction; this is a '
            '{} format, which carries no J or Lande-factor columns.'
            .format(infile, layout))

    header = lines[:header_offset]
    energy_in_cm = detect_energy_in_cm(header)
    if energy_in_cm:
        print('Energies given in cm^-1; converting to eV', file=sys.stderr)
    wl_vac = detect_vacuum(header)

    iw, ie, ig, ir, is_, iv = _COLUMNS[(kind, layout)]
    ijlo, ijup, ilo, iup = _ZEEMAN_COLUMNS
    group = 4 if layout == 'long' else 1

    stats = Stats()
    stats.processed = nlin
    records = []
    molecules = set()

    for k in range(nlin):
        base = header_offset + k * group
        data = lines[base]
        vals = data.split(',')
        try:
            species = vals[0].replace("'", '').strip()
            tokens = species.split()
            if len(tokens) != 2:
                raise ValueError
            el = tokens[0]
            ion = int(tokens[1])
            w = float(vals[iw])
            e = float(vals[ie])
            lggf = float(vals[ig])
            grad = float(vals[ir])
            gqst = float(vals[is_])
            gvdw = float(vals[iv])
            if zeeman:
                jlo = float(vals[ijlo])
                jup = float(vals[ijup])
                lande_lo = float(vals[ilo])
                lande_up = float(vals[iup])
        except (ValueError, IndexError):
            raise VALDFormatError(
                '{}: line {}: cannot parse as an extract-{} {}-format '
                'transition line: {!r}'
                .format(infile, base + 1, kind, layout, data.strip()[:70]))

        if el not in ELEMENTS:
            stats.n_mol += 1
            molecules.add(species)
            continue
        if el == 'H' and not hlist:
            stats.n_h += 1
            continue
        if ion > MAXION:
            stats.n_ion += 1
            continue

        if energy_in_cm:
            e /= CM_PER_EV
        if zeeman:
            glo, gup, category, schemes_used = classify_lande(
                jlo, jup, lande_lo, lande_up,
                clean_config(lines[base + 1]), clean_config(lines[base + 2]))
            if category == 'vald':
                stats.vald_lande += 1
            elif category == 'computed':
                stats.computed_lande += 1
                if 'LS' in schemes_used:
                    stats.computed_ls += 1
                if 'JJ' in schemes_used:
                    stats.computed_jj += 1
                if 'JK' in schemes_used:
                    stats.computed_jk += 1
            elif category == 'forbidden':
                stats.n_forbidden += 1
            else:  # 'failed'
                stats.n_failed += 1
            records.append((el, ion, w, e, lggf, grad, gqst, gvdw,
                            jlo, jup, glo, gup))
        else:
            records.append((el, ion, w, e, lggf, grad, gqst, gvdw))

    if stats.n_mol:
        print('Excluded {} molecular lines ({})'
              .format(stats.n_mol, ', '.join(sorted(molecules))),
              file=sys.stderr)
    else:
        print('Excluded 0 molecular lines', file=sys.stderr)
    print('Excluded {} hydrogen lines'.format(stats.n_h)
          + (' (kept: hydrogen list)' if hlist else ''), file=sys.stderr)
    print('Excluded {} lines with ion > {}'.format(stats.n_ion, MAXION),
          file=sys.stderr)

    if zeeman:
        print('{} lines with VALD Lande factors'.format(stats.vald_lande),
              file=sys.stderr)
        print('{} lines with computed Lande factors ({}/{}/{} for LS/JJ/JK)'
              .format(stats.computed_lande, stats.computed_ls,
                      stats.computed_jj, stats.computed_jk), file=sys.stderr)
        print('{} forbidden lines retained as non-magnetic'
              .format(stats.n_forbidden), file=sys.stderr)
        print('{} lines with no computable Lande factor retained as '
              'non-magnetic'.format(stats.n_failed), file=sys.stderr)

    return records, stats, wl_vac


# =============================================================================
#  HDF5 output
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


def _lines_np_dtype(strsize, zeeman=False):
    """Packed numpy dtype for the /Lines columns.  The atomic list has eight
    columns; the --zeeman list appends JLO, JUP, GLO, GUP (all f4)."""
    fields = [('EL', 'S{}'.format(strsize)), ('ION', '<i2'),
              ('W', '<f8'), ('E', '<f4'), ('LGGF', '<f4'),
              ('GRAD', '<f4'), ('GQST', '<f4'), ('GVDW', '<f4')]
    if zeeman:
        fields += [('JLO', '<f4'), ('JUP', '<f4'),
                   ('GLO', '<f4'), ('GUP', '<f4')]
    return np.dtype(fields)


# Physical units / conventions for each /Lines column, written as per-column
# <FIELD>_UNIT attributes so the file is self-describing.  The damping columns
# follow the VALD convention (stored as log10; GVDW packs ABO when positive).
COLUMN_UNITS = {
    'EL':   'element symbol',
    'ION':  'ionization stage (1 = neutral)',
    'W':    'Angstrom (air or vacuum per the WL_VAC attribute)',
    'E':    'eV (lower-level excitation energy)',
    'LGGF': 'log10(gf) (dimensionless)',
    'GRAD': 'log10(radiative damping /s) (VALD Rad.)',
    'GQST': 'log10(Stark damping /s per electron at 1e4 K) (VALD Stark)',
    'GVDW': 'VALD Waals: packed ABO sigma.alpha if >0, else '
            'log10(vdW damping /s per H atom at 1e4 K)',
    'JLO':  'J of the lower level (dimensionless)',
    'JUP':  'J of the upper level (dimensionless)',
    'GLO':  'Lande factor of the lower level (99.0 = unknown or forbidden)',
    'GUP':  'Lande factor of the upper level (99.0 = unknown or forbidden)',
}


def write_hdf(outfile, arr, strsize, wl_vac, hflag, base, provenance,
              zeeman=False):
    """Write the FastSpec atomic-list HDF5 file (columnar layout: one
    compressed dataset per column under the /Lines group).  With ``zeeman``
    the file carries the extra JLO/JUP/GLO/GUP columns and the z-1.0 format
    version."""
    n = len(arr)
    with h5py.File(outfile, 'w') as f:
        g = f.create_group('Lines')
        # Per-column datasets, each with its UNIT attribute (W uncompressed
        # for a fast bound scan; payload columns compressed).  The set of
        # columns follows the array's dtype (eight, or twelve with --zeeman).
        for field in arr.dtype.names:
            _write_column(g, field, arr[field], COLUMN_UNITS[field],
                          compress=(field != 'W'))
        # File-level metadata on the /Lines group.
        _write_scalar_attr(g, 'NL', np.int32(n))
        _write_scalar_attr(g, 'NAME', 'atom Zeeman' if zeeman else 'atom')
        _write_scalar_attr(g, 'WL_VAC', np.int16(1 if wl_vac else 0))
        _write_scalar_attr(g, 'HFLAG', np.int16(1 if hflag else 0))
        _write_scalar_attr(g, 'FILE', base)
        _write_scalar_attr(g, 'DATE', datetime.date.today().isoformat())
        _write_scalar_attr(g, 'PROVENANCE', provenance)
        _write_scalar_attr(g, 'FORMAT_VERSION',
                           ZEEMAN_FORMAT_VERSION if zeeman else ATOMIC_FORMAT_VERSION)
        # Wavelength range (sorted by W) for cheap file-level culling.
        _write_scalar_attr(g, 'WMIN', np.float64(arr['W'][0] if n > 0 else 0.0))
        _write_scalar_attr(g, 'WMAX', np.float64(arr['W'][-1] if n > 0 else 0.0))


# =============================================================================
#  Conversion driver
# =============================================================================
def convert(infile, outfile=None, wmin=None, wmax=None, hlist=False,
            provenance=None, zeeman=False):
    """Convert a (possibly compressed) VALD extraction to the FastSpec
    atomic-line HDF5 format.

    Returns (outfile, stats).  ``wmin``/``wmax`` restrict the output to
    lines with wavelength >= wmin and/or <= wmax (Angstroms, inclusive).
    ``hlist`` marks a dedicated hydrogen line list: H lines are kept and
    the HFLAG attribute is set to 1.  ``provenance`` is stored in the
    PROVENANCE attribute; it defaults to the input file name.  ``zeeman``
    writes the Zeeman-ready format (z-1.0) with the extra JLO/JUP/GLO/GUP
    columns; it requires a VALD long-format extraction.
    """
    base = strip_extensions(infile)
    if provenance is None:
        provenance = os.path.basename(infile)
    if outfile is None:
        suffix = '.zh5' if zeeman else '.h5'
        outfile = os.path.join(os.path.dirname(infile), base + suffix)
    print('{} -> {} ({} format)'
          .format(infile, outfile, 'Zeeman' if zeeman else 'atomic'),
          file=sys.stderr)

    records, stats, wl_vac = read_vald(infile, hlist=hlist, zeeman=zeeman)
    if not records:
        raise VALDFormatError(
            '{}: no atomic lines retained; nothing to write.'.format(infile))

    strsize = max(len(r[0]) for r in records) + 1
    arr = np.array([(r[0].encode('ascii'),) + r[1:] for r in records],
                   dtype=_lines_np_dtype(strsize, zeeman=zeeman))

    # Sort in wavelength (stable: lines of equal wavelength keep file order).
    arr = arr[np.argsort(arr['W'], kind='stable')]

    # Limit lines to a specific wavelength interval.
    if wmin is not None or wmax is not None:
        keep = np.ones(len(arr), dtype=bool)
        if wmin is not None:
            keep &= arr['W'] >= wmin
        if wmax is not None:
            keep &= arr['W'] <= wmax
        nkeep = int(keep.sum())
        if nkeep > 0:
            print('Initial number of lines {}'.format(len(arr)),
                  file=sys.stderr)
            print('Final number of lines {}'.format(nkeep), file=sys.stderr)
            arr = arr[keep]

    write_hdf(outfile, arr, strsize, wl_vac, hlist, base, provenance,
              zeeman=zeeman)
    stats.written = len(arr)
    return outfile, stats


# =============================================================================
#  Command-line interface
# =============================================================================
def main(argv=None):
    p = argparse.ArgumentParser(
        description='Convert a VALD3 extraction (extract stellar or extract '
                    'all, short or long format) into the atomic-line HDF5 '
                    'or Zeeman atomic-line HDF5 file read by FastSpec. '
                    'Molecular, hydrogen, and Ion > 5 lines are skipped.')
    p.add_argument('infile', help='input VALD extraction '
                                  '(plain or gzip/bzip2/xz compressed)')
    p.add_argument('outfile', nargs='?',
                   help='output file (default: <infile stem>.h5, or .zh5 '
                        'with --zeeman)')
    p.add_argument('-wmin', dest='wmin', type=float, default=None,
                   metavar='WL',
                   help='output only lines with wavelength >= WL (Angstroms)')
    p.add_argument('-wmax', dest='wmax', type=float, default=None,
                   metavar='WL',
                   help='output only lines with wavelength <= WL (Angstroms)')
    p.add_argument('--hlist', action='store_true',
                   help='the input is a dedicated hydrogen line list: keep '
                        'H lines and set the HFLAG attribute to 1')
    p.add_argument('--zeeman', action='store_true',
                   help='write the Zeeman-ready format (FORMAT_VERSION '
                        'z-1.0): the atomic list plus per-level JLO/JUP and '
                        'Lande factors GLO/GUP (missing factors estimated as '
                        'in vald2z.py, 99.0 where uncomputable or forbidden). '
                        'Requires a VALD long-format extraction.')
    p.add_argument('--provenance', dest='provenance', default=None,
                   metavar='TEXT',
                   help='text stored in the PROVENANCE attribute '
                        '(default: the input file name)')
    args = p.parse_args(argv)

    if args.outfile is not None and args.outfile == args.infile:
        p.error('Same input and output file names. Aborting.')

    try:
        outfile, stats = convert(args.infile, args.outfile,
                                 wmin=args.wmin, wmax=args.wmax,
                                 hlist=args.hlist, provenance=args.provenance,
                                 zeeman=args.zeeman)
    except (VALDFormatError, OSError) as e:
        print('Error: {}'.format(e), file=sys.stderr)
        return 1

    print('{} lines processed, {} lines written'
          .format(stats.processed, stats.written), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
