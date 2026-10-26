#!/usr/bin/env python3.12
"""
vald2z.py -- convert a VALD3 "extract stellar" or "extract all" long-format
line list into a Synthmag/Synmast/FastSpec input list with explicit Zeeman
splitting patterns. The two formats are told apart automatically: "extract
all" starts with the column headers (no extraction-summary line) and carries
no central-depth column (depth is set to 0.5) and no microturbulence (0 km/s
is assumed). The input may be gzip-, bzip2- or xz-compressed; compression is
detected from the file content, not from the file name.

------------------------------------------------------------------------------
Physics notes
------------------------------------------------------------------------------
Zeeman pattern
    In the weak-field / anomalous-Zeeman regime a spectral line splits into
    pi components (delta m = 0) and sigma components (delta m = +/-1). The
    wavelength shift of a component (in units of the normal Lorentz unit) is

        pi  :  m * (g_up - g_lo)                 (m_up = m_lo = m)
        sig :  m*g_up - (m+1)*g_lo               (m_up = m, m_lo = m+1)

    The relative strengths of the components are the exact squared
    Clebsch-Gordan / Wigner-3j coefficients, which for electric-dipole
    transitions reduce to the closed-form polynomials in J and m used below
    (see Landi Degl'Innocenti & Landolfi 2004).

Effective Lande factor
    g_eff = 0.5*(g_up + g_lo) + 0.25*(g_up - g_lo)*(Jup(Jup+1) - Jlo(Jlo+1))
    (Shenstone & Blair 1929; Landi Degl'Innocenti & Landolfi 2004 eq. 3.44).
    This is written to the output "Lande factor" column, encoded as
    ``-90 - g_eff`` (the Synthmag convention that flags a line carrying an
    explicit split pattern), or 99.000 for non-magnetic / uncomputable lines.

Lande factor estimation (for values VALD flags as missing, 99.000)
    A level's Lande factor is estimated from its J and its VALD term/config
    string according to the level's coupling scheme (LS, JJ or JK), by
    projecting the constituent orbital and spin moments onto J with g_s = 2.
    For JK the projection follows the K = J1 + L2, J = K + S2 coupling chain
    and holds for any S2.  A level whose quantum numbers cannot be parsed, or
    whose (K, J) violate the coupling triangle condition, keeps the 99.000
    sentinel.
------------------------------------------------------------------------------
"""

import argparse
import os
import re
import sys

# eV -> cm^-1 factor
CM_PER_EV = 8065.543937

# Letters used by VALD for orbital angular-momentum quantum numbers.
LVAL = ('S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
        'N', 'O', 'Q', 'R', 'T', 'U', 'V')
LVAL_LOWER = tuple(s.lower() for s in LVAL)

# All chemical element symbols, padded to width 2.
ALL_ELEMENTS = [
    'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
    'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
]

# Light elements for which LS coupling is expected to be a good approximation
# (used only when --LS-light is requested).
LIGHT_ELEMENTS = [
    'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na',
    'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc',
]

TOL = 1e-5  # tolerance for zeroing tiny shifts / dropping null components

# Fixed name of the optional --save-lande report (written to the current
# directory), holding the J values and Lande factors behind every output line.
LANDE_LOG = 'vald2z_lande.log'


# =============================================================================
#  Lande-factor estimation (LS / JJ / JK coupling)
# =============================================================================
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


def effective_lande(lande_lo, lande_up, jlo, jup):
    """Effective (centre-of-gravity) Lande factor of a transition."""
    return 0.5 * (lande_lo + lande_up) \
        + 0.25 * (lande_lo - lande_up) * (jlo * (jlo + 1.) - jup * (jup + 1.))


# =============================================================================
#  Zeeman splitting pattern
# =============================================================================
def _merge_coincident(shifts, intens):
    """Merge components that coincide in wavelength shift (within TOL) into a
    single component whose intensity is the sum of the merged intensities.

    Co-located delta-functions are physically indistinguishable, so this is an
    exact simplification.  When g_up == g_lo it turns an anomalous pattern into
    a normal Zeeman triplet automatically.  Order of first appearance is kept.
    """
    merged_g, merged_r = [], []
    for g, r in zip(shifts, intens):
        for k in range(len(merged_g)):
            if abs(g - merged_g[k]) < TOL:
                merged_r[k] += r
                break
        else:
            merged_g.append(g)
            merged_r.append(r)
    return merged_g, merged_r


def zeeman_pattern(jlo, jup, lande_lo, lande_up):
    """Compute the Zeeman splitting pattern of a transition.

    Returns ``(gpi, rpi, gsig, rsig)`` where the ``g*`` arrays hold the
    wavelength shifts (in Lorentz units) of the pi and sigma components and
    the ``r*`` arrays hold their normalised relative intensities.

    The shifts and strengths use the exact analytic (Condon-Shortley)
    expressions.  Components with vanishing intensity are dropped, near-zero
    shifts are snapped to exactly zero, and components that coincide in shift
    are merged (which turns a g_up == g_lo pattern into a normal triplet).
    """
    gpi, gsig, rpi, rsig = [], [], [], []

    # m = m_up runs from +Jup down to -Jup in unit steps (2*Jup+1 values).
    i = 0
    m = jup
    while abs(m) <= jup + TOL:
        gpi.append(m * (lande_up - lande_lo))
        gsig.append(m * lande_up - (m + 1.0) * lande_lo)
        if jlo == jup:
            rpi.append(m * m)
            rsig.append(0.25 * (jup - m) * (jup + m + 1.0))
        elif jlo > jup:  # Jup = Jlo - 1
            rpi.append((jup + 1.0) ** 2 - m * m)
            rsig.append(0.25 * (jup + m + 1.0) * (jup + m + 2.0))
        else:            # Jup = Jlo + 1
            rpi.append(jup ** 2 - m * m)
            rsig.append(0.25 * (jup - m) * (jup - m - 1.0))
        i += 1
        m = jup - i

    # Snap near-zero shifts to exactly zero (avoids "-0.0000" in the output).
    gpi = [0.0 if abs(g) < TOL else g for g in gpi]
    gsig = [0.0 if abs(g) < TOL else g for g in gsig]

    # Drop components with (essentially) zero relative intensity.
    keep_pi = [k for k in range(len(rpi)) if abs(rpi[k]) > TOL]
    gpi = [gpi[k] for k in keep_pi]
    rpi = [rpi[k] for k in keep_pi]
    keep_sig = [k for k in range(len(rsig)) if abs(rsig[k]) > TOL]
    gsig = [gsig[k] for k in keep_sig]
    rsig = [rsig[k] for k in keep_sig]

    # Wavelength-shift sign convention for the sigma components.
    gsig = [-g for g in gsig]

    # Merge components that coincide in wavelength shift, per branch.  For
    # g_up == g_lo all pi land at 0 and all sigma at +/-g, so this yields a
    # normal Zeeman triplet; otherwise shifts are distinct and nothing merges.
    gpi, rpi = _merge_coincident(gpi, rpi)
    gsig, rsig = _merge_coincident(gsig, rsig)

    # Normalise the relative intensities.
    spi = sum(rpi)
    ssig = sum(rsig)
    rpi = [r / spi for r in rpi]
    rsig = [r / ssig for r in rsig]
    return gpi, rpi, gsig, rsig


# =============================================================================
#  Input handling
# =============================================================================
class VALDFormatError(ValueError):
    """The input file is not a readable VALD long-format extraction."""


def open_vald(filename):
    """Open a VALD text file, transparently decompressing gzip/bzip2/xz.

    Compression is detected from the file's magic bytes rather than its name,
    so a mis-named compressed file still opens.  Text is decoded as Latin-1
    (a superset of the ASCII that VALD emits) so that stray bytes cannot
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


def read_vald_lines(infile):
    """Read all lines of a (possibly compressed) VALD extraction."""
    try:
        with open_vald(infile) as f:
            return [ln.rstrip('\r\n') for ln in f]
    except EOFError:
        raise VALDFormatError(
            '{}: compressed file is truncated or corrupt.'.format(infile))


# =============================================================================
#  VALD parsing helpers
# =============================================================================
_ELOW_UNIT_RE = re.compile(r'E[_ ]?low\s*\(([^)]*)\)', re.IGNORECASE)


def detect_energy_in_cm(header_lines):
    """Decide whether the excitation energies are in cm^-1 (vs eV).

    VALD writes the energy unit into the column header, e.g. ``E_low(eV)`` or
    ``E_low(cm-1)``.  Returns True if the unit indicates cm^-1 (so the values
    must be converted to eV), False if it is eV.  If no unit string is found,
    assumes eV (VALD's usual default) and returns False.
    """
    for line in header_lines:
        m = _ELOW_UNIT_RE.search(line)
        if m:
            unit = m.group(1).lower()
            if 'cm' in unit:      # cm-1, cm^-1, 1/cm, ...
                return True
            if 'ev' in unit:      # eV
                return False
    return False


def clean_config(raw):
    """Strip surrounding quotes/whitespace from a VALD config/term line."""
    return raw.strip().strip("'").strip()


def species_symbol(ion_field):
    """Extract the 2-character, space-padded element symbol from the first
    (quoted) field of a VALD data line, e.g. "'V 1'" -> "V ", "'Fe 1'" -> "Fe",
    "'CH 1'" -> "CH".
    """
    inner = ion_field.strip().strip("'").strip()  # e.g. "V 1"
    cut = inner.find(' ')
    if cut < 2:
        cut = 2
    return inner[:cut]


# =============================================================================
#  Output formatting
# =============================================================================
def fmt_pattern_line(count, values):
    """Format a shift/count line: ' N,  v1,  v2, ...' (I2 + F8.4 fields)."""
    return '{:2d},'.format(count) + ''.join('{:8.4f},'.format(v) for v in values)


def fmt_intensity_line(values):
    """Format an intensity line: '   i1,  i2, ...' (3x lead + F8.4 fields)."""
    return '   ' + ''.join('{:8.4f},'.format(v) for v in values)


def lande_log_lines(rows):
    """Build the contents of the --save-lande report from collected records.

    ``rows`` holds one ``(ion_field, ellen, wl, elo, jlo, jup, g_lo, g_up,
    magnetic)`` tuple per line written to the .z list, in the same order.  The
    first output line is the bare row count (so a reader can pick it up
    directly); every explanatory line, the column names included, starts with
    '#'.  Fields are comma separated and the species name is quoted exactly as
    in the .z list, so a Fortran list-directed read handles a row as-is.

    The Lande factors are written with five decimals rather than VALD's three:
    the computed (LS/JJ/JK) values are exact rationals such as 9/7, and three
    decimals would lose enough of them that a g_eff rebuilt downstream could
    drift by a few 1e-3 for high-J transitions.
    """
    out = ['{:d}'.format(len(rows))]
    #out.append('# J values and Lande factors used to compute the Zeeman patterns: one row per')
    #out.append('# line of the accompanying .z list, in the same order.')
    #out.append('# g_lo/g_up = 99 marks a Lande factor that could not be determined.')
    #out.append('# Mag = 1: line carries an explicit Zeeman pattern; 0: treated as non-magnetic.')
    ellen = max([6] + [len(ion_field) for ion_field, *_ in rows])
    out.append('#' + ' ' + 'Ion'.ljust(ellen - 2) + ' ' + 'Wl'.rjust(11)
               + ' ' + 'Elo'.rjust(8) + ' ' + 'Jlo'.rjust(5)
               + ' ' + 'Jup'.rjust(5) + ' ' + 'glo'.rjust(9)
               + ' ' + 'gup'.rjust(9) + ' ' + 'mag'.rjust(3))
    for ion_field, _ellen, wl, elo, jlo, jup, glo, gup, magnetic in rows:
        out.append('{ion:>{ellen}s},{wl:11.4f},{elo:8.4f},{jlo:5.1f},'
                   '{jup:5.1f},{glo:9.5f},{gup:9.5f},{mag:3d}'
                   .format(ion=ion_field, ellen=ellen, wl=wl, elo=elo,
                           jlo=jlo, jup=jup, glo=glo, gup=gup,
                           mag=1 if magnetic else 0))
    return out


# =============================================================================
#  Core conversion
# =============================================================================
class Stats:
    def __init__(self):
        self.processed = 0      # total input lines processed (X)
        self.written = 0        # lines written to the output
        self.vald_lande = 0     # magnetic lines using VALD Lande factors (Y)
        self.computed_lande = 0  # magnetic lines with >=1 computed factor (Z)
        self.computed_ls = 0    # of those, lines that used LS coupling (Z1)
        self.computed_jj = 0    # ... JJ coupling (Z2)
        self.computed_jk = 0    # ... JK coupling (Z3)
        self.n_h = 0            # hydrogen lines (M)
        self.n_mol = 0          # molecular lines (N)
        self.n_forbidden = 0    # E1-forbidden lines (K)
        self.n_failed = 0       # atomic lines with no computable Lande factor
        self.n_outside = 0      # lines outside the [wmin, wmax] range


def convert(infile, outfile, schemes=('LS', 'JJ', 'JK'),
            skip_failed=False, default_lande=None, skip_mol=False,
            skip_h=False, save_lin=False, ls_light=False,
            wmin=None, wmax=None, save_lande=False):
    """Convert ``infile`` (VALD extract-stellar long format) to ``outfile``
    (Synthmag .z format with explicit Zeeman patterns).

    The excitation-energy unit (eV or cm^-1) is detected automatically from
    the VALD column header; cm^-1 values are converted to eV.  The
    microturbulent velocity is taken from the VALD header (or defaults to
    0 km/s if absent).  ``wmin`` / ``wmax`` restrict the output to lines
    with wavelength >= wmin and/or <= wmax (Angstroms).

    The extraction-summary line of the output header (wl0, wl1, Nselected,
    Nprocessed, Vmicro) is recomputed for both input formats: Nselected is
    the number of lines actually written (after -wmin/-wmax and --skip-*),
    so downstream codes can rely on it.

    ``save_lande`` additionally writes LANDE_LOG, reporting the J values and
    Lande factors (VALD's or computed) behind every line of the output list.
    """

    lines = read_vald_lines(infile)

    if not lines:
        raise VALDFormatError('{}: input file is empty.'.format(infile))

    # -- Determine format: "extract all" has 'lande' in its first line. -------
    extract_all = 'lande' in lines[0].lower()

    stats = Stats()
    out_lines = []
    out_lines.append(' 0, 0, 0,')

    if extract_all:
        header_offset, nlin = _scan_extract_all(lines)
        if nlin == 0:
            raise VALDFormatError(
                '{}: no transition data found after the extract-all column '
                'headers.'.format(infile))
        vturb_file = None
        wl0_file = None            # no extraction summary in the input:
        wl1_file = None            # region taken from the output lines
        nproc = nlin
    else:
        # extract stellar: "wl0, wl1, Nlines, Nprocessed, Vmicro ..."
        fields = lines[0].split(',')
        try:
            wl0_file = float(fields[0])
            wl1_file = float(fields[1])
            nlin = int(fields[2])
        except (ValueError, IndexError):
            raise VALDFormatError(
                '{}: not a VALD long-format extraction: the first line is '
                'neither an extract-all column header nor an extract-stellar '
                'extraction summary: {!r}'
                .format(infile, lines[0].strip()[:70]))
        if nlin <= 0:
            raise VALDFormatError(
                '{}: the extraction summary reports no selected lines.'
                .format(infile))
        vturb_file = float(fields[4].split()[0]) if len(fields) > 4 else None
        header_offset = 3          # 3 header lines precede the data
        try:
            nproc = int(fields[3])
        except (ValueError, IndexError):
            nproc = nlin
        _validate_stellar_long(infile, lines, header_offset, nlin)

    stats.processed = nlin

    if vturb_file is not None:
        vturb = vturb_file
    else:
        vturb = 0.0
        print('Setting default microturbulence 0 km/s', file=sys.stderr)

    # -- Detect the excitation-energy unit from the VALD column header. -------
    energy_in_cm = detect_energy_in_cm(lines[:header_offset])
    if energy_in_cm:
        print('Energies given in cm^-1; converting to eV', file=sys.stderr)

    default_lande = 99.0 if default_lande is None else float(default_lande)

    # -- Standard short-VALD header.  The extraction-summary line is only a
    #    placeholder here; it is rewritten after the loop, once the number of
    #    lines actually retained (after -wmin/-wmax and --skip-*) is known.
    summary_idx = len(out_lines)
    out_lines.append('')
    out_lines.append('                                         '
                     'Damping parameters  Lande  Central')
    out_lines.append('Elm Ion  WL(A)    Excit(eV) Vmic log(gf) '
                     'Rad.   Stark  Waals factor  depth')

    # -- Optional --save-lin copy: a Zeeman-free list (.lin) that pySME and
    #    other short-VALD legacy readers accept.
    lin_lines = []
    lin_summary_idx = None
    if save_lin:
        lin_summary_idx = len(lin_lines)
        lin_lines.append('')
        lin_lines.append('                                         '
                         'Damping parameters  Lande  Central')
        lin_lines.append('Elm Ion  WL_air(A) Excit(eV) Vmic log(gf) '
                         'Rad.   Stark  Waals factor  depth')

    # Wavelengths of the first/last retained line (input is sorted in
    # wavelength, so these bound the output region when none is requested).
    wl_out_first = None
    wl_out_last = None

    # Per-line Lande/J records for the optional --save-lande report, kept in
    # step with the output list (one record per line actually written).
    lande_rows = []

    def emit(entry, wl, lin_entry=None, lande=None):
        """Append a retained line entry; track count and output region.

        ``lin_entry`` is the corresponding record for the --save-lin copy when
        it differs from the .z record (magnetic lines carry the effective
        Lande factor instead of the -90-g encoding); it defaults to ``entry``
        (identical record) for hydrogen, molecular and non-magnetic lines.
        ``lande`` is the record tuple for the --save-lande report.
        """
        nonlocal wl_out_first, wl_out_last
        out_lines.append(entry)
        if save_lin:
            lin_lines.append(entry if lin_entry is None else lin_entry)
        if save_lande and lande is not None:
            lande_rows.append(lande)
        stats.written += 1
        if wl_out_first is None:
            wl_out_first = wl
        wl_out_last = wl

    # -- Process every transition (4 file-lines each) ------------------------
    for k in range(nlin):
        base = header_offset + k * 4
        data = lines[base]
        cfg_lo_raw = lines[base + 1]
        cfg_up_raw = lines[base + 2]
        ref_raw = lines[base + 3]

        vals = data.split(',')
        try:
            wl = float(vals[1])
        except (ValueError, IndexError):
            raise _bad_data_line(infile, base, data)
        if (wmin is not None and wl < wmin) or (wmax is not None and wl > wmax):
            stats.n_outside += 1
            continue
        try:
            ion_field = vals[0].strip()
            loggf = float(vals[2])
            elo = float(vals[3])
            jlo = float(vals[4])
            eup = float(vals[5])
            jup = float(vals[6])
            lande_lo = float(vals[7])
            lande_up = float(vals[8])
            # vals[9] is VALD's mean Lande factor -- recomputed below, not used.
            rad = float(vals[10])
            stark = float(vals[11])
            waals = float(vals[12])
            depth = 0.5 if extract_all else float(vals[13])
        except (ValueError, IndexError):
            raise _bad_data_line(infile, base, data)

        if energy_in_cm:
            elo /= CM_PER_EV
            eup /= CM_PER_EV

        cfg_lo = clean_config(cfg_lo_raw)
        cfg_up = clean_config(cfg_up_raw)

        # Reference string, re-quoted: drop the leading "'_ " marker and re-wrap.
        ref = "' " + ref_raw.strip().lstrip("'_ ")

        # Field width for the (quoted) species name: at least 6.
        ellen = max(len(ion_field), 6)

        symbol = species_symbol(ion_field)
        is_element = symbol in ALL_ELEMENTS

        # -- Hydrogen and molecular lines: written directly, g = 99 (unless
        #    suppressed by --skip-H / --skip-mol).  These are never affected
        #    by --skip-failed or --default-Lande.
        if symbol == 'H ':
            stats.n_h += 1
            if not skip_h:
                emit(_fmt_data_line(ion_field, ellen, wl, elo,
                                    vturb, loggf, rad, stark,
                                    waals, 99.0, depth, ref,
                                    width=7), wl,
                     lande=(ion_field, ellen, wl, elo, jlo, jup,
                            lande_lo, lande_up, False))
            continue
        if not is_element:
            stats.n_mol += 1
            if not skip_mol:
                emit(_fmt_data_line(ion_field, ellen, wl, elo,
                                    vturb, loggf, rad, stark,
                                    waals, 99.0, depth, ref,
                                    width=7), wl,
                     lande=(ion_field, ellen, wl, elo, jlo, jup,
                            lande_lo, lande_up, False))
            continue

        # -- Electric-dipole (E1) selection rules: flag forbidden lines. -----
        zflag = True
        forbidden = False
        schemes_used = set()
        if (jlo == 0.0 and jup == 0.0) or abs(jlo - jup) > 1:
            forbidden = True
            zflag = False
        elif lande_up == 99.0 or lande_lo == 99.0:
            # -- Estimate missing Lande factor(s) from the coupling scheme. --
            # By default every scheme applies to every element.  With
            # --LS-light, LS is trusted only for light elements (LS is a poor
            # approximation for heavy species); JJ/JK still apply everywhere.
            use_schemes = schemes
            if ls_light and symbol not in LIGHT_ELEMENTS:
                use_schemes = tuple(s for s in schemes if s != 'LS')

            # Keep whatever could be determined (a level that fails retains the
            # 99.000 sentinel), so --save-lande can report a partial result.
            lo_missing = lande_lo == 99.0
            up_missing = lande_up == 99.0
            if lo_missing:
                lande_lo = estimate_lande(jlo, cfg_lo, use_schemes)
            if up_missing:
                lande_up = estimate_lande(jup, cfg_up, use_schemes)

            if lande_lo == 99.0 or lande_up == 99.0:
                zflag = False
            else:
                if lo_missing:
                    schemes_used.add(cfg_lo.split()[0])
                if up_missing:
                    schemes_used.add(cfg_up.split()[0])
                # A pattern with both g = 0 does not split -> non-magnetic.
                if lande_lo == 0.0 and lande_up == 0.0:
                    zflag = False

        # -- Non-magnetic outcome: forbidden / uncomputable / null split. ----
        if not zflag:
            if forbidden:
                stats.n_forbidden += 1
            else:
                stats.n_failed += 1
            if not skip_failed:
                emit(_fmt_data_line(ion_field, ellen, wl, elo,
                                    vturb, loggf, rad, stark,
                                    waals, default_lande, depth,
                                    ref, width=8), wl,
                     lande=(ion_field, ellen, wl, elo, jlo, jup,
                            lande_lo, lande_up, False))
            continue

        # -- Magnetic line: count VALD vs computed, then write. --------------
        if schemes_used:
            stats.computed_lande += 1
            if 'LS' in schemes_used:
                stats.computed_ls += 1
            if 'JJ' in schemes_used:
                stats.computed_jj += 1
            if 'JK' in schemes_used:
                stats.computed_jk += 1
        else:
            stats.vald_lande += 1

        geff = effective_lande(lande_lo, lande_up, jlo, jup)
        gcode = -90.0 - geff
        z_entry = _fmt_data_line(ion_field, ellen, wl, elo, vturb,
                                 loggf, rad, stark, waals, gcode, depth,
                                 ref, width=8, loggf_width=7)
        # .lin copy carries the actual effective Lande factor (same record
        # otherwise); the Zeeman pattern lines below are written only to .z.
        lin_entry = _fmt_data_line(ion_field, ellen, wl, elo, vturb,
                                   loggf, rad, stark, waals, geff, depth,
                                   ref, width=8, loggf_width=7)
        emit(z_entry, wl, lin_entry,
             lande=(ion_field, ellen, wl, elo, jlo, jup,
                    lande_lo, lande_up, True))

        gpi, rpi, gsig, rsig = zeeman_pattern(jlo, jup, lande_lo, lande_up)
        out_lines.append(fmt_pattern_line(len(gpi), gpi))
        out_lines.append(fmt_intensity_line(rpi))
        out_lines.append(fmt_pattern_line(len(gsig), gsig))
        out_lines.append(fmt_intensity_line(rsig))

    # -- Extraction-summary line: the retained-line count reflects the
    #    -wmin/-wmax window and the --skip-* options.  The region is the 
    #    requested window where given; the bound(s) not given default to 
    #    the input extraction region (extract stellar) or to the 
    #    first/last retained wavelength (extract all).
    wl0 = wmin if wmin is not None else \
        (wl0_file if wl0_file is not None else wl_out_first)
    wl1 = wmax if wmax is not None else \
        (wl1_file if wl1_file is not None else wl_out_last)
    summary = (
        ' {:.5f}, {:.5f}, {:d}, {:d}, {:.1f}, Wavelength region, '
        'lines selected, lines processed, Vmicro'
        .format(wl0 if wl0 is not None else 0.0,
                wl1 if wl1 is not None else 0.0,
                stats.written, nproc, vturb))
    out_lines[summary_idx] = summary
    if save_lin:
        lin_lines[lin_summary_idx] = summary

    # -- Copy the trailing reference/abundance section verbatim. -------------
    tail_start = header_offset + nlin * 4
    for ln in lines[tail_start:]:
        out_lines.append(ln)
        if save_lin:
            lin_lines.append(ln)
        if ln.strip() == '':
            break

    # Write line by line.  The .z gets the magnetic list; the .lin
    # (when requested) gets the Zeeman-free, effective-Lande copy.
    def _write(target, buf):
        with open(target, 'w') as f:
            for ln in buf:
                f.write(ln)
                f.write('\n')

    _write(outfile, out_lines)
    if save_lin:
        linpath = (outfile.rsplit('.z', 1)[0] + '.lin'
                   if outfile.endswith('.z') else outfile + '.lin')
        _write(linpath, lin_lines)
    if save_lande:
        _write(LANDE_LOG, lande_log_lines(lande_rows))

    return stats


def _fmt_data_line(ion_field, ellen, wl, elo, vturb, loggf, rad, stark, waals,
                   gcode, depth, ref, width=8, loggf_width=None):
    """Format a single transition data line.

    ``width`` is the field width (7 for H/molecular lines, 8 for atomic).
    ``loggf_width`` overrides the log(gf) width (the magnetic-line format uses
    F7.3 for log(gf) but F8.3 for the remaining numeric columns).
    """
    if loggf_width is None:
        loggf_width = width
    fmt = ('{ion:>{ellen}s},{wl:10.4f},{elo:8.4f},{vturb:4.1f},'
           '{loggf:{lw}.3f},{rad:{w}.3f},{stark:{w}.3f},{waals:{w}.3f},'
           '{g:{w}.3f},{depth:{w}.3f},{ref}')
    return fmt.format(ion=ion_field, ellen=ellen, wl=wl, elo=elo, vturb=vturb,
                      loggf=loggf, rad=rad, stark=stark, waals=waals, g=gcode,
                      depth=depth, ref=ref, w=width, lw=loggf_width)


def _bad_data_line(infile, base, data):
    """VALDFormatError for an unparseable transition data line."""
    return VALDFormatError(
        '{}: line {}: cannot parse as a long-format transition data line: '
        '{!r}'.format(infile, base + 1, data.strip()[:70]))


def _validate_stellar_long(infile, lines, header_offset, nlin):
    """Sanity-check an "extract stellar" input claiming ``nlin`` transitions.

    A VALD *short*-format extraction starts with the same summary line as the
    long format, so it passes the format detection; it is recognised here (by
    its Vmic column and its one-line-per-transition data records with only 11
    fields) and reported explicitly.  Also catches other malformed data lines
    and files shorter than the declared number of 4-line transition groups.
    """
    colhdr = lines[2].lower() if len(lines) > 2 else ''
    first_data = lines[header_offset] if len(lines) > header_offset else ''
    nfields = len(first_data.split(','))
    if nfields < 14:
        if 'vmic' in colhdr and 'j lo' not in colhdr:
            raise VALDFormatError(
                '{}: this is a VALD SHORT-format extraction. '
                'vald2z needs the LONG format.'.format(infile))
        raise VALDFormatError(
            '{}: line {}: a long-format transition line has 14 '
            'comma-separated fields, found {}: {!r}'
            .format(infile, header_offset + 1, nfields,
                    first_data.strip()[:70]))
    need = header_offset + 4 * nlin
    if len(lines) < need:
        raise VALDFormatError(
            '{}: file truncated? The extraction summary declares {} '
            'transitions ({} lines including the header), but the file has '
            'only {} lines.'.format(infile, nlin, need, len(lines)))


def _scan_extract_all(lines):
    """Handle the VALD "extract all" long format.

    Returns (header_offset, nlin).  The line count is obtained by scanning
    groups of 4 until the first group-boundary line that is not a transition
    data line.  Data lines start with a quoted species name ("'Fe 2', ..."),
    which distinguishes them from the footnote and 'References' lines that
    end the list (and from the "'_ ..." reference line inside each group,
    in case of misalignment).
    """
    # Two header lines precede the data in extract-all format.
    header_offset = 2
    nlin = 0
    i = header_offset
    while i + 3 < len(lines):
        ln = lines[i].lstrip()
        if not ln.startswith("'") or ln.startswith("'_"):
            break
        nlin += 1
        i += 4
    return header_offset, nlin


# =============================================================================
#  Command-line interface
# =============================================================================
def main(argv=None):
    p = argparse.ArgumentParser(
        description='Convert a VALD extract-stellar or extract-all '
                    'long-format line list into a Synthmag/Synmast/FastSpec '
                    'list with explicit Zeeman patterns.')
    p.add_argument('infile', help='input VALD long-format extraction '
                                  '(plain or gzip/bzip2/xz compressed)')
    p.add_argument('outfile', nargs='?',
                   help='output file (default: <infile stem>.z)')
    p.add_argument('-wmin', dest='wmin', type=float, default=None,
                   metavar='WL',
                   help='output only lines with wavelength >= WL (Angstroms)')
    p.add_argument('-wmax', dest='wmax', type=float, default=None,
                   metavar='WL',
                   help='output only lines with wavelength <= WL (Angstroms)')
    p.add_argument('--skip-failed', dest='skip_failed', action='store_true',
                   help='drop atomic lines for which no splitting pattern '
                        'could be calculated')
    p.add_argument('--default-Lande', dest='default_lande', type=float,
                   default=None, metavar='G',
                   help='effective Lande factor assigned to atomic lines that '
                        'keep no splitting pattern, instead of 99.0 '
                        '(H and molecular lines are not affected)')
    p.add_argument('--skip-mol', dest='skip_mol', action='store_true',
                   help='skip molecular lines')
    p.add_argument('--skip-H', dest='skip_h', action='store_true',
                   help='skip hydrogen lines')
    p.add_argument('--save-lin', dest='save_lin', action='store_true',
                   help='also write a Zeeman-free copy with a .lin extension, '
                        'readable by pySME and other short-VALD legacy codes')
    p.add_argument('--save-lande', dest='save_lande', action='store_true',
                   help='also write ' + LANDE_LOG + ', reporting the J values '
                        'and Lande factors (VALD or computed) behind every '
                        'output line, for use by other codes')
    p.add_argument('--LS-light', dest='ls_light', action='store_true',
                   help='trust LS-coupling Lande estimates only for light '
                        'elements (JJ/JK still applied to all elements; by '
                        'default LS is applied to all elements too)')
    args = p.parse_args(argv)

    outfile = args.outfile
    if outfile is None:
        # Strip a compression suffix first, then the format extension, so
        # e.g. list.dat.gz -> list.z (and list.dat -> list.z as before).
        stem, ext = os.path.splitext(args.infile)
        if ext.lower() in ('.gz', '.bz2', '.xz'):
            stem = os.path.splitext(stem)[0]
        outfile = (stem if stem else args.infile) + '.z'
        if outfile == args.infile:
            p.error('Same input and output file names. Aborting.')

    try:
        stats = convert(args.infile, outfile,
                        skip_failed=args.skip_failed,
                        default_lande=args.default_lande,
                        skip_mol=args.skip_mol, skip_h=args.skip_h,
                        save_lin=args.save_lin, ls_light=args.ls_light,
                        wmin=args.wmin, wmax=args.wmax,
                        save_lande=args.save_lande)
    except (VALDFormatError, OSError) as e:
        print('Error: {}'.format(e), file=sys.stderr)
        return 1

    def fate(n, noun, skipped):
        # H and molecular lines always keep g = 99 (non-magnetic); they are
        # not affected by --default-Lande.
        verb = 'skipped' if skipped else 'retained as non-magnetic'
        return '{:d} {:s} {:s}'.format(n, noun, verb)

    def fate_atomic(n, noun):
        # Atomic lines without a pattern (forbidden / uncomputable) are dropped
        # by --skip-failed, else assigned the --default-Lande value (99 = the
        # non-magnetic flag when the option is not given).
        if args.skip_failed:
            verb = 'skipped'
        elif args.default_lande is not None:
            verb = 'retained with {:.3f} Lande factor'.format(args.default_lande)
        else:
            verb = 'retained as non-magnetic'
        return '{:d} {:s} {:s}'.format(n, noun, verb)

    print('{:d} lines processed'.format(stats.processed), file=sys.stderr)
    if args.wmin is not None or args.wmax is not None:
        print('{:d} lines outside the wavelength range skipped'
              .format(stats.n_outside), file=sys.stderr)
    print('{:d} lines with VALD Lande factors'.format(stats.vald_lande),
          file=sys.stderr)
    print('{:d} lines with computed Lande factors ({:d}/{:d}/{:d} for LS/JJ/JK)'
          .format(stats.computed_lande, stats.computed_ls, stats.computed_jj,
                  stats.computed_jk), file=sys.stderr)
    print(fate(stats.n_h, 'H lines', args.skip_h), file=sys.stderr)
    print(fate(stats.n_mol, 'molecular lines', args.skip_mol), file=sys.stderr)
    print(fate_atomic(stats.n_forbidden, 'forbidden lines'), file=sys.stderr)
    print(fate_atomic(stats.n_failed, 'lines with no computable Lande factor'),
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
