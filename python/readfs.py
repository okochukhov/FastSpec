#!/usr/bin/env python3.12
"""
readfs.py -- read FastSpec binary output files (Python counterpart of the
IDL readfs.pro).  All outputs are unformatted Fortran stream files written
in native (little-endian) byte order.

  read_linst : line-strength list          (fs_linst.bin)
  read_spec  : intensity / Stokes spectra  (fs_spec.bin)
  read_flux  : disk-integrated flux        (fs_flux.bin)

The same functions handle the non-magnetic and magnetic scenarios: the
magnetic variant of fs_spec.bin (Stokes IQUV + continuum instead of
line + continuum intensity) is detected automatically, while fs_linst.bin
and fs_flux.bin have identical layouts in both cases.

Wavelengths are in Angstroms.  Intensities and fluxes are in the internal
FastSpec units; divide by the corresponding continuum spectrum to obtain
continuum-normalised spectra (this also applies to Stokes Q, U, V, which
are stored in the same units as Stokes I).

Typical usage:

    import readfs
    s = readfs.read_spec('fs_spec.bin')
    plot(s['wave'], s['spl'][0] / s['spc'][0])   # normalised I at mu = xmu[0]
"""

import numpy as np


def read_linst(file='fs_linst.bin'):
    """Read a FastSpec line-strength list (fs_linst.bin).

    The file contains one block per processed line list; the blocks are
    concatenated here.  Returns a dict of arrays, one element per retained
    spectral line:
      species  : species name, e.g. 'Fe' or 'OH'
      ion      : ionisation stage (1 = neutral)
      wave     : central wavelength [A]
      strength : central line depth 1 - I/Ic at the line centre
                 (disk-integrated when several mu angles are computed)
    """
    species, ion, wave, strength = [], [], [], []
    with open(file, 'rb') as f:
        while True:
            nrec = np.fromfile(f, '<i4', 1)
            if nrec.size == 0:
                break
            nrec = nrec[0]
            species.append(np.char.strip(np.fromfile(f, 'S8', nrec).astype('U8')))
            ion.append(np.fromfile(f, '<i4', nrec))
            wave.append(np.fromfile(f, '<f8', nrec))
            strength.append(np.fromfile(f, '<f4', nrec))
    return {'species': np.concatenate(species), 'ion': np.concatenate(ion),
            'wave': np.concatenate(wave), 'strength': np.concatenate(strength)}


def read_spec(file='fs_spec.bin'):
    """Read FastSpec intensity spectra (fs_spec.bin).

    Returns a dict:
      wave     : wavelength [A], (Nwl)
      xmu      : mu angles, (Nmu), ordered disk centre -> limb; for spherical
                 models mu refers to the outermost model layer
      spl      : line intensity, (Nmu, Nwl); Stokes I in the magnetic case
      spc      : continuum intensity, (Nmu, Nwl)
      spq, spu, spv : Stokes Q, U, V, (Nmu, Nwl); magnetic output only
      magnetic : True if the file contains Stokes spectra

    spl[i] is the spectrum at limb angle xmu[i].  The magnetic variant is
    recognised by the number of arrays in the file (5 instead of 2).
    """
    with open(file, 'rb') as f:
        nwl, nmu = np.fromfile(f, '<i4', 2)
        xmu = np.fromfile(f, '<f8', nmu)
        wave = np.fromfile(f, '<f8', nwl)
        # remaining data: 2 (I, Ic) or 5 (I, Ic, Q, U, V) arrays,
        # each written mu-first (Fortran order)
        sp = np.fromfile(f, '<f8').reshape(-1, nwl, nmu).transpose(0, 2, 1)
    out = {'wave': wave, 'xmu': xmu, 'spl': sp[0], 'spc': sp[1],
           'magnetic': len(sp) == 5}
    if out['magnetic']:
        out['spq'], out['spu'], out['spv'] = sp[2], sp[3], sp[4]
    return out


def read_flux(file='fs_flux.bin'):
    """Read a FastSpec flux spectrum (fs_flux.bin).

    Returns a dict:
      wave : wavelength [A], (Nwl)
      sfl  : disk-integrated line flux (Stokes I in the magnetic case)
      sfc  : disk-integrated continuum flux
    """
    with open(file, 'rb') as f:
        nwl = np.fromfile(f, '<i4', 1)[0]
        wave = np.fromfile(f, '<f8', nwl)
        sfl = np.fromfile(f, '<f8', nwl)
        sfc = np.fromfile(f, '<f8', nwl)
    return {'wave': wave, 'sfl': sfl, 'sfc': sfc}
