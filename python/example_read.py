#!/usr/bin/env python3.12
"""
example_read.py -- read and plot FastSpec output files with readfs.py.

Opens three figure windows, one per example under ../examples/:

  1. molecular-synthesis : FeH flux spectrum (unnormalised and normalised)
  2. atomic-synthesis    : solar flux around the Mg I b triplet, with the
                           strongest lines identified from the line-strength file
  3. magnetic-synthesis  : Stokes I and V versus limb angle around the same region

Paths are resolved relative to this file, so it can be run from any working
directory:

    python3 example_read.py
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                    # make "import readfs" work anywhere
import readfs

EXROOT = os.path.normpath(os.path.join(HERE, os.pardir, 'examples'))

# Mg I b triplet region (5167.3, 5172.7, 5183.6 A), used for the zoomed plots
MG_LO, MG_HI = 5160.0, 5190.0
# number of strongest lines to identify in the atomic zoom window (20-50)
N_STRONGEST = 30

_ROMAN = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI'}
def ion_label(species, ion):
    return f"{species} {_ROMAN.get(int(ion), str(int(ion)))}"


# =====================================================================
# 1. Molecular synthesis (FeH): disk-integrated flux
# =====================================================================
f = readfs.read_flux(os.path.join(EXROOT, 'molecular-synthesis', 'fs_flux.bin'))

fig1, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
ax[0].plot(f['wave'], f['sfl'], 'k', lw=0.5, label='line flux')
ax[0].plot(f['wave'], f['sfc'], 'r', lw=0.8, label='continuum')
ax[0].set_ylabel('Flux')
ax[0].legend(fontsize=8, loc='lower right')
ax[1].plot(f['wave'], f['sfl'] / f['sfc'], 'k', lw=0.5)
ax[1].set_ylabel('Normalised flux')
ax[1].set_xlabel('Wavelength (Å)')
ax[1].set_xlim(f['wave'][0], f['wave'][-1])
fig1.suptitle('1. Molecular synthesis (FeH): disk-integrated flux')
fig1.tight_layout()


# =====================================================================
# 2. Atomic synthesis (Sun): flux + line IDs, Mg I b triplet region
# =====================================================================
f = readfs.read_flux(os.path.join(EXROOT, 'atomic-synthesis', 'fs_flux.bin'))
lin = readfs.read_linst(os.path.join(EXROOT, 'atomic-synthesis', 'fs_linst.bin'))

fig2, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
ax[0].plot(f['wave'], f['sfl'], 'k', lw=0.6, label='line flux')
ax[0].plot(f['wave'], f['sfc'], 'r', lw=0.8, label='continuum')
ax[0].set_ylabel('Flux')
ax[0].legend(fontsize=8, loc='lower right')

norm = f['sfl'] / f['sfc']
ax[1].plot(f['wave'], norm, 'k', lw=0.6)
ax[1].set_ylabel('Normalised flux')
ax[1].set_xlabel('Wavelength (Å)')
ax[1].set_ylim(0.0, 1.75)

# identify the strongest lines in the zoom window.  Each ID is a bar rising from
# the line core (on the spectrum); its length ABOVE the continuum is proportional
# to the line strength, with the label at the top of the bar.
inwin = np.where((lin['wave'] >= MG_LO) & (lin['wave'] <= MG_HI))[0]
strongest = inwin[np.argsort(lin['strength'][inwin])[::-1][:N_STRONGEST]]
strongest = strongest[np.argsort(lin['wave'][strongest])]
BAR_H = 0.40                                # bar height above the continuum at strength = 1
print(f"{len(strongest)} strongest lines identified in {MG_LO:.0f}-{MG_HI:.0f} Å:")
for i in strongest:
    w, s = lin['wave'][i], lin['strength'][i]
    lab = ion_label(lin['species'][i], lin['ion'][i])
    core = np.interp(w, f['wave'], norm)
    top = 1.0 + BAR_H * s               
    ax[1].plot([w, w], [core, top], color='tab:blue', lw=0.9)
    ax[1].text(w, top + 0.02, lab, rotation=90, va='bottom', ha='center',
               fontsize=6, color='red')
    print(f"  {lab:6s} {w:9.3f} Å  depth {s:.3f}")

ax[0].set_xlim(MG_LO, MG_HI)                # sharex zooms both panels
fig2.suptitle('2. Atomic synthesis (Sun): flux and line IDs, Mg I b triplet')
fig2.tight_layout()


# =====================================================================
# 3. Magnetic synthesis: Stokes I and V versus limb angle, same region
# =====================================================================
m = readfs.read_spec(os.path.join(EXROOT, 'magnetic-synthesis', 'fs_spec.bin'))
if not m['magnetic']:
    sys.exit('magnetic-synthesis/fs_spec.bin does not contain Stokes spectra')

nmu = len(m['xmu'])
sel = (m['wave'] >= MG_LO) & (m['wave'] <= MG_HI)
I = m['spl'] / m['spc']                     # Stokes I / Ic
V = m['spv'] / m['spc']                     # Stokes V / Ic

off_I = 1.15 * np.ptp(I[:, sel])
off_V = 1.15 * np.ptp(V[:, sel])
colors = plt.cm.tab10(np.arange(nmu))

fig3, ax = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
for imu in range(nmu):
    ax[0].plot(m['wave'], I[imu] + imu * off_I, color=colors[imu], lw=0.8)
    ax[1].plot(m['wave'], V[imu] + imu * off_V, color=colors[imu], lw=0.8)
    for a, base in ((ax[0], 1.0 + imu * off_I), (ax[1], imu * off_V)):
        a.text(1.01, base, f"μ = {m['xmu'][imu]:.2f}", color=colors[imu],
               transform=a.get_yaxis_transform(), fontsize=7, va='center',
               ha='left', clip_on=False)
ax[0].set_ylabel('Stokes I / Ic  (offset by μ)')
ax[1].set_ylabel('Stokes V / Ic  (offset by μ)')
ax[1].set_xlabel('Wavelength (Å)')
ax[0].set_xlim(MG_LO, MG_HI)
fig3.suptitle('3. Magnetic synthesis: Stokes I and V versus limb angle')
fig3.tight_layout()
fig3.subplots_adjust(right=0.9)

plt.show()
