"""
Fetch a research-grade spectral library subset via spextra and export to
plain 2-column ASCII (wavelength_A, flam) in spectral_library/.

Libraries (as used by the ESO / Gemini / HST ETCs):
  * Pickles 1998, PASP 110, 863      - stellar flux library (subset)
  * Kinney+ 1996, ApJ 467, 38        - galaxy templates (all 12)
  * Calzetti+ 1994                    - starburst internal-reddening series
  * Francis+ 1991 / spextra "agn"    - QSO composite, Seyferts, LINER
  * Dobos+ 2012, MNRAS 420, 1217     - SDSS composite galaxy spectra (subset)
"""

import numpy as np
import pathlib
from spextra import Spextrum
from synphot import units

OUT = pathlib.Path("spectral_library")
OUT.mkdir(exist_ok=True)

PICK = ["o5v", "b0v", "b5v", "a0v", "a5v", "f0v", "f5v", "g0v", "g2v",
        "g5v", "k0v", "k2v", "k5v", "m0v", "m2v", "m5v",
        "g8iii", "k3iii", "m0iii", "m5iii", "a0i", "g2i", "m2i"]
KC96 = ["elliptical", "bulge", "s0", "sa", "sb", "sc",
        "starb1", "starb2", "starb3", "starb4", "starb5", "starb6"]
AGN  = ["qso", "seyfert1", "seyfert2", "liner", "ngc1068"]
DOBO = ["SF1", "SF2", "SF3", "SF4", "RED0", "RED2", "RED4", "BG", "RG"]

REFS = {
    "pickles": "Pickles 1998, PASP 110, 863",
    "kc96":    "Kinney et al. 1996, ApJ 467, 38 (starb: Calzetti et al. 1994)",
    "agn":     "Francis et al. 1991, ApJ 373, 465 (qso composite) + spextra AGN set",
    "dobos":   "Dobos et al. 2012, MNRAS 420, 1217 (SDSS composites)",
}

n_ok = n_fail = 0
for lib, items in [("pickles", PICK), ("kc96", KC96), ("agn", AGN),
                   ("dobos", DOBO)]:
    for item in items:
        name = f"{lib}/{item}"
        fout = OUT / f"{lib}_{item.replace('.', 'p')}.txt"
        try:
            sp = Spextrum(name)
            wav = sp.waveset.to_value("Angstrom")
            flx = units.convert_flux(sp.waveset, sp(sp.waveset),
                                     units.FLAM).value
            # No UV trim: keep everything down to the atlas blue limit so
            # templates can be redshifted into the LRS band (KC96 reaches
            # ~1200 A, the Francis QSO composite ~800 A).
            sel = (wav >= 800) & (wav <= 10500)
            if sel.sum() < 50:
                raise ValueError("too few points in LRS range")
            hdr = (f"{name}  |  {REFS[lib]}\n"
                   f"columns: wavelength_Angstrom  flam_erg_s_cm2_A "
                   f"(absolute scale arbitrary - ETC renormalizes)")
            np.savetxt(fout, np.column_stack([wav[sel], flx[sel]]),
                       fmt="%10.3f %12.5e", header=hdr)
            n_ok += 1
            print(f"  ok  {name:22s} -> {fout.name}  ({sel.sum()} pts)")
        except Exception as e:
            n_fail += 1
            print(f"FAIL  {name:22s}  {type(e).__name__}: {str(e)[:80]}")

print(f"\n{n_ok} exported, {n_fail} failed -> {OUT}/")

# ---------------------------------------------------------------------------
# Full-coverage extensions (fix red-truncated CDBS/Francis templates):
#   * ESO ETC spliced QSO (Francis+91 UV/blue + Turler+99 red/NIR)
#   * Brown+2014 atlas galaxies/AGN with complete 0.09-30 um SEDs
# ---------------------------------------------------------------------------
BROWN = [
    ("brown/NGC5033", "brown_NGC5033_sy1.5",
     "Brown+14 atlas (ApJS 212, 18): NGC 5033, Seyfert 1.5"),
    ("brown/NGC4579", "brown_NGC4579_liner",
     "Brown+14 atlas: NGC 4579, LINER/Sy1.9"),
    ("brown/NGC3379", "brown_NGC3379_elliptical",
     "Brown+14 atlas: NGC 3379, E1 elliptical"),
    ("brown/NGC4450", "brown_NGC4450_sab",
     "Brown+14 atlas: NGC 4450, Sab bulge-dominated"),
    ("brown/NGC0628", "brown_NGC0628_sc",
     "Brown+14 atlas: NGC 628, Sc grand-design spiral"),
]
for name, stem, ref in BROWN:
    try:
        sp = Spextrum(name)
        wav = sp.waveset.to_value("Angstrom")
        flx = units.convert_flux(sp.waveset, sp(sp.waveset), units.FLAM).value
        sel = (wav >= 800) & (wav <= 10500)
        hdr = (f"{name}  |  {ref}\n"
               "columns: wavelength_Angstrom  flam_erg_s_cm2_A "
               "(absolute scale arbitrary - ETC renormalizes)")
        np.savetxt(OUT / f"{stem}.txt", np.column_stack([wav[sel], flx[sel]]),
                   fmt="%10.3f %12.5e", header=hdr)
        print(f"  ok  {name:16s} -> {stem}.txt")
    except Exception as e:
        print(f"FAIL  {name}: {type(e).__name__} {str(e)[:80]}")

# ESO spliced QSO is served as a plain .dat (nm, FLAM), not a Spextrum item
try:
    import pooch
    f = pooch.retrieve("https://scopesim.univie.ac.at/spextra/database/"
                       "libraries/etc/misc/qso.dat", known_hash=None,
                       path=OUT, fname="_qso_raw.dat")
    d = np.genfromtxt(f, comments="#", invalid_raise=False)
    d = d[~np.isnan(d).any(axis=1)]
    wav_aa = d[:, 0] * 10.0
    sel = (wav_aa >= 800) & (wav_aa <= 10500)
    hdr = ("etc/misc/qso  |  ESO ETC QSO composite: Francis+91 (800-6000 A) "
           "spliced with Turler+99 (6000 A - 2 um)\n"
           "columns: wavelength_Angstrom  flam_erg_s_cm2_A "
           "(absolute scale arbitrary - ETC renormalizes)")
    np.savetxt(OUT / "agn_qso_ext.txt",
               np.column_stack([wav_aa[sel], d[sel, 1]]),
               fmt="%10.3f %12.5e", header=hdr)
    print("  ok  etc/misc/qso     -> agn_qso_ext.txt")
except Exception as e:
    print(f"FAIL  etc/misc/qso: {type(e).__name__} {str(e)[:80]}")
