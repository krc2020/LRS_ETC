"""
Export the real SDSS u/g/r/i/z filter responses (Doi et al. 2010, via the
speclite package's bundled sdss2010 curves, which include the reference
atmosphere at airmass 1.3) to plain ASCII in filters/.

    pip install speclite
    python3 fetch_sdss_filters.py
"""

import numpy as np
import pathlib
import speclite.filters as sf

OUT = pathlib.Path("filters")
OUT.mkdir(exist_ok=True)

for band in "ugriz":
    f = sf.load_filter(f"sdss2010-{band}")
    wav = np.asarray(f.wavelength, float)
    resp = np.asarray(f.response, float)
    hdr = (f"SDSS {band}-band filter response (sdss2010, Doi et al. 2010, "
           f"AJ 139, 1628)\n"
           "via speclite; includes atmosphere at airmass 1.3 at APO\n"
           "columns: wavelength_Angstrom  response(photon-counting)")
    np.savetxt(OUT / f"sdss_{band}.txt", np.column_stack([wav, resp]),
               fmt="%9.1f %12.6e", header=hdr)
    print(f"wrote filters/sdss_{band}.txt ({len(wav)} rows)")
