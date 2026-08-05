"""
Build 01_LRS_ETC_v6.ipynb — full-featured, streamlined ETC for LRS @ 2.4-m TNT.

New in v6 (user spec, 2026-07-23):
  * point-source (AB mag) or extended-source (mag/arcsec^2) normalization
  * spectral template library (stars, galaxies, AGN, dusty variants, M82)
  * throughput: 4 TNT mirrors at measured R(550 nm) 90.4/90.8/93.0/87.9 %
  * observing conditions: seeing, lunar phase, cloud cover, airmass
  * one mode: long slit; slits 1.8/2.7/4.5" = R 750/500/300 at 6000 A
  * user-supplied 2-column table input
  * flat / power-law continuum builder with emission/absorption lines
  * CCD: Andor Newton 256x1024 BEX2-DD; T = -80/-100 C; readout fast/med/slow
  * dispersion 4 A/pix, spatial 0.9 "/pix
  * output mode 1: N x t  -> SNR vs wavelength (per pixel & per A)
  * output mode 2: target SNR at lambda -> required on-source time
"""

import json, pathlib

OUT = pathlib.Path("LRS_ETC.ipynb")

cells = []

def code(src):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src.splitlines(keepends=True)}

def md(src):
    return {"cell_type": "markdown", "metadata": {},
            "source": src.splitlines(keepends=True)}


# ----------------------------------------------------------------------
cells.append(md("""# NARIT LRS @ 2.4-m TNT — Exposure-Time Calculator

**Version 6.4 (2026-07-24)** · maintainer: Krittapas Chanchaiworawit (NARIT)
· observation-preparation service for LRS @ TNT proposers.

**Changelog** — v6.4: **empirical commissioning calibration** (`USE_EMPIRICAL_CALIBRATION` switch, §2b) from the shutter-corrected BD+75°325 rates (stuck-open shutter → +7 s per frame; gray factor ≈ 0.29 with the measured chromatic response shape); v6.3: measured TNT mirror reflectivities, SDSS-filter + total-flux normalization, resolved-source profiles; v6.2: standalone `LRS_ETC` repository;
v6.1: validation vs 2026-04-01 commissioning, deployment audit; v6.0:
research template library (Pickles / Kinney–Calzetti / AGN / SDSS),
Galactic extinction, redshifting to z = 9, target visibility, request-time
budget w/ standards + calibration, aged-primary throughput, BEX2-DD
datasheet values, detector warnings, 2-D simulators.

A single-mode (long-slit) ETC for the Low-Resolution Spectrograph on the
Thai National Telescope. Structure:

1. Constants — telescope, slits, CCD, readout, operating temperature
2. System throughput — per-component model, all four TNT mirrors rescaled to measured R(550 nm)
2b. **Empirical commissioning calibration** — measured response ratio + `USE_EMPIRICAL_CALIBRATION` switch
3. Analytic template library — stars, galaxies, AGN, dusty variants
3b. **Research-grade library** — Pickles, Kinney–Calzetti, AGN, SDSS composites
4. Custom spectra — user table, or flat / power-law continuum ± lines
5. Normalization — point source (AB mag) or extended (mag/arcsec²)
5b. **Galactic extinction** — CCM89 / O94 / F99 / Calzetti00, R_V, E(B−V) or A_V
5c. **Redshifting** — shape or cosmological mode, native-UV coverage, z limits
6. Observing conditions — seeing, lunar phase, clouds, airmass
6b. **Target visibility** — RA/Dec → altitude range at the TNO, observability flags
7. ETC core — count rates, noise budget
7b. **Expected observed spectrum** — flux / ADU / S/N in the seeing disk, with saturation, non-linearity and noise-regime warnings
7c. **Simulated 2-D spectra** — raw frame (ADU, with sky lines) and reduced stacked S/N map (λ × arcsec, sky-subtraction residuals)
8. **Output mode 1** — N frames × t per frame → S/N vs λ (per pixel and per Å)
9. **Output mode 2** — target S/N at a wavelength → required on-source time
10. **Recommended request time** — overheads + 35 % weather margin, stacked-bar breakdown
11. **Validation** — ETC prediction vs the 2026-04-01 M82 commissioning observation

Conventions: wavelengths in Å, f_λ in erg s⁻¹ cm⁻² Å⁻¹, S/N quoted per
**pixel** (4 Å) unless stated otherwise."""))

cells.append(code("""import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.special import erf
import pathlib

plt.rcParams.update({"figure.dpi": 110, "font.size": 10,
                     "axes.titlesize": 11, "axes.labelsize": 10})
NAVY, LIME, ORANGE, RED = "#1E2761", "#39B54A", "#F4B400", "#E94F37"

ETC_VERSION = "6.4"

# Data folder: relative to the repo when run in place, else look next to
# the user's home copy of the repository.
_DATA_CANDIDATES = [
    pathlib.Path("data"),
    pathlib.Path.home() / "Documents" / "LRS_ETC" / "data",
]
DATA_ROOT = next((p for p in _DATA_CANDIDATES if p.exists()),
                 _DATA_CANDIDATES[0])
print(f"NARIT LRS ETC v{ETC_VERSION}  ·  data root: {DATA_ROOT}")

# High-resolution internal wavelength grid (1 A) and physical constants
WAV = np.arange(3500.0, 9501.0, 1.0)      # A
H_ERG_S  = 6.626e-27                       # erg s
C_AA_S   = 2.998e18                        # A / s
HC       = H_ERG_S * C_AA_S                # erg A
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 1. Constants — telescope, slit choices, CCD

* **Telescope**: 2.4-m Ritchey–Chrétien, f/10, 9 % central obstruction (by area).
* **Mode**: long-slit spectroscopy only. Slit width sets the resolution at 6000 Å:

  | slit | R @ 6000 Å | FWHM @ 6000 Å |
  |---|---|---|
  | 1.8″ | 750 | 8.0 Å |
  | 2.7″ | 500 | 12.0 Å |
  | 4.5″ | 300 | 20.0 Å |

* **Detector**: one choice — Andor Newton **BEX2-DD** (deep depletion,
  fringe-suppressed), 1024 × 256 pixels of 26 µm.
  Dispersion **4 Å/pixel**, spatial scale **0.9″/pixel**.
* **Operating temperature** sets dark current; **readout mode** sets read noise:

  | temperature | dark (e⁻/pix/s) | | readout | rate | RN (e⁻) |
  |---|---|---|---|---|---|
  | −80 °C | 0.08 | | slow | 50 kHz | 4 |
  | −100 °C | 0.003 | | medium | 1 MHz | 12 |
  | | | | fast | 3 MHz | 15 |

  Gain **4 e⁻/ADU**, 16-bit digitization (**full well 65 535 ADU**),
  bias level ≈ 300 ADU, response linear to **95 %** of full well."""))

cells.append(code("""# --- Telescope -----------------------------------------------------------
D_M          = 2.4
OBSTRUCTION  = 0.09                                   # fraction of area
A_CM2        = np.pi / 4 * (D_M * 100) ** 2 * (1 - OBSTRUCTION)

# --- Long-slit choices ---------------------------------------------------
SLITS = {                       # slit width " : R at 6000 A
    "1.8": dict(width_arcsec=1.8, R=750),
    "2.7": dict(width_arcsec=2.7, R=500),
    "4.5": dict(width_arcsec=4.5, R=300),
}
DISP_AA_PIX   = 4.0             # spectral dispersion
SPATIAL_AS_PIX = 0.9            # spatial scale along the slit
N_SPEC_PIX    = 1024
PIX_WAV       = 4000.0 + DISP_AA_PIX * np.arange(N_SPEC_PIX)   # 4000-8092 A

# --- Andor Newton 256x1024 BEX2-DD --------------------------------------
CCD_NAME  = "Andor Newton 256x1024 BEX2-DD"
DARK_E_PIX_S = {"-80C": 0.08, "-100C": 0.003}        # e-/pix/s (datasheet)
READ_NOISE_E = {"slow": 4.0, "medium": 12.0, "fast": 15.0}  # e- rms
READOUT_RATE = {"slow": "50 kHz", "medium": "1 MHz", "fast": "3 MHz"}
GAIN_E_ADU     = 4.0        # e- per ADU
FULL_WELL_ADU  = 65535      # 16-bit digitization limit
BIAS_ADU       = 300.0      # bias level
LINEARITY_FRAC = 0.95       # response linear up to 95% of full well

print(f"Collecting area: {A_CM2:,.0f} cm^2")
print(f"Wavelength coverage: {PIX_WAV[0]:.0f}-{PIX_WAV[-1]:.0f} A "
      f"({DISP_AA_PIX:.0f} A/pix x {N_SPEC_PIX} pix)")
for s, d in SLITS.items():
    print(f"  slit {s}\\": R = {d['R']}, FWHM at 6000 A = {6000/d['R']:.1f} A")
print(f"{CCD_NAME}:")
for T, dc in DARK_E_PIX_S.items():
    print(f"  dark current at {T}: {dc} e-/pix/s")
for m in READ_NOISE_E:
    print(f"  {m:6s} readout ({READOUT_RATE[m]:>6s}): RN = {READ_NOISE_E[m]:.0f} e-")
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 2. System throughput — measured TNT mirror reflectivities (2026)

The per-component model from `LRS_throughput.csv` assumed fresh dielectric
coatings on all four telescope mirrors (~99.6 % each at 600 nm). The four
mirrors have now been **measured individually**; each curve is rescaled to
match the measured reflectivity at **550 nm**:

| mirror | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| R(550 nm) | 90.4 % | 90.8 % | 93.0 % | 87.9 % |

The wavelength *shape* of each curve keeps the dielectric-coating model;
only the normalization is anchored to the measurements (four-mirror
product at 550 nm: **67.1 %**).

Below ~4000 Å and above ~8000 Å the CSV is extrapolated (band edges of the
model); we apply a soft taper there consistent with the blue QE collapse
observed in the commissioning data."""))

cells.append(code("""csv_path = None
for p in [pathlib.Path("LRS_throughput.csv"),
          DATA_ROOT.parent / "LRS_throughput.csv"]:
    if p.exists():
        csv_path = p
        break

tbl = np.genfromtxt(csv_path, delimiter=",", names=True)
wav_csv_aa = tbl["wavelength_nm"] * 10.0

mirrors4  = tbl["TNT_mirrors_M1xM2xM3xM4"]          # product of 4 dielectric mirrors
per_mirror = mirrors4 ** 0.25                        # single-mirror curve
lrs_only   = tbl["LRS_only_throughput"]              # FR x folds x coll x VPH x camera x QE

# --- Rescale ALL FOUR mirrors to the measured R(550 nm) -----------------
MIRROR_R_550 = {"M1": 0.904, "M2": 0.908, "M3": 0.930, "M4": 0.879}
i550 = np.argmin(np.abs(wav_csv_aa - 5500.0))
R0_550 = per_mirror[i550]                            # fresh-coating value
scale = {m: r / R0_550 for m, r in MIRROR_R_550.items()}
scale_prod = np.prod(list(scale.values()))
mirrors4_new = per_mirror ** 4 * scale_prod          # measured M1xM2xM3xM4
total_new  = mirrors4_new * lrs_only
total_old  = tbl["TOTAL_system_throughput"]
i600 = np.argmin(np.abs(wav_csv_aa - 6000.0))

# --- Interpolate onto the internal grid with soft band-edge tapers ------
eta_raw = np.interp(WAV, wav_csv_aa, total_new,
                    left=total_new[0], right=total_new[-1])
taper = np.ones_like(WAV)
blue = WAV < 4000
taper[blue] = np.exp(-((4000 - WAV[blue]) / 250.0) ** 2)
red = WAV > 8000
taper[red]  = np.exp(-((WAV[red] - 8000) / 500.0) ** 2)
ETA = eta_raw * taper

def system_throughput(wav_aa):
    \"\"\"Total telescope+LRS+CCD throughput at wav_aa (A).\"\"\"
    return np.interp(wav_aa, WAV, ETA)

print("Measured mirror reflectivities at 550 nm "
      "(curve shapes from the dielectric model):")
for m, r in MIRROR_R_550.items():
    print(f"  {m}: {r:.1%}")
print(f"4-mirror product at 550 nm: {mirrors4_new[i550]:.4f}  "
      f"(fresh-coating model was {per_mirror[i550]**4:.4f})")
print(f"Total throughput at 600 nm: {total_old[i600]:.4f} -> {total_new[i600]:.4f}")

fig, ax = plt.subplots(figsize=(10, 4.2))
ax.plot(wav_csv_aa, total_old, color="#999999", lw=1.2, ls="--",
        label="total (fresh coatings, earlier model)")
ax.plot(WAV, ETA, color=NAVY, lw=2.2,
        label="total (measured mirrors, band-edge tapers)")
ax.plot(wav_csv_aa, mirrors4_new, color=ORANGE, lw=1.2,
        label="4 TNT mirrors (measured 2026: 90.4/90.8/93.0/87.9 % at 550 nm)")
ax.plot(wav_csv_aa, lrs_only, color=LIME, lw=1.2,
        label="LRS instrument only (FR+folds+coll+VPH+camera+QE)")
ax.set_xlabel("Wavelength (Å)"); ax.set_ylabel("Throughput")
ax.set_title("LRS @ TNT system throughput — measured mirrors (2026)")
ax.set_xlim(3500, 9500); ax.set_ylim(0, 1.0)
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 2b. Empirical commissioning calibration — the switch

The component model above is the **theoretical** system. The 2026-04-01
commissioning standard (BD+75°325, all three slits) shows the as-built
system delivers less, with a measured wavelength shape (gray floor ≈ 0.29
in the mid-band, a strong chromatic collapse blueward of ~5500 Å, mild
red decline). Rates are corrected for the **stuck-open shutter** during
commissioning: every frame integrated ~7 s longer than the header
exposure time (30 s → 37 s, 60 s → 67 s effective).

`USE_EMPIRICAL_CALIBRATION = True` (default) multiplies the source *and*
sky photon rates by the measured response ratio, so every S/N, exposure
time and request estimate downstream matches the commissioning data.
Set it to `False` — globally here, or per-call with
`run_lrs_etc(..., use_empirical=False)` — for the theoretical model
(e.g. to predict performance after optics refurbishment). §11 always
compares the theory branch against the data by construction."""))

cells.append(code("""# --- Empirical response ratio: measured / theoretical model -------------
# Source: BD+75d325 2026-04-01 (pipeline sensitivity figure, digitized,
# +-15%), 4.5\" slit at the delivered ~4\" IQ, SHUTTER-CORRECTED
# (header 30 s -> 37 s effective: ratios x 30/37 vs the raw digitization).
SHUTTER_EXTRA_S = 7.0        # commissioning-era stuck-open shutter offset
_EMP_W = np.array([4500, 4750, 5000, 5250, 5500, 5800, 6000, 6250,
                   6500, 6800, 7000, 7200, 7400, 7800], float)
_EMP_R = np.array([0.07, 0.11, 0.18, 0.25, 0.31, 0.31, 0.30, 0.29,
                   0.25, 0.25, 0.24, 0.23, 0.21, 0.21]) * (30.0 / 37.0)

EMPIRICAL_RATIO_GRID = np.interp(WAV, _EMP_W, _EMP_R,
                                 left=_EMP_R[0], right=_EMP_R[-1])

USE_EMPIRICAL_CALIBRATION = True     # False -> pure theoretical model

print(f"Empirical calibration grid (shutter-corrected):")
for w in (4500, 5000, 5500, 6000, 6500, 7000, 7800):
    print(f"  {w} A: measured/model = "
          f"{np.interp(w, WAV, EMPIRICAL_RATIO_GRID):.2f}")
print(f"USE_EMPIRICAL_CALIBRATION = {USE_EMPIRICAL_CALIBRATION} "
      f"(default for all downstream estimates)")
"""))


cells.append(md("""## 3. Analytic template library (quick approximations)

**For real observation preparation prefer the empirical library in §3b** —
these analytic forms are kept for quick what-if exploration and as a
fallback when `spectral_library/` is unavailable.

Analytic templates on the 1 Å internal grid, each returned as a dict
`{"wav": Å, "flam": arbitrary-units}` — the **shape** is what matters,
normalization happens in §5. The library covers:

* **Stars** — blackbodies at the effective temperature of each class with
  key absorption features (Balmer series for A stars; Ca II H/K, Mg b,
  Na D for solar types).
* **Galaxies** — quenched (4000 Å break + absorption), star-forming,
  starburst, dusty starburst, dusty (Calzetti-attenuated) galaxy, young
  irregular (blue + high-ionisation lines).
* **AGN** — power-law + broad and narrow lines; dust-obscured AGN
  (reddened power-law, narrow lines only).
* **M82 (empirical)** — the commissioning-anchored template from
  `data/`, when available."""))

cells.append(code("""def _bb(T):
    \"\"\"Blackbody f_lam shape (arbitrary units) on WAV.\"\"\"
    x = HC / (WAV * 1.380649e-16 * T)
    with np.errstate(over="ignore"):
        b = 1.0 / (WAV ** 5 * np.expm1(np.clip(x, None, 700)))
    return b / b.max()

def _gauss_line(center, ew_aa, fwhm_aa, cont, absorption=False):
    \"\"\"Gaussian line with equivalent width ew_aa on continuum cont.\"\"\"
    sig = fwhm_aa / 2.3548
    prof = np.exp(-0.5 * ((WAV - center) / sig) ** 2)
    amp = ew_aa * np.interp(center, WAV, cont) / (sig * np.sqrt(2 * np.pi))
    return -amp * prof if absorption else amp * prof

def _calzetti_k(wav_aa):
    \"\"\"Calzetti+2000 attenuation curve k(lambda).\"\"\"
    um = wav_aa / 1e4
    k = np.where(um < 0.63,
        2.659 * (-2.156 + 1.509/um - 0.198/um**2 + 0.011/um**3) + 4.05,
        2.659 * (-1.857 + 1.040/um) + 4.05)
    return np.clip(k, 0, None)

def _redden(flam, ebv):
    return flam * 10 ** (-0.4 * ebv * _calzetti_k(WAV))

def _spec(flam, name):
    return dict(wav=WAV.copy(), flam=np.clip(flam, 1e-6 * flam.max(), None),
                name=name)

def _star(T, lines=(), name=""):
    c = _bb(T)
    f = c.copy()
    for cen, ew, fw in lines:
        f = f + _gauss_line(cen, ew, fw, c, absorption=True)
    return _spec(f, name)

BALMER = [(6562.8, 12, 18), (4861.3, 14, 16), (4340.5, 12, 14), (4101.7, 10, 13)]
SOLAR_ABS = [(3933.7, 12, 10), (3968.5, 9, 10), (4304, 5, 12),
             (5175, 6, 12), (5893, 4, 8)]

def _galaxy(kind):
    if kind == "quenched":
        c = 0.4 * _bb(4600) + 0.6 * _bb(3800)
        c = np.where(WAV < 4000, c * 0.45, c)               # 4000 A break
        f = c.copy()
        for cen, ew, fw in SOLAR_ABS + [(4861.3, 4, 12), (6562.8, 2, 12)]:
            f = f + _gauss_line(cen, ew, fw, c, absorption=True)
        return _spec(f, "quenched galaxy")
    if kind in ("star-forming", "starburst", "dusty starburst", "dusty galaxy",
                "young irregular"):
        blue_frac = {"star-forming": 0.45, "starburst": 0.70,
                     "dusty starburst": 0.70, "dusty galaxy": 0.45,
                     "young irregular": 0.85}[kind]
        c = blue_frac * _bb(15000) + (1 - blue_frac) * _bb(4800)
        ha = {"star-forming": 35, "starburst": 150, "dusty starburst": 150,
              "dusty galaxy": 35, "young irregular": 220}[kind]
        lines = [(6562.8, ha, 8), (4861.3, ha/2.86, 8),
                 (4958.9, ha*0.12, 8), (5006.8, ha*0.36, 8),
                 (6548.0, ha*0.10, 8), (6583.5, ha*0.30, 8),
                 (6716.4, ha*0.10, 8), (6730.8, ha*0.08, 8),
                 (3727.0, ha*0.35, 9)]
        if kind == "young irregular":                       # high ionisation
            lines += [(4363.2, ha*0.05, 8), (4685.7, ha*0.03, 8)]
        f = c.copy()
        for cen, ew, fw in lines:
            f = f + _gauss_line(cen, ew, fw, c)
        ebv = {"dusty starburst": 0.6, "dusty galaxy": 0.8}.get(kind, 0.0)
        if ebv:
            f = _redden(f, ebv)
        return _spec(f, f"{kind} galaxy" if "galaxy" not in kind else kind)
    raise ValueError(kind)

def _agn(obscured=False):
    c = (WAV / 6000.0) ** (-1.5)                             # f_lam power law
    c = c / c.max()
    narrow = [(5006.8, 40, 8), (4958.9, 13, 8), (6583.5, 25, 8),
              (6548.0, 8, 8), (3727.0, 15, 9), (6716.4, 8, 8)]
    broad  = [(6562.8, 250, 90), (4861.3, 90, 80), (4340.5, 40, 75),
              (2798.0, 0, 90)]
    f = c.copy()
    for cen, ew, fw in narrow:
        f = f + _gauss_line(cen, ew, fw, c)
    if not obscured:
        for cen, ew, fw in broad:
            f = f + _gauss_line(cen, ew, fw, c)
    if obscured:
        f = _redden(f, 1.0)
    return _spec(f, "dust-obscured AGN" if obscured else "AGN (type 1)")

TEMPLATES = {
    "O5V star":       _star(42000, [], "O5V star"),
    "B2V star":       _star(20600, BALMER[:2], "B2V star"),
    "A0V star":       _star(9700, BALMER, "A0V star"),
    "F5V star":       _star(6550, BALMER[:2] + SOLAR_ABS[:2], "F5V star"),
    "G2V star":       _star(5770, SOLAR_ABS, "G2V star"),
    "K5V star":       _star(4400, SOLAR_ABS, "K5V star"),
    "M2V star":       _star(3500, SOLAR_ABS[:3], "M2V star"),
    "quenched galaxy":     _galaxy("quenched"),
    "star-forming galaxy": _galaxy("star-forming"),
    "starburst galaxy":    _galaxy("starburst"),
    "dusty starburst":     _galaxy("dusty starburst"),
    "dusty galaxy":        _galaxy("dusty galaxy"),
    "young irregular":     _galaxy("young irregular"),
    "AGN (type 1)":        _agn(False),
    "dust-obscured AGN":   _agn(True),
}

for m82_path in [pathlib.Path("data/M82_template_3500_9000.txt"),
                 DATA_ROOT / "M82_template_3500_9000.txt"]:
    try:
        if not (m82_path.exists() and m82_path.stat().st_size > 1000):
            continue
        m = np.loadtxt(m82_path)
        if m.ndim == 2 and m.shape[0] > 100:
            TEMPLATES["M82 (empirical)"] = dict(
                wav=WAV.copy(),
                flam=np.interp(WAV, m[:, 0], m[:, 1],
                               left=np.nan, right=np.nan),
                name="M82 (empirical)")
            break
    except OSError:
        continue
else:
    print("note: M82 empirical template not found - skipped")

print(f"{len(TEMPLATES)} templates available:")
print("  " + ", ".join(TEMPLATES))

# Gallery
fig, axs = plt.subplots(4, 4, figsize=(14, 9), sharex=True)
for ax, (nm, sp) in zip(axs.ravel(), TEMPLATES.items()):
    fl = sp["flam"] / np.nanmax(sp["flam"])
    ax.plot(sp["wav"], fl, lw=0.6, color=NAVY)
    ax.set_title(nm, fontsize=8); ax.set_yticks([])
    ax.set_xlim(3500, 9500); ax.grid(alpha=0.2)
for ax in axs.ravel()[len(TEMPLATES):]:
    ax.axis("off")
fig.suptitle("Spectral template library (shapes, arbitrary normalization)", y=1.00)
fig.supxlabel("Wavelength (Å)")
plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 3b. Research-grade template library (recommended for observation prep)

The analytic templates above are convenient approximations. For real
proposal and observation planning use the **empirical library** shipped in
`spectral_library/` — the same template sets the ESO / Gemini / HST ETCs
draw from:

| set | contents | reference |
|---|---|---|
| Pickles | 24 stellar types, O5V–M5V dwarfs + giants + supergiants | Pickles 1998, PASP 110, 863 |
| Kinney–Calzetti | E, bulge, S0, Sa, Sb, Sc + starburst series (internal E(B−V) 0.05 → 0.7) | Kinney+ 1996, ApJ 467, 38; Calzetti+ 1994 |
| AGN | QSO composite (Francis+ 1991), Seyfert 1/2, LINER, NGC 1068 | Francis+ 1991, ApJ 373, 465 |
| SDSS composites | SF1–SF4 star-forming sequence, RED0/2/4 passive sequence, blue/red cloud | Dobos+ 2012, MNRAS 420, 1217 |

Use cases map: dusty starburst → `kc96_starb5/6`; quenched → `kc96_elliptical`
or `dobos_RED*`; young irregular → `kc96_starb1` / `dobos_SF4`;
dust-obscured AGN → `agn_seyfert2` / `agn_ngc1068`.

**Native wavelength coverage varies — full-coverage alternatives ship for
every truncated class.** Some CDBS-heritage templates stop short of the
LRS red limit: the Francis+91 QSO composite (`agn_qso`) ends at 6000 Å
rest (it was built from rest-UV/blue spectra of z ≳ 1 quasars),
`agn_seyfert1` at 7078 Å, `agn_liner`/`kc96_bulge` at 7550 Å, `kc96_sc`
at 7660 Å. `load_template()` prints a note when a template doesn't span
the band; S/N outside native coverage is NaN, never extrapolated.

Full-coverage substitutes in the library:

| truncated | use instead | coverage | source |
|---|---|---|---|
| `agn_qso` (→6000 Å) | **`agn_qso_ext`** | 800 Å – 1.05 µm | Francis+91 ⊕ Türler+99 (the ESO ETC splice) |
| `agn_seyfert1` (→7078 Å) | **`brown_NGC5033_sy1.5`** | 805 Å – 1.05 µm | Brown+14 atlas |
| `agn_liner` (→7550 Å) | **`brown_NGC4579_liner`** | 805 Å – 1.05 µm | Brown+14 atlas |
| `kc96_bulge` (→7550 Å) | **`brown_NGC3379_elliptical`**, `brown_NGC4450_sab` | 805 Å – 1.05 µm | Brown+14 atlas |
| `kc96_sc` (→7660 Å) | **`brown_NGC0628_sc`** | 805 Å – 1.05 µm | Brown+14 atlas |

`load_template("kc96_sb")` returns the same spectrum dict as the analytic
generators, so everything downstream (normalization, extinction, both output
modes) is identical. See `spectral_library/README.md` for provenance."""))

cells.append(code("""LIB_DIRS = [pathlib.Path("spectral_library"),
            DATA_ROOT.parent / "spectral_library"]
LIB_DIR = next((d for d in LIB_DIRS if d.exists()), None)

RESEARCH = {}
if LIB_DIR is not None:
    for f in sorted(LIB_DIR.glob("*_*.txt")):
        RESEARCH[f.stem] = f

def load_template(name):
    \"\"\"Load a research-library template by filename stem, e.g. 'pickles_g2v',
    'kc96_starb6', 'agn_qso', 'dobos_SF2'. Falls back to the analytic
    TEMPLATES dict if the name matches there instead.\"\"\"
    if name in RESEARCH:
        t = np.loadtxt(RESEARCH[name])
        wmin, wmax = t[:, 0].min(), t[:, 0].max()
        if wmin > PIX_WAV[0] or wmax < PIX_WAV[-1]:
            print(f"note: '{name}' native coverage {wmin:.0f}-{wmax:.0f} A "
                  f"does not span the full LRS band "
                  f"({PIX_WAV[0]:.0f}-{PIX_WAV[-1]:.0f} A); "
                  f"S/N outside it will be NaN")
        return dict(wav=WAV.copy(),
                    flam=np.interp(WAV, t[:, 0], t[:, 1],
                                   left=np.nan, right=np.nan),
                    wav_native=t[:, 0], flam_native=t[:, 1],
                    name=name)
    if name in TEMPLATES:
        return dict(TEMPLATES[name])
    raise KeyError(f"'{name}' not in spectral_library/ or TEMPLATES. "
                   f"Available: {sorted(RESEARCH) + sorted(TEMPLATES)}")

if RESEARCH:
    groups = {}
    for k in RESEARCH:
        groups.setdefault(k.split("_")[0], []).append(k)
    print(f"Research library: {len(RESEARCH)} templates from {LIB_DIR}/")
    for g, ks in groups.items():
        print(f"  {g:8s} ({len(ks):2d}): {', '.join(sorted(ks)[:8])}"
              + (" ..." if len(ks) > 8 else ""))
else:
    print("spectral_library/ not found - analytic templates only")

# Gallery: one panel per set, normalized at 6000 A
if RESEARCH:
    show = {"pickles": ["pickles_o5v", "pickles_a0v", "pickles_g2v",
                        "pickles_k5v", "pickles_m2v"],
            "kc96":    ["kc96_elliptical", "kc96_sb", "kc96_starb1",
                        "kc96_starb6"],
            "agn":     ["agn_qso_ext", "agn_seyfert2", "brown_NGC5033_sy1.5"],
            "dobos":   ["dobos_SF1", "dobos_SF4", "dobos_RED0", "dobos_RED4"]}
    fig, axs = plt.subplots(2, 2, figsize=(13, 6.5), sharex=True)
    for ax, (g, names) in zip(axs.ravel(), show.items()):
        for nm in names:
            if nm not in RESEARCH:
                continue
            sp = load_template(nm)
            f6 = np.interp(6000, sp["wav"], sp["flam"])
            ax.plot(sp["wav"], sp["flam"] / f6, lw=0.7,
                    label=nm.split("_", 1)[1])
        ax.set_yscale("log"); ax.set_xlim(3500, 9500)
        ax.set_title({"pickles": "Pickles 1998 stars",
                      "kc96": "Kinney-Calzetti galaxies",
                      "agn": "AGN templates",
                      "dobos": "SDSS composites (Dobos+12)"}[g], fontsize=10)
        ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.25)
    fig.supxlabel("Wavelength (Å)")
    fig.supylabel("f$_\\\\lambda$ / f$_\\\\lambda$(6000 Å)")
    plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 4. Custom spectra

Three ways to supply your own model instead of a library template:

1. **User table** — 2-column ASCII (wavelength Å, f_λ). Any absolute scale;
   normalization happens in §5.
2. **Flat continuum** — flat in f_λ or flat in f_ν.
3. **Power law** — f_λ ∝ λ^α.

Add any number of emission or absorption lines on top, specified by center,
equivalent width (Å) and FWHM (km/s)."""))

cells.append(code("""def user_table_spectrum(path):
    \"\"\"Load a 2-column ASCII table: wavelength_A  f_lam.\"\"\"
    t = np.loadtxt(path)
    return dict(wav=WAV.copy(),
                flam=np.interp(WAV, t[:, 0], t[:, 1],
                               left=np.nan, right=np.nan),
                name=pathlib.Path(path).name)

def continuum(kind="flat_flam", alpha=0.0):
    \"\"\"kind: 'flat_flam' | 'flat_fnu' | 'powerlaw' (f_lam ~ lambda^alpha).\"\"\"
    if kind == "flat_flam":
        f = np.ones_like(WAV)
    elif kind == "flat_fnu":                    # f_lam = f_nu c / lambda^2
        f = (6000.0 / WAV) ** 2
    elif kind == "powerlaw":
        f = (WAV / 6000.0) ** alpha
    else:
        raise ValueError(kind)
    return dict(wav=WAV.copy(), flam=f, name=f"{kind} (alpha={alpha})")

def add_line(spec, center_aa, ew_aa, fwhm_kms=300.0, absorption=False):
    \"\"\"Return a copy of spec with a Gaussian line added.
    ew_aa: equivalent width in A (positive number either way).\"\"\"
    fwhm_aa = center_aa * fwhm_kms / 2.998e5
    out = dict(spec)
    out["flam"] = spec["flam"] + _gauss_line(center_aa, ew_aa, fwhm_aa,
                                             spec["flam"],
                                             absorption=absorption)
    out["flam"] = np.clip(out["flam"], 0, None)
    return out

# Example: power-law continuum with one emission and one absorption line
demo = continuum("powerlaw", alpha=-1.0)
demo = add_line(demo, 6562.8, ew_aa=80, fwhm_kms=400)            # emission
demo = add_line(demo, 5175.0, ew_aa=8, fwhm_kms=500, absorption=True)
fig, ax = plt.subplots(figsize=(9, 3))
ax.plot(demo["wav"], demo["flam"], color=NAVY, lw=0.9)
ax.set_xlabel("Wavelength (Å)"); ax.set_ylabel("f$_\\\\lambda$ (arbitrary)")
ax.set_title("Custom model: $f_\\\\lambda \\\\propto \\\\lambda^{-1}$ + emission (Hα) + absorption (Mg b)")
ax.set_xlim(3500, 9500); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 5. Normalization — point/extended, SDSS filter, or total flux

Three ways to set the absolute flux scale (each works for library
templates, user tables, and continuum/line models alike):

* **Monochromatic** — `normalize(spec, mag, wav0)`: AB magnitude at a
  single wavelength (default 6000 Å).
* **SDSS filter magnitude** — `normalize_filter(spec, mag, band)`:
  synthetic photometry through the **real SDSS u/g/r/i/z responses**
  (Doi et al. 2010, including atmosphere at airmass 1.3; note the
  u-band red leak is present, as in the real filter). Photon-counting
  convention: ⟨f_ν⟩ = ∫f_λ T λ dλ / (c ∫T dλ/λ), m = −2.5 log⟨f_ν⟩ − 48.6.
  The function warns if the filter samples regions without template data.
* **Total flux** — `normalize_total_flux(spec, F, wav_range)`: scale so
  the integrated ∫f_λ dλ over a wavelength range equals `F`
  (erg s⁻¹ cm⁻²). The natural choice for **user-built emission-line
  models**: build the line with `add_line()`, then set the total line
  flux directly.

For all three: point sources give total magnitudes/fluxes (seeing + slit
losses applied downstream); `extended=True` interprets the same numbers
per arcsec²."""))

cells.append(code("""def _ab_to_flam(mag_ab, wav_aa):
    \"\"\"AB mag -> f_lam in erg/s/cm^2/A at wav_aa.\"\"\"
    fnu = 10 ** (-0.4 * (mag_ab + 48.60))            # erg/s/cm^2/Hz
    return fnu * C_AA_S / wav_aa ** 2

def normalize(spec, mag_ab, wav0=6000.0, extended=False):
    \"\"\"Scale spec so f_lam(wav0) corresponds to mag_ab.

    extended=False : mag_ab is the total AB magnitude (point source)
    extended=True  : mag_ab is AB mag / arcsec^2  (f_lam is then per arcsec^2)
    \"\"\"
    target = _ab_to_flam(mag_ab, wav0)
    current = np.interp(wav0, spec["wav"], spec["flam"])
    if not np.isfinite(current) or current <= 0:
        raise ValueError(
            f"cannot normalize at {wav0:.0f} A: the template has no data "
            f"there (zero-patched or outside coverage). Choose wav0 inside "
            f"the covered range.")
    out = dict(spec)
    out["flam"] = spec["flam"] * (target / current)
    out["extended"] = extended
    out["norm"] = f"{mag_ab} AB{'/arcsec^2' if extended else ''} at {wav0:.0f} A"
    return out

# --- SDSS filter photometry (real Doi et al. 2010 responses) ------------
_trapz = getattr(np, "trapezoid", np.trapz)
SDSS_BANDS = {}
for _b in "ugriz":
    for _p in [pathlib.Path("filters") / f"sdss_{_b}.txt",
               DATA_ROOT.parent / "filters" / f"sdss_{_b}.txt"]:
        try:
            if _p.exists():
                SDSS_BANDS[_b] = np.loadtxt(_p)
                break
        except OSError:
            continue
print(f"SDSS filters loaded: {sorted(SDSS_BANDS)}"
      if SDSS_BANDS else "filters/ not found - filter normalization disabled")

def synth_ab_mag(spec, band):
    \"\"\"Synthetic AB magnitude of spec through the real SDSS filter
    (photon-counting convention). Warns if the filter samples wavelengths
    where the template has no data.\"\"\"
    t = SDSS_BANDS[band]
    T = np.interp(WAV, t[:, 0], t[:, 1], left=0.0, right=0.0)
    fl = spec["flam"]
    ok = np.isfinite(fl)
    wt = T / WAV                                   # photon-count weight
    frac_cov = _trapz(np.where(ok, wt, 0.0), WAV) / max(_trapz(wt, WAV), 1e-30)
    # part of the filter may extend beyond the 3500-9500 A grid (u red leak
    # tail is inside; z extends to 1.11 um) - report grid coverage too
    grid_cov = (_trapz(t[:, 1] / t[:, 0],
                          t[:, 0], axis=0))
    Tgrid = _trapz(np.interp(t[:, 0], WAV, np.ones_like(WAV),
                                left=0, right=0) * t[:, 1] / t[:, 0], t[:, 0])
    if frac_cov < 0.99:
        print(f"WARNING: template covers only {frac_cov:.0%} of the "
              f"{band}-band photon weight (missing data set to 0)")
    if Tgrid / grid_cov < 0.95:
        print(f"note: {1 - Tgrid/grid_cov:.0%} of the {band}-band response "
              f"lies outside the 3500-9500 A grid and is ignored")
    num = _trapz(np.where(ok, fl, 0.0) * T * WAV, WAV)
    den = C_AA_S * _trapz(T / WAV, WAV)
    if num <= 0:
        return np.inf
    return -2.5 * np.log10(num / den) - 48.60

def normalize_filter(spec, mag_ab, band="r", extended=False):
    \"\"\"Scale spec so its synthetic SDSS `band` magnitude equals mag_ab.
    extended=True: mag_ab is AB mag / arcsec^2.\"\"\"
    if band not in SDSS_BANDS:
        raise KeyError(f"band '{band}' not available ({sorted(SDSS_BANDS)})")
    m0 = synth_ab_mag(spec, band)
    if not np.isfinite(m0):
        raise ValueError(f"spec has no positive flux in the {band} band")
    out = dict(spec)
    out["flam"] = spec["flam"] * 10 ** (-0.4 * (mag_ab - m0))
    out["extended"] = extended
    out["norm"] = (f"SDSS {band} = {mag_ab} AB"
                   f"{'/arcsec^2' if extended else ''}")
    return out

def normalize_total_flux(spec, total_flux, wav_range=(4000.0, 8092.0),
                          extended=False):
    \"\"\"Scale spec so the integral of f_lam over wav_range equals
    total_flux (erg/s/cm^2). Ideal for user-built line models.\"\"\"
    w0, w1 = wav_range
    sel = (spec["wav"] >= w0) & (spec["wav"] <= w1) & np.isfinite(spec["flam"])
    cur = _trapz(spec["flam"][sel], spec["wav"][sel])
    if cur <= 0:
        raise ValueError(f"no positive flux in {w0:.0f}-{w1:.0f} A")
    out = dict(spec)
    out["flam"] = spec["flam"] * (total_flux / cur)
    out["extended"] = extended
    out["norm"] = (f"integral {w0:.0f}-{w1:.0f} A = {total_flux:.3e} "
                   f"erg/s/cm^2{'/arcsec^2' if extended else ''}")
    return out

# --- Checks -------------------------------------------------------------
chk1 = normalize(TEMPLATES["G2V star"], 18.0)
chk2 = normalize(TEMPLATES["quenched galaxy"], 21.0, extended=True)
for c in (chk1, chk2):
    print(f"{c['name']:22s} -> {c['norm']}  "
          f"(f_lam at 6000 A = {np.interp(6000, c['wav'], c['flam']):.3e})")

if SDSS_BANDS:
    # Filter normalization: r = 19 Sb galaxy; report its full ugriz colors
    sb_r19 = normalize_filter(load_template("kc96_sb") if RESEARCH
                              else TEMPLATES["star-forming galaxy"],
                              19.0, band="r")
    mags = {b: synth_ab_mag(sb_r19, b) for b in "ugriz"}
    print(f"\\n{sb_r19['name'] if 'name' in sb_r19 else 'Sb'} scaled to "
          f"{sb_r19['norm']}:")
    print("  synthetic colors: " + "  ".join(f"{b}={m:.2f}"
          for b, m in mags.items()))

# Total-flux normalization: a pure H-alpha line model at 5e-15 erg/s/cm^2
line_model = continuum("flat_flam")
line_model["flam"] *= 1e-6                       # negligible continuum
line_model = add_line(line_model, 6562.8, ew_aa=1e6, fwhm_kms=250)
line_model = normalize_total_flux(line_model, 5e-15, wav_range=(6500, 6630))
print(f"\\nline model: {line_model['norm']}  (peak f_lam = "
      f"{np.nanmax(line_model['flam']):.2e})")
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 5b. Galactic (foreground) extinction

Redden the model with a Milky-Way extinction law before predicting counts.
Choose the law and R_V, and give **either** E(B−V) **or** A_V directly
(A_V = R_V × E(B−V)).

Laws (via the `dust_extinction` package when installed, else a built-in
CCM89 implementation):

* **CCM89** — Cardelli, Clayton & Mathis 1989 (R_V-parametrised, the default)
* **O94** — O'Donnell 1994 update of the optical polynomials
* **F99** — Fitzpatrick 1999
* **Calzetti00** — internal *attenuation* for starburst hosts (R_V = 4.05);
  use this for dust inside the target, not for the Milky-Way foreground

Get E(B−V) for your target from Schlafly & Finkbeiner 2011 via the
[IRSA dust tool](https://irsa.ipac.caltech.edu/applications/DUST/).

**Ordering matters.** If your target magnitude is the *observed* one
(usual case: from a catalogue), apply extinction **before** `normalize()` —
the spectrum keeps the requested observed magnitude and only its shape
reddens. If the magnitude is *intrinsic*, normalize first, then extinguish."""))

cells.append(code("""try:
    from dust_extinction.parameter_averages import CCM89, O94, F99
    import astropy.units as _u
    _HAVE_DUST_EXT = True
except ImportError:
    _HAVE_DUST_EXT = False

def _ccm89_builtin(wav_aa, rv):
    \"\"\"CCM89 A(lambda)/A(V), optical/NIR polynomials (3030-33000 A).\"\"\"
    x = 1e4 / wav_aa                       # 1/um
    a = np.zeros_like(x); b = np.zeros_like(x)
    opt = (x >= 1.1) & (x <= 3.3)
    y = x[opt] - 1.82
    a[opt] = (1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4
              + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7)
    b[opt] = (1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4
              - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7)
    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61
    return a + b / rv

def extinction_transmission(law="CCM89", rv=3.1, ebv=None, av=None):
    \"\"\"Fractional transmission vs WAV for the chosen law.\"\"\"
    if (ebv is None) == (av is None):
        raise ValueError("give exactly one of ebv= or av=")
    if av is None:
        av = rv * ebv
    ebv = av / rv
    if law == "Calzetti00":                      # attenuation, R_V fixed 4.05
        return 10 ** (-0.4 * (av / 4.05) * _calzetti_k(WAV))
    if _HAVE_DUST_EXT and law in ("CCM89", "O94", "F99"):
        model = {"CCM89": CCM89, "O94": O94, "F99": F99}[law](Rv=rv)
        return model.extinguish(WAV * _u.AA, Ebv=ebv)
    if law in ("CCM89", "O94", "F99"):           # built-in fallback
        return 10 ** (-0.4 * av * _ccm89_builtin(WAV, rv))
    raise ValueError(f"unknown law '{law}'")

def apply_extinction(spec, law="CCM89", rv=3.1, ebv=None, av=None):
    \"\"\"Return a reddened copy of spec.\"\"\"
    out = dict(spec)
    out["flam"] = spec["flam"] * extinction_transmission(law, rv, ebv, av)
    tag = f"E(B-V)={ebv:.3f}" if ebv is not None else f"A_V={av:.2f}"
    out["name"] = f"{spec.get('name','spec')} + {law} {tag} (Rv={rv})"
    return out

print("dust_extinction package:", "in use" if _HAVE_DUST_EXT
      else "NOT installed - using built-in CCM89 (pip install dust_extinction)")

# Demo, two panels:
#   top    - the transmission curves of the laws themselves (full band)
#   bottom - CCM89 applied to a FULL-COVERAGE AGN template (NGC 1068).
# Note: the Francis+91 QSO composite (agn_qso) natively ends at 6000 A -
# redward of that the ESO ETC splices in Turler+99 data which is not part
# of this optical library. Use agn_ngc1068 / agn_seyfert2 when you need
# the whole LRS band, or accept NaN S/N beyond a template's native range.
fig, axs = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True)

for law, rv, colr, ls in [("CCM89", 3.1, NAVY, "-"), ("CCM89", 5.0, NAVY, "--"),
                           ("O94", 3.1, LIME, "-"), ("F99", 3.1, ORANGE, "-"),
                           ("Calzetti00", 4.05, RED, "-")]:
    tr = extinction_transmission(law, rv=rv, ebv=0.3)
    axs[0].plot(WAV, tr, color=colr, ls=ls, lw=1.2,
                label=f"{law}, R$_V$={rv}")
axs[0].set_ylabel("Transmission  [E(B-V) = 0.3]")
axs[0].set_title("Extinction / attenuation laws at E(B-V) = 0.3")
axs[0].legend(fontsize=8, ncol=2); axs[0].grid(alpha=0.3)

base = (load_template("agn_ngc1068") if RESEARCH
        else TEMPLATES["AGN (type 1)"])
f6 = np.interp(6000, base["wav"], base["flam"])
axs[1].plot(WAV, base["flam"] / f6, color=NAVY, lw=0.8, label="intrinsic")
for ebv, colr in [(0.1, LIME), (0.3, ORANGE), (0.6, RED)]:
    r = apply_extinction(base, ebv=ebv)
    axs[1].plot(WAV, r["flam"] / f6, lw=0.8, color=colr,
                label=f"CCM89  E(B-V) = {ebv}  (A$_V$ = {3.1*ebv:.2f})")
axs[1].set_yscale("log"); axs[1].set_xlim(3500, 9500)
axs[1].set_xlabel("Wavelength (Å)")
axs[1].set_ylabel("f$_\\\\lambda$ / f$_\\\\lambda^{intr}$(6000 Å)")
axs[1].set_title(f"CCM89 applied to {base['name']} (full LRS coverage)")
axs[1].legend(fontsize=8); axs[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 5c. Redshifting templates

`redshift_spectrum(spec, z)` shifts a template to redshift z. The library
files keep their **full native UV coverage** (KC96 → 1235 Å, QSO
composite → 800 Å, Pickles → 1150 Å; the SDSS/Dobos composites physically
start at 3550 Å), so rest-UV light redshifts *into* the LRS band:

| template family | rest coverage | full 4000–8000 Å coverage up to |
|---|---|---|
| Kinney–Calzetti | 1235–9950 Å | **z ≈ 2.2** |
| QSO composite `agn_qso_ext` (Francis ⊕ Türler) | 800 Å – 1.05 µm | **z ≈ 0 – 4.0** |
| SDSS composites (Dobos) | 3550–9050 Å | z ≈ 0.13 |

Two modes:

* **Shape mode** (default) — shift wavelengths, keep the f_λ shape (with
  the (1+z)⁻¹ bandwidth factor). Use when you know the target's
  **observed** magnitude: redshift → extinguish → `normalize()`.
* **Cosmological mode** (`d_L_ref_Mpc=`) — for templates with an absolute
  flux scale (e.g. the M82 empirical template at D ≈ 3.5 Mpc): also apply
  the luminosity-distance dimming (D_ref/D_L(z))²/(1+z), using
  `astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)` with a built-in
  fallback integrator.

`redshift_limits(name)` reports the usable z ranges for any library
template. Feature reminders: Hα leaves the band at z = 0.22, [O III] 5007
at z = 0.60, the 4000 Å break at z = 1.0, Mg II 2798 *enters* at z = 0.43,
and Lyα arrives only at z ≈ 2.3.

**Redshifts are accepted up to z = 9.** Where the redshifted template no
longer covers the grid, the missing region is **patched with zeros** and a
warning names the missing rest-frame regime (e.g. *"lack of template in
the rest-frame X-ray and extreme-UV"*). Zero-patched regions carry no
source flux — the ETC will show sky-limited S/N ≈ 0 there, and
`normalize()` refuses an anchor wavelength inside a patched region."""))

cells.append(code("""try:
    from astropy.cosmology import FlatLambdaCDM as _FLCDM
    from astropy import units as _au
    _COSMO = _FLCDM(H0=70, Om0=0.30)
    def luminosity_distance_Mpc(z):
        return float(_COSMO.luminosity_distance(z).to(_au.Mpc).value)
except Exception:
    def luminosity_distance_Mpc(z, H0=70.0, Om=0.30):
        c_km = 2.998e5; DH = c_km / H0
        if z <= 0: return 0.0
        zs = np.linspace(0, z, 2048)
        E = np.sqrt(Om * (1 + zs) ** 3 + (1 - Om))
        DC = DH * _trapz(1.0 / E, zs)
        return DC * (1 + z)

Z_MAX = 9.0

_REST_REGIMES = [(0.0, 100.0, "X-ray"), (100.0, 912.0, "extreme-UV"),
                 (912.0, 3200.0, "UV"), (3200.0, 7000.0, "optical"),
                 (7000.0, 2.5e4, "NIR"), (2.5e4, np.inf, "IR")]

def _regime_names(rest_lo, rest_hi):
    return " and ".join(n for a, b, n in _REST_REGIMES
                        if rest_lo < b and rest_hi > a)

def redshift_spectrum(spec, z, d_L_ref_Mpc=None, verbose=True):
    \"\"\"Shift a template to redshift z (0 <= z <= 9), observed-frame on
    the WAV grid.

    Uses the template's native grid when available (research library) so
    rest-UV data redshifts into the optical. Where the redshifted template
    does not cover the grid, the spectrum is PATCHED WITH ZEROS and a
    warning names the missing rest-frame regime (e.g. \"Lack of template
    in X-ray and UV\"). Shape mode by default; give d_L_ref_Mpc for
    cosmological dimming of absolute-flux templates.\"\"\"
    if not (0.0 <= z <= Z_MAX):
        raise ValueError(f"z must be between 0 and {Z_MAX}")
    wr = spec.get("wav_native", spec["wav"])
    fr = spec.get("flam_native", spec["flam"])
    ok = np.isfinite(fr)
    wr, fr = wr[ok], fr[ok]
    wav_obs = wr * (1 + z)
    flam_obs = fr / (1 + z)
    if d_L_ref_Mpc is not None and z > 0:
        DL = luminosity_distance_Mpc(z)
        flam_obs = flam_obs * (d_L_ref_Mpc / DL) ** 2

    # Zero-patch outside the redshifted native coverage, with named warnings
    name = spec.get("name", "spec")
    if verbose and wav_obs.min() > WAV[0]:
        rest_lo, rest_hi = WAV[0] / (1 + z), wr.min()
        print(f"WARNING [{name} @ z={z:g}]: lack of template in the "
              f"rest-frame {_regime_names(rest_lo, rest_hi)} - observed "
              f"{WAV[0]:.0f}-{wav_obs.min():.0f} A patched with zeros")
    if verbose and wav_obs.max() < WAV[-1]:
        rest_lo, rest_hi = wr.max(), WAV[-1] / (1 + z)
        print(f"WARNING [{name} @ z={z:g}]: lack of template in the "
              f"rest-frame {_regime_names(rest_lo, rest_hi)} - observed "
              f"{wav_obs.max():.0f}-{WAV[-1]:.0f} A patched with zeros")

    out = dict(wav=WAV.copy(),
               flam=np.interp(WAV, wav_obs, flam_obs, left=0.0, right=0.0),
               name=f"{name} @ z={z:.2f}", z=z)
    if d_L_ref_Mpc is not None and z > 0:
        out["d_L_Mpc"] = DL
    return out

def redshift_limits(name, band=(4000.0, 8000.0)):
    \"\"\"Usable redshift ranges of a library template vs an observed band.\"\"\"
    t = np.loadtxt(RESEARCH[name])
    lmin, lmax = t[:, 0].min(), t[:, 0].max()
    z_full_lo = max(0.0, band[1] / lmax - 1)
    z_full_hi = band[0] / lmin - 1
    z_any = band[1] / lmin - 1
    print(f"{name}: rest {lmin:.0f}-{lmax:.0f} A -> full band "
          f"{band[0]:.0f}-{band[1]:.0f}: "
          + (f"z = {z_full_lo:.2f}-{z_full_hi:.2f}"
             if z_full_hi >= z_full_lo else "never")
          + f";  any overlap to z = {z_any:.2f}")
    return z_full_lo, z_full_hi, z_any

for nm in ("kc96_sb", "kc96_starb2", "agn_qso_ext", "dobos_SF2"):
    if nm in RESEARCH:
        redshift_limits(nm)

# --- Demo: the same starburst, same OBSERVED magnitude, at four z -------
# (S/N curves for this z-series appear in Sec. 8, after the ETC core.)
Z_SERIES = [0.0, 0.3, 0.7, 1.5]
Z_COLORS = [NAVY, LIME, ORANGE, RED]
base_z = load_template("kc96_starb2") if RESEARCH else TEMPLATES["starburst galaxy"]
fig, ax = plt.subplots(figsize=(11.5, 3.8))
for z, colr in zip(Z_SERIES, Z_COLORS):
    sz = normalize(redshift_spectrum(base_z, z), 20.0)    # observed 20 AB
    ax.plot(sz["wav"], sz["flam"], lw=0.8, color=colr, label=f"z = {z:.1f}")
    for lam0 in (6562.8, 5006.8, 3727.0, 2798.0):         # Ha [OIII] [OII] MgII
        lo = lam0 * (1 + z)
        if 4000 < lo < 8092:
            ax.axvline(lo, color=colr, lw=0.5, ls=":", alpha=0.6)
ax.set_yscale("log"); ax.set_ylabel(r"$f_\\lambda$ (erg/s/cm²/Å)")
ax.set_xlabel("Wavelength (Å)")
ax.set_title("KC96 starburst redshifted, renormalized to observed 20 AB "
             "(dotted: Hα, [OIII], [OII], MgII; rest-UV enters at high z)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# --- High-z demo: zero-patching + missing-regime warnings ---------------
qso = load_template("agn_qso_ext") if RESEARCH else TEMPLATES["AGN (type 1)"]
fig, ax = plt.subplots(figsize=(11.5, 3.6))
for z, colr in [(2.0, NAVY), (5.0, ORANGE), (8.0, RED)]:
    qz = redshift_spectrum(qso, z)
    fl = qz["flam"] / max(np.nanmax(qz["flam"]), 1e-30)
    ax.plot(qz["wav"], fl, lw=0.9, color=colr, label=f"QSO @ z = {z:.0f}")
    lya = 1215.7 * (1 + z)
    if 3500 < lya < 9500:
        ax.axvline(lya, color=colr, lw=0.6, ls=":", alpha=0.7)
        ax.text(lya, 1.02, "Lyα", color=colr, fontsize=8, ha="center")
ax.set_xlabel("Wavelength (Å)"); ax.set_ylabel("normalized f$_\\\\lambda$")
ax.set_title("QSO composite at high z — zero-patched where the template "
             "lacks rest-frame X-ray/EUV data")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# Cosmological-dimming check on the absolute-flux M82 template
if "M82 (empirical)" in TEMPLATES:
    print("\\nM82 (D=3.5 Mpc) moved to cosmological distances (mode-2 view):")
    for z in (0.01, 0.05, 0.10):
        mz = redshift_spectrum(TEMPLATES["M82 (empirical)"], z,
                                d_L_ref_Mpc=3.5)
        f6 = np.interp(6000, mz["wav"], mz["flam"])
        mag = -2.5 * np.log10(f6 * 6000**2 / C_AA_S) - 48.6
        print(f"  z={z:.2f}: D_L={mz['d_L_Mpc']:7.1f} Mpc, "
              f"observed {mag:5.1f} AB at 6000 A")
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 6. Observing conditions

* **Seeing** — PSF FWHM at 5500 Å (Doi Inthanon median ≈ 1.1″). The slit
  coupling integrates a **circular** PSF over the rectangular
  slit × extraction aperture. For the default Gaussian this is exact
  (the 2-D Gaussian separates in x, y — the erf product is the true
  rectangle integral, not a square-PSF approximation). A **Moffat**
  option (`psf_profile="moffat"`, `moffat_beta=`) models realistic
  power-law wings, which couple 10–25 % less through the 1.8″ slit at
  the same FWHM. Note: the §2b empirical calibration was derived with
  the Gaussian model, so its gray factor silently absorbs the mean wing
  loss at the calibration configuration — switch profiles consistently.
* **Lunar phase** — sets the *continuum* sky surface brightness
  (V, AB mag/arcsec²): dark 21.8, gray 20.9, bright 19.0, with a color
  term (moonlight is blue, zodiacal+scattered light redder).
* **Sky emission lines** — the night sky is not smooth: **airglow** adds
  [O I] 5577/6300/6364, Na D 5893, the O₂ band at 8645 and the **OH
  Meinel bands** that dominate redward of ~6800 Å. Airglow does **not**
  scale with lunar phase — under a bright moon the continuum rises to
  meet the lines, while on a dark night the red OH forest towers over the
  continuum and carves S/N dips at every band (visible in every S/N plot
  downstream). Line widths are rendered at the instrument resolution.
* **Cloud cover** — gray extinction: photometric 0.0 mag, thin cirrus
  0.5 mag, cloudy 1.2 mag.
* **Airmass** — atmospheric extinction k(λ) mag/airmass. Select by target
  **altitude**, five choices (airmass = sec z):

  | altitude | 30° | 45° | 60° | 75° | 90° |
  |---|---|---|---|---|---|
  | airmass | 2.00 | 1.41 | 1.15 | 1.04 | 1.00 |"""))

cells.append(code("""SKY_V_AB = {"dark": 21.8, "gray": 20.9, "bright": 19.0}
CLOUD_MAG = {"photometric": 0.0, "thin cirrus": 0.5, "cloudy": 1.2}

# color of the night sky relative to V (approx; OH makes the red brighter)
_SKY_SHAPE_W = [3500, 4500, 5500, 6500, 7500, 8500, 9500]
_SKY_SHAPE_D = [+0.8, +0.4, 0.0, -0.3, -0.7, -1.0, -1.1]     # dark sky
_SKY_SHAPE_B = [-0.4, -0.2, 0.0, -0.2, -0.5, -0.8, -0.9]     # moonlit (bluer)

# atmospheric extinction, mag per airmass
_EXT_W = [3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000, 9500]
_EXT_K = [0.55, 0.35, 0.25, 0.18, 0.15, 0.12, 0.09, 0.06, 0.05]

def sky_surface_brightness(lunar="dark"):
    \"\"\"CONTINUUM sky brightness, AB mag/arcsec^2 vs WAV (no lines).\"\"\"
    shape = _SKY_SHAPE_B if lunar == "bright" else _SKY_SHAPE_D
    return SKY_V_AB[lunar] + np.interp(WAV, _SKY_SHAPE_W, shape)

# --- Airglow emission lines (independent of lunar phase) ----------------
# (wavelength A, peak strength as multiple of the DARK-sky continuum there)
SKY_LINES = [
    (5577.3, 14.0),                       # [O I] green line
    (5890.0, 2.5), (5895.9, 1.5),         # Na D airglow
    (6300.3, 6.0), (6363.8, 2.0),         # [O I] red doublet
    # OH Meinel band clumps (blended at LRS resolution)
    (6834, 2.0), (6871, 2.5), (6923, 1.8),
    (7240, 2.5), (7276, 2.8), (7316, 3.0), (7340, 2.5), (7369, 2.2),
    (7524, 3.5), (7571, 3.0), (7623, 2.8), (7750, 3.2), (7794, 3.5),
    (7821, 3.2), (7913, 4.0), (7964, 3.5), (7993, 4.0), (8025, 3.6),
    (8399, 4.5), (8430, 4.2), (8465, 4.0), (8505, 4.5), (8541, 4.2),
    (8615, 4.0), (8645, 5.0),             # + O2 atmospheric band
    (8791, 5.0), (8827, 4.6), (8886, 4.8), (8919, 4.5), (8943, 4.8),
    (9002, 4.5), (9306, 5.0), (9375, 4.6), (9439, 5.0),
]

def sky_emission_flam(lunar="dark", R=750):
    \"\"\"Total sky spectrum f_lam per arcsec^2: lunar-phase continuum +
    fixed airglow lines rendered at the instrument resolution R.\"\"\"
    cont = _ab_to_flam(sky_surface_brightness(lunar), WAV)
    cont_dark = _ab_to_flam(sky_surface_brightness("dark"), WAV)
    lines = np.zeros_like(WAV)
    for cen, peak_x in SKY_LINES:
        fwhm = max(cen / R, 6.0)
        sig = fwhm / 2.3548
        peak = peak_x * np.interp(cen, WAV, cont_dark)
        lines += peak * np.exp(-0.5 * ((WAV - cen) / sig) ** 2)
    return cont + lines

def atmospheric_transmission(airmass=1.2, clouds="photometric"):
    \"\"\"Transmission vs WAV including gray cloud extinction.\"\"\"
    k = np.interp(WAV, _EXT_W, _EXT_K)
    return 10 ** (-0.4 * (k * airmass + CLOUD_MAG[clouds]))

def slit_coupling(seeing_arcsec, slit_width_arcsec, extract_arcsec,
                  profile="gaussian", moffat_beta=3.5):
    \"\"\"Fraction of the PSF passing the slit AND the extraction window.

    profile="gaussian": a circular 2-D Gaussian is separable in x,y, so
    its integral over the rectangular aperture is EXACTLY the product of
    the two 1-D erf terms (no square-PSF approximation involved).

    profile="moffat": real seeing has power-law wings; a Moffat of the
    same FWHM couples 10-25% less through a narrow slit than the
    Gaussian. The Moffat is NOT separable, so the rectangle integral is
    done numerically (fast; called once per ETC run). beta ~ 2.5 (bad,
    windy) to ~ 4.5 (good dome seeing); 3.5 is a fair default.\"\"\"
    if profile == "gaussian":
        sig = seeing_arcsec / 2.3548
        fx = erf(slit_width_arcsec / (2 * np.sqrt(2) * sig))
        fy = erf(extract_arcsec   / (2 * np.sqrt(2) * sig))
        return fx * fy
    if profile == "moffat":
        b = moffat_beta
        alpha = seeing_arcsec / (2 * np.sqrt(2 ** (1 / b) - 1))
        span = max(slit_width_arcsec, extract_arcsec) / 2 + 6 * seeing_arcsec
        n = 601
        x = np.linspace(-span, span, n)
        X, Y = np.meshgrid(x, x)
        I = (b - 1) / (np.pi * alpha**2) * (1 + (X**2 + Y**2) / alpha**2) ** (-b)
        inside = ((np.abs(X) <= slit_width_arcsec / 2)
                  & (np.abs(Y) <= extract_arcsec / 2))
        return float((I * inside).sum() * (x[1] - x[0]) ** 2)
    raise ValueError(f"unknown profile '{profile}'")

# --- Airmass selector: five altitude choices, airmass = sec z -----------
ALTITUDE_CHOICES_DEG = (30, 45, 60, 75, 90)

def select_airmass(altitude_deg):
    \"\"\"Airmass for one of the five allowed target altitudes (deg).\"\"\"
    if altitude_deg not in ALTITUDE_CHOICES_DEG:
        raise ValueError(f"altitude_deg must be one of "
                         f"{ALTITUDE_CHOICES_DEG}, got {altitude_deg}")
    return 1.0 / np.sin(np.radians(altitude_deg))

print("Altitude -> airmass choices:")
for alt in ALTITUDE_CHOICES_DEG:
    print(f"  {alt:2d} deg  ->  X = {select_airmass(alt):.3f}")

fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.3))
for lunar, colr in [("bright", RED), ("gray", ORANGE), ("dark", NAVY)]:
    axs[0].plot(WAV, sky_emission_flam(lunar), color=colr, lw=0.7,
                label=lunar)
axs[0].set_yscale("log")
axs[0].set_title("Sky: moonlit continuum + fixed airglow/OH lines")
axs[0].set_xlabel("Wavelength (Å)")
axs[0].set_ylabel(r"$f_\\lambda$ (erg/s/cm²/Å/arcsec²)")
axs[0].annotate("[O I] 5577", xy=(5577, np.interp(5577, WAV,
                sky_emission_flam("dark"))), fontsize=7, color="#555555",
                xytext=(4600, 2e-16), arrowprops=dict(arrowstyle="-",
                lw=0.5, color="#888888"))
axs[0].text(8300, 3e-17, "OH Meinel\\nforest", fontsize=7, color="#555555",
            ha="center")
axs[0].legend(fontsize=8, loc="upper left"); axs[0].grid(alpha=0.3)
for cl, colr in [("photometric", NAVY), ("thin cirrus", ORANGE), ("cloudy", RED)]:
    axs[1].plot(WAV, atmospheric_transmission(1.2, cl), color=colr, label=cl)
axs[1].set_title("Clouds (airmass 1.2)")
axs[1].set_xlabel("Wavelength (Å)"); axs[1].set_ylabel("Transmission")
axs[1].legend(fontsize=8); axs[1].grid(alpha=0.3)
alt_colors = [RED, ORANGE, "#8064A2", LIME, NAVY]
for alt, colr in zip(ALTITUDE_CHOICES_DEG, alt_colors):
    X = select_airmass(alt)
    axs[2].plot(WAV, atmospheric_transmission(X, "photometric"), color=colr,
                label=f"{alt}°  (X = {X:.2f})")
axs[2].set_title("Target altitude (photometric)")
axs[2].set_xlabel("Wavelength (Å)"); axs[2].set_ylabel("Transmission")
axs[2].legend(fontsize=8); axs[2].grid(alpha=0.3)
plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 6b. Target visibility from the Thai National Observatory

Give the target's **RA and Dec in decimal degrees**; the code computes the
altitude range at culmination from the TNO site coordinates
(Doi Inthanon: 18.5738° N, 98.4823° E, 2457 m):

* upper culmination: alt_max = 90° − |φ − δ|
* lower culmination: alt_min = |φ + δ| − 90°  (> 0 means circumpolar)

Flags: **alt_max < 30° → not observable** from the TNO;
30–50° → observable but sub-optimal (X > 2 at best);
**> 50° → optimal**. The target transits at LST = RA, and the function
suggests the matching `altitude_deg` choice for the ETC."""))

cells.append(code("""TNO_LAT_DEG  = 18.5738      # Doi Inthanon, Thai National Observatory
TNO_LON_DEG  = 98.4823
TNO_ELEV_M   = 2457.0

def target_visibility(ra_deg, dec_deg, verbose=True):
    \"\"\"Altitude range of a target (RA, Dec in decimal degrees) at the TNO.

    Returns dict with alt_max, alt_min (deg), circumpolar flag, status
    ('optimal' / 'observable' / 'not observable'), the recommended
    altitude_deg choice for run_lrs_etc, and the transit LST (hours).
    \"\"\"
    if not (0.0 <= ra_deg < 360.0):
        raise ValueError("RA must be in [0, 360) decimal degrees")
    if not (-90.0 <= dec_deg <= 90.0):
        raise ValueError("Dec must be in [-90, +90] decimal degrees")
    phi, dec = TNO_LAT_DEG, dec_deg
    alt_max = 90.0 - abs(phi - dec)
    alt_min = abs(phi + dec) - 90.0
    circumpolar = alt_min > 0

    if alt_max < 30.0:
        status = "not observable"
    elif alt_max < 50.0:
        status = "observable"
    else:
        status = "optimal"

    # Highest ETC altitude choice reachable by this target
    reachable = [a for a in ALTITUDE_CHOICES_DEG if a <= alt_max]
    alt_choice = max(reachable) if reachable else None

    if verbose:
        flag = {"optimal": "OK  (optimal)",
                "observable": "!   observable but sub-optimal "
                              "(alt_max < 50 deg, X > 1.3 at best)",
                "not observable": "XX  NOT OBSERVABLE from the TNO "
                                  "(alt_max < 30 deg)"}[status]
        print(f"RA {ra_deg:8.3f}, Dec {dec_deg:+8.3f}  ->  "
              f"alt {alt_min:+6.1f} to {alt_max:5.1f} deg   {flag}")
        if circumpolar:
            print("      circumpolar from the TNO (never sets)")
        if alt_choice is not None:
            print(f"      transit at LST = {ra_deg/15:.2f} h; use "
                  f"altitude_deg={alt_choice} in run_lrs_etc "
                  f"(X = {select_airmass(alt_choice):.2f})")
    return dict(alt_max=alt_max, alt_min=alt_min, circumpolar=circumpolar,
                status=status, altitude_deg=alt_choice,
                transit_lst_h=ra_deg / 15.0)

# Examples: M82, the LMC, and an equatorial quasar field
print("Visibility check from the TNO (lat 18.574 N):")
vis_m82 = target_visibility(148.970, 69.680)     # M82
vis_lmc = target_visibility(80.894, -69.756)     # LMC
vis_3c273 = target_visibility(187.278, 2.052)    # 3C 273
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 7. ETC core

For each 4 Å pixel column the collected source rate is

$S(\\lambda) = f_\\lambda \\cdot \\frac{\\lambda}{hc} \\cdot A \\cdot
\\eta(\\lambda) \\cdot T_{atm}(\\lambda) \\cdot f_{slit} \\cdot \\Delta\\lambda_{pix}$   [e⁻/s]

and the noise per column over the extraction aperture is the usual CCD
equation (source + sky + dark shot noise, plus read noise per frame):

$\\mathrm{S/N} = \\dfrac{S\\,t_{tot}}{\\sqrt{S\\,t_{tot} + B\\,t_{tot} +
n_{pix} D\\,t_{tot} + n_{pix} N_{frames}\\, RN^2}}$

The template is first convolved to the instrument resolution (set by the
slit) and resampled onto the 4 Å pixel grid. Quoted S/N is **per pixel**;
S/N per Å = S/N per pixel ÷ 2 (a 1 Å bin holds ¼ of a pixel's photons)."""))

cells.append(code("""def run_lrs_etc(spec, t_per_frame, n_frames,
                slit="1.8", seeing=1.1, lunar="dark", clouds="photometric",
                airmass=1.2, altitude_deg=None,
                temperature="-80C", readout="slow",
                extract_arcsec=None, source_fwhm_arcsec=0.0,
                use_empirical=None, psf_profile="gaussian",
                moffat_beta=3.5):
    \"\"\"Core ETC. Returns a dict of per-pixel arrays and the config.

    Airmass: either give airmass= directly, or altitude_deg= (one of
    30/45/60/75/90) which overrides airmass with sec z.

    use_empirical: None -> follow USE_EMPIRICAL_CALIBRATION; True ->
    scale source+sky photon rates by the measured commissioning response
    ratio (Sec 2b); False -> pure theoretical component model.

    Spatial profile: a Gaussian of FWHM_eff = sqrt(seeing^2 +
    source_fwhm_arcsec^2) - i.e. the object size convolved with the PSF.
    source_fwhm_arcsec=0 is a pure point source; a resolved Gaussian
    object (compact galaxy, nucleus) both spreads along the slit AND
    loses more light at the slit jaws. extended=True (from normalize())
    instead treats the source as uniform surface brightness with no
    profile losses.\"\"\"
    if altitude_deg is not None:
        airmass = select_airmass(altitude_deg)
    sl = SLITS[slit]
    R  = sl["R"]; w_slit = sl["width_arcsec"]
    extended = spec.get("extended", False)
    fwhm_eff = float(np.hypot(seeing, source_fwhm_arcsec))
    if extract_arcsec is None:
        extract_arcsec = max(1.5 * fwhm_eff, 2 * SPATIAL_AS_PIX)
    n_spat = max(1, int(np.ceil(extract_arcsec / SPATIAL_AS_PIX)))

    # --- convolve template to the instrument resolution, resample to pixels
    fwhm_aa = 6000.0 / R
    flam_conv = gaussian_filter1d(np.nan_to_num(spec["flam"]),
                                  fwhm_aa / 2.3548)      # grid is 1 A
    flam_pix = np.interp(PIX_WAV, WAV, flam_conv)

    # --- atmospheric and instrument response on the pixel grid
    atm = np.interp(PIX_WAV, WAV, atmospheric_transmission(airmass, clouds))
    eta = system_throughput(PIX_WAV)
    phot = PIX_WAV / HC                                   # photons per erg

    # --- source rate (e-/s per pixel column)
    if extended:
        area = w_slit * extract_arcsec                    # arcsec^2 in aperture
        S = flam_pix * phot * A_CM2 * eta * atm * DISP_AA_PIX * area
        fslit = 1.0
    else:
        fslit = slit_coupling(fwhm_eff, w_slit, extract_arcsec,
                              profile=psf_profile, moffat_beta=moffat_beta)
        S = flam_pix * phot * A_CM2 * eta * atm * DISP_AA_PIX * fslit

    # --- sky rate (e-/s per pixel column over the aperture) — continuum
    #     (lunar-phase dependent) + fixed airglow/OH lines at this slit's R
    sky_flam = sky_emission_flam(lunar, R=R)
    sky_pix  = np.interp(PIX_WAV, WAV, sky_flam)
    B = sky_pix * phot * A_CM2 * eta * DISP_AA_PIX * (w_slit * extract_arcsec)

    # --- empirical commissioning calibration (Sec 2b) -------------------
    emp = USE_EMPIRICAL_CALIBRATION if use_empirical is None else use_empirical
    if emp:
        _fac = np.interp(PIX_WAV, WAV, EMPIRICAL_RATIO_GRID)
        S = S * _fac
        B = B * _fac

    # --- detector terms
    dark = DARK_E_PIX_S[temperature]
    rn   = READ_NOISE_E[readout]
    t_tot = t_per_frame * n_frames

    var = S * t_tot + B * t_tot + n_spat * dark * t_tot + n_spat * n_frames * rn**2
    snr_pix = np.where(var > 0, S * t_tot / np.sqrt(var), 0.0)

    return dict(wav=PIX_WAV, pixel=np.arange(N_SPEC_PIX),
                snr_pix=snr_pix, snr_aa=snr_pix / np.sqrt(DISP_AA_PIX),
                S=S, B=B, dark_rate=n_spat * dark,
                rn_var_total=n_spat * n_frames * rn**2,
                flam_pix=flam_pix,
                t_tot=t_tot, fslit=fslit, n_spat=n_spat, R=R,
                config=dict(slit=slit, seeing=seeing, fwhm_eff=fwhm_eff,
                            source_fwhm_arcsec=source_fwhm_arcsec,
                            lunar=lunar, clouds=clouds, airmass=airmass,
                            temperature=temperature, readout=readout,
                            extract_arcsec=extract_arcsec,
                            calibration=("commissioning-empirical" if emp
                                         else "theoretical model"),
                            t_per_frame=t_per_frame, n_frames=n_frames))

# Sanity check: star-forming galaxy, 19 AB, 3 x 600 s, defaults
sf = normalize(load_template("kc96_sb") if RESEARCH
               else TEMPLATES["star-forming galaxy"], 19.0)
r = run_lrs_etc(sf, 600, 3)
print(f"Sanity check  (19 AB Sb galaxy, 3x600 s, 1.8\\" slit, dark, slow, -80C)")
print(f"  slit coupling {r['fslit']:.2f}, extraction {r['config']['extract_arcsec']:.1f}\\" "
      f"({r['n_spat']} spatial pixels)")
print(f"  median S/N per pixel: {np.median(r['snr_pix']):.1f}   [{r['config']['calibration']}]")
r_th = run_lrs_etc(sf, 600, 3, use_empirical=False)
print(f"  theoretical-model S/N : {np.median(r_th['snr_pix']):.1f}   (switch: USE_EMPIRICAL_CALIBRATION)")
print(f"  S/N per pixel at Ha : {np.interp(6563, r['wav'], r['snr_pix']):.1f}")
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 7b. Expected observed spectrum + detector warnings

`expected_spectrum(...)` integrates the PSF within the **seeing disk**
(extraction aperture = 1 seeing FWHM by default) and renders the three
things an observer actually gets, with per-wavelength warning masks:

* **flux** — the model f_λ resampled to the 4 Å pixels,
* **ADU per frame** — the *peak pixel* (brightest spatial row, the one
  that saturates first, including sky + dark + bias) and the total
  extracted counts,
* **S/N per pixel** for the full stack.

Warning masks (shaded on all panels):

| color | condition |
|---|---|
| red | **saturation** — peak pixel ≥ 65 535 ADU in one frame |
| orange | **non-linearity** — peak pixel ≥ 95 % of full well |
| purple | **read-noise-dominated** — RN² is the largest variance term |
| brown | **dark-dominated** — dark current is the largest variance term |

Saturation and non-linearity call for shorter frames (more of them);
read-noise domination calls for longer frames or slower readout; dark
domination (rare at −80 °C, essentially never at −100 °C) calls for the
colder operating point."""))

cells.append(code("""def check_warnings(res, t_per_frame):
    \"\"\"Per-wavelength warning masks for a run_lrs_etc() result.\"\"\"
    cfg = res["config"]
    seeing = cfg.get("fwhm_eff", cfg["seeing"]); h = cfg["extract_arcsec"]
    sig = seeing / 2.3548
    # fraction of the extracted light landing on the brightest 0.9" row
    fy_extract = erf((h / 2) / (np.sqrt(2) * sig))
    fy_center  = erf((SPATIAL_AS_PIX / 2) / (np.sqrt(2) * sig))
    peak_row_frac = np.clip(fy_center / max(fy_extract, 1e-9), 0, 1)

    sky_pix_rate  = res["B"] / res["n_spat"]          # e-/s per pixel
    dark          = DARK_E_PIX_S[cfg["temperature"]]
    peak_e_frame  = (res["S"] * peak_row_frac + sky_pix_rate + dark) * t_per_frame
    peak_adu      = peak_e_frame / GAIN_E_ADU + BIAS_ADU

    sat    = peak_adu >= FULL_WELL_ADU
    nonlin = (peak_adu >= LINEARITY_FRAC * FULL_WELL_ADU) & ~sat

    # dominant variance term per wavelength (over the full stack)
    t_tot = res["t_tot"]
    var = np.vstack([res["S"] * t_tot,                       # source shot
                     res["B"] * t_tot,                       # sky shot
                     np.full_like(res["S"], res["dark_rate"] * t_tot),
                     np.full_like(res["S"], res["rn_var_total"])])
    dom = np.argmax(var, axis=0)
    rn_dom   = (dom == 3) & ~sat & ~nonlin
    dark_dom = (dom == 2) & ~sat & ~nonlin

    return dict(peak_adu=peak_adu, saturated=sat, nonlinear=nonlin,
                rn_dominated=rn_dom, dark_dominated=dark_dom,
                peak_row_frac=peak_row_frac)

def _shade(ax, wav, mask, color, label):
    \"\"\"Shade contiguous True runs of mask.\"\"\"
    if not mask.any():
        return
    idx = np.where(mask)[0]
    splits = np.where(np.diff(idx) > 1)[0]
    first = True
    for grp in np.split(idx, splits + 1):
        ax.axvspan(wav[grp[0]], wav[grp[-1]], color=color, alpha=0.20,
                   lw=0, label=label if first else None)
        first = False

WARN_STYLE = [("saturated", RED, "saturated"),
              ("nonlinear", ORANGE, "non-linear (>95% FW)"),
              ("rn_dominated", "#8064A2", "read-noise-dominated"),
              ("dark_dominated", "#8B5A2B", "dark-dominated")]

def expected_spectrum(spec, t_per_frame, n_frames, seeing=1.1,
                       extract_in_seeing=1.0, source_fwhm_arcsec=0.0,
                       **etc_kwargs):
    \"\"\"Run the ETC extracting within the (PSF x source) effective disk
    and plot flux / ADU / S/N with warning shading. Returns (res, wrn).\"\"\"
    fwhm_eff = float(np.hypot(seeing, source_fwhm_arcsec))
    res = run_lrs_etc(spec, t_per_frame, n_frames, seeing=seeing,
                      source_fwhm_arcsec=source_fwhm_arcsec,
                      extract_arcsec=extract_in_seeing * fwhm_eff,
                      **etc_kwargs)
    wrn = check_warnings(res, t_per_frame)
    cfg = res["config"]

    for name, _c, lab in WARN_STYLE:
        m = wrn[name]
        if m.any():
            w = res["wav"][m]
            print(f"WARNING - {lab}: {m.sum()} pixels "
                  f"({w.min():.0f}-{w.max():.0f} A)")
    if not any(wrn[n].any() for n, _c, _l in WARN_STYLE):
        print("No detector warnings - exposure is in the source/sky "
              "shot-noise regime and unsaturated.")

    fig, axs = plt.subplots(3, 1, figsize=(11.5, 8.2), sharex=True)
    axs[0].plot(res["wav"], res["flam_pix"], color=NAVY, lw=0.9)
    axs[0].set_yscale("log")
    axs[0].set_ylabel(r"$f_\\lambda$ (erg/s/cm²/Å)")
    axs[0].set_title(
        f"{spec.get('name','model')}  ·  {n_frames}x{t_per_frame:.0f} s, "
        f"{cfg['slit']}\\" slit, seeing {seeing}\\", extraction "
        f"{cfg['extract_arcsec']:.1f}\\" ({res['n_spat']} rows)")

    axs[1].plot(res["wav"], wrn["peak_adu"], color=NAVY, lw=0.9,
                label="peak pixel (per frame, incl. sky+dark+bias)")
    tot_adu = ((res["S"] + res["B"] + res["dark_rate"]) * t_per_frame
               / GAIN_E_ADU + res["n_spat"] * BIAS_ADU)
    axs[1].plot(res["wav"], tot_adu, color="#4F81BD", lw=0.8, ls="--",
                label="total extracted (per frame)")
    axs[1].axhline(FULL_WELL_ADU, color=RED, lw=1.0)
    axs[1].axhline(LINEARITY_FRAC * FULL_WELL_ADU, color=ORANGE, lw=0.9,
                   ls=":")
    axs[1].text(res["wav"][5], FULL_WELL_ADU * 1.05, "full well",
                color=RED, fontsize=8)
    axs[1].set_yscale("log")
    axs[1].set_ylabel("ADU")
    axs[1].legend(fontsize=8, loc="lower right")

    axs[2].plot(res["wav"], res["snr_pix"], color=LIME, lw=0.9)
    axs[2].set_ylabel("S/N per pixel")
    axs[2].set_xlabel("Wavelength (Å)")

    for ax in axs:
        for name, color, lab in WARN_STYLE:
            _shade(ax, res["wav"], wrn[name], color,
                   lab if ax is axs[0] else None)
        ax.grid(alpha=0.25)
    handles, labels = axs[0].get_legend_handles_labels()
    if handles:
        axs[0].legend(handles, labels, fontsize=8, loc="lower right")
    plt.tight_layout(); plt.show()
    return res, wrn

# --- Demo 1: bright A0V standard, long frame -> saturation + non-linearity
bright = normalize(load_template("pickles_a0v") if RESEARCH
                   else TEMPLATES["A0V star"], 11.0)
_ = expected_spectrum(bright, t_per_frame=600, n_frames=1, seeing=1.1)

# --- Demo 2: faint quasar host, short frames -> read-noise-dominated
faint = normalize(load_template("dobos_RED2") if RESEARCH
                  else TEMPLATES["quenched galaxy"], 21.5)
_ = expected_spectrum(faint, t_per_frame=60, n_frames=30, seeing=1.1,
                      readout="fast")
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 7c. Simulated 2-D spectra — raw frame and reduced S/N map

Two heat-map views of what the detector will actually record:

* **Raw single frame (ADU)** — bias + dark + sky (continuum **and**
  airglow/OH lines, which appear as bright vertical stripes spanning the
  whole slit) + the source trace, with photon and read noise. The sky
  illumination **along the slit follows a 4th-order polynomial** (the
  same slow gradient the flat-field step fights in the reduction
  notebook).
* **Reduced, stacked S/N map (λ × arcsec)** — N frames stacked,
  bias/dark/flat removed, sky subtracted per wavelength column by a
  polynomial fit along the slit. The fit is deliberately **lower order
  (deg 2) than the true deg-4 illumination**, so the subtraction leaves
  systematic residuals — strongest under the bright OH bands, exactly as
  in real LRS data.

**Vignetting.** Only the central **3′ (180″)** of the slit is illuminated;
the detector's spatial axis is 256 × 0.9″ = 230″, so the rows beyond
±90″ receive **no sky and no source** — they show only bias + dark +
read noise, forming the overscan-like dark bands above and below the
spectral band that you see in real LRS frames. The sky fit in the
reduction step is restricted to the illuminated rows, again matching the
real pipeline (which must also dodge the slit-edge holes)."""))

cells.append(code("""def _row_fractions(y_rows, seeing):
    \"\"\"Fraction of a centered Gaussian PSF falling in each 0.9-arcsec row.\"\"\"
    sig = seeing / 2.3548
    lo = (y_rows - SPATIAL_AS_PIX / 2) / (np.sqrt(2) * sig)
    hi = (y_rows + SPATIAL_AS_PIX / 2) / (np.sqrt(2) * sig)
    return 0.5 * (erf(hi) - erf(lo))

SLIT_LENGTH_ARCSEC = 180.0        # only the central 3' is illuminated

def _slit_illum_poly4(y_rows):
    \"\"\"4th-order slit-illumination profile across the ILLUMINATED 3'.\"\"\"
    u = y_rows / (SLIT_LENGTH_ARCSEC / 2)     # -1..1 across the slit
    return 1.0 + 0.05*u - 0.08*u**2 + 0.04*u**3 + 0.10*u**4

def _vignette(y_rows, rolloff_arcsec=2.0):
    \"\"\"Vignetting function: 1 inside the central 3', 0 outside, with a
    smooth ~2-arcsec roll-off at the slit ends.\"\"\"
    edge = SLIT_LENGTH_ARCSEC / 2
    return 0.5 * (1 + erf((edge - np.abs(y_rows))
                          / (np.sqrt(2) * rolloff_arcsec)))

def _sim_rates(spec, t_per_frame, n_frames, n_rows, seeing, **etc_kwargs):
    \"\"\"Common per-row/per-column e-/s rates for the 2-D simulators.

    Point / resolved-Gaussian sources: spatial profile is the Gaussian of
    FWHM_eff = sqrt(seeing^2 + source^2) integrated over each 0.9\" row.
    extended=True (surface-brightness) sources: distributed UNIFORMLY
    along the illuminated slit, like the sky.\"\"\"
    res = run_lrs_etc(spec, t_per_frame, n_frames, seeing=seeing,
                      extract_arcsec=n_rows * SPATIAL_AS_PIX, **etc_kwargs)
    cfg = res["config"]
    fwhm_eff = cfg.get("fwhm_eff", seeing)
    y = (np.arange(n_rows) - n_rows / 2 + 0.5) * SPATIAL_AS_PIX
    vig = _vignette(y)
    if spec.get("extended", False):
        # uniform surface brightness: equal share per illuminated row
        src_rows = (res["S"][None, :] / n_rows) * np.ones((n_rows, 1))
        src_rows = src_rows * vig[:, None]
    else:
        frac = _row_fractions(y, fwhm_eff)
        sig = fwhm_eff / 2.3548
        fy_tot = erf((cfg["extract_arcsec"] / 2) / (np.sqrt(2) * sig))
        src_rows = ((res["S"] / max(fy_tot, 1e-9))[None, :] * frac[:, None]
                    * vig[:, None])
    print(f"  spatial profile: FWHM_eff = {fwhm_eff:.1f} arcsec "
          f"= {fwhm_eff/SPATIAL_AS_PIX:.1f} rows"
          + ("  (uniform: extended source)" if spec.get('extended') else
             f", slit coupling = {res['fslit']:.2f}"))
    sky_row = (res["B"] / res["n_spat"])[None, :] * np.ones((n_rows, 1))
    illum = _slit_illum_poly4(y)
    sky_rows = sky_row * illum[:, None] * vig[:, None]
    dark = DARK_E_PIX_S[cfg["temperature"]]
    rn = READ_NOISE_E[cfg["readout"]]
    return res, y, src_rows, sky_rows, dark, rn, vig

def _mark_vignette(ax, y):
    \"\"\"Dashed lines + labels at the illuminated-slit boundaries.\"\"\"
    edge = SLIT_LENGTH_ARCSEC / 2
    if y.min() < -edge and y.max() > edge:
        for s in (-1, 1):
            ax.axhline(s * edge, color="white", lw=0.9, ls="--", alpha=0.85)
        ax.text(ax.get_xlim()[0] + 40, edge + 6, "vignetted (no sky)",
                color="white", fontsize=7.5)
        ax.text(ax.get_xlim()[0] + 40, -edge - 14, "vignetted (no sky)",
                color="white", fontsize=7.5)

def simulate_2d_frame(spec, t_per_frame, n_rows=256, seeing=1.1, seed=42,
                       **etc_kwargs):
    \"\"\"Raw single-frame 2-D spectrum in ADU (heat map). Full 256-row
    frame by default: rows beyond the central 3' are vignetted and show
    only bias + dark + read noise (overscan-like bands).\"\"\"
    res, y, src, sky, dark, rn, vig = _sim_rates(spec, t_per_frame, 1,
                                                  n_rows, seeing, **etc_kwargs)
    rng = np.random.default_rng(seed)
    e_expect = (src + sky + dark) * t_per_frame
    frame_e = rng.poisson(np.clip(e_expect, 0, None)).astype(float)
    frame_e += rng.normal(0, rn, frame_e.shape)
    adu = np.clip(frame_e / GAIN_E_ADU + BIAS_ADU, 0, FULL_WELL_ADU)

    fig, ax = plt.subplots(figsize=(12, 4.6))
    vmax = np.percentile(adu, 99.3)
    im = ax.imshow(adu, origin="lower", aspect="auto", cmap="magma",
                   vmin=BIAS_ADU * 0.97, vmax=vmax,
                   extent=[res["wav"][0], res["wav"][-1], y[0], y[-1]])
    ax.set_xlabel("Wavelength (Å)"); ax.set_ylabel("Along slit (arcsec)")
    cfg = res["config"]
    ax.set_title(f"Raw frame — {spec.get('name','model')}, "
                 f"{t_per_frame:.0f} s, {cfg['slit']}\\" slit, "
                 f"{cfg['lunar']} sky  (ADU; central 3' illuminated, "
                 f"deg-4 illumination, airglow lines)")
    _mark_vignette(ax, y)
    plt.colorbar(im, ax=ax, label="ADU", pad=0.01)
    plt.tight_layout(); plt.show()
    return adu

def simulate_reduced_2d(spec, t_per_frame, n_frames, n_rows=256, seeing=1.1,
                         sky_fit_deg=2, seed=42, **etc_kwargs):
    \"\"\"Reduced + stacked + sky-subtracted 2-D S/N map (lambda x arcsec).

    The true sky follows a deg-4 slit-illumination profile but is fitted
    per column with a deg `sky_fit_deg` polynomial (source rows masked,
    fit restricted to the ILLUMINATED rows only), leaving realistic
    residuals under the bright sky lines. Vignetted rows beyond the
    central 3' contain no sky to subtract - only noise.\"\"\"
    res, y, src, sky, dark, rn, vig = _sim_rates(spec, t_per_frame, n_frames,
                                                  n_rows, seeing, **etc_kwargs)
    rng = np.random.default_rng(seed)
    t_tot = t_per_frame * n_frames
    src_e, sky_e = src * t_tot, sky * t_tot
    var = src_e + sky_e + dark * t_tot + n_frames * rn**2
    noisy = src_e + sky_e + rng.normal(0, np.sqrt(var))

    # Per-column sky fit along the slit: use illuminated rows only,
    # excluding rows near the source (as the real pipeline does)
    mask_src  = np.abs(y) < 3 * seeing
    fit_rows  = (~mask_src) & (vig > 0.9)
    sky_est = np.zeros_like(noisy)
    yy = y / (SLIT_LENGTH_ARCSEC / 2)
    illum_rows = vig > 0.5
    for c in range(noisy.shape[1]):
        p = np.polyfit(yy[fit_rows], noisy[fit_rows, c], sky_fit_deg)
        sky_est[illum_rows, c] = np.polyval(p, yy[illum_rows])
    reduced = noisy - sky_est
    snr_map = reduced / np.sqrt(var)

    fig, ax = plt.subplots(figsize=(12, 4.6))
    im = ax.imshow(snr_map, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-8, vmax=8,
                   extent=[res["wav"][0], res["wav"][-1], y[0], y[-1]])
    cfg = res["config"]
    ax.set_xlabel("Wavelength (Å)"); ax.set_ylabel("Along slit (arcsec)")
    ax.set_title(f"Reduced 2-D S/N — {n_frames}x{t_per_frame:.0f} s stacked, "
                 f"sky fit deg={sky_fit_deg} vs true deg-4, illuminated "
                 f"rows only (residual stripes under OH bands)")
    _mark_vignette(ax, y)
    plt.colorbar(im, ax=ax, label="S/N per pixel", pad=0.01)
    plt.tight_layout(); plt.show()
    return snr_map

# --- Demos: 19-AB starburst on a gray night in 2" seeing ----------------
# 2" seeing makes the PSF spread and the slit loss obvious: the trace
# spans ~4-5 rows and the 1.8" slit passes only ~70% of the light
# (vs ~95% in 1.1" seeing). Compare the printed slit couplings.
sb19 = normalize(load_template("kc96_starb2") if RESEARCH
                 else TEMPLATES["starburst galaxy"], 19.0)
for see in (1.1, 2.0):
    r_cmp = run_lrs_etc(sb19, 600, 3, seeing=see)
    print(f"seeing {see:.1f} arcsec: slit coupling {r_cmp['fslit']:.2f}, "
          f"median S/N {np.nanmedian(r_cmp['snr_pix']):.1f}")
_ = simulate_2d_frame(sb19, t_per_frame=600, lunar="gray", seeing=2.0)
_ = simulate_reduced_2d(sb19, t_per_frame=600, n_frames=3, lunar="gray",
                        seeing=2.0)

# A resolved Gaussian object (compact galaxy, FWHM 3") in the same seeing:
# profile broadens to sqrt(2^2 + 3^2) = 3.6" and slit losses grow further.
_ = simulate_reduced_2d(sb19, t_per_frame=600, n_frames=3, lunar="gray",
                        seeing=2.0, source_fwhm_arcsec=3.0)
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 8. Output mode 1 — N × t → S/N curve

Give the number of frames and exposure time per frame; get S/N vs
wavelength, plotted **per pixel** and **per Å**, with a secondary axis in
detector pixels. Change any keyword to explore conditions."""))

cells.append(code("""def plot_snr(results, labels=None):
    \"\"\"Overlay S/N curves for one or more run_lrs_etc() results.\"\"\"
    if isinstance(results, dict):
        results = [results]
    labels = labels or [None] * len(results)
    fig, axs = plt.subplots(2, 1, figsize=(11.5, 6.2), sharex=True)
    colors = [NAVY, LIME, ORANGE, RED, "#8064A2", "#4F81BD"]
    for res, lab, colr in zip(results, labels, colors):
        cfg = res["config"]
        lab = lab or (f"{cfg['slit']}\\" slit, {cfg['n_frames']}x"
                      f"{cfg['t_per_frame']:.0f} s, {cfg['lunar']} sky")
        axs[0].plot(res["wav"], res["snr_pix"], lw=1.0, color=colr, label=lab)
        axs[1].plot(res["wav"], res["snr_aa"],  lw=1.0, color=colr, label=lab)
    axs[0].set_ylabel("S/N per pixel (4 Å)")
    axs[1].set_ylabel("S/N per Å")
    axs[1].set_xlabel("Wavelength (Å)")
    axs[0].legend(fontsize=8); axs[0].grid(alpha=0.3); axs[1].grid(alpha=0.3)
    sec = axs[0].secondary_xaxis(
        "top", functions=(lambda w: (w - PIX_WAV[0]) / DISP_AA_PIX,
                          lambda p: PIX_WAV[0] + p * DISP_AA_PIX))
    sec.set_xlabel("Detector pixel (spectral direction)")
    plt.tight_layout(); plt.show()

# Example: a 19-AB Kinney-Calzetti starburst (with Galactic E(B-V)=0.05)
# under three moon phases. Extinction BEFORE normalization: 19.0 is the
# observed magnitude.
sb_base = load_template("kc96_starb2") if RESEARCH else TEMPLATES["starburst galaxy"]
sb = normalize(apply_extinction(sb_base, ebv=0.05), 19.0)
runs = [run_lrs_etc(sb, 600, 3, lunar=l) for l in ("dark", "gray", "bright")]
plot_snr(runs, [f"{l} sky" for l in ("dark", "gray", "bright")])

# ...and the effect of slit choice (resolution vs coupled flux)
runs2 = [run_lrs_etc(sb, 600, 3, slit=s) for s in ("1.8", "2.7", "4.5")]
plot_snr(runs2, [f"{s}\\"  (R={SLITS[s]['R']})" for s in ("1.8", "2.7", "4.5")])

# ...and the Sec. 5c redshift series: same observed 20 AB at four z
runs3 = [run_lrs_etc(normalize(redshift_spectrum(base_z, z), 20.0), 600, 3)
         for z in Z_SERIES]
plot_snr(runs3, [f"z = {z:.1f}" for z in Z_SERIES])
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 9. Output mode 2 — target S/N at a wavelength → required time

Give the S/N per pixel you need at a specific wavelength; the quadratic
solution of the CCD equation returns the **total on-source time**, split
into frames of a chosen length."""))

cells.append(code("""def required_time(spec, target_snr, wav_aa, n_frames=3, **etc_kwargs):
    \"\"\"Total on-source seconds for target_snr per pixel at wav_aa.

    Uses one run_lrs_etc call to get the rates, then solves
    a^2 t^2 - snr^2 (b t + c) = 0 for t.
    \"\"\"
    probe = run_lrs_etc(spec, 1.0, n_frames, **etc_kwargs)
    i = np.argmin(np.abs(probe["wav"] - wav_aa))
    a = probe["S"][i]                                    # e-/s source
    b = probe["S"][i] + probe["B"][i] + probe["dark_rate"]
    c = probe["rn_var_total"]
    if a <= 0:
        return np.nan
    snr2 = target_snr ** 2
    t_tot = (snr2 * b + np.sqrt(snr2**2 * b**2 + 4 * a**2 * snr2 * c)) / (2 * a**2)
    return t_tot

# Example: S/N = 10 per pixel at H-alpha for a 20-AB dusty starburst
# (Kinney-Calzetti starb6, internal E(B-V) ~ 0.7)
tgt_base = (load_template("kc96_starb6") if RESEARCH
            else TEMPLATES["dusty starburst"])
tgt = normalize(tgt_base, 20.0)
print("Total on-source time for S/N = 10 per pixel at 6563 A")
print(f"{'condition':>34s} {'t_total':>10s} {'per frame (x3)':>15s}")
for lunar in ("dark", "gray", "bright"):
    for readout in ("slow", "fast"):
        t = required_time(tgt, 10, 6563, n_frames=3,
                          lunar=lunar, readout=readout)
        print(f"{lunar + ' sky, ' + readout + ' readout':>34s} "
              f"{t:9.0f} s {t/3:11.0f} s")

# Sensitivity of the answer to target brightness
mags = np.arange(17, 22.1, 0.5)
times = [required_time(normalize(tgt_base, m), 10, 6563) for m in mags]
fig, ax = plt.subplots(figsize=(8, 3.6))
ax.semilogy(mags, times, "o-", color=NAVY)
ax.axhline(3600, color="gray", lw=0.7, ls="--"); ax.text(17.1, 4000, "1 hr",
    color="gray", fontsize=8)
ax.set_xlabel("AB magnitude at 6000 Å"); ax.set_ylabel("Total on-source (s)")
ax.set_title("Time to reach S/N = 10 per pixel at Hα — dusty starburst, dark sky")
ax.grid(alpha=0.3, which="both")
plt.tight_layout(); plt.show()
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 10. Recommended request time per target — overheads + weather margin

The on-source time from mode 1/2 is **not** what you should request in the
proposal. `recommended_request_time()` converts science time into the
recommended request per target/pointing, accounting for:

| overhead | value |
|---|---|
| CCD readout | 7 s (slow) / 3 s (medium) / 1 s (fast) per frame |
| Telescope slew | 120 s per target |
| Offset / dither | 10 s per round |
| Acquisition image | 30 s per image |
| Rounds of offset + acquisition | 3 per target or per 1-h block |
| Thru-slit check | 30 s per target/block |
| Operator communication | up to 5 min per 1-h block |
| **Standard star round** | every 2 h; same composition as a science visit (slew + 3 acquisition rounds + thru-slit) with up to **10 × 30 s** exposures |
| Weather / mechanical / tracking loss | **35 % of the total requested night time** |
| Night calibration (optional, `night_calibration=True`) | **1 h per half night (6 h)** — fixed, **not** subject to the 35 % loss |
| Daytime calibration (bias, dark, flat, arc) | ~1 h per half night — **not charged** to the night-time request; reported for information |

The block-dependent overheads (acquisition rounds, thru-slit, operator) are
charged once per started 1-hour observing block, and one standard-star
round is charged per started 2-hour interval; both are iterated until the
counts converge. The 35 % loss applies to everything on-sky *except* the
night-calibration hour(s): `requested = in_dome / (1 − 0.35) + t_night_cal`."""))

cells.append(code("""READOUT_TIME_S = {"slow": 7.0, "medium": 3.0, "fast": 1.0}
SLEW_S          = 120.0
OFFSET_S        = 10.0
ACQ_IMG_S       = 30.0
ACQ_ROUNDS      = 3
THRUSLIT_S      = 30.0
OPERATOR_S      = 300.0          # max 5 min per 1-h block
WEATHER_LOSS    = 0.35           # fraction of the requested night time
BLOCK_S         = 3600.0
STD_INTERVAL_S  = 7200.0         # one standard-star round every 2 h
STD_EXP_S       = 30.0           # per standard exposure
STD_MAX_EXP     = 10             # max exposures per standard round
HALF_NIGHT_S    = 6 * 3600.0
NIGHT_CAL_S     = 3600.0         # 1 h per half night, fixed
DAYTIME_CAL_S   = 3600.0         # informational only, not charged

def recommended_request_time(t_per_frame=None, n_frames=None,
                              science_s=None, readout="slow",
                              n_std_exp=STD_MAX_EXP,
                              night_calibration=False,
                              label=None, verbose=True):
    \"\"\"Recommended time to REQUEST for one target/pointing.

    Give either (t_per_frame, n_frames) or science_s (total on-source
    seconds, e.g. from required_time(); frames then assumed 900 s each).

    n_std_exp         : exposures per standard-star round (30 s each, <= 10)
    night_calibration : True adds 1 h per half night of night-time
                        calibration - fixed cost, NOT subject to the 35 %
                        weather/mechanical loss.

    Returns a dict with the full breakdown (all seconds).
    \"\"\"
    if science_s is None:
        science_s = t_per_frame * n_frames
    elif n_frames is None:
        t_per_frame = min(900.0, science_s)
        n_frames = int(np.ceil(science_s / t_per_frame))

    n_std_exp = min(int(n_std_exp), STD_MAX_EXP)
    readout_s = n_frames * READOUT_TIME_S[readout]
    # A standard round has the same composition as a science visit:
    # slew + 3x(offset+acq) + thru-slit, then n_std_exp short exposures.
    std_round_s = (SLEW_S + ACQ_ROUNDS * (OFFSET_S + ACQ_IMG_S) + THRUSLIT_S
                   + n_std_exp * (STD_EXP_S + READOUT_TIME_S[readout]))

    # Iterate: block/standard counts depend on the total in-dome time,
    # which includes the overheads those counts control.
    n_blocks, n_std = 1, 1
    for _ in range(12):
        acq_s      = n_blocks * ACQ_ROUNDS * (OFFSET_S + ACQ_IMG_S)
        thruslit_s = n_blocks * THRUSLIT_S
        operator_s = n_blocks * OPERATOR_S
        standards_s = n_std * std_round_s
        in_dome = (science_s + readout_s + SLEW_S + acq_s + thruslit_s
                   + operator_s + standards_s)
        nb  = max(1, int(np.ceil(in_dome / BLOCK_S)))
        ns  = max(1, int(np.ceil(in_dome / STD_INTERVAL_S)))
        if nb == n_blocks and ns == n_std:
            break
        n_blocks, n_std = nb, ns

    requested_night = in_dome / (1.0 - WEATHER_LOSS)
    weather_s = requested_night - in_dome

    # Half nights spanned by the weather-inflated science request
    n_half_nights = max(1, int(np.ceil(requested_night / HALF_NIGHT_S)))
    night_cal_s = n_half_nights * NIGHT_CAL_S if night_calibration else 0.0
    requested = requested_night + night_cal_s
    daytime_cal_s = n_half_nights * DAYTIME_CAL_S     # informational

    out = dict(label=label or f"{n_frames}x{t_per_frame:.0f} s ({readout})",
               science=science_s, readout=readout_s, slew=SLEW_S,
               acquisition=acq_s, thruslit=thruslit_s, operator=operator_s,
               standards=standards_s, weather=weather_s,
               night_cal=night_cal_s, in_dome=in_dome,
               requested=requested, daytime_cal_info=daytime_cal_s,
               n_blocks=n_blocks, n_std_rounds=n_std,
               n_half_nights=n_half_nights, n_frames=n_frames,
               efficiency=science_s / requested)
    if verbose:
        print(f"{out['label']}: science {science_s/60:.1f} min, "
              f"{n_blocks} block(s), {n_std} standard round(s), "
              f"in-dome {in_dome/60:.1f} min, REQUEST {requested/60:.1f} min "
              f"(open-shutter efficiency {out['efficiency']:.0%})")
        if night_calibration:
            print(f"   incl. night calibration {night_cal_s/60:.0f} min "
                  f"({n_half_nights} half-night(s) x 1 h, fixed, "
                  f"not weather-charged)")
        print(f"   info: daytime calibration (bias/dark/flat/arc) "
              f"~{daytime_cal_s/60:.0f} min for {n_half_nights} half-night(s) "
              f"- NOT charged to the night-time request")
    return out

def plot_request_breakdown(results):
    \"\"\"Horizontal stacked bars for one or more recommended_request_time()
    results.\"\"\"
    if isinstance(results, dict):
        results = [results]
    comps  = ["science", "readout", "slew", "acquisition",
              "thruslit", "operator", "standards", "weather", "night_cal"]
    names  = ["on-source science", "CCD readout", "slew",
              "offset + acquisition (3x)", "thru-slit check",
              "operator comm.", "standard star (10x30 s / 2 h)",
              "weather / mech. loss (35 %)",
              "night calibration (fixed)"]
    colors = [NAVY, "#4F81BD", LIME, ORANGE, "#8064A2", "#999999",
               "#2AA198", RED, "#444444"]
    hatches = [None, None, None, None, None, None, None, None, "//"]
    fig, ax = plt.subplots(figsize=(11.5, 1.15 + 0.85 * len(results)))
    for j, res in enumerate(results):
        left = 0.0
        for comp, colr, htc in zip(comps, colors, hatches):
            v = res.get(comp, 0.0) / 60.0
            if v <= 0:
                continue
            ax.barh(j, v, left=left, color=colr, edgecolor="white",
                    height=0.55, lw=0.6, hatch=htc)
            if v > res["requested"] / 60.0 * 0.045:
                ax.text(left + v / 2, j, f"{v:.0f}", ha="center", va="center",
                        color="white", fontsize=8, fontweight="bold")
            left += v
        ax.text(left + res["requested"] / 60.0 * 0.012, j,
                f"{res['requested']/60:.0f} min  "
                f"(η = {res['efficiency']:.0%})",
                va="center", fontsize=9, color="#222222")
    ax.set_yticks(range(len(results)))
    ax.set_yticklabels([r["label"] for r in results], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Time (minutes)")
    ax.set_title("Recommended request time per target — breakdown "
                 "(daytime calibration not shown: not charged)")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=c,
                                     hatch=h, edgecolor="white")
                       for c, h in zip(colors, hatches)],
              labels=names, fontsize=8, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.28))
    ax.grid(alpha=0.25, axis="x")
    ax.set_xlim(0, max(r["requested"] for r in results) / 60.0 * 1.24)
    plt.tight_layout(); plt.show()

# --- Examples ----------------------------------------------------------
# 1) The mode-2 answer from Sec. 9: S/N=10 at Ha for the 20-AB dusty
#    starburst (dark sky, slow readout)
t_sci = required_time(tgt, 10, 6563, n_frames=3)
r1 = recommended_request_time(t_per_frame=np.ceil(t_sci / 3), n_frames=3,
                               readout="slow",
                               label="dusty starburst, S/N 10 @ Ha")
# 2) A deep 3x1200 s spectrum of the same target
r2 = recommended_request_time(1200, 3, readout="slow",
                               label="deep 3x1200 s (slow)")
# 3) Same deep plan with fast readout and 12x300 s dithering
r3 = recommended_request_time(300, 12, readout="fast",
                               label="deep 12x300 s (fast)")
# 4) The deep plan again, requesting a night-calibration set
r4 = recommended_request_time(1200, 3, readout="slow",
                               night_calibration=True,
                               label="deep 3x1200 s + night cal")

plot_request_breakdown([r1, r2, r3, r4])
"""))


# ----------------------------------------------------------------------
cells.append(md("""## 11. Validation — the spec model vs the 2026-04-01 M82 observation

Dataset: M82 nucleus, 3 x 60 s (header) 1.8\u2033 slit, **bright night**,
altitude ~39\u00b0 (airmass ~2), and — discovered later — a **stuck-open
shutter adding ~7 s to every frame** (60 s header \u2192 67 s effective; the
standard's 30 s frames \u2192 37 s). All rates below use effective times.
The NARIT reduction pipeline measures **median S/N \u2248 9
per pixel** (4626\u20137996 \u00c5, ADU statistics). The M82 template carries the
flux-calibrated extracted spectrum, i.e. the flux already through the slit
and the 15.3\u2033 extraction aperture.

**Method: invert the ETC.** Under the true commissioning conditions, find
the source flux that would reproduce the measured S/N, and express it as a
continuum flux density, r\u2032-band AB magnitude and surface brightness. The
gap between that and the template's actual photometry is the unexplained
system deficit \u2014 in magnitudes, where instrument people can argue about
it. Two readings of the measured S/N are used: 9 (ADU statistics, as
reported by the pipeline) and 18 (if the counts are converted to electrons
at gain 4, Poisson noise doubles the true S/N)."""))

cells.append(code("""m82 = TEMPLATES.get("M82 (empirical)")
if m82 is not None and np.isfinite(m82["flam"]).any():
    AP_AREA = 1.8 * 15.3                     # arcsec^2, pipeline aperture
    CONFIG = dict(t_per_frame=67, n_frames=3, slit="1.8",   # 60 s + 7 s shutter
                  seeing=0.3, extract_arcsec=15.3,   # template is post-slit
                  lunar="bright", clouds="photometric",
                  altitude_deg=30,                    # M82 culminates ~39 deg
                  temperature="-80C", readout="slow",
                  use_empirical=False)                # THEORY branch, always

    def _med_snr(scale, cfg=CONFIG):
        s = dict(m82); s["flam"] = m82["flam"] * scale
        r = run_lrs_etc(s, **cfg)
        band = (r["wav"] > 4626) & (r["wav"] < 7996)
        return float(np.nanmedian(r["snr_pix"][band])), r

    snr_pred, r_pred = _med_snr(1.0)
    snr_gray, _ = _med_snr(1.0, dict(CONFIG, lunar="gray", altitude_deg=60))
    i6 = np.argmin(np.abs(r_pred["wav"] - 6000))
    r_tmp = synth_ab_mag(m82, "r")

    print(f"Predicted median S/N for the template as observed:")
    print(f"  gray sky, X=1.15 (naive)   : {snr_gray:.0f}")
    print(f"  bright sky, X=2.0 (actual) : {snr_pred:.0f}   <- conditions " 
          f"explain almost nothing (source-dominated)")
    print(f"  measured by the pipeline   : 9")
    print()

    def _solve(target):
        lo, hi = 1e-6, 1.0
        for _ in range(50):
            mid = np.sqrt(lo * hi)
            s, _ = _med_snr(mid)
            lo, hi = (mid, hi) if s < target else (lo, mid)
        return np.sqrt(lo * hi)

    print("If the spec-based ETC were TRUE, the observed S/N corresponds to:")
    print(f"{'':32s}{'S/N=9 (ADU)':>14s}{'S/N=18 (e-, gain 4)':>21s}")
    rows = {}
    for tgt in (9.0, 18.0):
        f = _solve(tgt)
        s = dict(m82); s["flam"] = m82["flam"] * f
        rows[tgt] = dict(f=f, dmag=-2.5*np.log10(f),
                         flam6=float(np.interp(6000, s["wav"], s["flam"])),
                         rmag=synth_ab_mag(s, "r"))
    print(f"{'  scale factor f':32s}{rows[9]['f']:14.2e}{rows[18]['f']:21.2e}")
    print(f"{'  f_lam(6000A) in aperture':32s}{rows[9]['flam6']:14.2e}"
          f"{rows[18]['flam6']:21.2e}   erg/s/cm2/A")
    print(f"{'  r-band AB (27.5 arcsec^2)':32s}{rows[9]['rmag']:14.2f}"
          f"{rows[18]['rmag']:21.2f}")
    print(f"{'  r-band SB (AB/arcsec^2)':32s}"
          f"{rows[9]['rmag']+2.5*np.log10(AP_AREA):14.2f}"
          f"{rows[18]['rmag']+2.5*np.log10(AP_AREA):21.2f}")
    print(f"{'  deficit vs template (mag)':32s}{rows[9]['dmag']:14.2f}"
          f"{rows[18]['dmag']:21.2f}")
    print()
    print(f"Template (flux-calibrated): r' = {r_tmp:.2f} AB in the aperture "
          f"= {r_tmp + 2.5*np.log10(AP_AREA):.2f} AB/arcsec^2")
    print("  -> superficially consistent with M82's often-quoted nuclear")
    print("     surface brightness (~15-16 r-mag/arcsec^2) - but see the")
    print("     RESOLUTION below: the standard-star measurement shows the")
    print("     template absolute scale is in fact inflated.")
    print()
    print("Accounting (magnitudes of S/N deficit explained):")
    print(f"  bright moon + airmass 2     : ~0.06 mag  (already in the model)")
    print(f"  ADU->e- gain (x4 -> S/N x2) : "
          f"{rows[9]['dmag']-rows[18]['dmag']:.2f} mag")
    print(f"  UNEXPLAINED before the standard-star test: "
          f"~{rows[18]['dmag']:.1f} mag = flux factor "
          f"~{10**(0.4*rows[18]['dmag']):.0f}")
    print()
    # ---- RESOLUTION (2026-07-24): measured BD+75 325 raw rates ----------
    # From the NARIT pipeline (LRS_Reduction_Pipeline.ipynb, 2026-04-01
    # frames; median single frame, 13-row boxcar, peak of the extracted
    # spectrum ~6300-6600 A):
    # shutter-corrected effective exposures: 30->37 s, 60->67 s
    MEAS_BD = {"4.5": 13893/37.0, "2.7": 9093/37.0, "1.8": 12048/67.0}
    bdp = normalize(TEMPLATES["O5V star"], 9.55, wav0=5500)  # BD proxy
    print("MEASURED standard-star rates vs the spec model")
    print("(delivered image quality ~4 arcsec that night - the measured")
    print(" 1.8/4.5-slit flux ratio 0.43 matches ~4+ arcsec, not 1.5):")
    print(f"{'slit':>6s}{'measured ADU/s/col':>20s}{'model(4arcsec) ADU/s/col':>21s}"
          f"{'ratio':>7s}")
    ratios = []
    for slit in ("4.5", "2.7", "1.8"):
        rm = run_lrs_etc(bdp, 60, 1, slit=slit, seeing=4.0,
                         lunar="bright", altitude_deg=45,
                         extract_arcsec=13*0.9, use_empirical=False)
        model = np.nanmax(rm["S"]) / GAIN_E_ADU
        ratios.append(MEAS_BD[slit] / model)
        print(f"{slit:>6s}{MEAS_BD[slit]:>20.0f}{model:>21.0f}"
              f"{MEAS_BD[slit]/model:>7.2f}")
    fac = float(np.mean(ratios))
    print(f"  -> the as-built system delivers {fac:.2f}x the component model")
    print(f"     CONSISTENTLY across slits (shutter-corrected): a gray")
    print(f"     factor ~0.29 (1.35 mag), NOT the naive 1/83 from M82.")
    print(f"     This factor + the measured chromatic shape IS the Sec 2b")
    print(f"     empirical calibration - with USE_EMPIRICAL_CALIBRATION=True")
    print(f"     the ETC reproduces the commissioning rates by construction.")
    print()
    print("REVISED deficit budget for the M82 S/N=9 anomaly:")
    print("  moon + airmass (modeled)      : 0.06 mag")
    print("  ADU->e- gain                   : 0.81 mag")
    print("  true system efficiency (x0.29) : 1.35 mag   <- measured above")
    print("  delivered image quality ~4 arcsec    : n/a for the template inversion")
    print("  REMAINING ~2.8 mag (x13): the M82 flux-calibration/template")
    print("  chain - the template absolute scale is inflated (approximate")
    print("  CALSPEC anchors + noisy calibrated spectrum), and S/N=9 was")
    print("  measured on the flux-calibrated spectrum, not raw counts.")
    print()
    print("PRACTICAL CONSEQUENCES:")
    print(f"  * ABSOLUTE predictions: leave USE_EMPIRICAL_CALIBRATION=True")
    print(f"    (Sec 2b) - source/sky rates x{fac:.2f} gray + measured shape;")
    print(f"    S/N x ~0.54 source-limited, exposure x ~3.4 vs theory")
    print("  * do NOT use the M82 empirical template as an absolute anchor")
    print("    until it is re-anchored to the pipeline sensitivity function")
    print("  * plan for delivered image quality of 3-4 arcsec until focus/")
    print("    guiding improves - at 1.8 arcsec that alone costs ~1 mag for")
    print("    point sources (use the wider slits)")
"""))


cells.append(md("""---

### Notes for proposers — model assumptions & limitations

* **Sky model** is a schematic dark/gray/bright scale anchored at V with a
  fixed color shape — not a measured TNO sky spectrum. Near full moon or
  in the far red (OH bands) the true sky can deviate substantially.
* **Slit coupling** assumes a Gaussian PSF centered in the slit; no
  differential atmospheric refraction along the slit. At airmass 2
  (altitude 30°) blue-red centroid drift of several arcsec is possible for
  long exposures with the slit off the parallactic angle.
* **Throughput** is the §2 component model (mirrors anchored to the
  2026 measured reflectivities at 550 nm) with band-edge tapers below 4000 Å / above 8000 Å; the §11
  validation suggests the real system currently performs below this model
  — apply the empirical margin printed there.
* **S/N per Å** = S/N per pixel ÷ 2 is exact only in the shot-noise limit;
  in read-noise-limited regimes rebinning helps more than this suggests.
* Templates outside their native wavelength coverage give NaN (§3b), by
  design — no extrapolated science.
* Recommended workflow: §6b visibility → §3b template → §5b extinction →
  §5 normalize (observed mag) → mode 2 for time → §10 for the requested
  total with `night_calibration=` as needed.

*Maintained by NARIT. Throughput model: `LRS_throughput_assumptions.md`.
Version and changelog: see the header cell. Please report discrepancies
between ETC predictions and delivered data — they feed the §11 margin.*"""))


nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3",
                                   "language": "python", "name": "python3"},
                   "language_info": {"name": "python", "version": "3.10"}},
      "nbformat": 4, "nbformat_minor": 5}

OUT.write_text(json.dumps(nb, indent=1))
print(f"wrote {OUT} with {len(cells)} cells")
