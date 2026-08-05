<p align="center">
  <img src="assets/monet_calc_logo.png" alt="MONET Calc logo" width="260">
</p>

<h1 align="center">MONET Calc</h1>
<p align="center"><b>Module of Numerical Exposure Time Calculator</b><br>
NARIT LRS @ 2.4-m Thai National Telescope</p>

<p align="center">
  <a href="https://monet-calc.streamlit.app"><img src="https://static.streamlit.io/badges/streamlit_badge_black_white.svg" alt="Streamlit app"></a>
</p>

Observation-preparation tool for the **Low-Resolution Spectrograph (LRS)**
on the 2.4-m Thai National Telescope, Doi Inthanon. Single mode: long-slit
spectroscopy, 4000–8092 Å at 4 Å/pixel.

**Maintainer:** Krittapas Chanchaiworawit (NARIT) · **Version 6.4** (2026-07-24)

---

## Quick start — three ways to use it

**1. Python package** (see `COOKBOOK.md` for the full parameter manual):

```bash
git clone https://github.com/krc2020/LRS_ETC.git && cd LRS_ETC
pip install -e ".[full]"
```
```python
import lrs_etc as etc
spec = etc.normalize_filter(etc.load_template("kc96_sb"), 19.0, "r")
res  = etc.run_lrs_etc(spec, t_per_frame=600, n_frames=3, slit="1.8")
t    = etc.required_time(spec, target_snr=10, wav_aa=6563)
```

**2. Web app** (Streamlit — three calculators: S/N ↔ time, simulated
1-D/2-D spectra, recommended request time; all with warnings):

```bash
pip install -e ".[full,app]"
streamlit run streamlit_app.py
```
To publish: push this repo to GitHub → [share.streamlit.io](https://share.streamlit.io)
→ new app → pick `streamlit_app.py`. Data files ship in the repo, so no
extra setup is needed.

**3. Notebook** (the annotated derivation of everything):

```bash
pip install numpy scipy matplotlib astropy dust_extinction
jupyter lab LRS_ETC.ipynb
```

Everything runs top-to-bottom in ~1 minute. Recommended workflow for a
proposal:

1. **§6b** — check your target's visibility from the TNO (RA/Dec → altitude flags)
2. **§3b** — pick a template (`load_template("kc96_sb")`) or supply your own (§4)
3. **§5c** — redshift it (`redshift_spectrum`, z ≤ 9)
4. **§5b** — apply Galactic extinction (CCM89/O94/F99, R_V, E(B−V) or A_V)
5. **§5** — normalize to the observed magnitude (point AB or extended AB/arcsec²)
6. **§8/§9** — S/N for a given N × t, or required time for a target S/N
7. **§7b** — check detector warnings (saturation, non-linearity, noise regime)
8. **§10** — convert on-source time to the **requested** time (overheads,
   standard stars, optional night calibration, 35 % weather margin)

## Repository layout

| path | contents |
|---|---|
| `lrs_etc/` | **the Python package** (config · core physics · figure builders) |
| `streamlit_app.py` | **the web app** (3 calculators, sidebar target/conditions setup) |
| `COOKBOOK.md` | parameter manual for every module |
| `LRS_ETC.ipynb` | the calculator notebook (executed, with all example outputs) |
| `spectral_library/` | 49 research-grade templates: Pickles 1998 stars, Kinney–Calzetti 1996 galaxies, AGN (Francis+91 QSO composite, Seyferts, LINER), Dobos+12 SDSS composites — full native UV coverage, see its README for provenance |
| `data/M82_template_3500_9000.txt` | empirical M82 template (absolute flux, from LRS commissioning) |
| `LRS_throughput.csv` + `LRS_throughput_assumptions.md` | per-component system throughput model |
| `ETC_flowchart.{pdf,png}` | one-page visual summary of the calculator |
| `build_scripts/` | generators: notebook builder, library fetcher (spextra), flowchart |

## Instrument summary

Slits 1.8″/2.7″/4.5″ → R = 750/500/300 at 6000 Å. Detector: Andor Newton
BEX2-DD 1024×256, 0.9″/pixel spatial; gain 4 e⁻/ADU, full well 65 535 ADU,
bias ≈ 300 ADU, linear to 95 % FW; dark 0.08 (−80 °C) / 0.003 (−100 °C)
e⁻/pix/s; read noise 4/12/15 e⁻ at 50 kHz/1 MHz/3 MHz. Only the central
3′ of the slit is illuminated.

## Absolute calibration — commissioning-anchored (v6.4)

The ETC ships **commissioning-calibrated by default**
(`USE_EMPIRICAL_CALIBRATION = True`, §2b): source and sky rates are scaled
by the response ratio measured from the 2026-04-01 BD+75°325 standard
(shutter-corrected; gray factor ≈ 0.29 with a measured chromatic shape —
strong blue collapse below 5500 Å). Set the switch to False for the
theoretical component model (e.g. post-refurbishment forecasts). §11
documents the calibration and its decomposition. Please report
prediction-vs-delivery discrepancies — they refine this calibration.

## Citation

> Chanchaiworawit, K. (2026). *NARIT LRS Exposure-Time Calculator*,
> National Astronomical Research Institute of Thailand.

Template provenance: Pickles 1998, PASP 110, 863 · Kinney et al. 1996,
ApJ 467, 38 · Calzetti et al. 1994 · Francis et al. 1991, ApJ 373, 465 ·
Dobos et al. 2012, MNRAS 420, 1217 · retrieved via speXtra
(Leschinski 2021).

Contact: **krittapas [at] narit.or.th**
