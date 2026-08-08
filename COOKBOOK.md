# lrs_etc cookbook

Parameter reference for the `lrs_etc` Python package (v6.5). Everything
mirrors `LRS_ETC.ipynb`; the notebook remains the annotated derivation.

```python
import lrs_etc as etc
```

## 1. Build a spectrum

Spectra are dicts `{"wav", "flam", ...}` on a 3500–9500 Å grid.

| function | parameters | notes |
|---|---|---|
| `list_templates()` | — | 55 names: `pickles_*` (stars), `kc96_*` (galaxies), `agn_*`, `dobos_*` (SDSS composites), `brown_*` (full-coverage SEDs) |
| `load_template(name)` | `name` from the list | prints a note if the native coverage misses part of 4000–8092 Å |
| `user_table_spectrum(path)` | 2-column ASCII: λ(Å), f_λ | absolute scale irrelevant (renormalized later) |
| `continuum(kind, alpha)` | `kind`: `"flat_flam"`, `"flat_fnu"`, `"powerlaw"`; `alpha`: f_λ ∝ λ^α | |
| `add_line(spec, center_aa, ew_aa, fwhm_kms, absorption)` | EW in Å, FWHM in km/s | emission by default |

Recommended full-coverage substitutes for red-truncated templates:
`agn_qso_ext` (not `agn_qso`), `brown_NGC5033_sy1.5` (Seyfert 1),
`brown_NGC4579_liner`, `brown_NGC3379_elliptical`, `brown_NGC0628_sc`.

## 2. Redshift and extinction (order matters)

Workflow for an observed magnitude: **redshift → extinguish → normalize**.

| function | parameters |
|---|---|
| `redshift_spectrum(spec, z, d_L_ref_Mpc=None)` | `0 ≤ z ≤ 9`; missing coverage zero-patched with a named warning; `d_L_ref_Mpc` adds cosmological dimming (absolute-flux templates) |
| `redshift_limits(name, band=(4000, 8000))` | → `(z_full_lo, z_full_hi, z_any)` |
| `apply_extinction(spec, law, rv, ebv=… or av=…)` | `law`: `"CCM89"` (default), `"O94"`, `"F99"`, `"Calzetti00"` (internal attenuation); give exactly one of `ebv`/`av` |

## 3. Normalize (sets the absolute flux)

| function | parameters |
|---|---|
| `normalize(spec, mag_ab, wav0=6000, extended=False)` | monochromatic AB at λ₀ |
| `normalize_filter(spec, mag_ab, band="r", extended=False)` | synthetic photometry through the real SDSS `u g r i z` curves (Doi+2010) |
| `normalize_total_flux(spec, F, wav_range, extended=False)` | ∫f_λ dλ = F (erg s⁻¹ cm⁻²); ideal for line-only models |
| `synth_ab_mag(spec, band)` | read a synthetic magnitude/color off any spectrum |

`extended=True` interprets magnitudes/fluxes **per arcsec²** (uniform
surface brightness; no PSF/slit losses; aperture = slit × extraction).
Give `source_extent_arcsec=` in `run_lrs_etc` for an object of finite
extent along the slit — the extraction window matches the extent and the
2-D simulators keep sky rows outside the object for a valid subtraction.

## 4. Run the ETC

```python
res = etc.run_lrs_etc(spec, t_per_frame, n_frames, **options)
```

| option | values (default first) | meaning |
|---|---|---|
| `slit` | `"1.8"`, `"2.7"`, `"4.5"` | width ↔ R = 750/500/300 at 6000 Å |
| `seeing` | `1.1` | PSF FWHM (″) |
| `source_fwhm_arcsec` | `0.0` | resolved Gaussian source; profile FWHM_eff = √(seeing² + source²) |
| `psf_profile`, `moffat_beta` | `"gaussian"`, `3.5` | `"moffat"` = realistic wings (10–25 % more slit loss) |
| `lunar` | `"dark"`, `"gray"`, `"bright"` | sky continuum (airglow/OH lines always included) |
| `clouds` | `"photometric"`, `"thin cirrus"`, `"cloudy"` | gray 0 / 0.5 / 1.2 mag |
| `altitude_deg` or `airmass` | 30/45/60/75/90 → X = sec z | altitude overrides airmass |
| `temperature` | `"-80C"`, `"-100C"` | dark 0.08 / 0.003 e⁻/pix/s |
| `readout` | `"slow"`, `"medium"`, `"fast"` | RN 4/12/15 e⁻; readout 7/3/1 s |
| `source_extent_arcsec` | `None` | extended sources only: object extent along the slit (″). Sets the extraction window, confines the object band in the 2-D simulators, and reserves sky rows outside it for the sky fit. Omit → object fills the slit, no on-slit sky subtraction (offset-sky warning) |
| `extract_arcsec` | 1.5 × FWHM_eff (point) or the source extent (extended) | extraction window along the slit |
| `use_empirical` | `None` (→ config default `True`) | **True** = commissioning-calibrated (×0.29 gray + measured chromatic shape, shutter-corrected); **False** = theoretical model |

Result keys: `wav, snr_pix, snr_aa, S, B` (e⁻/s), `flam_pix, fslit,
n_spat, config`.

Helpers: `required_time(spec, target_snr, wav_aa, n_frames, **options)` →
total on-source seconds; `check_warnings(res, t_per_frame)` +
`warning_summary(...)` → saturation / non-linearity / read-noise- /
dark-dominated masks; `target_visibility(ra_deg, dec_deg)` → altitude
flags from the TNO; `select_airmass(altitude_deg)`.

## 5. Request-time budget

```python
b = etc.recommended_request_time(t_per_frame=600, n_frames=3,
                                 readout="slow", night_calibration=False)
```
Charges readout, slew 120 s, 3×(offset 10 s + acquisition 30 s),
thru-slit 30 s, operator ≤5 min per 1-h block, a standard-star round
(≤10 × 30 s) every 2 h, then ÷0.65 (35 % weather/mechanical margin).
`night_calibration=True` adds a fixed 1 h per 6-h half night *outside*
the margin. Daytime calibration is reported (`daytime_cal_info`) but never
charged. Alternative input: `science_s=` (e.g. from `required_time`).

## 6. Figures (all return `matplotlib` Figures)

| function | shows |
|---|---|
| `snr_figure(res_or_list, labels)` | S/N vs λ, per pixel + per Å |
| `expected_spectrum_figure(res, t_per_frame)` | flux / ADU (peak pixel + full well lines) / S/N with warning shading; returns `(fig, warnings)` |
| `raw_frame_figure(spec, t, ...)` | single raw 2-D frame in ADU: trace, airglow lines, deg-4 slit illumination, 3′ vignetting, noise |
| `reduced_2d_figure(spec, t, N, sky_fit_deg=2, ...)` | stacked, reduced, sky-subtracted 2-D S/N map (λ × arcsec) with OH residual stripes |
| `request_breakdown_figure(budget_or_list)` | horizontal stacked time-budget bars |

## 7. Calibration switch

`lrs_etc.config.USE_EMPIRICAL_CALIBRATION` (default `True`) applies the
response ratio measured from the 2026-04-01 BD+75°325 standard
(shutter-corrected: +7 s per frame) to source and sky rates. Gray factor
≈ 0.29 with a strong blue collapse below 5500 Å. Use `False` to forecast
post-refurbishment performance. The empirical grid was derived with the
Gaussian PSF model — if you adopt `psf_profile="moffat"` as standard,
re-derive the grid consistently (`config._EMP_R`).

## 8. Worked example

```python
import lrs_etc as etc

# a z = 0.3 starburst, observed r' = 20.5, E(B-V) = 0.08 foreground
s = etc.load_template("kc96_starb2")
s = etc.redshift_spectrum(s, 0.3)
s = etc.apply_extinction(s, ebv=0.08)
s = etc.normalize_filter(s, 20.5, "r")

etc.target_visibility(150.1, 2.2)              # observable?
t = etc.required_time(s, 10, 5007 * 1.3,       # S/N 10 at [O III] observed
                      slit="2.7", seeing=1.5, lunar="gray", altitude_deg=60)
budget = etc.recommended_request_time(science_s=t, night_calibration=False)
print(f"on-source {t/60:.0f} min -> request {budget['requested']/60:.0f} min")
```
