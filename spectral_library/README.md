# Research-grade spectral library for the LRS ETC

Plain 2-column ASCII (wavelength Å, f_λ erg s⁻¹ cm⁻² Å⁻¹; absolute scale
arbitrary — the ETC renormalizes). Trimmed to 3000–10 500 Å. Retrieved via
the `spextra` database (Leschinski 2021), which packages the standard
template sets used by the ESO, Gemini and HST exposure-time calculators.

| prefix | source | reference |
|---|---|---|
| `pickles_*` | Stellar flux library, 24 types: O5V–M5V main sequence, G8III–M5III giants, A0I/G2I/M2I supergiants (`b57v` = B5–7V composite) | Pickles 1998, PASP 110, 863 |
| `kc96_*` | Galaxy templates: elliptical, bulge, S0, Sa, Sb, Sc + starburst series `starb1–6` ordered by internal reddening E(B−V) ≈ 0.05 → 0.7 | Kinney et al. 1996, ApJ 467, 38; Calzetti et al. 1994, ApJ 429, 582 |
| `agn_qso` | QSO composite | Francis et al. 1991, ApJ 373, 465 |
| `agn_seyfert1/2`, `agn_liner`, `agn_ngc1068` | AGN class templates | spextra AGN set (STScI CDBS heritage) |
| `dobos_*` | SDSS DR7 composite galaxy spectra: SF1–SF4 (star-forming sequence), RED0/2/4 (passive sequence), BG/RG (blue/red cloud) | Dobos et al. 2012, MNRAS 420, 1217 |
| `agn_qso_ext` | ESO ETC spliced QSO composite: Francis+91 (800–6000 Å) ⊕ Türler+99 3C 273 SED (6000 Å–2 µm) — **use this instead of `agn_qso` for full LRS coverage** | Francis et al. 1991; Türler et al. 1999, A&AS 134, 89 |
| `brown_*` | Full-coverage (0.09–30 µm) real-galaxy SEDs replacing red-truncated CDBS/KC96 classes: NGC 5033 (Sy 1.5), NGC 4579 (LINER), NGC 3379 (E1), NGC 4450 (Sab), NGC 628 (Sc) | Brown et al. 2014, ApJS 212, 18 |

Retrieved 2026-07-23 from `https://scopesim.univie.ac.at/spextra/database/`.

If you add files, keep the naming pattern `<library>_<item>.txt` — the ETC
notebook scans this folder and labels templates from the filename.
