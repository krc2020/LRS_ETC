# TNT + LRS System Throughput Model — assumptions

Wavelength grid: **400 – 800 nm** at **1 nm** spacing (401 points).

The full system throughput is the product of seven element groups:

η_sys(λ) = R_tel(λ)⁴ × T_L1(λ) × R_fold(λ)² × T_coll(λ) × η_grating(λ) × T_cam(λ) × QE_CCD(λ)

> **v6.3 note (LRS_ETC repo):** the four TNT mirrors were **measured
> individually in 2026**; the ETC rescales each fresh-coating curve to the
> measured reflectivity at 550 nm — M1 90.4 %, M2 90.8 %, M3 93.0 %,
> M4 87.9 % (four-mirror product 67.1 % at 550 nm), lowering the 600 nm
> total from 0.52 to 0.36. See notebook §2. The §11 validation suggests
> the as-built system still delivers less than this measured-mirror model
> — see the open-issue note in `CLAUDE.md` before quoting absolute
> sensitivities.

## Per-component assumptions

| Component | Spec / model |
|---|---|
| **TNT M1·M2·M3·M4 mirrors** | Four reflective surfaces (primary, secondary, tertiary, and cube M4 fold) with a broadband-visible **dielectric coating peaking at 600 nm**: R_single ≈ 99.6 % at 600 nm, ~98 % at 500 & 700 nm, ~93 % at the band edges (400 and 800 nm). Net: R_single⁴, which doubly penalises off-peak wavelengths. |
| **Focal reducer L1** | Front-half doublet of the TNT focal reducer (Prasit et al. 2019, SPIE 11116). Two air-glass surfaces with a broadband VIS AR coating designed for 400 – 800 nm — modeled with per-surface R ≈ 0.3-1.0 % across the band. |
| **LRS fold mirrors (× 2)** | 3-inch flats with Thorlabs broadband dielectric (E02-class) coating; R_avg ≈ 99.5 % from 420 – 750 nm, dropping to ≈ 97 % at 800 nm and ≈ 98.5 % at 400 nm. |
| **LRS collimator (250 mm)** | Thorlabs achromatic doublet with `-A` AR coating (400 – 700 nm). Per-surface R from spec: < 0.5 % in 450 – 700 nm, rising to ≈ 5 % at 800 nm. Two air-glass surfaces. |
| **VPH grating** | Peak η = 90 % at 650 nm; η ≈ 40 % at 450 nm and 800 nm. Asymmetric Gaussian (σ_blue ≈ 175 nm, σ_red ≈ 132 nm). |
| **Nikon camera lens (58 mm)** | Eight air-glass surfaces with **no AR coating**. Fresnel loss ≈ 4 % per surface (n ≈ 1.515) → T ≈ (0.96)⁸ ≈ 72 %, ~λ-independent. |
| **CCD QE** | Andor Newton **BEX2-DD** back-illuminated deep-depletion CCD; anchors read directly from Andor's QE chart. Nearly flat 88 – 94 % across 400 – 820 nm, shallow dip near 540 nm, peak ~94 % around 750 – 800 nm. |

## Headline numbers (full TNT + L1 + LRS — fresh coatings, pre-ageing)

| Quantity | LRS-only | Full system (× telescope × L1) |
|---|---|---|
| Peak throughput  | **57.2%** @ 655 nm | **54.9%** @ 649 nm |
| η at 450 nm      | 25.5% | 21.8% |
| η at 550 nm      | 45.2% | 43.7% |
| η at 650 nm      | 57.2% | 54.9% |
| η at 750 nm      | 39.0% | 32.2% |
| η at 800 nm      | 22.3% | 16.3% |
| Mean 400 – 800 nm | 41.1% | 37.6% |

## Sources / vendor pages

- **Phetra+ 2016 (Optics Express 24, 1416)** — TNT focal reducer optical design.
  https://opg.optica.org/oe/fulltext.cfm?uri=oe-24-2-1416
- **Prasit+ 2019 (SPIE 11116, 111161A)** — TNT focal reducer installation & test.
  https://spie.org/Publications/Proceedings/Paper/10.1117/12.2529833
- **NARIT TNT factsheet** — 2.4 m Ritchey-Chrétien, f/10 final.
  https://www.narit.or.th/index.php/thai-national-observatory-tno-location-menu-2
- **Andor Newton BEX2-DD QE chart** — Andor / Oxford Instruments datasheet.
- **Thorlabs broadband dielectric (E02) coating** — BB1-E02 and 2-inch variants.
- **Thorlabs achromatic doublet, -A AR coating** — AC254-250-A reflectance graphs.
- **VPH grating profile** — user-supplied anchors (90 % @ 650 nm; 40 % @ 450 & 800 nm).

> Replace any anchor array at the top of `make_lrs_throughput.py` with measured
> numbers as they become available; the rest of the model and the outputs
> regenerate in seconds.
