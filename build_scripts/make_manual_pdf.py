"""
MONET Calc — user manual PDF (brief), branded, ~4 pages.
"""

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, PageBreak)

NAVY = colors.HexColor("#1E2761")
LIME = colors.HexColor("#39B54A")
ORANGE = colors.HexColor("#F4B400")
RED = colors.HexColor("#E94F37")
LIGHT = colors.HexColor("#EEF0F7")

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], textColor=NAVY, fontSize=16,
                    spaceAfter=6)
H2 = ParagraphStyle("H2", parent=ss["Heading2"], textColor=NAVY, fontSize=12,
                    spaceBefore=10, spaceAfter=4)
BODY = ParagraphStyle("BODY", parent=ss["Normal"], fontSize=9.2, leading=12.5)
SMALL = ParagraphStyle("SMALL", parent=ss["Normal"], fontSize=8.2, leading=11,
                       textColor=colors.HexColor("#444444"))
CELL = ParagraphStyle("CELL", parent=ss["Normal"], fontSize=8.4, leading=11)
CELLB = ParagraphStyle("CELLB", parent=CELL, fontName="Helvetica-Bold",
                       textColor=colors.white)

def T(data, widths, header=True):
    rows = [[Paragraph(c, CELLB if (header and i == 0) else CELL)
             for c in row] for i, row in enumerate(data)]
    t = Table(rows, colWidths=widths, repeatRows=1 if header else 0)
    style = [("VALIGN", (0, 0), (-1, -1), "TOP"),
             ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#BBBBBB")),
             ("LEFTPADDING", (0, 0), (-1, -1), 4),
             ("RIGHTPADDING", (0, 0), (-1, -1), 4),
             ("TOPPADDING", (0, 0), (-1, -1), 3),
             ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
    if header:
        style += [("BACKGROUND", (0, 0), (-1, 0), NAVY),
                  ("TEXTCOLOR", (0, 0), (-1, 0), colors.white)]
    for r in range(2, len(data), 2):
        style.append(("BACKGROUND", (0, r), (-1, r), LIGHT))
    t.setStyle(TableStyle(style))
    return t

doc = SimpleDocTemplate("MONET_Calc_manual.pdf", pagesize=letter,
                        leftMargin=0.75*inch, rightMargin=0.75*inch,
                        topMargin=0.6*inch, bottomMargin=0.6*inch,
                        title="MONET Calc User Manual",
                        author="Krittapas Chanchaiworawit / NARIT")
S = []

# ---------------- header -------------------------------------------------
logo = Image("spectral_library/../assets_logo.png") if False else None
try:
    logo = Image("monet_logo_local.png", width=1.15*inch, height=1.15*inch)
    hdr = Table([[logo,
                  [Paragraph("MONET Calc", H1),
                   Paragraph("Module of Numerical Exposure Time Calculator — "
                             "NARIT LRS @ 2.4-m Thai National Telescope",
                             BODY),
                   Paragraph("User manual v6.5 · "
                             "https://monet-calc.streamlit.app · "
                             "krittapas [at] narit.or.th", SMALL)]]],
                colWidths=[1.3*inch, 5.7*inch])
    hdr.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE")]))
    S.append(hdr)
except Exception:
    S.append(Paragraph("MONET Calc — User Manual v6.5", H1))
S.append(Spacer(1, 8))

# ---------------- quick start -------------------------------------------
S.append(Paragraph("Quick start (4 steps)", H2))
S.append(Paragraph(
    "<b>1.</b> In the <b>sidebar</b>, define your target: pick a spectral "
    "template (or upload a 2-column ASCII table of wavelength [Å] and "
    "f<sub>λ</sub>, or build a continuum ± line), set redshift and "
    "Galactic E(B−V), then set the observed brightness (SDSS filter "
    "magnitude, monochromatic AB, or total flux). "
    "<b>2.</b> Set conditions and instrument: seeing, moon, clouds, target "
    "altitude, slit, CCD temperature, readout. "
    "<b>3.</b> Pick one of the three tabs and press its button. "
    "<b>4.</b> Read the plots and the color-coded warnings.", BODY))
S.append(Spacer(1, 4))
S.append(Paragraph(
    "<b>Calibration toggle (important):</b> ON (default) scales predictions "
    "to match the measured 2026-04-01 commissioning performance of the "
    "as-built system (×0.29 with a measured color-dependent shape; "
    "conservative, realistic <i>now</i>). OFF gives the theoretical "
    "component model (use for forecasts after optics refurbishment).", BODY))

# ---------------- three tabs --------------------------------------------
S.append(Paragraph("The three calculators", H2))
S.append(T([
    ["Tab", "Give", "Get"],
    ["① S/N ↔ exposure time",
     "N frames × exposure per frame — or — a target S/N at a chosen "
     "wavelength",
     "S/N vs wavelength (per 4-Å pixel and per Å), the expected "
     "spectrum in flux / ADU / S/N with warning shading — or — the total "
     "on-source time and per-frame split"],
    ["② Simulated 1-D & 2-D spectra",
     "N frames × exposure per frame",
     "Expected 1-D spectrum (flux, ADU incl. peak pixel vs full well, S/N); "
     "raw 2-D frame in ADU (trace, airglow lines, slit vignetting, noise); "
     "reduced, stacked, sky-subtracted 2-D S/N map (λ × arcsec) "
     "with realistic OH-band residuals; extended sources appear as a "
     "band of the extent you set, with sky fitted outside it"],
    ["③ Recommended request time",
     "N × t and whether nighttime calibration is needed",
     "The time to REQUEST in the proposal: on-source + readout, slew, "
     "acquisition, thru-slit, operator, standard-star rounds (every 2 h), "
     "all ÷0.65 for the 35% weather margin; night calibration (1 h per "
     "half night) added outside the margin; daytime calibration shown as "
     "info, never charged"],
], [1.45*inch, 1.85*inch, 3.7*inch]))
S.append(Spacer(1, 6))
S.append(Paragraph("Warnings (shaded bands on the plots)", H2))
S.append(T([
    ["Color", "Meaning", "What to do"],
    ["Red — saturation", "peak pixel ≥ 65,535 ADU in one frame",
     "shorten frames, take more of them"],
    ["Orange — non-linear", "peak pixel above 95% of full well",
     "shorten frames slightly"],
    ["Purple — read-noise-dominated", "RN² is the largest noise term",
     "longer frames and/or slow readout"],
    ["Brown — dark-dominated", "dark current is the largest noise term",
     "use the −100 °C operating point"],
], [1.7*inch, 2.8*inch, 2.5*inch]))
S.append(PageBreak())

# ---------------- template guide ----------------------------------------
S.append(Paragraph("Template guide — which template for which object", H2))
S.append(Paragraph(
    "Templates are research-grade spectra from the standard ETC libraries "
    "(also used by ESO/Gemini/HST calculators). Names are "
    "<b>family_item</b>. All are renormalized to your requested brightness; "
    "only the <i>shape</i> matters. If a template does not cover part of "
    "4000–8092 Å the app reports it and S/N there is undefined "
    "rather than extrapolated.", BODY))
S.append(Spacer(1, 4))

S.append(Paragraph("Stars — <b>pickles_*</b> (Pickles 1998 library)", H2))
S.append(T([
    ["Template(s)", "Object type"],
    ["pickles_o5v, b0v, b57v", "hot O/B main-sequence stars (blue, "
     "Balmer + He lines); flux standards, young clusters"],
    ["pickles_a0v, a5v", "A dwarfs (strong Balmer absorption); Vega-like "
     "standards"],
    ["pickles_f0v, f5v, g0v, g2v, g5v", "F/G dwarfs; g2v = solar analog"],
    ["pickles_k0v, k2v, k5v, m0v, m2v, m5v", "K/M dwarfs (red, molecular "
     "bands); cool stars, exoplanet hosts"],
    ["pickles_g8iii, k3iii, m0iii, m5iii", "G/K/M giants; evolved stars, "
     "stellar populations of early-type galaxies"],
    ["pickles_a0i, g2i, m2i", "supergiants; luminous evolved massive stars"],
], [2.5*inch, 4.5*inch]))
S.append(Spacer(1, 5))

S.append(Paragraph("Galaxies — <b>kc96_*</b> (Kinney–Calzetti atlas)",
                   H2))
S.append(T([
    ["Template(s)", "Object type"],
    ["kc96_elliptical", "quenched / passive elliptical: red continuum, "
     "4000-Å break, absorption lines only"],
    ["kc96_bulge, kc96_s0", "bulge-dominated and lenticular galaxies "
     "(note: red coverage ends ~7550 Å; see brown_* substitutes)"],
    ["kc96_sa, kc96_sb, kc96_sc", "spirals from early (Sa) to late (Sc) "
     "type: growing blue continuum and emission lines"],
    ["kc96_starb1 … starb6", "starburst galaxies ordered by internal "
     "dust: starb1 = blue, nearly unobscured (young irregular-like); "
     "starb5/starb6 = dusty starbursts, E(B−V) up to ~0.7"],
], [2.5*inch, 4.5*inch]))
S.append(Spacer(1, 5))

S.append(Paragraph("Active nuclei — <b>agn_*</b>", H2))
S.append(T([
    ["Template(s)", "Object type"],
    ["agn_qso_ext  (recommended)", "luminous unobscured quasar; Francis+91 "
     "composite spliced with Türler+99 — covers the full band at "
     "any redshift 0–4"],
    ["agn_qso", "same quasar composite, native version: rest coverage ends "
     "at 6000 Å — use only when redshift brings the rest-UV into "
     "the band (z ≈ 0.33–4)"],
    ["agn_seyfert1", "broad-line (type 1) Seyfert nucleus + host"],
    ["agn_seyfert2, agn_ngc1068", "obscured (type 2) AGN, narrow lines "
     "only; ngc1068 is the archetype with full-band coverage"],
    ["agn_liner", "low-ionization nuclear emission (LINER)"],
], [2.5*inch, 4.5*inch]))
S.append(Spacer(1, 5))

S.append(Paragraph("SDSS composites — <b>dobos_*</b> (Dobos+2012) "
                   "and real-galaxy SEDs — <b>brown_*</b> (Brown+2014)", H2))
S.append(T([
    ["Template(s)", "Object type"],
    ["dobos_SF1 … SF4", "average SDSS star-forming galaxies, "
     "increasing star-formation activity; best for typical low-z emission-"
     "line galaxies (full band only to z ≈ 0.13)"],
    ["dobos_RED0, RED2, RED4", "average SDSS passive (red-sequence) "
     "galaxies"],
    ["dobos_BG, dobos_RG", "blue-cloud and red-sequence population means"],
    ["brown_NGC5033_sy1.5", "real Seyfert 1.5 galaxy — full-band "
     "substitute for agn_seyfert1"],
    ["brown_NGC4579_liner", "real LINER — full-band substitute for "
     "agn_liner"],
    ["brown_NGC3379_elliptical", "real elliptical — full-band "
     "substitute for kc96_elliptical/bulge"],
    ["brown_NGC4450_sab", "real early-type (Sab) spiral"],
    ["brown_NGC0628_sc", "real late-type (Sc) spiral — full-band "
     "substitute for kc96_sc"],
], [2.5*inch, 4.5*inch]))
S.append(PageBreak())

# ---------------- practical notes ---------------------------------------
S.append(Paragraph("Practical notes", H2))
S.append(Paragraph(
    "<b>Order of operations.</b> The app applies your settings as: "
    "redshift → Galactic extinction → normalization. The "
    "magnitude you enter is therefore the <i>observed</i> brightness — "
    "what a catalog would list.", BODY))
S.append(Spacer(1, 3))
S.append(Paragraph(
    "<b>Redshift.</b> Accepted to z = 9. Where the redshifted template has "
    "no rest-frame data (e.g., no X-ray/extreme-UV), the spectrum is padded "
    "with zeros and the app tells you which rest-frame regime is missing; "
    "S/N there is zero, and normalization in a padded region is refused. "
    "Kinney–Calzetti templates cover the full band to z ≈ 2.2; "
    "the extended QSO to z ≈ 4; SDSS composites only to z ≈ 0.13.",
    BODY))
S.append(Spacer(1, 3))
S.append(Paragraph(
    "<b>Extended sources.</b> Check “Extended source” to interpret "
    "brightness per arcsec² (uniform surface brightness; no slit "
    "losses; flux scales with slit × extraction area). Also set "
    "“Object extent along the slit” (new in v6.5): the extraction "
    "window matches the extent, the 2-D maps show the object as a "
    "seeing-blurred band of that size, and the sky is fitted only "
    "outside it — so the reduced map keeps your object instead of "
    "subtracting it as sky. If the object fills the whole slit, the "
    "app warns you and shows the map without sky subtraction: plan "
    "offset-sky exposures for such targets. For a compact "
    "but resolved object, instead keep a total magnitude and set "
    "“Source FWHM” — the profile widens to "
    "√(seeing² + source²), increasing slit losses "
    "realistically.", BODY))
S.append(Spacer(1, 3))
S.append(Paragraph(
    "<b>Slit choice.</b> 1.8″ gives R ≈ 750, 2.7″ → 500, "
    "4.5″ → 300 at 6000 Å. Under poor image quality "
    "(seeing worse than 2.5″) the wider slits recover large slit losses at modest "
    "resolution cost. The “moffat” PSF option models realistic "
    "seeing wings (10–25% more slit loss than Gaussian).", BODY))
S.append(Spacer(1, 3))
S.append(Paragraph(
    "<b>What to expect (calibrated mode, 1.8″ slit, 1.5″ seeing, "
    "dark sky, 30 min).</b> Point source r′ ≈ 18: S/N ≈ 30 "
    "per pixel; r′ ≈ 20: S/N ≈ 7; practical continuum limit "
    "r′ ≈ 21.5–22 for long stacks. Divide exposure times by "
    "~3 when forecasting with the theoretical model instead.", BODY))
S.append(Spacer(1, 3))
S.append(Paragraph(
    "<b>Requesting time.</b> Always quote tab ③’s "
    "“REQUEST” number in proposals — it includes all "
    "overheads, standard-star visits, and the 35% weather/mechanical "
    "margin. Daytime calibrations are free; nighttime calibration costs a "
    "fixed 1 h per half night and must be stated.", BODY))
S.append(Spacer(1, 10))
S.append(Paragraph(
    "MONET Calc v6.5 · NARIT · source: github.com/krc2020/LRS_ETC "
    "· The empirical calibration is anchored to the 2026-04-01 "
    "BD+75°325 commissioning standard (shutter-corrected). Template "
    "provenance: Pickles 1998; Kinney et al. 1996; Calzetti et al. 1994; "
    "Francis et al. 1991; Türler et al. 1999; Dobos et al. 2012; "
    "Brown et al. 2014.", SMALL))

doc.build(S)
print("wrote MONET_Calc_manual.pdf")
