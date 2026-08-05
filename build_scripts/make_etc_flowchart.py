"""
Infographic-style flow chart of the LRS ETC v6 (01_LRS_ETC_v6.ipynb).
Course palette, matplotlib patches. Saves PNG + PDF.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

NAVY, LIME, ORANGE, RED = "#1E2761", "#39B54A", "#F4B400", "#E94F37"
GRAY_BG = "#F4F5F9"

fig, ax = plt.subplots(figsize=(13.5, 9.0))
ax.set_xlim(0, 13.5); ax.set_ylim(0, 9.0)
ax.axis("off")

# ----------------------------------------------------------------- helpers
def card(x, y, w, h, title, lines, color, title_fs=10.5, body_fs=8.6,
         body_color="#222222"):
    """Card with colored header band and white body."""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.03,rounding_size=0.10",
                                fc="white", ec=color, lw=1.8, zorder=2))
    head_h = 0.42
    ax.add_patch(FancyBboxPatch((x, y + h - head_h), w, head_h,
                                boxstyle="round,pad=0.03,rounding_size=0.10",
                                fc=color, ec=color, lw=1.8, zorder=3))
    ax.text(x + w / 2, y + h - head_h / 2, title, ha="center", va="center",
            fontsize=title_fs, fontweight="bold", color="white", zorder=4)
    n = len(lines)
    body_h = h - head_h - 0.14
    for i, ln in enumerate(lines):
        yy = y + body_h - (i + 0.5) * body_h / max(n, 1) + 0.05
        ax.text(x + 0.16, yy, ln, ha="left", va="center",
                fontsize=body_fs, color=body_color, zorder=4)

def step(x, y, w, h, text, color, fs=8.8, fc="white", tc="#222222", lw=1.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.03,rounding_size=0.10",
                                fc=fc, ec=color, lw=lw, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=tc, zorder=4, linespacing=1.45)

def arrow(p0, p1, color="#666666", lw=1.8, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style,
                                 mutation_scale=14, lw=lw, color=color,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1))

# ------------------------------------------------------------------ title
ax.text(6.75, 8.68, "LRS @ 2.4-m TNT  —  Exposure-Time Calculator",
        ha="center", fontsize=17, fontweight="bold", color=NAVY)
ax.text(6.75, 8.32, "long-slit spectroscopy  ·  4000–8092 Å  ·  4 Å/pixel  ·  0.9″/pixel",
        ha="center", fontsize=10, color="#555555")

# ------------------------------------------------------------ input cards
y_in, h_in = 6.10, 1.85
card(0.25, y_in, 3.15, h_in, "1 · SOURCE MODEL",
     ["Research library (49):",
      "   Pickles stars · KC96 galaxies",
      "   AGN/QSO · SDSS composites",
      "16 analytic templates",
      "User table (λ, f$_\\lambda$)",
      "Flat / power-law ± lines"], NAVY)
card(3.65, y_in, 3.05, h_in, "2 · NORM + EXTINCTION",
     ["Point: AB mag at λ$_0$",
      "Extended: AB mag / arcsec²",
      "   flux = SB × (slit × aperture)",
      "Galactic extinction:",
      "   CCM89 / O94 / F99 / Cal00",
      "   choose R$_V$; E(B−V) or A$_V$"], LIME)
card(6.95, y_in, 3.10, h_in, "3 · CONDITIONS",
     ["Seeing PSF (FWHM, ″)",
      "Moon: dark 21.8 / gray 20.9 /",
      "   bright 19.0  AB/arcsec²",
      "Clouds: photometric / cirrus /",
      "   cloudy  (0 / 0.5 / 1.2 mag)",
      "Airmass × extinction k(λ)"], ORANGE)
card(10.35, y_in, 2.90, h_in, "4 · INSTRUMENT",
     ["Slit 1.8″ / 2.7″ / 4.5″",
      "   R 750 / 500 / 300 @6000 Å",
      "BEX2-DD 256×1024 CCD",
      "   dark: 0.08 (−80 °C),",
      "           0.003 (−100 °C) e⁻/s",
      "   RN: 4 / 12 / 15 e⁻"], RED)

# ------------------------------------------------------------- core panel
yc, hc = 3.30, 2.10
ax.add_patch(FancyBboxPatch((0.25, yc), 13.0, hc,
                            boxstyle="round,pad=0.03,rounding_size=0.12",
                            fc=GRAY_BG, ec=NAVY, lw=2.2, zorder=1))
ax.text(0.55, yc + hc - 0.30, "ETC CORE", fontsize=11.5,
        fontweight="bold", color=NAVY, zorder=4)

ys, hs = yc + 0.28, 1.18
step(0.55, ys, 2.30, hs,
     "convolve to slit R\n(FWHM = 6000 Å / R)\nresample → 4 Å pixels", NAVY)
step(3.15, ys, 2.65, hs,
     "throughput η(λ)\nmirrors @ measured R(550 nm)\n× LRS optics × QE × T$_{atm}$(λ)", NAVY)
step(6.10, ys, 2.30, hs,
     "slit coupling\npoint: erf(w / 2√2 σ)\nextended: w × h area", NAVY)
step(8.70, ys, 2.30, hs,
     "e⁻ rates per pixel\nsource S(λ) · sky B(λ)\ndark n$_{pix}$D · read RN²", NAVY)
step(11.30, ys, 1.85, hs,
     "S/N =\n$\\frac{S\\,t}{\\sqrt{S t + B t + n D t + n N RN^2}}$",
     NAVY, fs=9.0, fc="white", tc=NAVY, lw=2.4)

for x0 in (2.85, 5.80, 8.40, 11.00):
    arrow((x0, ys + hs / 2), (x0 + 0.32, ys + hs / 2), color=NAVY, lw=2.0)

# arrows from inputs into the core
for xin in (1.82, 5.17, 8.50, 11.80):
    arrow((xin, y_in - 0.02), (xin, yc + hc + 0.03), color="#888888", lw=1.7)

# ---------------------------------------------------------- output cards
y_out, h_out = 0.55, 1.95
card(1.60, y_out, 4.80, h_out, "OUTPUT MODE 1  ·  forward",
     ["Give:  N frames × t per frame",
      "Get:   S/N vs wavelength",
      "          per pixel (4 Å)  and  per Å",
      "          + detector-pixel axis"], NAVY, body_fs=9.2)
card(7.10, y_out, 4.80, h_out, "OUTPUT MODE 2  ·  inverse",
     ["Give:  target S/N at wavelength λ",
      "Get:   total on-source time",
      "          closed-form quadratic of the",
      "          CCD equation → t, split into frames"], LIME, body_fs=9.2)

arrow((4.60, yc - 0.02), (4.00, y_out + h_out + 0.03), color=NAVY, lw=2.0, rad=-0.08)
arrow((8.90, yc - 0.02), (9.50, y_out + h_out + 0.03), color=LIME, lw=2.0, rad=0.08)

ax.text(6.75, 0.16,
        "NARIT LRS Exposure-Time Calculator  ·  LRS_ETC.ipynb  ·  "
        "2.4-m Thai National Telescope",
        ha="center", fontsize=8, color="#888888")

plt.tight_layout(pad=0.4)
for ext in ("png", "pdf"):
    fig.savefig(f"ETC_flowchart.{ext}", dpi=200, bbox_inches="tight",
                facecolor="white")
print("saved ETC_flowchart.png / .pdf")
