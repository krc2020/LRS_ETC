"""
Standalone figure: S/N (at 6000 A, per 4-A pixel) vs total integration time
for point sources of r'(AB) = 12...22, under dark / gray / bright moon.

Physics is taken DIRECTLY from LRS_ETC.ipynb (v6.2): this script executes
the notebook's definition cells up to the ETC core, so the figure can never
drift from the calculator. Single exposure, slow readout, -80 C, 1.8" slit,
seeing 1.1", altitude 60 deg. Saturation is ignored (idealised curves).
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Load the ETC's own functions by executing its cells ----------------
nb = json.load(open("LRS_ETC.ipynb"))
g = {"__name__": "etc_exec"}
_show = plt.show
plt.show = lambda *a, **k: plt.close("all")
for cell in nb["cells"]:
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"])
    exec(compile(src, "<etc-cell>", "exec"), g)
    if "def run_lrs_etc" in src:
        break
plt.show = _show
plt.close("all")

continuum, normalize = g["continuum"], g["normalize"]
run_lrs_etc = g["run_lrs_etc"]
NAVY, LIME, ORANGE, RED = g["NAVY"], g["LIME"], g["ORANGE"], g["RED"]

# --- Compute ------------------------------------------------------------
MAGS   = [12, 14, 16, 18, 20, 22]
PHASES = [("dark", "-"), ("gray", "--"), ("bright", ":")]
TIMES  = np.logspace(0, 4.3, 36)          # 1 s ... ~5.5 h
LAM0   = 6000.0

mag_colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(MAGS)))

snr = {}
for m in MAGS:
    spec = normalize(continuum("flat_fnu"), m, wav0=LAM0)   # flat-fnu = same AB everywhere
    for lunar, _ls in PHASES:
        vals = []
        for t in TIMES:
            r = run_lrs_etc(spec, t_per_frame=t, n_frames=1,
                            slit="1.8", seeing=1.1, lunar=lunar,
                            clouds="photometric", altitude_deg=60,
                            temperature="-80C", readout="slow")
            vals.append(np.interp(LAM0, r["wav"], r["snr_pix"]))
        snr[(m, lunar)] = np.array(vals)

# --- Plot ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.5, 7.0))
for m, colr in zip(MAGS, mag_colors):
    for lunar, ls in PHASES:
        ax.plot(TIMES, snr[(m, lunar)], ls=ls, color=colr, lw=1.6,
                alpha=0.95)
    ax.text(TIMES[-1] * 1.12, snr[(m, "dark")][-1], f"r' = {m}",
            color=colr, fontsize=10, va="center", fontweight="bold")

# sqrt(t) guide
ax.plot([30, 3000], [3 * (30 / 30) ** 0.5, 3 * (3000 / 30) ** 0.5],
        color="#999999", lw=0.9)
ax.text(700, 22, r"S/N $\propto \sqrt{t}$", color="#777777", fontsize=9,
        rotation=20)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 4e4); ax.set_ylim(0.03, 3e3)
ax.set_xlabel("Total integration time (s)", fontsize=11)
ax.set_ylabel("S/N per pixel at 6000 Å", fontsize=11)
cal = ("commissioning-calibrated (×0.29 gray + measured shape)"
       if g.get("USE_EMPIRICAL_CALIBRATION", False) else "theoretical model")
ax.set_title("LRS @ 2.4-m TNT — point source S/N vs integration time\n"
             "1.8″ slit · seeing 1.1″ · altitude 60° · photometric · "
             f"slow readout · −80 °C · {cal}", fontsize=10.5)
ax.grid(alpha=0.3, which="both")

from matplotlib.lines import Line2D
style_handles = [Line2D([0], [0], color="#333333", ls=ls, lw=1.6,
                        label=f"{lunar} night")
                 for lunar, ls in PHASES]
ax.legend(handles=style_handles, loc="upper left", fontsize=10,
          framealpha=0.95, title="lunar phase")

for hline, lab in [(3, "S/N = 3"), (10, "S/N = 10")]:
    ax.axhline(hline, color="#bbbbbb", lw=0.7, ls="-.")
    ax.text(1.15, hline * 1.1, lab, color="#999999", fontsize=8)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"LRS_SNR_vs_time.{ext}", dpi=200, bbox_inches="tight",
                facecolor="white")
print("saved LRS_SNR_vs_time.png / .pdf")
for m in MAGS:
    t10 = np.interp(10.0, snr[(m, "dark")], TIMES)
    print(f"  r'={m}: S/N=10 (dark) at t = {t10:8.1f} s")
