"""
lrs_etc.plots — figure builders. Every function RETURNS a matplotlib
Figure (no plt.show), so they drop straight into Streamlit via st.pyplot.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

from .config import (NAVY, LIME, ORANGE, RED, DISP_AA_PIX, SPATIAL_AS_PIX,
                     GAIN_E_ADU, FULL_WELL_ADU, BIAS_ADU, LINEARITY_FRAC,
                     DARK_E_PIX_S, READ_NOISE_E, SLIT_LENGTH_ARCSEC)
from .core import run_lrs_etc, check_warnings

WARN_STYLE = [("saturated", RED, "saturated"),
              ("nonlinear", ORANGE, "non-linear (>95% FW)"),
              ("rn_dominated", "#8064A2", "read-noise-dominated"),
              ("dark_dominated", "#8B5A2B", "dark-dominated")]

def _shade(ax, wav, mask, color, label=None):
    if not mask.any():
        return
    idx = np.where(mask)[0]
    for j, grp in enumerate(np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)):
        ax.axvspan(wav[grp[0]], wav[grp[-1]], color=color, alpha=0.20,
                   lw=0, label=label if j == 0 else None)

def snr_figure(results, labels=None):
    """S/N vs wavelength (per pixel + per A) for one or more ETC results."""
    if isinstance(results, dict):
        results = [results]
    labels = labels or [None] * len(results)
    fig, axs = plt.subplots(2, 1, figsize=(10.5, 6), sharex=True)
    colors = [NAVY, LIME, ORANGE, RED, "#8064A2", "#4F81BD"]
    for res, lab, colr in zip(results, labels, colors):
        cfg = res["config"]
        lab = lab or (f"{cfg['slit']}\" slit, {cfg['n_frames']}x"
                      f"{cfg['t_per_frame']:.0f} s, {cfg['lunar']}")
        axs[0].plot(res["wav"], res["snr_pix"], lw=1.0, color=colr, label=lab)
        axs[1].plot(res["wav"], res["snr_aa"], lw=1.0, color=colr)
    axs[0].set_ylabel("S/N per pixel (4 Å)")
    axs[1].set_ylabel("S/N per Å")
    axs[1].set_xlabel("Wavelength (Å)")
    axs[0].legend(fontsize=8)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig

def expected_spectrum_figure(res, t_per_frame, spec_name=""):
    """3-panel flux / ADU / S/N view with warning shading. Returns
    (fig, warnings_dict)."""
    wrn = check_warnings(res, t_per_frame)
    cfg = res["config"]
    fig, axs = plt.subplots(3, 1, figsize=(10.5, 8), sharex=True)
    axs[0].plot(res["wav"], res["flam_pix"], color=NAVY, lw=0.9)
    axs[0].set_yscale("log")
    axs[0].set_ylabel(r"$f_\lambda$ (erg/s/cm²/Å)")
    axs[0].set_title(f"{spec_name}  ·  {cfg['n_frames']}x"
                     f"{cfg['t_per_frame']:.0f} s, {cfg['slit']}\" slit, "
                     f"extraction {cfg['extract_arcsec']:.1f}\"  "
                     f"[{cfg['calibration']}]", fontsize=10)
    tot_adu = ((res["S"] + res["B"] + res["dark_rate"]) * t_per_frame
               / GAIN_E_ADU + res["n_spat"] * BIAS_ADU)
    axs[1].plot(res["wav"], wrn["peak_adu"], color=NAVY, lw=0.9,
                label="peak pixel / frame")
    axs[1].plot(res["wav"], tot_adu, color="#4F81BD", lw=0.8, ls="--",
                label="total extracted / frame")
    axs[1].axhline(FULL_WELL_ADU, color=RED, lw=1.0)
    axs[1].axhline(LINEARITY_FRAC * FULL_WELL_ADU, color=ORANGE, lw=0.9,
                   ls=":")
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
    h, l = axs[0].get_legend_handles_labels()
    if h:
        axs[0].legend(h, l, fontsize=8, loc="lower right")
    fig.tight_layout()
    return fig, wrn

# ------------------------- 2-D simulators -------------------------------

def _row_fractions(y, fwhm):
    sig = fwhm / 2.3548
    lo = (y - SPATIAL_AS_PIX / 2) / (np.sqrt(2) * sig)
    hi = (y + SPATIAL_AS_PIX / 2) / (np.sqrt(2) * sig)
    return 0.5 * (erf(hi) - erf(lo))

def _slit_illum_poly4(y):
    u = y / (SLIT_LENGTH_ARCSEC / 2)
    return 1.0 + 0.05*u - 0.08*u**2 + 0.04*u**3 + 0.10*u**4

def _vignette(y, roll=2.0):
    edge = SLIT_LENGTH_ARCSEC / 2
    return 0.5 * (1 + erf((edge - np.abs(y)) / (np.sqrt(2) * roll)))

def _extent_profile(y, extent_arcsec, seeing):
    """Seeing-smoothed top-hat of the given along-slit extent: ~1 inside
    the object, rolling off over the seeing scale at its edges."""
    sig = max(seeing, 0.3) / 2.3548
    half = extent_arcsec / 2.0
    return 0.5 * (erf((half - y) / (np.sqrt(2) * sig))
                  + erf((half + y) / (np.sqrt(2) * sig)))

def _sim_rates(spec, t_per_frame, n_frames, n_rows, seeing, **etc_kwargs):
    res = run_lrs_etc(spec, t_per_frame, n_frames, seeing=seeing,
                      extract_arcsec=n_rows * SPATIAL_AS_PIX, **etc_kwargs)
    cfg = res["config"]
    fwhm = cfg.get("fwhm_eff", seeing)
    y = (np.arange(n_rows) - n_rows / 2 + 0.5) * SPATIAL_AS_PIX
    vig = _vignette(y)
    if spec.get("extended", False):
        extent = cfg.get("source_extent_arcsec")
        if extent:
            # finite object: distribute its total through-slit flux over a
            # seeing-smoothed top-hat of the given extent
            prof = _extent_profile(y, extent, seeing) * vig
            norm = max(prof.sum(), 1e-9)
            src = res["S"][None, :] * (prof / norm)[:, None]
        else:
            # legacy behavior: uniform SB filling the whole illuminated slit
            src = ((res["S"][None, :] / n_rows) * np.ones((n_rows, 1))
                   * vig[:, None])
    else:
        frac = _row_fractions(y, fwhm)
        sig = fwhm / 2.3548
        fy = erf((cfg["extract_arcsec"] / 2) / (np.sqrt(2) * sig))
        src = ((res["S"] / max(fy, 1e-9))[None, :] * frac[:, None]
               * vig[:, None])
    sky = ((res["B"] / res["n_spat"])[None, :] * _slit_illum_poly4(y)[:, None]
           * vig[:, None])
    dark = DARK_E_PIX_S[cfg["temperature"]]
    rn = READ_NOISE_E[cfg["readout"]]
    return res, y, src, sky, dark, rn, vig

def _mark_vignette(ax, y):
    edge = SLIT_LENGTH_ARCSEC / 2
    if y.min() < -edge and y.max() > edge:
        for s in (-1, 1):
            ax.axhline(s * edge, color="white", lw=0.9, ls="--", alpha=0.85)

def raw_frame_figure(spec, t_per_frame, n_rows=256, seeing=1.1, seed=42,
                     **etc_kwargs):
    """Raw single-frame 2-D spectrum in ADU (heat map). Returns fig."""
    res, y, src, sky, dark, rn, vig = _sim_rates(spec, t_per_frame, 1,
                                                 n_rows, seeing, **etc_kwargs)
    rng = np.random.default_rng(seed)
    e = rng.poisson(np.clip((src + sky + dark) * t_per_frame, 0, None)
                    ).astype(float)
    e += rng.normal(0, rn, e.shape)
    adu = np.clip(e / GAIN_E_ADU + BIAS_ADU, 0, FULL_WELL_ADU)
    fig, ax = plt.subplots(figsize=(11, 4.2))
    im = ax.imshow(adu, origin="lower", aspect="auto", cmap="magma",
                   vmin=BIAS_ADU * 0.97, vmax=np.percentile(adu, 99.3),
                   extent=[res["wav"][0], res["wav"][-1], y[0], y[-1]])
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("Along slit (arcsec)")
    cfg = res["config"]
    ax.set_title(f"Raw frame — {t_per_frame:.0f} s, {cfg['slit']}\" slit, "
                 f"{cfg['lunar']} sky (ADU)", fontsize=10)
    _mark_vignette(ax, y)
    fig.colorbar(im, ax=ax, label="ADU", pad=0.01)
    fig.tight_layout()
    return fig

def reduced_2d_figure(spec, t_per_frame, n_frames, n_rows=256, seeing=1.1,
                      sky_fit_deg=2, seed=42, **etc_kwargs):
    """Reduced + stacked + sky-subtracted 2-D S/N map. Returns fig."""
    res, y, src, sky, dark, rn, vig = _sim_rates(spec, t_per_frame, n_frames,
                                                 n_rows, seeing, **etc_kwargs)
    rng = np.random.default_rng(seed)
    t_tot = t_per_frame * n_frames
    src_e, sky_e = src * t_tot, sky * t_tot
    var = src_e + sky_e + dark * t_tot + n_frames * rn ** 2
    noisy = src_e + sky_e + rng.normal(0, np.sqrt(var))
    cfg = res["config"]
    fwhm = cfg.get("fwhm_eff", seeing)
    # Sky-fit rows must EXCLUDE the object. For an extended source of
    # known extent, mask its extent plus a seeing margin; a point/compact
    # source is masked over 3x the effective FWHM. If the object (or an
    # extended source with no stated extent) leaves too few clean rows,
    # on-slit sky subtraction is impossible - warn and skip it, exactly
    # as at the telescope you would need offset-sky exposures.
    if spec.get("extended", False):
        extent = cfg.get("source_extent_arcsec")
        half_excl = (extent / 2 + max(seeing, 2.0)) if extent else np.inf
    else:
        half_excl = 3 * fwhm
    fit_rows = (np.abs(y) >= half_excl) & (vig > 0.9)
    illum = vig > 0.5
    yy = y / (SLIT_LENGTH_ARCSEC / 2)
    sky_est = np.zeros_like(noisy)
    sky_subtracted = fit_rows.sum() >= 8
    if sky_subtracted:
        for c in range(noisy.shape[1]):
            p = np.polyfit(yy[fit_rows], noisy[fit_rows, c], sky_fit_deg)
            sky_est[illum, c] = np.polyval(p, yy[illum])
    else:
        print("WARNING: the object occupies (nearly) the whole slit - no "
              "clean sky rows for on-slit subtraction. Showing the map "
              "WITHOUT sky subtraction; plan offset-sky exposures for "
              "such targets (or set source_extent_arcsec).")
    snr_map = (noisy - sky_est) / np.sqrt(var)
    fig, ax = plt.subplots(figsize=(11, 4.2))
    im = ax.imshow(snr_map, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-8, vmax=8,
                   extent=[res["wav"][0], res["wav"][-1], y[0], y[-1]])
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("Along slit (arcsec)")
    ax.set_title(f"Reduced 2-D S/N — {n_frames}x{t_per_frame:.0f} s, "
                 + (f"sky fit deg={sky_fit_deg} vs true deg-4"
                    if sky_subtracted else
                    "NO sky subtraction (object fills the slit)"),
                 fontsize=10)
    _mark_vignette(ax, y)
    fig.colorbar(im, ax=ax, label="S/N per pixel", pad=0.01)
    fig.tight_layout()
    return fig

# ------------------------- request-time bars ----------------------------

def request_breakdown_figure(results):
    """Horizontal stacked bars for recommended_request_time() results."""
    if isinstance(results, dict):
        results = [results]
    comps = ["science", "readout", "slew", "acquisition", "thruslit",
             "operator", "standards", "weather", "night_cal"]
    names = ["on-source science", "CCD readout", "slew",
             "offset + acquisition", "thru-slit", "operator",
             "standard star", "weather/mech. (35%)", "night cal (fixed)"]
    colors = [NAVY, "#4F81BD", LIME, ORANGE, "#8064A2", "#999999",
               "#2AA198", RED, "#444444"]
    hatches = [None] * 8 + ["//"]
    fig, ax = plt.subplots(figsize=(10.5, 1.2 + 0.9 * len(results)))
    for j, res in enumerate(results):
        left = 0.0
        for comp, colr, htc in zip(comps, colors, hatches):
            v = res.get(comp, 0.0) / 60.0
            if v <= 0:
                continue
            ax.barh(j, v, left=left, color=colr, edgecolor="white",
                    height=0.55, lw=0.6, hatch=htc)
            if v > res["requested"] / 60.0 * 0.05:
                ax.text(left + v / 2, j, f"{v:.0f}", ha="center",
                        va="center", color="white", fontsize=8,
                        fontweight="bold")
            left += v
        ax.text(left + res["requested"] / 60.0 * 0.012, j,
                f"{res['requested']/60:.0f} min (η={res['efficiency']:.0%})",
                va="center", fontsize=9)
    ax.set_yticks(range(len(results)))
    ax.set_yticklabels([r["label"] for r in results], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Time (minutes)")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=c, hatch=h,
                                     edgecolor="white")
                       for c, h in zip(colors, hatches)],
              labels=names, fontsize=7.5, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.30))
    ax.grid(alpha=0.25, axis="x")
    ax.set_xlim(0, max(r["requested"] for r in results) / 60.0 * 1.28)
    fig.tight_layout()
    return fig
