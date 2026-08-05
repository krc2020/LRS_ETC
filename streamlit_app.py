"""
MONET Calc (Module of Numerical Exposure Time Calculator) —
Streamlit front-end for the NARIT LRS ETC.

Run locally:   streamlit run streamlit_app.py
Deploy:        push the repo to GitHub -> share.streamlit.io -> pick this file.
"""

from pathlib import Path

import numpy as np
import streamlit as st

import lrs_etc as etc
from lrs_etc import config

_LOGO = Path(__file__).parent / "assets" / "monet_calc_logo.png"

st.set_page_config(page_title="MONET Calc — NARIT LRS ETC",
                   page_icon=str(_LOGO) if _LOGO.exists() else "🔭",
                   layout="wide")

NAVY = "#1E2761"
c_logo, c_title = st.columns([1, 5])
if _LOGO.exists():
    c_logo.image(str(_LOGO), width=130)
c_title.markdown(
    f"<h1 style='color:{NAVY};margin-bottom:0'>MONET Calc</h1>"
    f"<h4 style='color:{NAVY};margin-top:0'>Module of Numerical Exposure "
    f"Time Calculator — LRS @ 2.4-m TNT</h4>", unsafe_allow_html=True)
c_title.caption(f"v{etc.__version__} · long-slit 4000–8092 Å · 4 Å/pixel · "
                f"empirical commissioning calibration available")

# ======================= sidebar: configuration =========================
with st.sidebar:
    if _LOGO.exists():
        st.image(str(_LOGO), width=110)
    st.header("Target")
    src_mode = st.radio("Spectrum source",
                        ["Library template", "Upload 2-column table",
                         "Continuum ± line builder"])
    if src_mode == "Library template":
        tmpl = st.selectbox("Template", etc.list_templates(),
                            index=etc.list_templates().index("kc96_sb")
                            if "kc96_sb" in etc.list_templates() else 0)
        spec0 = etc.load_template(tmpl, verbose=False)
    elif src_mode == "Upload 2-column table":
        up = st.file_uploader("ASCII: wavelength_Å  f_λ", type=["txt", "dat"])
        if up is None:
            st.stop()
        t = np.loadtxt(up)
        spec0 = dict(wav=etc.WAV.copy(),
                     flam=np.interp(etc.WAV, t[:, 0], t[:, 1],
                                    left=np.nan, right=np.nan),
                     wav_native=t[:, 0], flam_native=t[:, 1],
                     name=up.name)
    else:
        kind = st.selectbox("Continuum", ["flat_flam", "flat_fnu", "powerlaw"])
        alpha = st.slider("Power-law α (f_λ ∝ λ^α)", -3.0, 2.0, -1.0, 0.1) \
            if kind == "powerlaw" else 0.0
        spec0 = etc.continuum(kind, alpha)
        if st.checkbox("Add a line"):
            c1, c2, c3 = st.columns(3)
            lam_l = c1.number_input("λ₀ (Å)", 4000.0, 8092.0, 6562.8)
            ew = c2.number_input("EW (Å)", 0.1, 5000.0, 50.0)
            fw = c3.number_input("FWHM (km/s)", 50.0, 20000.0, 300.0)
            absorb = st.checkbox("Absorption (else emission)")
            spec0 = etc.add_line(spec0, lam_l, ew, fw, absorption=absorb)

    z = st.number_input("Redshift z", 0.0, 9.0, 0.0, 0.01)
    ebv = st.number_input("Galactic E(B−V)", 0.0, 2.0, 0.0, 0.01)
    rv = st.number_input("R_V", 2.0, 6.0, 3.1, 0.1)

    st.header("Normalization (observed)")
    norm_mode = st.radio("Mode", ["SDSS filter mag", "AB mag at λ₀",
                                  "Total flux in range"])
    extended = st.checkbox("Extended source (per arcsec²)")
    if norm_mode == "SDSS filter mag":
        band = st.selectbox("Band", list("ugriz"), index=2)
        mag = st.number_input("AB magnitude", 5.0, 28.0, 19.0, 0.1)
    elif norm_mode == "AB mag at λ₀":
        lam0 = st.number_input("λ₀ (Å)", 4000.0, 8092.0, 6000.0)
        mag = st.number_input("AB magnitude", 5.0, 28.0, 19.0, 0.1)
    else:
        fx = st.number_input("Total flux (erg/s/cm²)", value=5e-15,
                             format="%.3e")
        w0 = st.number_input("Range start (Å)", 4000.0, 8092.0, 6500.0)
        w1 = st.number_input("Range end (Å)", 4000.0, 8092.0, 6630.0)

    st.header("Conditions & instrument")
    slit = st.selectbox("Slit (R at 6000 Å)",
                        ["1.8", "2.7", "4.5"],
                        format_func=lambda s: f'{s}\"  (R={etc.SLITS[s]["R"]})')
    seeing = st.slider("Seeing FWHM (″)", 0.7, 5.0, 1.5, 0.1)
    src_fwhm = st.slider("Source FWHM (″, 0 = point)", 0.0, 10.0, 0.0, 0.5)
    psf = st.selectbox("PSF profile", ["gaussian", "moffat"])
    beta = st.slider("Moffat β", 2.0, 5.0, 3.5, 0.5) if psf == "moffat" else 3.5
    lunar = st.selectbox("Moon", ["dark", "gray", "bright"])
    clouds = st.selectbox("Clouds", ["photometric", "thin cirrus", "cloudy"])
    alt = st.selectbox("Target altitude (airmass)", [90, 75, 60, 45, 30],
                       index=2, format_func=lambda a:
                       f"{a}°  (X={1/np.sin(np.radians(a)):.2f})")
    temp = st.selectbox("CCD temperature", ["-80C", "-100C"])
    readout = st.selectbox("Readout", ["slow", "medium", "fast"],
                           format_func=lambda m:
                           f"{m} ({config.READOUT_RATE[m]}, "
                           f"RN {config.READ_NOISE_E[m]:.0f} e⁻)")
    use_emp = st.toggle("Commissioning-calibrated (×0.29 + measured shape)",
                        value=True,
                        help="Off = theoretical component model")

# ----------------------- build the observed spectrum --------------------
try:
    spec = spec0
    if z > 0:
        spec = etc.redshift_spectrum(spec, z, verbose=False)
    if ebv > 0:
        spec = etc.apply_extinction(spec, ebv=ebv, rv=rv)
    if norm_mode == "SDSS filter mag":
        spec = etc.normalize_filter(spec, mag, band, extended=extended)
    elif norm_mode == "AB mag at λ₀":
        spec = etc.normalize(spec, mag, lam0, extended=extended)
    else:
        spec = etc.normalize_total_flux(spec, fx, (w0, w1), extended=extended)
except (ValueError, KeyError) as e:
    st.error(f"Spectrum setup failed: {e}")
    st.stop()

ETC_KW = dict(slit=slit, seeing=seeing, lunar=lunar, clouds=clouds,
              altitude_deg=alt, temperature=temp, readout=readout,
              source_fwhm_arcsec=src_fwhm, use_empirical=use_emp,
              psf_profile=psf, moffat_beta=beta)

def show_warnings(res, t_frame):
    wrn = etc.check_warnings(res, t_frame)
    msgs = etc.warning_summary(wrn, res["wav"])
    for m in msgs:
        (st.error if "saturat" in m else st.warning)(f"⚠ {m}")
    if not msgs:
        st.success("No detector warnings — source/sky shot-noise regime, "
                   "unsaturated.")

# ========================== main: three modes ===========================
tab1, tab2, tab3 = st.tabs(["① S/N ↔ exposure time",
                            "② Simulated 1-D & 2-D spectra",
                            "③ Recommended request time"])

with tab1:
    sub = st.radio("Calculation", ["S/N from N × t", "Time for target S/N"],
                   horizontal=True)
    if sub == "S/N from N × t":
        c1, c2 = st.columns(2)
        n_fr = c1.number_input("N frames", 1, 100, 3)
        t_fr = c2.number_input("Exposure / frame (s)", 1.0, 3600.0, 600.0)
        if st.button("Calculate S/N", type="primary"):
            res = etc.run_lrs_etc(spec, t_fr, n_fr, **ETC_KW)
            i6 = np.argmin(np.abs(res["wav"] - 6000))
            c1, c2, c3 = st.columns(3)
            c1.metric("Median S/N / pixel",
                      f"{np.nanmedian(res['snr_pix']):.1f}")
            c2.metric("S/N at 6000 Å", f"{res['snr_pix'][i6]:.1f}")
            c3.metric("Slit coupling", f"{res['fslit']:.2f}")
            st.pyplot(etc.snr_figure(res))
            fig, _ = etc.expected_spectrum_figure(res, t_fr,
                                                  spec.get("name", ""))
            st.pyplot(fig)
            show_warnings(res, t_fr)
            st.caption(f"Calibration: {res['config']['calibration']}")
    else:
        c1, c2, c3 = st.columns(3)
        snr_t = c1.number_input("Target S/N / pixel", 1.0, 500.0, 10.0)
        lam_t = c2.number_input("At wavelength (Å)", 4000.0, 8092.0, 6563.0)
        n_fr = c3.number_input("N frames", 1, 100, 3)
        if st.button("Calculate required time", type="primary"):
            t_tot = etc.required_time(spec, snr_t, lam_t, n_frames=n_fr,
                                      **ETC_KW)
            if not np.isfinite(t_tot):
                st.error("No solution — the template has no flux at that "
                         "wavelength (zero-patched or outside coverage).")
            else:
                c1, c2, c3 = st.columns(3)
                c1.metric("Total on-source", f"{t_tot/60:.1f} min")
                c2.metric("Per frame", f"{t_tot/n_fr:.0f} s")
                c3.metric("At wavelength", f"{lam_t:.0f} Å")
                res = etc.run_lrs_etc(spec, t_tot / n_fr, n_fr, **ETC_KW)
                st.pyplot(etc.snr_figure(res))
                show_warnings(res, t_tot / n_fr)
                st.caption(f"Calibration: {res['config']['calibration']}")

with tab2:
    c1, c2 = st.columns(2)
    n_fr2 = c1.number_input("N frames ", 1, 100, 3)
    t_fr2 = c2.number_input("Exposure / frame (s) ", 1.0, 3600.0, 600.0)
    if st.button("Simulate spectra", type="primary"):
        res = etc.run_lrs_etc(spec, t_fr2, n_fr2, **ETC_KW)
        st.subheader("Expected 1-D spectrum (flux · ADU · S/N)")
        fig, _ = etc.expected_spectrum_figure(res, t_fr2,
                                              spec.get("name", ""))
        st.pyplot(fig)
        show_warnings(res, t_fr2)
        st.subheader("Raw 2-D frame (single exposure, ADU)")
        st.pyplot(etc.raw_frame_figure(spec, t_fr2, seeing=seeing,
                                       **{k: v for k, v in ETC_KW.items()
                                          if k != "seeing"}))
        st.subheader("Reduced 2-D S/N map (stacked, sky-subtracted)")
        st.pyplot(etc.reduced_2d_figure(spec, t_fr2, n_fr2, seeing=seeing,
                                        **{k: v for k, v in ETC_KW.items()
                                           if k != "seeing"}))
        st.caption("Vignetting: only the central 3′ of the slit is "
                   "illuminated. Sky-subtraction residuals under the OH "
                   "bands are simulated with a deg-2 fit against the true "
                   "deg-4 slit illumination.")

with tab3:
    c1, c2, c3 = st.columns(3)
    n_fr3 = c1.number_input("N frames  ", 1, 100, 3)
    t_fr3 = c2.number_input("Exposure / frame (s)  ", 1.0, 3600.0, 600.0)
    night_cal = c3.checkbox("Night calibration (1 h / half night)")
    if st.button("Compute request", type="primary"):
        b = etc.recommended_request_time(t_per_frame=t_fr3, n_frames=n_fr3,
                                         readout=readout,
                                         night_calibration=night_cal,
                                         label="this program")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("On-source", f"{b['science']/60:.0f} min")
        c2.metric("In-dome", f"{b['in_dome']/60:.0f} min")
        c3.metric("REQUEST", f"{b['requested']/60:.0f} min")
        c4.metric("Open-shutter efficiency", f"{b['efficiency']:.0%}")
        st.pyplot(etc.request_breakdown_figure(b))
        st.info(f"ℹ Daytime calibration (bias/dark/flat/arc): "
                f"~{b['daytime_cal_info']/60:.0f} min for "
                f"{b['n_half_nights']} half-night(s) — NOT charged to the "
                f"request. Standard-star rounds included: "
                f"{b['n_std_rounds']}.")

st.divider()
st.caption("MONET Calc · NARIT LRS ETC · empirical calibration anchored to the "
           "2026-04-01 BD+75°325 commissioning standard (shutter-corrected)"
           " · theoretical mode available via the sidebar toggle · "
           "krittapas [at] narit.or.th")
