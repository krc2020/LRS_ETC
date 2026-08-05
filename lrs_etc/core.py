"""
lrs_etc.core — spectra, conditions, and the ETC engine.

Physics identical to LRS_ETC.ipynb v6.4. All spectra are dicts:
{"wav": Angstrom grid, "flam": erg/s/cm2/A, "name": str, ...}.
"""

import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter1d
from scipy.special import erf

from .config import (WAV, PIX_WAV, HC, C_AA_S, A_CM2, SLITS, DISP_AA_PIX,
                     SPATIAL_AS_PIX, SLIT_LENGTH_ARCSEC, DARK_E_PIX_S,
                     READ_NOISE_E, GAIN_E_ADU, FULL_WELL_ADU, BIAS_ADU,
                     LINEARITY_FRAC, MIRROR_R_550, SKY_V_AB, CLOUD_MAG,
                     ALTITUDE_CHOICES_DEG, TNO_LAT_DEG, DATA_ROOT,
                     EMPIRICAL_RATIO_GRID, Z_MAX, READOUT_TIME_S, SLEW_S,
                     OFFSET_S, ACQ_IMG_S, ACQ_ROUNDS, THRUSLIT_S,
                     OPERATOR_S, WEATHER_LOSS, BLOCK_S, STD_INTERVAL_S,
                     STD_EXP_S, STD_MAX_EXP, HALF_NIGHT_S, NIGHT_CAL_S,
                     DAYTIME_CAL_S)
from . import config

_trapz = getattr(np, "trapezoid", np.trapz)

# ========================================================================
# Throughput
# ========================================================================

def _load_throughput():
    tbl = np.genfromtxt(DATA_ROOT / "LRS_throughput.csv",
                        delimiter=",", names=True)
    wav_aa = tbl["wavelength_nm"] * 10.0
    per_mirror = tbl["TNT_mirrors_M1xM2xM3xM4"] ** 0.25
    lrs_only = tbl["LRS_only_throughput"]
    i550 = np.argmin(np.abs(wav_aa - 5500.0))
    scale_prod = np.prod([r / per_mirror[i550]
                          for r in MIRROR_R_550.values()])
    total = per_mirror ** 4 * scale_prod * lrs_only
    eta = np.interp(WAV, wav_aa, total, left=total[0], right=total[-1])
    taper = np.ones_like(WAV)
    blue = WAV < 4000
    taper[blue] = np.exp(-((4000 - WAV[blue]) / 250.0) ** 2)
    red = WAV > 8000
    taper[red] = np.exp(-((WAV[red] - 8000) / 500.0) ** 2)
    return eta * taper

ETA = _load_throughput()

def system_throughput(wav_aa):
    """Total telescope+LRS+CCD throughput (theoretical model) at wav_aa."""
    return np.interp(wav_aa, WAV, ETA)

# ========================================================================
# Spectral templates
# ========================================================================

LIB_DIR = DATA_ROOT / "spectral_library"
RESEARCH = ({f.stem: f for f in sorted(LIB_DIR.glob("*_*.txt"))}
            if LIB_DIR.exists() else {})

def list_templates():
    """Names accepted by load_template()."""
    return sorted(RESEARCH)

def load_template(name, verbose=True):
    """Load a research-library template by filename stem.

    Examples: 'pickles_g2v', 'kc96_starb2', 'agn_qso_ext', 'dobos_SF2',
    'brown_NGC5033_sy1.5'. Returns a spectrum dict with native-grid
    arrays retained for redshifting.
    """
    if name not in RESEARCH:
        raise KeyError(f"'{name}' not in spectral_library/. "
                       f"Available: {list_templates()}")
    t = np.loadtxt(RESEARCH[name])
    wmin, wmax = t[:, 0].min(), t[:, 0].max()
    if verbose and (wmin > PIX_WAV[0] or wmax < PIX_WAV[-1]):
        print(f"note: '{name}' native coverage {wmin:.0f}-{wmax:.0f} A does "
              f"not span the full LRS band; S/N outside it will be NaN")
    return dict(wav=WAV.copy(),
                flam=np.interp(WAV, t[:, 0], t[:, 1],
                               left=np.nan, right=np.nan),
                wav_native=t[:, 0], flam_native=t[:, 1], name=name)

def user_table_spectrum(path):
    """2-column ASCII: wavelength_A  f_lam (any absolute scale)."""
    t = np.loadtxt(path)
    return dict(wav=WAV.copy(),
                flam=np.interp(WAV, t[:, 0], t[:, 1],
                               left=np.nan, right=np.nan),
                wav_native=t[:, 0], flam_native=t[:, 1],
                name=Path(path).name)

def continuum(kind="flat_flam", alpha=0.0):
    """'flat_flam' | 'flat_fnu' | 'powerlaw' (f_lam ~ lambda^alpha)."""
    if kind == "flat_flam":
        f = np.ones_like(WAV)
    elif kind == "flat_fnu":
        f = (6000.0 / WAV) ** 2
    elif kind == "powerlaw":
        f = (WAV / 6000.0) ** alpha
    else:
        raise ValueError(kind)
    return dict(wav=WAV.copy(), flam=f, name=f"{kind} (alpha={alpha})")

def _gauss_line(center, ew_aa, fwhm_aa, cont, absorption=False):
    sig = fwhm_aa / 2.3548
    prof = np.exp(-0.5 * ((WAV - center) / sig) ** 2)
    amp = ew_aa * np.interp(center, WAV, cont) / (sig * np.sqrt(2 * np.pi))
    return -amp * prof if absorption else amp * prof

def add_line(spec, center_aa, ew_aa, fwhm_kms=300.0, absorption=False):
    """Add a Gaussian emission/absorption line (EW in A, FWHM in km/s)."""
    fwhm_aa = center_aa * fwhm_kms / 2.998e5
    out = dict(spec)
    out["flam"] = np.clip(spec["flam"] + _gauss_line(
        center_aa, ew_aa, fwhm_aa, spec["flam"], absorption), 0, None)
    return out

# ========================================================================
# Normalization & photometry
# ========================================================================

def _ab_to_flam(mag_ab, wav_aa):
    return 10 ** (-0.4 * (mag_ab + 48.60)) * C_AA_S / wav_aa ** 2

def normalize(spec, mag_ab, wav0=6000.0, extended=False):
    """Scale so f_lam(wav0) matches mag_ab (AB; per arcsec^2 if extended)."""
    current = np.interp(wav0, spec["wav"], spec["flam"])
    if not np.isfinite(current) or current <= 0:
        raise ValueError(f"no template data at {wav0:.0f} A - "
                         f"choose wav0 inside the covered range")
    out = dict(spec)
    out["flam"] = spec["flam"] * (_ab_to_flam(mag_ab, wav0) / current)
    out["extended"] = extended
    out["norm"] = f"{mag_ab} AB{'/arcsec^2' if extended else ''} at {wav0:.0f} A"
    return out

_FILT_DIR = DATA_ROOT / "filters"
SDSS_BANDS = ({b: np.loadtxt(_FILT_DIR / f"sdss_{b}.txt")
               for b in "ugriz"
               if (_FILT_DIR / f"sdss_{b}.txt").exists()}
              if _FILT_DIR.exists() else {})

def synth_ab_mag(spec, band, verbose=True):
    """Synthetic AB magnitude through the real SDSS filter (photon-counting)."""
    t = SDSS_BANDS[band]
    T = np.interp(WAV, t[:, 0], t[:, 1], left=0.0, right=0.0)
    fl = spec["flam"]
    ok = np.isfinite(fl)
    wt = T / WAV
    frac = _trapz(np.where(ok, wt, 0.0), WAV) / max(_trapz(wt, WAV), 1e-30)
    if verbose and frac < 0.99:
        print(f"WARNING: template covers only {frac:.0%} of the "
              f"{band}-band photon weight")
    num = _trapz(np.where(ok, fl, 0.0) * T * WAV, WAV)
    den = C_AA_S * _trapz(T / WAV, WAV)
    return np.inf if num <= 0 else -2.5 * np.log10(num / den) - 48.60

def normalize_filter(spec, mag_ab, band="r", extended=False):
    """Scale so the synthetic SDSS `band` magnitude equals mag_ab."""
    m0 = synth_ab_mag(spec, band)
    if not np.isfinite(m0):
        raise ValueError(f"no positive flux in the {band} band")
    out = dict(spec)
    out["flam"] = spec["flam"] * 10 ** (-0.4 * (mag_ab - m0))
    out["extended"] = extended
    out["norm"] = f"SDSS {band} = {mag_ab} AB{'/arcsec^2' if extended else ''}"
    return out

def normalize_total_flux(spec, total_flux, wav_range=(4000.0, 8092.0),
                         extended=False):
    """Scale so the integral of f_lam over wav_range equals total_flux."""
    w0, w1 = wav_range
    sel = (spec["wav"] >= w0) & (spec["wav"] <= w1) & np.isfinite(spec["flam"])
    cur = _trapz(spec["flam"][sel], spec["wav"][sel])
    if cur <= 0:
        raise ValueError(f"no positive flux in {w0:.0f}-{w1:.0f} A")
    out = dict(spec)
    out["flam"] = spec["flam"] * (total_flux / cur)
    out["extended"] = extended
    out["norm"] = f"integral {w0:.0f}-{w1:.0f} A = {total_flux:.3e}"
    return out

# ========================================================================
# Extinction & redshift
# ========================================================================

def _ccm89_builtin(wav_aa, rv):
    x = 1e4 / wav_aa
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

def _calzetti_k(wav_aa):
    um = wav_aa / 1e4
    k = np.where(um < 0.63,
                 2.659 * (-2.156 + 1.509/um - 0.198/um**2 + 0.011/um**3) + 4.05,
                 2.659 * (-1.857 + 1.040/um) + 4.05)
    return np.clip(k, 0, None)

def extinction_transmission(law="CCM89", rv=3.1, ebv=None, av=None):
    """Fractional transmission vs WAV. Laws: CCM89, O94, F99, Calzetti00."""
    if (ebv is None) == (av is None):
        raise ValueError("give exactly one of ebv= or av=")
    if av is None:
        av = rv * ebv
    ebv = av / rv
    if law == "Calzetti00":
        return 10 ** (-0.4 * (av / 4.05) * _calzetti_k(WAV))
    try:
        from dust_extinction.parameter_averages import CCM89, O94, F99
        import astropy.units as u
        model = {"CCM89": CCM89, "O94": O94, "F99": F99}[law](Rv=rv)
        return model.extinguish(WAV * u.AA, Ebv=ebv)
    except ImportError:
        if law in ("CCM89", "O94", "F99"):
            return 10 ** (-0.4 * av * _ccm89_builtin(WAV, rv))
        raise

def apply_extinction(spec, law="CCM89", rv=3.1, ebv=None, av=None):
    """Reddened copy of spec. Apply BEFORE normalize() for observed mags."""
    out = dict(spec)
    out["flam"] = spec["flam"] * extinction_transmission(law, rv, ebv, av)
    tag = f"E(B-V)={ebv}" if ebv is not None else f"A_V={av}"
    out["name"] = f"{spec.get('name','spec')} + {law} {tag}"
    return out

def luminosity_distance_Mpc(z, H0=70.0, Om=0.30):
    try:
        from astropy.cosmology import FlatLambdaCDM
        import astropy.units as u
        return float(FlatLambdaCDM(H0=H0, Om0=Om)
                     .luminosity_distance(z).to(u.Mpc).value)
    except ImportError:
        if z <= 0:
            return 0.0
        zs = np.linspace(0, z, 2048)
        E = np.sqrt(Om * (1 + zs) ** 3 + (1 - Om))
        return 2.998e5 / H0 * _trapz(1.0 / E, zs) * (1 + z)

_REST_REGIMES = [(0., 100., "X-ray"), (100., 912., "extreme-UV"),
                 (912., 3200., "UV"), (3200., 7000., "optical"),
                 (7000., 2.5e4, "NIR"), (2.5e4, np.inf, "IR")]

def _regimes(lo, hi):
    return " and ".join(n for a, b, n in _REST_REGIMES if lo < b and hi > a)

def redshift_spectrum(spec, z, d_L_ref_Mpc=None, verbose=True):
    """Shift to redshift z (0..9). Missing coverage is zero-patched with a
    named rest-frame-regime warning. d_L_ref_Mpc enables cosmological
    dimming for absolute-flux templates."""
    if not (0.0 <= z <= Z_MAX):
        raise ValueError(f"z must be between 0 and {Z_MAX}")
    wr = spec.get("wav_native", spec["wav"])
    fr = spec.get("flam_native", spec["flam"])
    ok = np.isfinite(fr)
    wr, fr = wr[ok], fr[ok]
    wav_obs, flam_obs = wr * (1 + z), fr / (1 + z)
    if d_L_ref_Mpc is not None and z > 0:
        DL = luminosity_distance_Mpc(z)
        flam_obs = flam_obs * (d_L_ref_Mpc / DL) ** 2
    name = spec.get("name", "spec")
    if verbose and wav_obs.min() > WAV[0]:
        print(f"WARNING [{name} @ z={z:g}]: lack of template in the "
              f"rest-frame {_regimes(WAV[0]/(1+z), wr.min())} - observed "
              f"{WAV[0]:.0f}-{wav_obs.min():.0f} A patched with zeros")
    if verbose and wav_obs.max() < WAV[-1]:
        print(f"WARNING [{name} @ z={z:g}]: lack of template in the "
              f"rest-frame {_regimes(wr.max(), WAV[-1]/(1+z))} - observed "
              f"{wav_obs.max():.0f}-{WAV[-1]:.0f} A patched with zeros")
    out = dict(wav=WAV.copy(),
               flam=np.interp(WAV, wav_obs, flam_obs, left=0.0, right=0.0),
               name=f"{name} @ z={z:.2f}", z=z)
    if d_L_ref_Mpc is not None and z > 0:
        out["d_L_Mpc"] = DL
    return out

def redshift_limits(name, band=(4000.0, 8000.0)):
    """(z_full_lo, z_full_hi, z_any) for full/partial band coverage."""
    t = np.loadtxt(RESEARCH[name])
    lmin, lmax = t[:, 0].min(), t[:, 0].max()
    return (max(0.0, band[1] / lmax - 1), band[0] / lmin - 1,
            band[1] / lmin - 1)

# ========================================================================
# Conditions
# ========================================================================

_SKY_SHAPE_W = [3500, 4500, 5500, 6500, 7500, 8500, 9500]
_SKY_SHAPE_D = [+0.8, +0.4, 0.0, -0.3, -0.7, -1.0, -1.1]
_SKY_SHAPE_B = [-0.4, -0.2, 0.0, -0.2, -0.5, -0.8, -0.9]
_EXT_W = [3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000, 9500]
_EXT_K = [0.55, 0.35, 0.25, 0.18, 0.15, 0.12, 0.09, 0.06, 0.05]

SKY_LINES = [
    (5577.3, 14.0), (5890.0, 2.5), (5895.9, 1.5), (6300.3, 6.0),
    (6363.8, 2.0), (6834, 2.0), (6871, 2.5), (6923, 1.8), (7240, 2.5),
    (7276, 2.8), (7316, 3.0), (7340, 2.5), (7369, 2.2), (7524, 3.5),
    (7571, 3.0), (7623, 2.8), (7750, 3.2), (7794, 3.5), (7821, 3.2),
    (7913, 4.0), (7964, 3.5), (7993, 4.0), (8025, 3.6), (8399, 4.5),
    (8430, 4.2), (8465, 4.0), (8505, 4.5), (8541, 4.2), (8615, 4.0),
    (8645, 5.0), (8791, 5.0), (8827, 4.6), (8886, 4.8), (8919, 4.5),
    (8943, 4.8), (9002, 4.5), (9306, 5.0), (9375, 4.6), (9439, 5.0)]

def sky_surface_brightness(lunar="dark"):
    """Continuum sky brightness, AB mag/arcsec^2 vs WAV."""
    shape = _SKY_SHAPE_B if lunar == "bright" else _SKY_SHAPE_D
    return SKY_V_AB[lunar] + np.interp(WAV, _SKY_SHAPE_W, shape)

def sky_emission_flam(lunar="dark", R=750):
    """Sky f_lam/arcsec^2: lunar continuum + fixed airglow/OH lines at R."""
    cont = _ab_to_flam(sky_surface_brightness(lunar), WAV)
    cont_dark = _ab_to_flam(sky_surface_brightness("dark"), WAV)
    lines = np.zeros_like(WAV)
    for cen, px in SKY_LINES:
        sig = max(cen / R, 6.0) / 2.3548
        lines += px * np.interp(cen, WAV, cont_dark) * np.exp(
            -0.5 * ((WAV - cen) / sig) ** 2)
    return cont + lines

def atmospheric_transmission(airmass=1.2, clouds="photometric"):
    k = np.interp(WAV, _EXT_W, _EXT_K)
    return 10 ** (-0.4 * (k * airmass + CLOUD_MAG[clouds]))

def select_airmass(altitude_deg):
    """Airmass (sec z) for one of the allowed altitudes: 30/45/60/75/90."""
    if altitude_deg not in ALTITUDE_CHOICES_DEG:
        raise ValueError(f"altitude_deg must be one of "
                         f"{ALTITUDE_CHOICES_DEG}")
    return 1.0 / np.sin(np.radians(altitude_deg))

def target_visibility(ra_deg, dec_deg, verbose=True):
    """Culmination altitude range at the TNO with observability flags."""
    if not (0.0 <= ra_deg < 360.0):
        raise ValueError("RA must be in [0, 360)")
    if not (-90.0 <= dec_deg <= 90.0):
        raise ValueError("Dec must be in [-90, 90]")
    alt_max = 90.0 - abs(TNO_LAT_DEG - dec_deg)
    alt_min = abs(TNO_LAT_DEG + dec_deg) - 90.0
    status = ("not observable" if alt_max < 30 else
              "observable" if alt_max < 50 else "optimal")
    reachable = [a for a in ALTITUDE_CHOICES_DEG if a <= alt_max]
    alt_choice = max(reachable) if reachable else None
    if verbose:
        print(f"RA {ra_deg:.3f}, Dec {dec_deg:+.3f} -> alt {alt_min:+.1f} "
              f"to {alt_max:.1f} deg: {status}"
              + (f"; use altitude_deg={alt_choice}" if alt_choice else ""))
    return dict(alt_max=alt_max, alt_min=alt_min, status=status,
                circumpolar=alt_min > 0, altitude_deg=alt_choice,
                transit_lst_h=ra_deg / 15.0)

# ========================================================================
# Slit coupling & the ETC engine
# ========================================================================

def slit_coupling(seeing_arcsec, slit_width_arcsec, extract_arcsec,
                  profile="gaussian", moffat_beta=3.5):
    """PSF fraction through the rectangular slit x extraction aperture.

    'gaussian': exact (circular Gaussian separates -> erf product).
    'moffat'  : numeric integral; realistic wings couple 10-25% less.
    """
    if profile == "gaussian":
        sig = seeing_arcsec / 2.3548
        return (erf(slit_width_arcsec / (2 * np.sqrt(2) * sig))
                * erf(extract_arcsec / (2 * np.sqrt(2) * sig)))
    if profile == "moffat":
        b = moffat_beta
        alpha = seeing_arcsec / (2 * np.sqrt(2 ** (1 / b) - 1))
        span = max(slit_width_arcsec, extract_arcsec) / 2 + 6 * seeing_arcsec
        x = np.linspace(-span, span, 601)
        X, Y = np.meshgrid(x, x)
        I = (b - 1) / (np.pi * alpha**2) * (1 + (X**2 + Y**2) / alpha**2) ** (-b)
        inside = ((np.abs(X) <= slit_width_arcsec / 2)
                  & (np.abs(Y) <= extract_arcsec / 2))
        return float((I * inside).sum() * (x[1] - x[0]) ** 2)
    raise ValueError(f"unknown profile '{profile}'")

def run_lrs_etc(spec, t_per_frame, n_frames, slit="1.8", seeing=1.1,
                lunar="dark", clouds="photometric", airmass=1.2,
                altitude_deg=None, temperature="-80C", readout="slow",
                extract_arcsec=None, source_fwhm_arcsec=0.0,
                use_empirical=None, psf_profile="gaussian",
                moffat_beta=3.5):
    """Core ETC: per-pixel S/N for `spec` under the given configuration.

    Returns a dict: wav, pixel, snr_pix, snr_aa, S, B (e-/s), dark_rate,
    rn_var_total, flam_pix, t_tot, fslit, n_spat, R, config.
    """
    if altitude_deg is not None:
        airmass = select_airmass(altitude_deg)
    sl = SLITS[slit]
    R, w_slit = sl["R"], sl["width_arcsec"]
    extended = spec.get("extended", False)
    fwhm_eff = float(np.hypot(seeing, source_fwhm_arcsec))
    if extract_arcsec is None:
        extract_arcsec = max(1.5 * fwhm_eff, 2 * SPATIAL_AS_PIX)
    n_spat = max(1, int(np.ceil(extract_arcsec / SPATIAL_AS_PIX)))

    flam_conv = gaussian_filter1d(np.nan_to_num(spec["flam"]),
                                  (6000.0 / R) / 2.3548)
    flam_pix = np.interp(PIX_WAV, WAV, flam_conv)
    atm = np.interp(PIX_WAV, WAV, atmospheric_transmission(airmass, clouds))
    eta = system_throughput(PIX_WAV)
    phot = PIX_WAV / HC

    if extended:
        S = (flam_pix * phot * A_CM2 * eta * atm * DISP_AA_PIX
             * w_slit * extract_arcsec)
        fslit = 1.0
    else:
        fslit = slit_coupling(fwhm_eff, w_slit, extract_arcsec,
                              profile=psf_profile, moffat_beta=moffat_beta)
        S = flam_pix * phot * A_CM2 * eta * atm * DISP_AA_PIX * fslit

    sky_pix = np.interp(PIX_WAV, WAV, sky_emission_flam(lunar, R=R))
    B = (sky_pix * phot * A_CM2 * eta * DISP_AA_PIX
         * w_slit * extract_arcsec)

    emp = (config.USE_EMPIRICAL_CALIBRATION if use_empirical is None
           else use_empirical)
    if emp:
        fac = np.interp(PIX_WAV, WAV, EMPIRICAL_RATIO_GRID)
        S, B = S * fac, B * fac

    dark = DARK_E_PIX_S[temperature]
    rn = READ_NOISE_E[readout]
    t_tot = t_per_frame * n_frames
    var = (S * t_tot + B * t_tot + n_spat * dark * t_tot
           + n_spat * n_frames * rn ** 2)
    snr_pix = np.where(var > 0, S * t_tot / np.sqrt(var), 0.0)

    return dict(wav=PIX_WAV, pixel=np.arange(len(PIX_WAV)),
                snr_pix=snr_pix, snr_aa=snr_pix / np.sqrt(DISP_AA_PIX),
                S=S, B=B, dark_rate=n_spat * dark,
                rn_var_total=n_spat * n_frames * rn ** 2,
                flam_pix=flam_pix, t_tot=t_tot, fslit=fslit,
                n_spat=n_spat, R=R,
                config=dict(slit=slit, seeing=seeing, fwhm_eff=fwhm_eff,
                            source_fwhm_arcsec=source_fwhm_arcsec,
                            lunar=lunar, clouds=clouds, airmass=airmass,
                            temperature=temperature, readout=readout,
                            extract_arcsec=extract_arcsec,
                            calibration=("commissioning-empirical" if emp
                                         else "theoretical model"),
                            t_per_frame=t_per_frame, n_frames=n_frames))

def required_time(spec, target_snr, wav_aa, n_frames=3, **etc_kwargs):
    """Total on-source seconds for target_snr per pixel at wav_aa."""
    probe = run_lrs_etc(spec, 1.0, n_frames, **etc_kwargs)
    i = np.argmin(np.abs(probe["wav"] - wav_aa))
    a = probe["S"][i]
    b = probe["S"][i] + probe["B"][i] + probe["dark_rate"]
    c = probe["rn_var_total"]
    if a <= 0:
        return np.nan
    s2 = target_snr ** 2
    return (s2 * b + np.sqrt(s2**2 * b**2 + 4 * a**2 * s2 * c)) / (2 * a**2)

def check_warnings(res, t_per_frame):
    """Saturation / non-linearity / noise-regime masks for an ETC result."""
    cfg = res["config"]
    fwhm = cfg.get("fwhm_eff", cfg["seeing"])
    h = cfg["extract_arcsec"]
    sig = fwhm / 2.3548
    fy_ext = erf((h / 2) / (np.sqrt(2) * sig))
    fy_cen = erf((SPATIAL_AS_PIX / 2) / (np.sqrt(2) * sig))
    peak_frac = np.clip(fy_cen / max(fy_ext, 1e-9), 0, 1)
    sky_pix = res["B"] / res["n_spat"]
    dark = DARK_E_PIX_S[cfg["temperature"]]
    peak_adu = ((res["S"] * peak_frac + sky_pix + dark) * t_per_frame
                / GAIN_E_ADU + BIAS_ADU)
    sat = peak_adu >= FULL_WELL_ADU
    nonlin = (peak_adu >= LINEARITY_FRAC * FULL_WELL_ADU) & ~sat
    t_tot = res["t_tot"]
    var = np.vstack([res["S"] * t_tot, res["B"] * t_tot,
                     np.full_like(res["S"], res["dark_rate"] * t_tot),
                     np.full_like(res["S"], res["rn_var_total"])])
    dom = np.argmax(var, axis=0)
    return dict(peak_adu=peak_adu, saturated=sat, nonlinear=nonlin,
                rn_dominated=(dom == 3) & ~sat & ~nonlin,
                dark_dominated=(dom == 2) & ~sat & ~nonlin,
                peak_row_frac=peak_frac)

def warning_summary(wrn, wav):
    """List of human-readable warning strings for check_warnings output."""
    out = []
    for key, lab in [("saturated", "saturation"),
                     ("nonlinear", "non-linearity (>95% full well)"),
                     ("rn_dominated", "read-noise-dominated"),
                     ("dark_dominated", "dark-current-dominated")]:
        m = wrn[key]
        if m.any():
            out.append(f"{lab}: {m.sum()} pixels "
                       f"({wav[m].min():.0f}-{wav[m].max():.0f} A)")
    return out

# ========================================================================
# Request-time budget
# ========================================================================

def recommended_request_time(t_per_frame=None, n_frames=None, science_s=None,
                             readout="slow", n_std_exp=STD_MAX_EXP,
                             night_calibration=False, label=None,
                             verbose=False):
    """Recommended REQUEST time per target: overheads, standards every 2 h,
    35% weather margin; night calibration (1 h/half-night) outside the
    margin. Returns a breakdown dict (all seconds)."""
    if science_s is None:
        science_s = t_per_frame * n_frames
    elif n_frames is None:
        t_per_frame = min(900.0, science_s)
        n_frames = int(np.ceil(science_s / t_per_frame))
    n_std_exp = min(int(n_std_exp), STD_MAX_EXP)
    readout_s = n_frames * READOUT_TIME_S[readout]
    std_round_s = (SLEW_S + ACQ_ROUNDS * (OFFSET_S + ACQ_IMG_S) + THRUSLIT_S
                   + n_std_exp * (STD_EXP_S + READOUT_TIME_S[readout]))
    n_blocks, n_std = 1, 1
    for _ in range(12):
        acq_s = n_blocks * ACQ_ROUNDS * (OFFSET_S + ACQ_IMG_S)
        thruslit_s = n_blocks * THRUSLIT_S
        operator_s = n_blocks * OPERATOR_S
        in_dome = (science_s + readout_s + SLEW_S + acq_s + thruslit_s
                   + operator_s + n_std * std_round_s)
        nb = max(1, int(np.ceil(in_dome / BLOCK_S)))
        ns = max(1, int(np.ceil(in_dome / STD_INTERVAL_S)))
        if nb == n_blocks and ns == n_std:
            break
        n_blocks, n_std = nb, ns
    requested_night = in_dome / (1.0 - WEATHER_LOSS)
    n_half = max(1, int(np.ceil(requested_night / HALF_NIGHT_S)))
    night_cal_s = n_half * NIGHT_CAL_S if night_calibration else 0.0
    requested = requested_night + night_cal_s
    out = dict(label=label or f"{n_frames}x{t_per_frame:.0f} s ({readout})",
               science=science_s, readout=readout_s, slew=SLEW_S,
               acquisition=acq_s, thruslit=thruslit_s, operator=operator_s,
               standards=n_std * std_round_s,
               weather=requested_night - in_dome, night_cal=night_cal_s,
               in_dome=in_dome, requested=requested,
               daytime_cal_info=n_half * DAYTIME_CAL_S,
               n_blocks=n_blocks, n_std_rounds=n_std,
               n_half_nights=n_half, n_frames=n_frames,
               efficiency=science_s / requested)
    if verbose:
        print(f"{out['label']}: science {science_s/60:.1f} min, "
              f"REQUEST {requested/60:.1f} min "
              f"(efficiency {out['efficiency']:.0%})")
    return out
