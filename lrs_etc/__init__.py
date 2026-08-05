"""
lrs_etc — NARIT LRS @ 2.4-m TNT Exposure-Time Calculator.

Quick start::

    import lrs_etc as etc

    spec = etc.normalize_filter(etc.load_template("kc96_sb"), 19.0, "r")
    res  = etc.run_lrs_etc(spec, t_per_frame=600, n_frames=3,
                           slit="1.8", seeing=1.5, lunar="dark")
    print(res["snr_pix"].max())

    t = etc.required_time(spec, target_snr=10, wav_aa=6563)
    budget = etc.recommended_request_time(science_s=t)

See COOKBOOK.md for the full parameter reference. The empirical
commissioning calibration is ON by default
(`lrs_etc.config.USE_EMPIRICAL_CALIBRATION`); pass
`use_empirical=False` to any run for the theoretical model.
"""

from .config import (__version__, WAV, PIX_WAV, SLITS, DISP_AA_PIX,
                     SPATIAL_AS_PIX, GAIN_E_ADU, FULL_WELL_ADU,
                     ALTITUDE_CHOICES_DEG, USE_EMPIRICAL_CALIBRATION,
                     EMPIRICAL_RATIO_GRID, DATA_ROOT)
from .core import (system_throughput, list_templates, load_template,
                   user_table_spectrum, continuum, add_line,
                   normalize, normalize_filter, normalize_total_flux,
                   synth_ab_mag, extinction_transmission, apply_extinction,
                   redshift_spectrum, redshift_limits,
                   luminosity_distance_Mpc,
                   sky_emission_flam, sky_surface_brightness,
                   atmospheric_transmission, select_airmass,
                   target_visibility, slit_coupling,
                   run_lrs_etc, required_time, check_warnings,
                   warning_summary, recommended_request_time)
from .plots import (snr_figure, expected_spectrum_figure, raw_frame_figure,
                    reduced_2d_figure, request_breakdown_figure)

__all__ = [
    "__version__", "WAV", "PIX_WAV", "SLITS", "DATA_ROOT",
    "USE_EMPIRICAL_CALIBRATION", "EMPIRICAL_RATIO_GRID",
    "ALTITUDE_CHOICES_DEG",
    "system_throughput", "list_templates", "load_template",
    "user_table_spectrum", "continuum", "add_line",
    "normalize", "normalize_filter", "normalize_total_flux", "synth_ab_mag",
    "extinction_transmission", "apply_extinction",
    "redshift_spectrum", "redshift_limits", "luminosity_distance_Mpc",
    "sky_emission_flam", "sky_surface_brightness",
    "atmospheric_transmission", "select_airmass", "target_visibility",
    "slit_coupling", "run_lrs_etc", "required_time",
    "check_warnings", "warning_summary", "recommended_request_time",
    "snr_figure", "expected_spectrum_figure", "raw_frame_figure",
    "reduced_2d_figure", "request_breakdown_figure",
]
