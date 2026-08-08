"""
lrs_etc.config — instrument constants, site, data paths, calibration switch.

All values mirror LRS_ETC.ipynb v6.4. Datasheet numbers are user-confirmed;
mirror reflectivities are the 2026 NARIT measurements.
"""

from pathlib import Path
import numpy as np

__version__ = "6.5.0"

# --- data root: repo checkout layout ------------------------------------
def _find_data_root():
    cands = [
        Path(__file__).resolve().parent.parent,        # repo root (editable)
        Path.cwd(),
        Path.home() / "Documents" / "LRS_ETC",
    ]
    for c in cands:
        if (c / "LRS_throughput.csv").exists():
            return c
    return cands[0]

DATA_ROOT = _find_data_root()

# --- internal wavelength grid & physical constants ----------------------
WAV = np.arange(3500.0, 9501.0, 1.0)         # Angstrom
H_ERG_S = 6.626e-27
C_AA_S  = 2.998e18
HC      = H_ERG_S * C_AA_S                    # erg Angstrom

# --- telescope -----------------------------------------------------------
D_M         = 2.4
OBSTRUCTION = 0.09
A_CM2       = np.pi / 4 * (D_M * 100) ** 2 * (1 - OBSTRUCTION)

# --- long-slit configuration --------------------------------------------
SLITS = {"1.8": dict(width_arcsec=1.8, R=750),
         "2.7": dict(width_arcsec=2.7, R=500),
         "4.5": dict(width_arcsec=4.5, R=300)}
DISP_AA_PIX    = 4.0
SPATIAL_AS_PIX = 0.9
N_SPEC_PIX     = 1024
PIX_WAV        = 4000.0 + DISP_AA_PIX * np.arange(N_SPEC_PIX)
SLIT_LENGTH_ARCSEC = 180.0                    # central 3' illuminated

# --- Andor Newton BEX2-DD (user-confirmed datasheet) --------------------
CCD_NAME       = "Andor Newton 256x1024 BEX2-DD"
DARK_E_PIX_S   = {"-80C": 0.08, "-100C": 0.003}
READ_NOISE_E   = {"slow": 4.0, "medium": 12.0, "fast": 15.0}
READOUT_RATE   = {"slow": "50 kHz", "medium": "1 MHz", "fast": "3 MHz"}
GAIN_E_ADU     = 4.0
FULL_WELL_ADU  = 65535
BIAS_ADU       = 300.0
LINEARITY_FRAC = 0.95

# --- measured TNT mirror reflectivities at 550 nm (2026) ----------------
MIRROR_R_550 = {"M1": 0.904, "M2": 0.908, "M3": 0.930, "M4": 0.879}

# --- observing conditions -----------------------------------------------
SKY_V_AB  = {"dark": 21.8, "gray": 20.9, "bright": 19.0}
CLOUD_MAG = {"photometric": 0.0, "thin cirrus": 0.5, "cloudy": 1.2}
ALTITUDE_CHOICES_DEG = (30, 45, 60, 75, 90)

# --- TNO site ------------------------------------------------------------
TNO_LAT_DEG, TNO_LON_DEG, TNO_ELEV_M = 18.5738, 98.4823, 2457.0

# --- overheads (Sec 10 of the notebook) ---------------------------------
READOUT_TIME_S = {"slow": 7.0, "medium": 3.0, "fast": 1.0}
SLEW_S, OFFSET_S, ACQ_IMG_S, ACQ_ROUNDS = 120.0, 10.0, 30.0, 3
THRUSLIT_S, OPERATOR_S = 30.0, 300.0
WEATHER_LOSS, BLOCK_S = 0.35, 3600.0
STD_INTERVAL_S, STD_EXP_S, STD_MAX_EXP = 7200.0, 30.0, 10
HALF_NIGHT_S, NIGHT_CAL_S, DAYTIME_CAL_S = 6 * 3600.0, 3600.0, 3600.0

# --- empirical commissioning calibration (v6.4, shutter-corrected) ------
# BD+75d325 2026-04-01: measured/model response ratio; +7 s stuck-shutter
# correction applied (header 30 s -> 37 s effective).
SHUTTER_EXTRA_S = 7.0
_EMP_W = np.array([4500, 4750, 5000, 5250, 5500, 5800, 6000, 6250,
                   6500, 6800, 7000, 7200, 7400, 7800], float)
_EMP_R = np.array([0.07, 0.11, 0.18, 0.25, 0.31, 0.31, 0.30, 0.29,
                   0.25, 0.25, 0.24, 0.23, 0.21, 0.21]) * (30.0 / 37.0)
EMPIRICAL_RATIO_GRID = np.interp(WAV, _EMP_W, _EMP_R,
                                 left=_EMP_R[0], right=_EMP_R[-1])

#: Default calibration mode: True = match 2026-04-01 commissioning data,
#: False = theoretical component model. Override per call with
#: run_lrs_etc(..., use_empirical=...).
USE_EMPIRICAL_CALIBRATION = True

# --- redshift limit ------------------------------------------------------
Z_MAX = 9.0

# --- palette -------------------------------------------------------------
NAVY, LIME, ORANGE, RED = "#1E2761", "#39B54A", "#F4B400", "#E94F37"
