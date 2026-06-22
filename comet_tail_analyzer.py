#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
  comet_tail_analyzer.py  —  Finson–Probstein Dust Tail Model
  Version 3.0   ·   Teerasak Thaluang (MPC O51/O58)
=============================================================================
  Changelog:
    v3.0  • compute_orbit_diagram() now captures the comet's heliocentric
            velocity from elem_to_state() (previously discarded with a
            throwaway `_`) and returns its unit direction as
            velocity_dir_xyz, plus r_max_plot for scaling. draw_orbit_
            diagram() uses these to draw an arrow at the comet's current
            position showing its direction of motion along the orbit.
    v3.0  • beta_to_size() gained an optional Qpr parameter (default 1.0,
            backward compatible) so the GUI's new standalone Dust particle
            radius… calculator can expose it as an editable value without
            duplicating the formula.
          • NEW: compute_dust_production_rate(afrho_cm, r_h, v_dust_kms,
            p_v) — factored out of generate_dust_analysis()'s inline Afρ
            section so the GUI's new standalone Dust production rate…
            calculator can share the exact same formula instead of
            duplicating it a third time.
          • generate_dust_analysis() no longer reports F-P grain radius
            (moved to the GUI's standalone calculator — that formula needs
            only β/ρ_d/Qpr, not a per-comet model) or accepts an rho_d
            parameter (now unused). The inline Afρ/dust-production section
            is unchanged.
    v3.0  • Physical parameters (grain density ρ_d, dust expansion velocity
            v_dust, geometric albedo p_v) are no longer hardcoded — they are
            now keyword arguments on generate_dust_analysis() and
            CometTailVisualizer (rho_d), with literature-referenced
            defaults (ρ_d=0.5 g/cm³, Fulle et al. 2016; p_v=0.04,
            SEPPCoN/Schambeau et al. 2021) so the GUI's new Calculation >
            Physical Parameters… dialog can let the user override them
            per-comet instead of editing source. beta_to_size()'s default
            ρ also changed 1.0 → 0.5 g/cm³ to match.
          • fetch_from_cobs()'s H₀/n least-squares fit now also returns
            n_err, the 1-σ standard error on n from the fit residuals.
            generate_dust_analysis() uses it (together with n_fit_pts) to
            gate display: fits with <10 points or σ_n>1.0 are labeled
            "⚠ PRELIMINARY" and the derived m_pred/outburst-flag claims
            that compound n's uncertainty are suppressed rather than
            shown as if solid (thresholds are a judgment call, documented
            inline — adjust if needed).
          • CRITICAL FIX: β→grain-radius coefficient was using an
            uncorrected two-parameter form (C_pr·Qpr/(ρr), C_pr=1.19e-3
            kg/m²) that omitted a factor of 2 relative to the convention
            that constant is normally paired with (cf. Schambeau et al.
            2021, ApJ — d=diameter form; Liu & Liu 2024 — explicit
            C_pr·Qpr/(2ρa) radius form. NOTE: Moreno 2025, A&A 695 A263
            cites the same C_pr=1.19e-3 kg/m² value but does not show the
            equation's exact form in the text — not independently
            verified as using the /2 convention; not relied upon here).
            Replaced with BETA_TO_RADIUS_UM = 0.574,
            derived directly from Burns, Lamy & Soter (1979) Icarus 40,
            1, Eq. 19, sidestepping the C_pr/2 ambiguity entirely. All
            v2.4 grain radii were too LARGE by a factor of 2.07 — see the
            BETA_TO_RADIUS_UM constant comment for the full derivation.
            Fixed in all 3 places it was duplicated (beta_to_size(),
            generate_dust_analysis()'s _sz_str(), and an inline copy in
            the summary line); the latter two now delegate to
            beta_to_size() so the formula can never drift out of sync
            between call sites again.
          • All "grain size" labels/output corrected to "grain radius"
            (the value was always a radius; the ambiguous label risked
            being misread as diameter — a 2× error stacking on top of
            the formula fix above).
          • compute_model() gained project_vector_sky() + an
            'antivel_dir' key: the negative of the comet's heliocentric
            velocity projected onto the sky, drawn as a second arrow
            alongside the existing Sun-direction arrow (same convention
            as Moreno 2025 Fig. 1 and Mariblanca-Escalona et al. 2026
            Fig. 7) — shows how far a real, non-zero-ejection-velocity
            tail is expected to lean from the pure antisolar line.
          • New compute_orbit_diagram() / draw_orbit_diagram(): a full
            3D heliocentric-ecliptic view of the comet's physical
            position on its orbit (Sun as a yellow circle, perihelion
            marked, Earth's orbit for scale, drop-lines showing
            out-of-plane height) — distinct from the existing sky-
            projected 'orbit' trail in compute_model(), which only shows
            tail-axis direction, not where the comet physically sits.
            Valid uniformly for elliptical/parabolic/hyperbolic orbits
            via the conic formula r(f)=q(1+e)/(1+e·cos f) directly in
            true anomaly.
    v2.4  • COMET_DB redesigned: orbital elements removed entirely.  The
            catalogue now stores only preferred observation date (obs) and
            editorial note.  Elements are always fetched live from JPL
            Horizons, eliminating stale/incorrect local elements.
          • C/2016 R2 (PanSTARRS) corrected: Omega 325.015°→80.569°,
            omega 5.538°→33.192° (previous values from preliminary MPC orbit).
          • _normalize_comet_desig: strips discoverer parenthetical suffix
            before querying Horizons ("C/2016 R2 (PanSTARRS)"→"C/2016 R2").
            Previously all long-period PRESET comets returned 403/no-result.
          • date_to_jd accepts decimal-day notation (2018-01-14.583),
            HH:MM, and HH:MM:SS in addition to plain YYYY-MM-DD.
          • jd_to_str outputs decimal-day format (2018-01-14.5830 UT).
          • MPC fetch removed from GUI source selector (MPC API returns
            HTTP 403; JPL Horizons is the sole source).
          • fetch_comet: simplified to Horizons-only; prefer parameter and
            MPC fallback path removed entirely.
          • --prefer CLI argument removed (was redundant; Horizons is sole source).
    v2.3  • matplotlib.use() is now conditional on no backend being set,
            preventing backend conflicts when imported by CometTailGUI.py.
          • print() calls in fetch_comet() replaced with logging.info/debug.
          • ax.spines[:] replaced with list(ax.spines.values()) for
            matplotlib < 3.4 compatibility.

    v2.2  • Increased default n_pts 80→200 to fix syndyne gaps for small β.

    v2.0  • Major cleanup: removed ~1140 lines of dead code.
          • Numerical output identical to v1.1 — verified on C/2020 F3.
=============================================================================
  References:
    Finson & Probstein (1968a), ApJ 154, 327  — model & equations
    Finson & Probstein (1968b), ApJ 154, 353  — application to Arend-Roland
    Bredichin (1884) — syndyne / synchrone concepts
=============================================================================
  Requirements:
    pip install numpy scipy matplotlib astropy astroquery Pillow requests

  Quick start:
    python comet_tail_analyzer.py                                     # interactive
    python comet_tail_analyzer.py --comet "C/2020 F3 (NEOWISE)" --date 2020-07-23
    python comet_tail_analyzer.py --fetch "C/2023 A3" --date 2024-10-10
    python comet_tail_analyzer.py --image photo.jpg                   # overlay mode
=============================================================================
"""

from __future__ import annotations

__version__ = "2.5"

import argparse
import logging
import os
import re
import sys
import time
import warnings
from datetime import datetime, timezone

import numpy as np

import matplotlib
# Only set the backend when no backend has been selected yet.
# When imported by CometTailGUI.py, QtAgg is already active; this guard
# prevents a backend-conflict warning and ensures the GUI renders correctly.
if not matplotlib.get_backend():
    matplotlib.use('TkAgg' if os.environ.get('DISPLAY') else 'Agg')
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
#  CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────
MU   = 2.9591220828e-4          # GM_sun  (AU³ day⁻²)
OBL  = np.radians(23.43929111)  # J2000 obliquity
J2K  = 2451545.0                # JD of J2000.0

# β → grain radius coefficient, derived DIRECTLY from Burns, Lamy & Soter
# (1979, Icarus 40, 1) Eq. 19:
#     β ≡ F_rad/F_grav = (3·L_sun·Qpr) / (16·π·c·G·M_sun·ρ·s) = 5.7×10⁻⁵ Qpr/(ρs)   [cgs: ρ in g/cm³, s in cm]
# Converting s from cm to µm (×10⁴):
#     β = 0.574 · Qpr / (ρ[g/cm³] · s[µm])         ⇒  s[µm] = 0.574 · Qpr / (β · ρ)
#
# NOTE (v3.0): earlier versions used a different two-parameter form,
#     β = C_pr·Qpr/(ρr) with C_pr = 1.19×10⁻³ kg/m²  (Finson & Probstein 1968a
#     convention, as quoted in e.g. Moreno 2025 A&A 695 A263 and
#     Mariblanca-Escalona et al. 2026 — both of whom pair this C_pr with an
#     explicit /2 in the denominator: β = C_pr·Qpr/(2ρr)).
# v2.4 omitted that factor of 2, giving grain radii 2.07× too LARGE.
# v3.0 sidesteps the C_pr/2 ambiguity entirely by deriving the coefficient
# straight from Eq. 19's own numerical result, verified independently from
# first principles against the primary source.
BETA_TO_RADIUS_UM = 0.574       # µm · g·cm⁻³,  Qpr = 1 (Burns et al. 1979, Eq. 19)

# ─────────────────────────────────────────────────────────────────────────────
#  COMET DATABASE  (30+ objects)
#  All elements: q(AU), e, i(°), Omega(°), omega(°), T(YYYY-MM-DD UT)
# ─────────────────────────────────────────────────────────────────────────────

# COMET_DB is a curated catalogue of interesting/historical comets.
# It stores ONLY editorial metadata:
#   obs  — the preferred observation date (peak brightness, best geometry, or
#           scientifically notable moment) — pre-filled as the default date
#           when the user picks this comet from the PRESET list.
#   note — short human-readable description shown in the GUI.
#
# Orbital elements are NOT stored here. They are always fetched live from
# JPL Horizons so they are never stale or wrong.
# The PRESET tab simply pre-fills the comet name + preferred date, then
# triggers the same Horizons fetch as the FETCH JPL tab.
COMET_DB = {
    # ── Great / recent naked-eye comets ──────────────────────────────────────
    "C/2025 R3 (PANSTARRS)":            dict(obs="2026-04-07", note="Best comet of 2026"),
    "C/2023 A3 (Tsuchinshan-ATLAS)":    dict(obs="2024-10-10", note="Spectacular Oct 2024 comet, mag ≈ −4 at elongation"),
    "C/2021 A1 (Leonard)":              dict(obs="2021-12-15", note="Best naked-eye comet of 2021"),
    "C/2020 F3 (NEOWISE)":              dict(obs="2020-07-23", note="Spectacular naked-eye comet, July 2020"),
    "C/2022 E3 (ZTF)":                  dict(obs="2023-02-01", note="Green comet, Jan–Feb 2023, slightly hyperbolic"),
    "C/2020 S3 (Erasmus)":              dict(obs="2020-11-25", note="Southern-hemisphere comet, Nov 2020"),
    "C/2019 Y4 (ATLAS)":                dict(obs="2020-04-06", note="Disintegrated before perihelion, Apr 2020"),
    "C/2019 U6 (Lemmon)":               dict(obs="2020-06-10", note="Southern comet, June 2020"),
    "C/2018 Y1 (Iwamoto)":              dict(obs="2019-02-12", note="Very close approach 0.28 AU from Earth"),
    "C/2017 T2 (PanSTARRS)":            dict(obs="2020-05-01", note="Northern circumpolar in spring 2020"),
    "C/2016 R2 (PanSTARRS)":            dict(obs="2018-01-14", note="Blue comet — CO-dominated coma, unusual spectrum"),
    "C/2015 V2 (Johnson)":              dict(obs="2017-05-10", note="Northern spring 2017"),
    "C/2015 ER61 (PanSTARRS)":          dict(obs="2017-04-15", note="Outburst April 2017"),
    "C/2014 Q2 (Lovejoy)":              dict(obs="2015-01-15", note="Naked-eye Jan 2015, alcohol detected"),
    "C/2013 A1 (Siding Spring)":        dict(obs="2014-10-19", note="Passed 140,000 km from Mars"),
    "C/2012 S1 (ISON)":                 dict(obs="2013-11-15", note="Sungrazer — disintegrated at perihelion"),
    "C/2011 L4 (PanSTARRS)":            dict(obs="2013-03-20", note="Bright southern comet, early 2013"),
    "C/2009 P1 (Garradd)":              dict(obs="2012-02-10", note="Displayed both gas and dust tails simultaneously"),
    "C/2007 N3 (Lulin)":                dict(obs="2009-02-24", note="Retrograde near-parabolic, dual opposing tails"),
    "C/2006 P1 (McNaught)":             dict(obs="2007-01-14", note="Brightest since 1965, mag ≈ −5.5, striated tail"),
    "C/2001 Q4 (NEAT)":                 dict(obs="2004-05-07", note="Northern naked-eye May 2004"),
    "C/1999 T1 (McNaught-Hartley)":     dict(obs="2001-11-20", note="Long-period, good dust tail"),
    # ── Classic / historical comets ──────────────────────────────────────────
    "C/1995 O1 (Hale-Bopp)":            dict(obs="1997-04-10", note="Great Comet of 1997, visible 18 months"),
    "C/1996 B2 (Hyakutake)":            dict(obs="1996-03-25", note="Closest comet approach (0.102 AU) since 1770"),
    "C/1993 A1 (Mueller)":              dict(obs="1993-11-01", note="Distant comet, q > 4 AU"),
    "C/1990 K1 (Levy)":                 dict(obs="1990-09-15", note="Best comet of 1990"),
    "C/1975 V1 (West)":                 dict(obs="1976-03-05", note="Fragmented nucleus, spectacular dust tail"),
    "C/1956 R1 (Arend-Roland)":         dict(obs="1957-04-25", note="Original F-P test case — famous anti-tail"),
    # ── Periodic comets ──────────────────────────────────────────────────────
    "1P/Halley":                         dict(obs="1986-03-15", note="Most famous periodic, P≈76yr"),
    "2P/Encke":                          dict(obs="2004-01-15", note="Shortest known period P≈3.3yr"),
    "17P/Holmes":                        dict(obs="2007-11-01", note="Dramatic million-fold outburst 2007"),
    "81P/Wild 2":                        dict(obs="2004-01-02", note="Stardust mission target"),
    "67P/Churyumov-Gerasimenko":         dict(obs="2015-08-01", note="Rosetta/Philae mission — rubber-duck nucleus"),
    "9P/Tempel 1":                       dict(obs="2005-07-04", note="Deep Impact mission target"),
}

# ─────────────────────────────────────────────────────────────────────────────
#  DATE / JD UTILITIES
# ─────────────────────────────────────────────────────────────────────────────
# Accepted date formats (all are converted to JD with full sub-day precision):
#   YYYY-MM-DD              plain calendar date (midnight UT)
#   YYYY-MM-DD.fff          decimal day  ← astronomical standard, e.g. 2018-01-14.583
#   YYYY-MM-DDTHH:MM        ISO-8601 with T separator
#   YYYY-MM-DD HH:MM        space separator
#   YYYY-MM-DD HH:MM:SS     with seconds
# The decimal-day and HH:MM forms are mutually exclusive; decimal day takes priority.
_DATE_RE = re.compile(
    r'(\d{4})-(\d{1,2})-(\d{1,2})'          # YYYY-MM-DD  (groups 1-3)
    r'(?:'
        r'(\.\d+)'                            # .fraction   (group 4)  — decimal day
        r'|[T\s](\d{1,2}):(\d{2})(?::(\d{2})(?:\.\d+)?)?'  # HH:MM[:SS] (groups 5-7)
    r')?'
)


def date_to_jd(s: str) -> float:
    """Convert a date string to Julian Date (UTC).

    Accepted formats::
        '2018-01-14'           → JD at 00:00 UT
        '2018-01-14.583'       → JD at fractional day (astronomical notation)
        '2018-01-14 13:55'     → JD at HH:MM UT
        '2018-01-14T13:55'     → ISO-8601
        '2018-01-14 13:55:12'  → HH:MM:SS UT
    """
    m = _DATE_RE.match(s.strip())
    if not m:
        raise ValueError(
            f"Cannot parse date: {s!r}\n"
            f"  Accepted: 'YYYY-MM-DD', 'YYYY-MM-DD.fff', 'YYYY-MM-DD HH:MM', "
            f"'YYYY-MM-DD HH:MM:SS'"
        )
    y, mo = int(m[1]), int(m[2])
    d_int = int(m[3])
    if m[4]:                                  # decimal day, e.g. .583
        day_frac = float(m[4])               # e.g. 0.583
    elif m[5]:                                # HH:MM[:SS]
        h  = int(m[5])
        mn = int(m[6])
        sc = float(m[7]) if m[7] else 0.0
        day_frac = h / 24.0 + mn / 1440.0 + sc / 86400.0
    else:
        day_frac = 0.0

    if mo <= 2:
        y  -= 1
        mo += 12
    A = y // 100
    B = 2 - A + A // 4
    return (int(365.25 * (y + 4716))
            + int(30.6001 * (mo + 1))
            + d_int + B - 1524.5
            + day_frac)


def jd_to_str(jd: float) -> str:
    """Convert JD to human-readable UTC string.

    Returns 'YYYY-MM-DD.fff UT' (decimal day) when sub-day precision is
    present, otherwise 'YYYY-MM-DD HH:MM UT'.
    """
    J = jd + 0.5
    Z = int(J)
    F = J - Z            # day fraction (0 = midnight, 0.5 = noon)
    A = Z
    if Z >= 2299161:
        a = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + a - a // 4
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)
    day   = B - D - int(30.6001 * E)
    month = E - 1 if E < 14 else E - 13
    year  = C - 4716 if month > 2 else C - 4715
    if F < 1e-5:          # exactly midnight — plain date
        return f"{year:04d}-{month:02d}-{day:02d} UT"
    # Decimal-day notation: e.g. 2018-01-14.583 UT
    # Round to 4 decimal places (≈ 8.6 s precision — sufficient for comets)
    frac_str = f"{F:.4f}"[1:]   # e.g. ".5830"
    return f"{year:04d}-{month:02d}-{day:02d}{frac_str} UT"


def today_jd() -> float:
    """Current JD in UTC."""
    dt = datetime.now(tz=timezone.utc)
    return date_to_jd(dt.strftime('%Y-%m-%d %H:%M'))

# ─────────────────────────────────────────────────────────────────────────────
#  JPL HORIZONS FETCHING
# ─────────────────────────────────────────────────────────────────────────────

_HTTP_HEADERS = {
    "User-Agent": "CometTailAnalyzer/2.5",
    "Accept":     "application/json",
}


def _normalize_comet_desig(designation: str):
    """
    Return a list of (id_string, id_type) pairs to try with JPL Horizons.

    Horizons accepts periodic comets as:
      - "88P"           with id_type='designation'
      - "88P/Howell"    with id_type='comet_name'
    Long-period comets MUST have the discoverer suffix stripped:
      - "C/2023 A3"     with id_type='designation'   ← Horizons packed form
      - "C/2016 R2 (PanSTARRS)" → strip → "C/2016 R2" before sending

    Horizons does NOT accept "(PanSTARRS)" / "(Hale-Bopp)" etc. as part of
    a designation query — it will return no results.  The full string is kept
    only as a comet_name fallback.
    """
    d = designation.strip()
    attempts: list[tuple[str, str | None]] = []

    if re.match(r'^\d+[PD]', d):
        # Periodic comet: starts with "NNP" or "NND"
        num_part = re.match(r'^(\d+[PD])', d).group(1)   # e.g. "88P"
        attempts.append((num_part, 'designation'))
        attempts.append((num_part, None))
        if '/' in d:
            full       = d                                # "88P/Howell"
            name_part  = d.split('/', 1)[1].strip()       # "Howell"
            attempts.append((full,      'comet_name'))
            attempts.append((name_part, 'comet_name'))
        else:
            attempts.append((d, 'comet_name'))
    else:
        # Long-period / other: strip discoverer name in parentheses
        # "C/2016 R2 (PanSTARRS)"         → "C/2016 R2"
        # "C/2023 A3 (Tsuchinshan-ATLAS)" → "C/2023 A3"
        # "C/2023 A3"                     → unchanged
        clean = re.sub(r'\s*\([^)]*\)\s*$', '', d).strip()
        if clean and clean != d:
            attempts.append((clean, 'designation'))   # correct packed form — try first
            attempts.append((clean, None))
        attempts.append((d, 'designation'))           # original as fallback
        attempts.append((d, 'comet_name'))
        attempts.append((d, None))
        if clean and clean != d:
            attempts.append((clean, 'comet_name'))

    return attempts


def _parse_ambiguous_record(error_text: str) -> str | None:
    """
    Extract the most-recent record ID from a Horizons 'Ambiguous target'
    error message. Returns None if no records are found.
    """
    matches = re.findall(r'(\d{6,10})\s+(\d{4})', error_text)
    if not matches:
        return None
    # Pick the record with the most recent epoch year
    return max(matches, key=lambda x: int(x[1]))[0]


def fetch_from_horizons(designation: str, date: str | None = None) -> dict:
    """
    Fetch comet orbital elements from JPL Horizons.

    Handles ambiguous periodic-comet targets (picks latest epoch) and both
    periodic (e.g. '88P') and long-period (e.g. 'C/2023 A3') comets.
    """
    try:
        from astroquery.jplhorizons import Horizons
    except ImportError:
        raise ImportError("astroquery not installed. Run: pip install astroquery")

    epoch_jd  = date_to_jd(date) if date else today_jd()
    last_exc: Exception | None = None

    def _query(id_str, id_type):
        kw = dict(id=id_str, location='@10', epochs=epoch_jd)
        if id_type:
            kw['id_type'] = id_type
        el = Horizons(**kw).elements()
        if len(el) == 0:
            raise RuntimeError("Empty result")
        return el

    def _extract(el, name=None) -> dict:
        Tp = float(el['Tp_jd'][0])
        return dict(
            q     = float(el['q'][0]),
            e     = float(el['e'][0]),
            i     = float(el['incl'][0]),
            Omega = float(el['Omega'][0]),
            omega = float(el['w'][0]),
            T     = jd_to_str(Tp)[:10],
            T_jd  = Tp,
            name  = name or designation,
            source= 'JPL Horizons',
        )

    for id_str, id_type in _normalize_comet_desig(designation):
        try:
            return _extract(_query(id_str, id_type))
        except Exception as ex:
            last_exc = ex
            if 'ambiguous' in str(ex).lower():
                rec_id = _parse_ambiguous_record(str(ex))
                if rec_id:
                    try:
                        return _extract(_query(rec_id, None), name=designation)
                    except Exception as ex2:
                        last_exc = ex2

    raise RuntimeError(
        f"Horizons could not find: {designation!r}\n"
        f"Last error: {last_exc}\n\n"
        f"Suggestions:\n"
        f"  • For periodic comets: try '88P/Howell' (include name after /)\n"
        f"  • Check internet connection\n"
        f"  • Verify at: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html"
    )


def fetch_from_cobs(designation: str, comet_el: dict | None = None) -> dict:
    """
    Fetch H₀ and n by combining:
      1. COBS obs_list.api  → date + magnitude for each observation
      2. Horizons ephemeris → r_helio, delta for those dates
         (or analytic propagation from comet_el as a fallback)
      3. Least-squares fit  → H₀, n  (active zone only, r ≤ r_max_fit)
    Final fallback: Horizons stored M1/k1 from elements().

    Returns a dict with keys:
        H0, n, H0_fit, n_fit, last_mag, last_date, last_r,
        obs_list, raw_obs, n_fit_pts, r_max_fit, source, comet_name
    """
    import requests

    comet_el  = comet_el or {}
    q_peri    = comet_el.get('q', 2.0)
    r_max_fit = min(3.5, max(q_peri * 2.5, 2.0))

    COBS_API = "https://cobs.si/api"

    # ── helpers ──────────────────────────────────────────────────────────
    def _clean_des(name: str) -> list[str]:
        d = name.strip()
        variants = [d.split('(')[0].strip()]
        if '/' in d:
            variants.append(d.split('/')[0].strip())
        return list(dict.fromkeys(variants))

    def _fit_H0_n(obs):
        """
        Least-squares fit of  m = H₀ + 5·log10(Δ) + 2.5·n·log10(r)
        on the 'active zone' (r ≤ r_max_fit, mag ≤ 19).
        Falls back to all valid observations if too few are in the active zone.
        Returns (H₀, n, n_err, points_used) — H₀, n, n_err are None if the fit fails.
        n_err is the 1-σ standard error on n from the fit's residuals
        (None if there are too few degrees of freedom to estimate it).
        """
        valid = [o for o in obs
                 if 1.0 <= o['mag'] <= 19.0
                 and o['r_helio'] > 0.1 and o['delta'] > 0.01
                 and o['r_helio'] <= r_max_fit]
        pts = (valid if len(valid) >= 5
               else [o for o in obs
                     if 1.0 <= o['mag'] <= 20.0 and o['r_helio'] > 0.1])
        if len(pts) < 5:
            return None, None, None, pts
        m  = np.array([o['mag']               for o in pts])
        lr = np.array([np.log10(o['r_helio']) for o in pts])
        ld = np.array([np.log10(o['delta'])   for o in pts])
        if np.std(lr) < 0.01:
            return None, None, None, pts
        y    = m - 5.0 * ld
        A    = np.column_stack([np.ones(len(y)), lr])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        n_val = float(coef[1]) / 2.5
        if not (0.5 <= n_val <= 25):
            return None, None, None, pts
        # 1-σ standard error on the slope (→ on n), from the fit residuals.
        # dof = N points − 2 free parameters (H0, slope).
        n_err = None
        dof = len(y) - 2
        if dof > 0:
            resid  = y - A @ coef
            sigma2 = float(np.sum(resid ** 2)) / dof
            try:
                cov = sigma2 * np.linalg.inv(A.T @ A)
                n_err = round(float(np.sqrt(cov[1, 1])) / 2.5, 2)
            except np.linalg.LinAlgError:
                n_err = None
        return round(float(coef[0]), 2), round(n_val, 2), n_err, pts


    # ── 1. COBS observations (date + magnitude) ──────────────────────────
    raw_obs: list[dict] = []
    for des in _clean_des(designation):
        try:
            for page in range(1, 12):
                r = requests.get(f"{COBS_API}/obs_list.api",
                                 params={"des": des, "format": "json", "page": page},
                                 headers=_HTTP_HEADERS, timeout=10)
                if r.status_code != 200:
                    break
                data = r.json()
                recs = data.get("objects", [])
                if not recs:
                    break
                for rec in recs:
                    try:
                        date_str = str(rec.get("obs_date", ""))[:10]
                        mag = float(rec.get("magnitude", rec.get("mag", 99)))
                        if date_str and 1.0 <= mag <= 22.0:
                            raw_obs.append({"date": date_str, "mag": mag})
                    except Exception:
                        pass
                info = data.get("info", {})
                if info.get("page", page) >= info.get("pages", page):
                    break
            if raw_obs:
                break
        except Exception:
            continue

    # ── 2. Get r_helio, delta for each observation ───────────────────────
    obs_with_eph: list[dict] = []
    eph_source = ""
    if raw_obs:
        try:
            from astroquery.jplhorizons import Horizons

            dates_sorted = sorted({o["date"] for o in raw_obs})
            eph_epochs = dict(start=dates_sorted[0], stop=dates_sorted[-1], step='5d')

            eph_table = None
            for id_str, id_type in _normalize_comet_desig(designation):
                try:
                    kw = dict(id=id_str, location='500', epochs=eph_epochs)
                    if id_type:
                        kw['id_type'] = id_type
                    eph_table = Horizons(**kw).ephemerides(quantities='19,20')
                    if len(eph_table) > 0:
                        break
                except Exception as ex:
                    if 'ambiguous' in str(ex).lower():
                        rec_id = _parse_ambiguous_record(str(ex))
                        if rec_id:
                            try:
                                eph_table = Horizons(
                                    id=rec_id, location='500',
                                    epochs=eph_epochs).ephemerides(quantities='19,20')
                                if len(eph_table) > 0:
                                    break
                            except Exception:
                                pass
                    eph_table = None

            if eph_table is not None and len(eph_table) > 0:
                eph_jd_map = {}
                for row in eph_table:
                    try:
                        jd  = float(row['datetime_jd'])
                        eph_jd_map[round(jd, 1)] = (float(row['r']), float(row['delta']))
                    except Exception:
                        pass
                # Match each obs to nearest ephemeris point (within ±3.5 d)
                for obs in raw_obs:
                    try:
                        obs_jd = round(date_to_jd(obs["date"]), 1)
                        best   = min(eph_jd_map.keys(), key=lambda jd: abs(jd - obs_jd))
                        if abs(best - obs_jd) <= 3.5:
                            r_h, dlt = eph_jd_map[best]
                            obs_with_eph.append(dict(
                                date=obs["date"], mag=obs["mag"],
                                r_helio=r_h, delta=dlt))
                    except Exception:
                        continue
                if obs_with_eph:
                    eph_source = "Horizons eph"
        except Exception:
            pass

    # ── Analytic fallback: propagate orbit from comet_el ────────────────
    if not obs_with_eph and raw_obs and comet_el.get('q'):
        try:
            for obs in raw_obs:
                try:
                    obs_jd = date_to_jd(obs["date"])
                    r_C, _ = elem_to_state(comet_el, obs_jd)
                    r_E    = earth_pos(obs_jd)
                    obs_with_eph.append(dict(
                        date=obs["date"], mag=obs["mag"],
                        r_helio=float(np.linalg.norm(r_C)),
                        delta  =float(np.linalg.norm(r_C - r_E))))
                except Exception:
                    continue
            if obs_with_eph:
                eph_source = "analytic eph"
        except Exception:
            pass

    # ── 3. Fit H₀, n ─────────────────────────────────────────────────────
    if obs_with_eph:
        H0_fit, n_fit, n_err, pts = _fit_H0_n(obs_with_eph)
        if H0_fit is not None:
            last = max(obs_with_eph, key=lambda o: o['date'])
            return dict(
                H0=H0_fit, n=n_fit, n_err=n_err,
                H0_fit=H0_fit, n_fit=n_fit,
                last_mag   = last.get("mag"),
                last_date  = last.get("date", ""),
                last_r     = last.get("r_helio"),
                obs_list   = obs_with_eph,
                raw_obs    = raw_obs,
                n_fit_pts  = len(pts),
                r_max_fit  = r_max_fit,
                source     = f"COBS {len(raw_obs)} obs → {eph_source} fit",
                comet_name = designation,
            )

    # Raw observations available, but fit not possible
    if raw_obs:
        last = max(raw_obs, key=lambda o: o['date'])
        return dict(
            H0=None, n=None, n_err=None, H0_fit=None, n_fit=None,
            last_mag=last.get("mag"), last_date=last.get("date", ""),
            last_r=None, obs_list=[], raw_obs=raw_obs,
            n_fit_pts=0, r_max_fit=r_max_fit,
            source=f"COBS {len(raw_obs)} obs (ephemeris/fit unavailable)",
            comet_name=designation,
        )

    # ── Final fallback: Horizons M1/k1 ──────────────────────────────────
    try:
        from astroquery.jplhorizons import Horizons
        from math import isnan

        obs_jd = today_jd()
        for id_str, id_type in _normalize_comet_desig(designation):
            try:
                kw = dict(id=id_str, location='@10', epochs=obs_jd)
                if id_type:
                    kw['id_type'] = id_type
                el = Horizons(**kw).elements()
            except Exception as ex:
                if 'ambiguous' in str(ex).lower():
                    rec_id = _parse_ambiguous_record(str(ex))
                    if rec_id:
                        try:
                            el = Horizons(id=rec_id, location='@10',
                                          epochs=obs_jd).elements()
                        except Exception:
                            continue
                    else:
                        continue
                else:
                    continue

            if 'M1' in el.colnames and 'k1' in el.colnames:
                M1v, k1v = float(el['M1'][0]), float(el['k1'][0])
                if not (isnan(M1v) or isnan(k1v)):
                    return dict(
                        H0=round(M1v, 2), n=round(k1v, 2), n_err=None,
                        H0_fit=None, n_fit=None,
                        last_mag=None, last_date='', last_r=None,
                        obs_list=[], raw_obs=[],
                        n_fit_pts=0, r_max_fit=r_max_fit,
                        source="Horizons M1/k1 (⚠ may differ from COBS)",
                        comet_name=designation,
                    )
    except Exception:
        pass

    raise RuntimeError(
        f"Light curve unavailable for {designation!r}.\n"
        f"Tip: Check internet connection or enter H₀/n manually."
    )


def compute_dust_production_rate(afrho_cm: float, r_h: float,
                                 v_dust_kms: float, p_v: float = 0.04) -> dict:
    """
    Rough Afρ-based dust production rate estimate — order-of-magnitude
    only, NOT a substitute for full coma photometry/Mie-scattering
    modeling (cf. Section 3.1.2-3.1.3 of Schambeau et al. 2021 for the
    rigorous version). Shared by generate_dust_analysis()'s inline Afρ
    section and the standalone Calculation > Dust production rate…
    calculator, so the two can never silently drift apart.

        Afρ_norm(1AU) = Afρ_obs · r_h²
        Q_d ≈ Afρ_norm · v_dust / (2·p_v)        [kg/s, v_dust in km/s
                                                    converted to cm/s via
                                                    the 1e5 factor below]

    Returns dict(afrho_norm, Qd_rough, activity).
    """
    afrho_norm = afrho_cm * r_h ** 2
    Qd_rough   = afrho_norm * v_dust_kms * 1e5 / (p_v * 2)   # kg/s
    if   Qd_rough > 1e4: act = "VERY HIGH (outburst / hyperactive candidate)"
    elif Qd_rough > 1e3: act = "HIGH (active near perihelion)"
    elif Qd_rough > 1e2: act = "MODERATE (typical active comet)"
    else:                act = "LOW (weakly active or distant)"
    return dict(afrho_norm=afrho_norm, Qd_rough=Qd_rough, activity=act)


def generate_dust_analysis(comet_el: dict, model_info: dict,
                           afrho_cm: float | None = None,
                           v_dust_kms: float | None = None,
                           p_v: float = 0.04) -> str:
    """
    Generate a plain-text preliminary dust analysis report.
    Combines F-P results and Afρ/dust-production interpretation.

    Two sections that used to live in this report were moved out to their
    own standalone tools/tabs in v3.0, since each is self-contained and
    duplicating them here risked them silently drifting apart from the
    place a user would actually go to refresh them:
      - F-P grain radius → Calculation > Dust particle radius… (needs
        only β/ρ_d/Qpr, not a per-comet model; beta_to_size() remains its
        single source of truth, shared with the β TABLE's radius column).
      - COBS light-curve fit (H₀/n, m_pred vs. m_obs outburst flagging)
        → the LIGHT CURVE tab/popup window, which already shows it.
    This function no longer takes a cobs_data argument as a result.

    Physical parameters (all user-adjustable directly in the relevant
    standalone calculator under the GUI's Calculation menu; defaults shown
    are literature values, not universal constants):
        v_dust_kms dust expansion velocity [km/s] used only for the rough
                   Q_d (Afρ-based) production-rate estimate. If None
                   (default), falls back to the empirical scaling law
                   v = 0.1·r_h^-0.5 km/s (uncited rule of thumb — override
                   with a measured value when available, e.g. quiescent
                   ~0.01-0.05 km/s vs. outburst ~0.1-0.3 km/s for 29P,
                   Schambeau et al. 2017, 2019).
        p_v        visible geometric albedo, default 0.04 (typical bare
                   cometary nucleus; SEPPCoN/Schambeau et al. 2021 NEATM
                   assumption).
    """
    r_h     = model_info.get("r_helio", 1.0)
    r_geo   = model_info.get("r_geo",   1.0)
    phase   = model_info.get("phase_angle", 0.0)
    T_jd    = comet_el.get("T_jd", 0)
    obs_jd  = model_info.get("obs_jd", T_jd)
    dt_peri = obs_jd - T_jd
    name    = comet_el.get("name", "Comet")
    q       = comet_el.get("q", r_h)
    e       = comet_el.get("e", 1.0)
    i_deg   = comet_el.get("i", 0.0)
    P_yr    = comet_el.get("P_yr")

    lines: list[str] = []
    lines.append(f"═══ DUST ANALYSIS  —  {name} ═══")
    lines.append(f"Date: {model_info.get('obs_str','')[:10]}  "
                 f"r☉={r_h:.4f} AU  Δ={r_geo:.4f} AU  Phase={phase:.1f}°")
    lines.append(f"T_peri offset: {dt_peri:+.0f} d  "
                 f"({'post' if dt_peri>0 else 'pre'}-perihelion)")
    lines.append("")

    # ── Orbital classification ──────────────────────────────────────────
    lines.append("── ORBIT ──────────────────────────────────")
    if e < 0.98:
        if P_yr is None:
            try:
                a = q / (1 - e)
                P_yr = a ** 1.5
            except Exception:
                P_yr = 1.0
        if P_yr < 20:
            family = "Jupiter Family Comet (JFC)"
        elif P_yr < 200:
            family = "Halley-type Comet (HTC)"
        else:
            family = "Long-period Comet (LPC)"
    elif e < 1.0:
        family = "Long-period Comet / Near-parabolic"
    elif e < 1.01:
        family = "Near-parabolic / Dynamically new"
    else:
        family = f"Hyperbolic / Interstellar candidate  (e={e:.4f})"
    lines.append(f"  Type    : {family}")
    lines.append(f"  q={q:.4f} AU  e={e:.6f}  i={i_deg:.2f}°")

    if r_h > 5:
        ice = "CO / CO₂ sublimation (r > 5 AU)"
    elif r_h > 3:
        ice = "CO₂ sublimation (3–5 AU)"
    else:
        ice = "H₂O sublimation (r < 3 AU)"
    lines.append(f"  Ice drv : {ice}")
    lines.append("")

    # NOTE (v3.0): the F-P grain-radius (β→r) section that used to live here
    # was moved out to its own standalone calculator — Calculation >
    # Dust particle radius… — since the formula needs only (β, ρ_d, Qpr),
    # not a computed model, and doesn't belong mixed into a per-comet
    # analysis report. beta_to_size() remains the single source of truth
    # for the formula; both that calculator and the β TABLE's radius
    # column call it directly.



    # ── Afρ / Dust production ───────────────────────────────────────────
    if afrho_cm and afrho_cm > 0:
        lines.append("── Afρ / DUST PRODUCTION ───────────────────")
        lines.append(f"  Afρ (obs)     = {afrho_cm:.0f} cm  at r={r_h:.3f} AU")
        # Rough Q_d estimate
        if v_dust_kms is not None:
            v_dust  = v_dust_kms
            v_label = f"{v_dust:.3f} km/s  (user-set)"
        else:
            v_dust  = 0.1 * r_h ** (-0.5)        # km/s, empirical r^-0.5 law
            v_label = f"{v_dust:.3f} km/s  (r^-0.5 law, no fixed reference — override if known)"
        qd = compute_dust_production_rate(afrho_cm, r_h, v_dust, p_v)
        lines.append(f"  Afρ norm(1AU) = {qd['afrho_norm']:.0f} cm")
        lines.append(f"  v_dust est.   ≈ {v_label}")
        lines.append(f"  p_v (albedo)  = {p_v:g}")
        lines.append(f"  Q_d (rough)   ≈ {qd['Qd_rough']:.0e} kg/s")
        lines.append(f"  Activity      : {qd['activity']}")
        lines.append("")

    # NOTE (v3.0): the COBS light-curve section (H₀/n fit, m_pred vs.
    # m_obs outburst-candidate flagging) that used to print here was
    # removed — that information belongs to, and is already shown in,
    # the dedicated LIGHT CURVE tab/popup window, so duplicating it in
    # the text report risked the two silently disagreeing if one was
    # refreshed and the other wasn't. This function no longer takes a
    # cobs_data argument at all as a result. fetch_from_cobs()'s
    # n_err/n_fit_pts robustness fields remain available on whatever the
    # caller's own cobs_data dict is, for the LIGHT CURVE tab itself to
    # surface a "preliminary fit" warning if it wants one.



    # ── Summary ─────────────────────────────────────────────────────────
    lines.append("── SUMMARY ────────────────────────────────")
    summary = f"{name} observed at r☉={r_h:.3f} AU"
    summary += (f", {dt_peri:.0f}d post-perihelion." if dt_peri > 0
                else f", {abs(dt_peri):.0f}d pre-perihelion.")
    if afrho_cm:
        summary += f" Afρ={afrho_cm:.0f} cm"
    if r_h > 5:
        summary += (" Activity at this distance is driven by supervolatile"
                    " ices (CO/CO₂), suggesting a pristine nucleus.")
    elif r_h > 3:
        summary += " CO₂ sublimation likely dominates at this distance."
    else:
        summary += " Classical H₂O sublimation active."
    lines.append(f"  {summary}")
    lines.append("")
    lines.append("  [Comet Tail Analyzer by Teerasak Thaluang]")
    return "\n".join(lines)


def fetch_comet(name_or_desig: str, date: str | None = None) -> dict:
    """
    Fetch live orbital elements from JPL Horizons for any comet.

    COMET_DB stores only editorial metadata (preferred date + note); there
    is no local element cache that could return stale or incorrect data.

    Args:
        name_or_desig: Any recognisable comet designation, e.g.:
                       "C/2023 A3", "C/2023 A3 (Tsuchinshan-ATLAS)",
                       "1P/Halley", "67P"
        date:          Observation date string (any format accepted by
                       date_to_jd).  Used to select the Horizons epoch.
                       Defaults to today if omitted.

    Returns:
        dict with keys: q, e, i, Omega, omega, T, T_jd, name, source

    Raises:
        RuntimeError if Horizons cannot find the object.
    """
    logging.info("  Fetching '%s' from JPL Horizons...", name_or_desig)
    try:
        el = fetch_from_horizons(name_or_desig, date)
        logging.info("  [Horizons] q=%.5f e=%.6f i=%.3f°", el['q'], el['e'], el['i'])
        return el
    except Exception as ex:
        is_periodic = bool(re.match(r'^\d+[PD]', name_or_desig.strip()))
        hint = ""
        if is_periodic:
            hint = (f"\n  Tip: For periodic comets, try including the name, e.g. "
                    f"'{name_or_desig.split('/')[0]}/Howell'")
        raise RuntimeError(
            f"Could not find orbital elements for: {name_or_desig!r}\n"
            f"Horizons error: {ex}{hint}"
        ) from ex


# ─────────────────────────────────────────────────────────────────────────────
#  VECTOR MATH
# ─────────────────────────────────────────────────────────────────────────────
def vmag(v):  return np.linalg.norm(v, axis=0) if v.ndim > 1 else np.linalg.norm(v)
def vhat(v):  m = vmag(v); return v / m if m > 1e-15 else v * 0
def vcrs(a, b): return np.cross(a, b)
def vdot(a, b): return np.dot(a, b)


# ─────────────────────────────────────────────────────────────────────────────
#  KEPLER SOLVERS
# ─────────────────────────────────────────────────────────────────────────────
def kepler_elliptic(M: float, e: float) -> float:
    """Solve Kepler's equation  M = E − e·sin E  (elliptic)."""
    M = ((M % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
    if M > np.pi:
        M -= 2 * np.pi
    E = M + e * np.sin(M)
    for _ in range(60):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E += dE
        if abs(dE) < 1e-13:
            break
    return E


def kepler_hyperbolic(M: float, e: float) -> float:
    """Solve hyperbolic Kepler's equation  M = e·sinh H − H."""
    H = np.sign(M or 1.0) * np.log(2 * abs(M) / e + 1.8)
    for _ in range(60):
        f  = e * np.sinh(H) - H - M
        fp = e * np.cosh(H) - 1.0
        if abs(fp) < 1e-15:
            break
        H -= f / fp
        if abs(f / fp) < 1e-13:
            break
    return H


# ─────────────────────────────────────────────────────────────────────────────
#  ORBITAL ELEMENTS  ↔  STATE VECTOR  (heliocentric ecliptic J2000)
# ─────────────────────────────────────────────────────────────────────────────
def elem_to_state(el: dict, jd: float, mu: float = MU):
    """
    Convert orbital elements (q, e, i, Omega, omega, T or T_jd)
    to a state vector (r, v) at the given JD, in heliocentric ecliptic J2000.
    """
    q, e = el['q'], el['e']
    iR, OR, oR = np.radians(el['i']), np.radians(el['Omega']), np.radians(el['omega'])
    T  = el.get('T_jd', date_to_jd(el['T']))

    cO, sO = np.cos(OR), np.sin(OR)
    cI, sI = np.cos(iR), np.sin(iR)
    co, so = np.cos(oR), np.sin(oR)

    # Perifocal axes in ecliptic
    P = np.array([cO*co - sO*so*cI,    sO*co + cO*so*cI,    so*sI])
    Q = np.array([-(cO*so + sO*co*cI), -(sO*so - cO*co*cI), co*sI])

    dt = jd - T
    eps_e = 1e-5

    if abs(e - 1.0) < eps_e:
        # Parabolic — Barker's equation
        W  = 3.0 * np.sqrt(mu / (2 * q ** 3)) * dt
        s  = np.cbrt(W + np.sqrt(W * W + 1))
        tf = s - 1.0 / s
        r  = q * (1.0 + tf * tf)
        h  = np.sqrt(2 * mu * q)
        f  = 2.0 * np.arctan(tf)
        xo, yo   = r * np.cos(f), r * np.sin(f)
        vxo, vyo = -mu / h * np.sin(f), mu / h * (1 + np.cos(f))
    elif e < 1.0:
        # Elliptic
        a   = q / (1 - e)
        n   = np.sqrt(mu / a ** 3)
        E   = kepler_elliptic(n * dt, e)
        cE, sE = np.cos(E), np.sin(E)
        rv  = a * (1 - e * cE)
        cf  = (cE - e) / (1 - e * cE)
        sf  = np.sqrt(1 - e * e) * sE / (1 - e * cE)
        h   = np.sqrt(mu * a * (1 - e * e))
        xo, yo   = rv * cf, rv * sf
        vxo, vyo = -mu / h * sf, mu / h * (e + cf)
    else:
        # Hyperbolic
        ah  = q / (e - 1)
        n   = np.sqrt(mu / ah ** 3)
        H   = kepler_hyperbolic(n * dt, e)
        cH, sH = np.cosh(H), np.sinh(H)
        rv  = ah * (e * cH - 1)
        cf  = (e - cH) / (e * cH - 1)
        sf  = np.sqrt(e * e - 1) * sH / (e * cH - 1)
        h   = np.sqrt(mu * ah * (e * e - 1))
        xo, yo   = rv * cf, rv * sf
        vxo, vyo = -mu / h * sf, mu / h * (e + cf)

    r3 = P * xo + Q * yo
    v3 = P * vxo + Q * vyo
    return r3, v3


def state_to_elem(r, v, mu: float = MU):
    """Convert state vector → orbital elements dict (heliocentric ecliptic)."""
    rm = vmag(r)
    vm = vmag(v)
    hv = vcrs(r, v)
    hm = vmag(hv)
    if hm < 1e-15 or rm < 1e-15:
        return None
    energy = 0.5 * vm ** 2 - mu / rm
    a      = -mu / (2 * energy) if abs(energy) > 1e-12 else 1e12
    eVec   = r * (vm ** 2 / mu - 1 / rm) - v * (vdot(r, v) / mu)
    e      = vmag(eVec)
    iR     = np.arccos(np.clip(hv[2] / hm, -1, 1))
    Nv     = np.array([-hv[1], hv[0], 0.0])
    Nm     = vmag(Nv)
    Omega  = np.degrees(np.arccos(np.clip(Nv[0] / Nm, -1, 1))) if Nm > 1e-12 else 0.0
    if Nm > 1e-12 and Nv[1] < 0:
        Omega = 360 - Omega
    omega = (np.degrees(np.arccos(np.clip(vdot(Nv, eVec) / (Nm * e), -1, 1)))
             if Nm > 1e-12 and e > 1e-8 else 0.0)
    if e > 1e-8 and eVec[2] < 0:
        omega = 360 - omega
    fv = (np.degrees(np.arccos(np.clip(vdot(eVec, r) / (e * rm), -1, 1)))
          if e > 1e-8 else 0.0)
    if vdot(r, v) < 0:
        fv = 360 - fv
    fR    = np.radians(fv)
    T_off = 0.0
    if e < 1 - 1e-5 and abs(a) < 1e10:
        cE2 = (e + np.cos(fR)) / (1 + e * np.cos(fR))
        sE2 = np.sqrt(1 - e ** 2) * np.sin(fR) / (1 + e * np.cos(fR))
        E2  = np.arctan2(sE2, cE2)
        T_off = (E2 - e * np.sin(E2)) / np.sqrt(mu / abs(a) ** 3)
    elif e > 1 + 1e-5 and abs(a) < 1e10:
        ah = abs(a)
        tH = np.sqrt((e - 1) / (e + 1)) * np.tan(fR / 2)
        if abs(tH) < 0.9999:
            H2    = 2 * np.arctanh(tH)
            T_off = (e * np.sinh(H2) - H2) / np.sqrt(mu / ah ** 3)
    q_out = abs(a) * (1 - e) if e < 1 else abs(a) * (e - 1)
    return dict(q=q_out, e=e, i=np.degrees(iR),
                Omega=Omega, omega=omega, T_off=T_off, a=a)


# ─────────────────────────────────────────────────────────────────────────────
#  EARTH POSITION  (astropy if available, analytic fallback)
# ─────────────────────────────────────────────────────────────────────────────
def earth_pos_analytic(jd: float) -> np.ndarray:
    """
    Earth's heliocentric ecliptic position (AU), J2000 — simplified VSOP87.

    Note: L is the *Sun's* geocentric ecliptic longitude.
          Earth's *heliocentric* longitude = L_sun + 180°.
    """
    T = (jd - J2K) / 36525.0
    M = np.radians((357.52910 + 35999.0503 * T) % 360)
    C = 1.9146 * np.sin(M) + 0.0200 * np.sin(2 * M) + 0.0003 * np.sin(3 * M)
    L_sun   = (280.46646 + 36000.76983 * T + C) % 360
    L_earth = np.radians((L_sun + 180.0) % 360)
    rv = 1.000140 - 0.016708 * np.cos(M) - 0.000141 * np.cos(2 * M)
    return np.array([rv * np.cos(L_earth), rv * np.sin(L_earth), 0.0])


def earth_pos(jd: float) -> np.ndarray:
    """Earth heliocentric ecliptic position (AU). Uses astropy when available."""
    try:
        from astropy.time import Time
        from astropy.coordinates import get_body_barycentric
        import astropy.units as u
        t = Time(jd, format='jd', scale='tdb')
        e = get_body_barycentric('earth', t)
        s = get_body_barycentric('sun',   t)
        r_eq = (e.xyz - s.xyz).to(u.AU).value   # heliocentric ICRS (≈equatorial J2000)
        # Rotate equatorial J2000 → ecliptic J2000
        cE, sE = np.cos(OBL), np.sin(OBL)
        return np.array([r_eq[0],
                          cE * r_eq[1] + sE * r_eq[2],
                         -sE * r_eq[1] + cE * r_eq[2]])
    except Exception:
        return earth_pos_analytic(jd)


# ─────────────────────────────────────────────────────────────────────────────
#  COORDINATE TRANSFORMS
# ─────────────────────────────────────────────────────────────────────────────
def ecl_to_eq(r):
    """Ecliptic J2000 → equatorial J2000."""
    cE, sE = np.cos(OBL), np.sin(OBL)
    return np.array([r[0], cE * r[1] - sE * r[2], sE * r[1] + cE * r[2]])


# ─────────────────────────────────────────────────────────────────────────────
#  RK4 INTEGRATOR  (fallback for β > 1 or degenerate orbits)
# ─────────────────────────────────────────────────────────────────────────────
def rk4_propagate(r0, v0, mu_eff, dt, steps: int = 200):
    """Integrate two-body motion under μ_eff for time dt with classical RK4."""
    h = dt / steps
    s = np.concatenate([r0, v0])

    def deriv(sv):
        r3 = np.dot(sv[:3], sv[:3]) ** 1.5
        if r3 < 1e-20:
            return np.concatenate([sv[3:], np.zeros(3)])
        return np.concatenate([sv[3:], -mu_eff / r3 * sv[:3]])

    for _ in range(steps):
        k1 = deriv(s) * h
        k2 = deriv(s + k1 / 2) * h
        k3 = deriv(s + k2 / 2) * h
        k4 = deriv(s + k3) * h
        s  = s + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if not np.all(np.isfinite(s)):
            return None
    rf = s[:3]
    rm = vmag(rf)
    return rf if (np.all(np.isfinite(rf)) and 1e-3 < rm < 1e4) else None


# ─────────────────────────────────────────────────────────────────────────────
#  DUST PARTICLE PROPAGATION  (core of Finson–Probstein)
#  β = F_radiation / F_gravity  (effective gravity reduction)
#  particle moves under μ_eff = μ_sun · (1 − β)
# ─────────────────────────────────────────────────────────────────────────────
def dust_position(beta: float, t_emit: float, t_obs: float, comet_el: dict):
    """
    Compute the position of a dust grain at t_obs, emitted at t_emit
    from the nucleus with zero relative velocity. Returns 3-vector
    [AU] in heliocentric ecliptic, or None on failure.
    """
    dt = t_obs - t_emit
    if dt < 0:
        return None
    if dt < 1e-8:
        try:
            r, _ = elem_to_state(comet_el, t_obs)
            return r
        except Exception:
            return None

    try:
        r0, v0 = elem_to_state(comet_el, t_emit)
    except Exception:
        return None

    # β = 0: particle stays on nucleus orbit
    if abs(beta) < 1e-9:
        try:
            r, _ = elem_to_state(comet_el, t_obs)
            return r
        except Exception:
            return None

    mu_eff = MU * (1.0 - beta)

    # β ≈ 1: radiation ≈ gravity → straight-line motion
    if abs(mu_eff) < 1e-10:
        return r0 + v0 * dt

    # β > 1: repulsive effective potential → numerical RK4
    if mu_eff < 0:
        return rk4_propagate(r0, v0, mu_eff, dt)

    # General case: analytic Keplerian propagation under μ_eff
    try:
        dust_el = state_to_elem(r0, v0, mu_eff)
        if dust_el is None:
            return rk4_propagate(r0, v0, mu_eff, dt)
        if (not np.isfinite(dust_el['a'])
                or abs(dust_el['a']) > 1e9
                or dust_el['q'] < 1e-5):
            return r0 + v0 * dt
        T_peri  = t_emit - dust_el['T_off']
        el_dust = {**dust_el, 'T_jd': T_peri, 'T': jd_to_str(T_peri)[:10]}
        r, _    = elem_to_state(el_dust, t_obs, mu=mu_eff)
        rm      = vmag(r)
        if not np.isfinite(rm) or rm > 1e4 or rm < 1e-3:
            return None
        return r
    except Exception:
        return rk4_propagate(r0, v0, mu_eff, dt)


# ─────────────────────────────────────────────────────────────────────────────
#  SKY PROJECTION
# ─────────────────────────────────────────────────────────────────────────────
def project_sky(r_particle, r_comet, r_earth):
    """
    Project a particle onto the sky plane relative to the comet nucleus.
    Returns (xi_AU, eta_AU, RA_deg, Dec_deg):
      xi  = East offset   (AU)
      eta = North offset  (AU)
    """
    rEC  = ecl_to_eq(r_comet - r_earth)
    rED  = ecl_to_eq(r_particle - r_earth)
    dist = vmag(rEC)
    if dist < 1e-10:
        return None
    RA   = np.arctan2(rEC[1], rEC[0])
    Dec  = np.arcsin(np.clip(rEC[2] / dist, -1, 1))
    east  = np.array([-np.sin(RA), np.cos(RA), 0.0])
    north = np.array([-np.sin(Dec) * np.cos(RA),
                      -np.sin(Dec) * np.sin(RA),
                       np.cos(Dec)])
    diff = rED - rEC
    return (float(vdot(diff, east)),
            float(vdot(diff, north)),
            np.degrees(RA) % 360,
            np.degrees(Dec))


def project_vector_sky(vec, r_comet, r_earth):
    """
    Project a DIRECTION vector (e.g. heliocentric velocity) onto the sky
    plane, using the same East/North basis that project_sky() derives from
    the comet's own position. Unlike project_sky(), `vec` is a displacement
    (not an absolute position) so no r_earth subtraction is applied to it.

    Used for the anti-velocity arrow (v3.0): the negative of the comet's
    heliocentric velocity vector, projected onto the sky, shows how much
    the dust tail is expected to lean away from the pure anti-solar
    direction due to non-zero ejection velocity / orbital motion — the
    same convention used by Moreno (2025, A&A 695 A263, Fig. 1) and
    Mariblanca-Escalona et al. (2026, MNRAS) to label their tail images.

    Returns (xi_component, eta_component) in the same units as `vec`,
    or None if the comet is at the origin (degenerate geometry).
    """
    rEC  = ecl_to_eq(r_comet - r_earth)
    dist = vmag(rEC)
    if dist < 1e-10:
        return None
    RA   = np.arctan2(rEC[1], rEC[0])
    Dec  = np.arcsin(np.clip(rEC[2] / dist, -1, 1))
    east  = np.array([-np.sin(RA), np.cos(RA), 0.0])
    north = np.array([-np.sin(Dec) * np.cos(RA),
                      -np.sin(Dec) * np.sin(RA),
                       np.cos(Dec)])
    v_eq = ecl_to_eq(vec)
    return float(vdot(v_eq, east)), float(vdot(v_eq, north))


def sky_to_pixel(xi, eta, nuc_x, nuc_y, au_per_px, north_pa_deg):
    """
    Convert sky offsets (AU East, AU North) → image pixel coordinates.
    north_pa_deg: clockwise angle from image-up to celestial North.
    """
    th = np.radians(north_pa_deg)
    dx = (xi * (-np.cos(th)) + eta * np.sin(th))    / au_per_px
    dy = (xi * (-np.sin(th)) + eta * (-np.cos(th))) / au_per_px
    return nuc_x + dx, nuc_y + dy


# ─────────────────────────────────────────────────────────────────────────────
#  MODEL COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────
def compute_model(comet_el: dict, obs_jd: float,
                  beta_values, sync_ages,
                  max_age: float = 200, n_pts: int = 200) -> dict:
    """
    Compute syndyne and synchrone curves for the given comet at obs_jd.
    Returns a dict with keys: syndynes, synchrones, sun_dir, orbit, info.
    """
    r_C, v_C = elem_to_state(comet_el, obs_jd)
    r_E      = earth_pos(obs_jd)

    ref = project_sky(r_C, r_C, r_E)
    if ref is None:
        raise RuntimeError("Sky projection failed — check orbital elements and date.")
    _, _, RA, Dec = ref

    # ── Syndynes (constant β, age 0 → max_age) ───────────────────────────
    syndynes = []
    for beta in beta_values:
        xi_arr  = np.empty(n_pts + 1)
        eta_arr = np.empty(n_pts + 1)
        age_arr = np.linspace(0, max_age, n_pts + 1)
        for k, age in enumerate(age_arr):
            r_D = dust_position(beta, obs_jd - age, obs_jd, comet_el)
            if r_D is None:
                xi_arr[k] = eta_arr[k] = np.nan
                continue
            proj = project_sky(r_D, r_C, r_E)
            if proj is None:
                xi_arr[k] = eta_arr[k] = np.nan
            else:
                xi_arr[k], eta_arr[k] = proj[0], proj[1]
        syndynes.append(dict(beta=beta, xi=xi_arr, eta=eta_arr, age=age_arr))

    # ── Synchrones (constant age, β log-spaced) ──────────────────────────
    b_min, b_max = 0.0001, max(max(beta_values), 1.0)
    synchrones = []
    for age in sync_ages:
        beta_arr = np.logspace(np.log10(b_min), np.log10(b_max), n_pts + 1)
        xi_arr   = np.empty(n_pts + 1)
        eta_arr  = np.empty(n_pts + 1)
        for k, beta in enumerate(beta_arr):
            r_D = dust_position(beta, obs_jd - age, obs_jd, comet_el)
            if r_D is None:
                xi_arr[k] = eta_arr[k] = np.nan
                continue
            proj = project_sky(r_D, r_C, r_E)
            if proj is None:
                xi_arr[k] = eta_arr[k] = np.nan
            else:
                xi_arr[k], eta_arr[k] = proj[0], proj[1]
        synchrones.append(dict(age=age, xi=xi_arr, eta=eta_arr, beta=beta_arr))

    # ── Orbital path (±40 days, every 2 days) ────────────────────────────
    orbit_pts = []
    for dt in range(-40, 41, 2):
        try:
            r_orb, _ = elem_to_state(comet_el, obs_jd + dt)
            proj     = project_sky(r_orb, r_C, r_E)
            if proj:
                orbit_pts.append((proj[0], proj[1], dt))
        except Exception:
            pass

    # ── Sun direction (points FROM comet TOWARD the Sun) ──────────────────
    r_sun_proj = project_sky(np.zeros(3), r_C, r_E)
    sun_xi  = r_sun_proj[0] if r_sun_proj else 0.0
    sun_eta = r_sun_proj[1] if r_sun_proj else 0.0

    # ── Anti-velocity direction (v3.0) ─────────────────────────────────────
    # Negative of the comet's heliocentric velocity, projected onto the sky.
    # A non-zero-ejection-velocity dust tail leans toward this direction
    # rather than the pure antisolar direction — see Moreno (2025) Fig. 1
    # and Mariblanca-Escalona et al. (2026) Fig. 7, who plot both vectors
    # on every tail image for exactly this reason.
    av_proj = project_vector_sky(-v_C, r_C, r_E)
    antivel_xi, antivel_eta = av_proj if av_proj else (0.0, 0.0)

    # ── Ephemeris info ──────────────────────────────────────────────────
    r_helio   = float(vmag(r_C))
    r_geo     = float(vmag(r_C - r_E))
    r_ES      = float(vmag(r_E))
    cos_phase = (r_helio ** 2 + r_geo ** 2 - r_ES ** 2) / (2 * r_helio * r_geo)
    phase_ang = float(np.degrees(np.arccos(np.clip(cos_phase, -1, 1))))

    info = dict(r_helio=r_helio, r_geo=r_geo, phase_angle=phase_ang,
                RA=RA, Dec=Dec, obs_str=jd_to_str(obs_jd), obs_jd=obs_jd)

    return dict(syndynes=syndynes, synchrones=synchrones,
                sun_dir=(sun_xi, sun_eta),
                antivel_dir=(antivel_xi, antivel_eta),
                orbit=orbit_pts, info=info)


# ─────────────────────────────────────────────────────────────────────────────
#  BETA → GRAIN RADIUS
# ─────────────────────────────────────────────────────────────────────────────
def compute_orbit_diagram(comet_el: dict, obs_jd: float, n_pts: int = 300) -> dict:
    """
    3D heliocentric-ecliptic diagram of the comet's position on its ACTUAL
    orbital ellipse/hyperbola, relative to the Sun, perihelion, and Earth's
    orbit (v3.0; full 3D as of this revision).

    This is physically distinct from compute_model()'s 'orbit' key, which
    is the comet's path PROJECTED ONTO THE SKY (±40 days around obs_jd,
    used to judge the tail-axis direction relative to the orbital trail).
    That sky-projection answers "which way does the tail point on this
    image"; this diagram answers "where does the comet physically sit in
    its orbit right now" — useful context when interpreting why dust
    production/Afρ changes at a given point in the orbit (e.g. post-
    perihelion fragmentation events, onset-of-activity distance, etc.),
    and the 3D tilt directly shows the orbital inclination relative to
    the ecliptic (= Earth's orbital plane).

    Uses the conic-section formula r(f) = q(1+e)/(1+e·cos f) directly in
    true anomaly, which is valid uniformly for elliptical, parabolic, and
    hyperbolic orbits alike — no period or Kepler-equation solving needed
    for the orbit shape itself (only the comet's *current* position still
    requires the usual time-based propagation via elem_to_state()).

    Returns a dict with keys (all positions are full 3D ecliptic AU,
    Z = height above/below the ecliptic plane = Earth's orbital plane):
        orbit_xyz        (n_pts,3)  the orbit's shape in 3D
        comet_xyz         (x,y,z)   comet's current position
        perihelion_xyz    (x,y,z)   perihelion point
        earth_orbit_xyz   (200,3)   Earth's orbit (assumed circular, 1 AU, z≈0)
        earth_xyz         (x,y,z)   Earth's current position
        sun_xyz           (0,0,0)
        q, e, r_helio, name
    """
    q, e = comet_el['q'], comet_el['e']
    iR  = np.radians(comet_el['i'])
    OR  = np.radians(comet_el['Omega'])
    oR  = np.radians(comet_el['omega'])
    cO, sO = np.cos(OR), np.sin(OR)
    cI, sI = np.cos(iR), np.sin(iR)
    co, so = np.cos(oR), np.sin(oR)
    P = np.array([cO*co - sO*so*cI,    sO*co + cO*so*cI,    so*sI])
    Q = np.array([-(cO*so + sO*co*cI), -(sO*so - cO*co*cI), co*sI])

    r_C, v_C    = elem_to_state(comet_el, obs_jd)
    r_helio_now = float(vmag(r_C))
    v_mag = float(vmag(v_C))
    # Unit direction of motion (heliocentric velocity) — (0,0,0) is the
    # sentinel for "no velocity data" (only possible if elem_to_state ever
    # returns a zero vector, which it shouldn't for a real orbit) so
    # draw_orbit_diagram can skip drawing the arrow rather than drawing a
    # zero-length one.
    v_hat = (v_C / v_mag) if v_mag > 1e-12 else np.zeros(3)

    # Cap the plotted curve by a physical DISTANCE (not a fixed fraction of
    # an angle), scaled to stay relevant to the comet's own perihelion and
    # current position — so Earth's 1 AU orbit is always a clearly visible
    # fraction of the frame without manual rescaling.
    r_max_plot = max(8.0 * q, 3.0 * r_helio_now, 4.0)

    # Separate, much larger absolute threshold for "is the aphelion itself
    # modest enough to just show the whole ellipse" (e.g. Halley's ~35 AU
    # is a perfectly fine, attractive full ellipse to render) — this must
    # NOT be compared against r_max_plot (~4-8 AU), which is the small
    # inner-system reference scale used only for capping genuinely huge
    # orbits; comparing against that instead would wrongly capt every
    # closed orbit larger than a few AU, including Halley's.
    APHELION_FULL_DISPLAY_LIMIT = 100.0   # AU

    show_full_ellipse = False
    if e < 1.0 - 1e-6:
        aphelion = q * (1.0 + e) / (1.0 - e)
        # Very-high-e bound orbits (e.g. e=0.999 with aphelion in the
        # hundreds of AU) are just as much a "dwarfs the inner system"
        # problem as a parabolic orbit, so above this limit they get the
        # same distance cap as the open-orbit branch below.
        show_full_ellipse = aphelion <= APHELION_FULL_DISPLAY_LIMIT

    if show_full_ellipse:
        f_arr = np.linspace(-np.pi, np.pi, n_pts)
    else:
        # Open orbit (e ≳ 1) or a closed orbit with a huge aphelion: r(f)
        # grows very fast near the (asymptotic, for e≥1) or maximum (for
        # e<1) angle — capping by a fixed fraction of that angle can still
        # reach hundreds of AU for e close to 1 (e.g. e=1.00036 reaches
        # ~157 AU at a 97%-of-asymptote cutoff). Solve the conic formula
        # directly for the true anomaly at which r = r_max_plot instead.
        cos_f_lim = np.clip((q * (1.0 + e) / r_max_plot - 1.0) / e, -1.0, 1.0)
        f_lim     = np.arccos(cos_f_lim)
        f_arr = np.linspace(-f_lim, f_lim, n_pts)

    r_arr = q * (1.0 + e) / (1.0 + e * np.cos(f_arr))
    xo, yo = r_arr * np.cos(f_arr), r_arr * np.sin(f_arr)
    orbit_xyz = np.outer(xo, P) + np.outer(yo, Q)   # (n_pts, 3), full 3D

    comet_xyz = r_C
    peri_xyz  = q * P

    theta = np.linspace(0, 2 * np.pi, 200)
    # Earth's orbit lies in the ecliptic plane by definition → Z = 0
    earth_orbit_xyz = np.column_stack([np.cos(theta), np.sin(theta),
                                        np.zeros_like(theta)])
    r_E       = earth_pos(obs_jd)
    earth_xyz = r_E

    return dict(orbit_xyz=orbit_xyz, comet_xyz=comet_xyz, perihelion_xyz=peri_xyz,
                earth_orbit_xyz=earth_orbit_xyz, earth_xyz=earth_xyz,
                sun_xyz=np.array([0.0, 0.0, 0.0]), q=q, e=e,
                r_helio=r_helio_now, name=comet_el.get('name', 'Comet'),
                obs_str=jd_to_str(obs_jd),
                velocity_dir_xyz=v_hat, r_max_plot=r_max_plot)


def _set_3d_equal_aspect(ax, points_list):
    """Force equal X/Y/Z scaling on a 3D Axes (mplot3d doesn't do this by default)."""
    all_pts  = np.vstack(points_list)
    mins, maxs = all_pts.min(axis=0), all_pts.max(axis=0)
    centers    = (mins + maxs) / 2.0
    half_range = max((maxs - mins).max() / 2.0, 1e-3)
    ax.set_xlim(centers[0] - half_range, centers[0] + half_range)
    ax.set_ylim(centers[1] - half_range, centers[1] + half_range)
    ax.set_zlim(centers[2] - half_range, centers[2] + half_range)
    try:
        ax.set_box_aspect((1, 1, 1))   # matplotlib ≥ 3.3
    except Exception:
        pass


def draw_orbit_diagram(ax, diagram: dict, dark: bool = True):
    """
    Render the 3D orbit diagram built by compute_orbit_diagram() onto a
    matplotlib 3D Axes (created with projection='3d'). Shared between the
    CLI/standalone path and the GUI's ORBIT VIEW window so the two never
    visually drift apart. Mouse-drag rotates the view interactively when
    embedded in a Qt canvas.

    `dark` selects between the dark-space and light-observatory color
    schemes (v3.0) — this module has no dependency on CometTailGUI's Qt
    theme system, so the caller just passes a bool (e.g.
    dark=(CURRENT_THEME == "dark")) rather than importing a theme dict.
    """
    if dark:
        c = dict(bg='#060b14', pane='#0a1220', pane_edge='#1a2a40',
                 tick='#2a4060', label='#3a6080', legend_bg='#060b14',
                 legend_edge='#1a2a40', legend_text='white', comet_fill='#ffffff')
    else:
        c = dict(bg='#ffffff', pane='#f3f5f8', pane_edge='#c6cdd4',
                 tick='#57606a', label='#3a6080', legend_bg='#ffffff',
                 legend_edge='#c6cdd4', legend_text='#24292f', comet_fill='#1c2128')

    ax.set_facecolor(c['bg'])
    try:
        ax.xaxis.pane.set_facecolor(c['pane'])
        ax.yaxis.pane.set_facecolor(c['pane'])
        ax.zaxis.pane.set_facecolor(c['pane'])
        ax.xaxis.pane.set_edgecolor(c['pane_edge'])
        ax.yaxis.pane.set_edgecolor(c['pane_edge'])
        ax.zaxis.pane.set_edgecolor(c['pane_edge'])
    except Exception:
        pass
    ax.tick_params(colors=c['tick'])

    op = diagram['orbit_xyz']
    ax.plot(op[:, 0], op[:, 1], op[:, 2], '-', color='#7ab8ff', lw=1.3, zorder=2,
            label=f"Orbit (q={diagram['q']:.3f} AU, e={diagram['e']:.4f})")

    eo = diagram['earth_orbit_xyz']
    ax.plot(eo[:, 0], eo[:, 1], eo[:, 2], '--', color='#2a7a4a', lw=0.8,
            alpha=0.7, zorder=1, label="Earth's orbit (1 AU, ecliptic plane)")

    # Sun — yellow circle (changed from star marker, v3.0)
    sx, sy, sz = diagram['sun_xyz']
    ax.scatter([sx], [sy], [sz], s=260, c='#ffe030', marker='o',
              edgecolors='#a06800', linewidths=1.2, depthshade=False,
              zorder=5, label='Sun')

    px, py, pz = diagram['perihelion_xyz']
    ax.scatter([px], [py], [pz], s=50, c='#ff8030', marker='x',
              linewidths=1.8, depthshade=False, zorder=4, label='Perihelion')

    ex, ey, ez = diagram['earth_xyz']
    ax.scatter([ex], [ey], [ez], s=42, c='#40c0ff', marker='o',
              depthshade=False, zorder=5,
              label=f"Earth (at {diagram['obs_str']})")

    cx, cy, cz = diagram['comet_xyz']
    # comet marker fill switches dark/light so it's never a white-on-white
    # (or near-invisible) blob — the pink edge stays constant either way.
    ax.scatter([cx], [cy], [cz], s=55, c=c['comet_fill'], marker='o',
              edgecolors='#ff5078', linewidths=1.4, depthshade=False,
              zorder=6,
              label=f"{diagram['name']} (at {diagram['obs_str']}, r☉={diagram['r_helio']:.3f} AU)")

    # Drop-lines to the ecliptic plane (Z=0) — makes the out-of-plane
    # height of the comet and perihelion immediately readable in 3D.
    ax.plot([cx, cx], [cy, cy], [0, cz], ':', color='#ff5078', lw=0.9, alpha=0.6, zorder=3)
    ax.plot([px, px], [py, py], [0, pz], ':', color='#ff8030', lw=0.9, alpha=0.5, zorder=3)

    # Direction-of-motion arrow (v3.0) — heliocentric velocity unit vector
    # from compute_orbit_diagram(), scaled to a fixed fraction of the same
    # r_max_plot used to size the whole diagram so it stays a sensibly
    # visible length regardless of how big/small this particular orbit is.
    vhx, vhy, vhz = diagram.get('velocity_dir_xyz', (0.0, 0.0, 0.0))
    if vhx or vhy or vhz:
        arrow_len = 0.18 * diagram.get('r_max_plot', max(diagram['q'], 1.0))
        ax.quiver(cx, cy, cz, vhx*arrow_len, vhy*arrow_len, vhz*arrow_len,
                  color='#00e0a0', linewidth=2.0, arrow_length_ratio=0.3,
                  zorder=7, label=f"{diagram['name']} direction of motion")

    _set_3d_equal_aspect(ax, [op, eo, np.array([[sx, sy, sz]]),
                              np.array([[px, py, pz]]), np.array([[ex, ey, ez]]),
                              np.array([[cx, cy, cz]])])

    ax.set_xlabel('Ecliptic X (AU)', color=c['label'], fontfamily='monospace', fontsize=8)
    ax.set_ylabel('Ecliptic Y (AU)', color=c['label'], fontfamily='monospace', fontsize=8)
    ax.set_zlabel('Ecliptic Z (AU)', color=c['label'], fontfamily='monospace', fontsize=8)
    ax.legend(loc='lower left', fontsize=6, framealpha=0.85,
              facecolor=c['legend_bg'], edgecolor=c['legend_edge'],
              labelcolor=c['legend_text'],
              prop={'family': 'monospace'}, markerscale=0.6,
              handletextpad=0.4, borderpad=0.5, labelspacing=0.35)


def beta_to_size(beta: float, rho_g_cm3: float = 0.5, Qpr: float = 1.0) -> str:
    """
    Grain RADIUS (not diameter) for given β and dust density ρ (g/cm³).

    From Burns, Lamy & Soter (1979), Icarus 40, 1, Eq. 19 — derived directly,
    without the C_pr/Q_pr two-parameter split (see BETA_TO_RADIUS_UM above
    for why this avoids a factor-of-2 ambiguity present in v2.4):
        r [µm] = 0.574 · Qpr / (ρ [g/cm³] · β)

    Default ρ = 0.5 g/cm³ (= 500 kg/m³), the Fulle et al. (2016) in-situ
    Rosetta/67P bulk-density measurement; default Qpr = 1 (appropriate for
    grains large relative to the wavelength — see Burns et al. 1979).
    Both are overridable from the GUI's Calculation > Dust particle
    radius… calculator (v3.0).

    Reference values at ρ = 1 g/cm³, Qpr = 1 (for comparison with older
    v3.0 output before the default ρ changed):
        β = 1.0   → r ≈ 0.57 µm   (radiation pressure ≈ gravity)
        β = 0.1   → r ≈ 5.7 µm
        β = 0.01  → r ≈ 57 µm
        β = 0.001 → r ≈ 0.57 mm   (near-orbital, debris trail)
        (at the ρ=0.5 default, every value above is exactly 2× larger)
    """
    if beta <= 0:
        return '∞'
    a_um = BETA_TO_RADIUS_UM * Qpr / (beta * rho_g_cm3)
    if   a_um >= 10000: return f'{a_um/1000:.1f} mm'
    elif a_um >= 1000:  return f'{a_um/1000:.2f} mm'
    elif a_um >= 1:     return f'{a_um:.1f} µm'
    else:               return f'{a_um:.3f} µm'


# ─────────────────────────────────────────────────────────────────────────────
#  VISUALIZER
# ─────────────────────────────────────────────────────────────────────────────
class CometTailVisualizer:
    """
    Interactive matplotlib visualizer for the Finson–Probstein model.
    Supports a standalone sky-plot mode and an image-overlay mode.
    """

    SYNDYNE_CMAP = plt.cm.autumn_r
    SYNC_CMAP    = plt.cm.YlOrBr_r

    DEFAULT_BETAS = [0.001, 0.01, 0.05, 0.1, 0.3, 0.6, 1.0]
    DEFAULT_AGES  = [10, 30, 60, 90, 120, 180]

    def __init__(self, comet_el: dict, obs_jd: float,
                 beta_values=None, sync_ages=None,
                 max_age: float = 200, n_pts: int = 80,
                 image_path: str | None = None,
                 au_per_px: float = 0.001,
                 north_pa: float = 0.0,
                 title: str = '',
                 rho_d: float = 0.5):

        self.comet_el    = comet_el
        self.obs_jd      = obs_jd
        self.beta_values = beta_values or self.DEFAULT_BETAS
        self.sync_ages   = sync_ages   or self.DEFAULT_AGES
        self.max_age     = max_age
        self.n_pts       = n_pts
        self.image_path  = image_path
        self.au_per_px   = au_per_px
        self.north_pa    = north_pa
        self.title       = title or comet_el.get('name', 'Comet')
        # Grain bulk density [g/cm³] for the β→radius annotation label.
        # Default = Fulle et al. (2016); settable directly (like nuc_x,
        # au_per_px, etc. below) from the GUI's Physical Parameters dialog.
        self.rho_d       = rho_d

        self.model      = None
        self.img_arr    = None
        self.nuc_x      = 0
        self.nuc_y      = 0
        self.show_synd  = True
        self.show_sync  = True
        self.show_orbit = True

        if image_path:
            try:
                from PIL import Image as PILImage
                img = PILImage.open(image_path)
                self.img_arr = np.array(img)
                h, w = self.img_arr.shape[:2]
                self.nuc_x, self.nuc_y = w // 2, h // 2
            except Exception as ex:
                print(f"  Warning: Could not load image: {ex}")
                self.img_arr = None

    # ── Compute / plot dispatch ──────────────────────────────────────────

    def _compute(self):
        print("\n  Computing Finson–Probstein model...")
        print(f"  β values: {self.beta_values}")
        print(f"  Synchrone ages: {self.sync_ages} days")
        t0 = time.time()
        self.model = compute_model(
            self.comet_el, self.obs_jd,
            self.beta_values, self.sync_ages,
            self.max_age, self.n_pts)
        print(f"  Done in {time.time()-t0:.1f}s")
        info = self.model['info']
        print(f"  r☉ = {info['r_helio']:.4f} AU  |  Δ = {info['r_geo']:.4f} AU"
              f"  |  Phase = {info['phase_angle']:.1f}°")
        print(f"  RA = {info['RA']:.4f}°  |  Dec = {info['Dec']:.4f}°")

    def plot(self, save_path: str | None = None):
        """Create the main visualization. Save if save_path is given."""
        self._compute()
        if self.img_arr is not None:
            self._plot_overlay(save_path)
        else:
            self._plot_standalone(save_path)

    # ── Standalone plot ──────────────────────────────────────────────────

    def _plot_standalone(self, save_path=None):
        fig = plt.figure(figsize=(14, 8), facecolor='#050810')
        fig.canvas.manager.set_window_title(f'Comet Tail Analyzer — {self.title}')

        gs_main = gs.GridSpec(1, 2, figure=fig, width_ratios=[3, 1],
                              left=0.05, right=0.98, top=0.93, bottom=0.08, wspace=0.02)
        ax      = fig.add_subplot(gs_main[0])
        ax_info = fig.add_subplot(gs_main[1])

        self._draw_standalone_ax(ax)
        self._draw_info_panel(ax_info)
        self._add_legend(ax)

        fig.suptitle(self.title, color='#7ab8ff', fontsize=13,
                     fontfamily='monospace', y=0.97)
        plt.figtext(0.5, 0.01,
                    f'Finson–Probstein model  |  Obs: {self.model["info"]["obs_str"]}',
                    ha='center', color='#2a4060', fontsize=9, fontfamily='monospace')

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight',
                        facecolor=fig.get_facecolor())
            print(f"  Saved → {save_path}")
        else:
            plt.show()

    def _draw_standalone_ax(self, ax):
        m = self.model

        # Auto-scale: collect all points including origin
        all_xi  = [0.0]
        all_eta = [0.0]
        for s in (*m['syndynes'], *m['synchrones']):
            xi, eta = s['xi'], s['eta']
            mask = np.isfinite(xi) & np.isfinite(eta)
            all_xi.extend(xi[mask])
            all_eta.extend(eta[mask])

        if len(all_xi) < 2:
            ax.text(0.5, 0.5, 'No valid points computed',
                    ha='center', va='center',
                    color='red', transform=ax.transAxes)
            return

        xm = np.percentile(all_xi,  [2, 98])
        ym = np.percentile(all_eta, [2, 98])
        pad = 0.15
        dx  = (xm[1] - xm[0]) * pad or 0.01
        dy  = (ym[1] - ym[0]) * pad or 0.01
        ax.set_xlim(xm[0] - dx, xm[1] + dx)
        ax.set_ylim(ym[0] - dy, ym[1] + dy)

        ax.set_facecolor('#060b14')
        ax.tick_params(colors='#2a4060')
        for sp in list(ax.spines.values()):
            sp.set_edgecolor('#1a2a40')
        ax.set_xlabel('Δ East (AU)  →',  color='#3a6080',
                      fontsize=9, fontfamily='monospace')
        ax.set_ylabel('↑  Δ North (AU)', color='#3a6080',
                      fontsize=9, fontfamily='monospace')
        ax.grid(color='#0d1a2e', linestyle='-', linewidth=0.4, alpha=0.8)

        # Star-field decoration
        rng = np.random.default_rng(42)
        xl, yl = ax.get_xlim(), ax.get_ylim()
        sx = rng.uniform(xl[0], xl[1], 200)
        sy = rng.uniform(yl[0], yl[1], 200)
        sz = rng.uniform(0.2, 1.0, 200)
        ax.scatter(sx, sy, s=sz * 0.6, c='white', alpha=0.3, zorder=0)

        # Orbital path
        if self.show_orbit and m['orbit']:
            op = np.array(m['orbit'])
            ax.plot(op[:, 0], op[:, 1], '--', color='#2050a0', lw=0.8,
                    alpha=0.5, zorder=1, label='Orbital path')
            if len(op) > 1:
                mid = len(op) // 2
                ax.annotate('', xy=(op[mid + 1, 0], op[mid + 1, 1]),
                            xytext=(op[mid, 0], op[mid, 1]),
                            arrowprops=dict(arrowstyle='->',
                                            color='#2050a0', lw=0.8),
                            zorder=2)

        self._draw_curves(ax, m, in_pixels=False)

        # Anti-solar / Sun direction arrow
        sx_dir, sy_dir = m['sun_dir']
        slen = np.hypot(sx_dir, sy_dir)
        if slen > 1e-10:
            scale = min(abs(ax.get_xlim()[1] - ax.get_xlim()[0]),
                        abs(ax.get_ylim()[1] - ax.get_ylim()[0])) * 0.18
            ex, ey = sx_dir / slen * scale, sy_dir / slen * scale
            ax.annotate('', xy=(ex, ey), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='->',
                                        color='#ffe030', lw=1.8),
                        zorder=4)
            ax.text(ex * 1.15, ey * 1.15, '☀',
                    color='#ffe030', fontsize=12,
                    ha='center', va='center', zorder=5)

        # Anti-velocity arrow (v3.0) — shows expected lean of a
        # non-zero-ejection-velocity tail away from the pure antisolar line
        avx_dir, avy_dir = m.get('antivel_dir', (0.0, 0.0))
        avlen = np.hypot(avx_dir, avy_dir)
        if avlen > 1e-10:
            scale = min(abs(ax.get_xlim()[1] - ax.get_xlim()[0]),
                        abs(ax.get_ylim()[1] - ax.get_ylim()[0])) * 0.18
            evx, evy = avx_dir / avlen * scale, avy_dir / avlen * scale
            ax.annotate('', xy=(evx, evy), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='->',
                                        color='#ff5078', lw=1.6),
                        zorder=4)
            ax.text(evx * 1.15, evy * 1.15, '−v',
                    color='#ff5078', fontsize=9, fontfamily='monospace',
                    ha='center', va='center', zorder=5)

        # Nucleus
        ax.plot(0, 0, '+', color='white', ms=12, mew=1.5, zorder=6)
        ax.plot(0, 0, 'o', color='white', ms=3,  zorder=7)

        # Compass
        cx = xl[0] + (xl[1] - xl[0]) * 0.88
        cy = yl[0] + (yl[1] - yl[0]) * 0.12
        aL = (xl[1] - xl[0]) * 0.06
        th = np.radians(self.north_pa)
        ax.annotate('N', xy=(cx + np.sin(th) * aL, cy - np.cos(th) * aL),
                    xytext=(cx, cy),
                    arrowprops=dict(arrowstyle='->', color='#60c8ff', lw=1.5),
                    color='#60c8ff', fontsize=8,
                    fontfamily='monospace', zorder=8)
        ax.annotate('E', xy=(cx - np.cos(th) * aL, cy - np.sin(th) * aL),
                    xytext=(cx, cy),
                    arrowprops=dict(arrowstyle='->', color='#4080b0', lw=1.2),
                    color='#4080b0', fontsize=8,
                    fontfamily='monospace', zorder=8)

        # Scale bar
        raw_scale = (xl[1] - xl[0]) * 0.15
        mag       = 10 ** np.floor(np.log10(raw_scale))
        nbar      = round(raw_scale / mag) * mag
        bx        = xl[0] + (xl[1] - xl[0]) * 0.05
        by        = yl[0] + (yl[1] - yl[0]) * 0.04
        ax.plot([bx, bx + nbar], [by, by],
                color='white', lw=2, alpha=0.7, zorder=8)
        ax.plot([bx] * 2, [by - dy * 0.2, by + dy * 0.2],
                color='white', lw=1, alpha=0.6)
        ax.plot([bx + nbar] * 2, [by - dy * 0.2, by + dy * 0.2],
                color='white', lw=1, alpha=0.6)
        label = f'{nbar:.3f} AU' if nbar < 0.1 else f'{nbar:.2f} AU'
        ax.text(bx + nbar / 2, by + dy * 0.35, label,
                ha='center', color='white', fontsize=8,
                alpha=0.8, fontfamily='monospace')

    # ── Shared curve drawer ──────────────────────────────────────────────

    def _draw_curves(self, ax, m, in_pixels: bool):
        """Draw syndyne/synchrone curves either in AU (standalone) or in pixels."""
        def _to_pixel(xi, eta):
            return sky_to_pixel(xi, eta, self.nuc_x, self.nuc_y,
                                self.au_per_px, self.north_pa)

        def _split_segments(xv, yv):
            """Yield (x_seg, y_seg) for each contiguous run of finite values."""
            mask = np.isfinite(xv) & np.isfinite(yv)
            if not mask.any():
                return
            idx = np.where(np.diff(mask.astype(int)) != 0)[0] + 1
            edges = [0, *idx, len(xv)]
            for i in range(len(edges) - 1):
                sl = slice(edges[i], edges[i + 1])
                if mask[sl].any():
                    yield xv[sl], yv[sl]

        # Syndynes
        if self.show_synd:
            n_s = len(m['syndynes'])
            for idx, synd in enumerate(m['syndynes']):
                c  = self.SYNDYNE_CMAP(idx / (n_s - 1) if n_s > 1 else 0.5)
                xi, eta = synd['xi'], synd['eta']
                if in_pixels:
                    px, py = _to_pixel(xi, eta)
                else:
                    px, py = xi, eta
                for x_seg, y_seg in _split_segments(px, py):
                    ax.plot(x_seg, y_seg, color=c, lw=1.4, alpha=0.85, zorder=3)
                # Tail-end label
                valid = np.where(np.isfinite(px) & np.isfinite(py))[0]
                if len(valid):
                    lp = valid[-1]
                    ax.text(px[lp] + (4 if in_pixels else 0),
                            py[lp],
                            f' β={synd["beta"]}',
                            color=c, fontsize=7.5, va='center',
                            fontfamily='monospace', zorder=5)

        # Synchrones
        if self.show_sync:
            n_sy = len(m['synchrones'])
            for idx, sync in enumerate(m['synchrones']):
                c = self.SYNC_CMAP(idx / (n_sy - 1) if n_sy > 1 else 0.5)
                xi, eta = sync['xi'], sync['eta']
                if in_pixels:
                    px, py = _to_pixel(xi, eta)
                else:
                    px, py = xi, eta
                for x_seg, y_seg in _split_segments(px, py):
                    ax.plot(x_seg, y_seg, color=c, lw=1.2, alpha=0.82,
                            linestyle='-.', zorder=3)
                valid = np.where(np.isfinite(px) & np.isfinite(py))[0]
                if len(valid):
                    lp = valid[-1]
                    ax.text(px[lp] + (4 if in_pixels else 0),
                            py[lp],
                            f' t−{sync["age"]}d',
                            color=c, fontsize=7.5, va='center',
                            fontfamily='monospace', zorder=5)

    # ── Info panel ───────────────────────────────────────────────────────

    def _draw_info_panel(self, ax):
        ax.set_facecolor('#060b14')
        ax.set_axis_off()
        info = self.model['info']

        y, dy = 0.97, 0.055
        kw = dict(transform=ax.transAxes, fontfamily='monospace', va='top')

        ax.text(0.05, y, 'EPHEMERIS', color='#3a6090',
                fontsize=8, fontweight='bold', **kw)
        y -= dy * 0.6
        for k, v in [
            ('r☉',    f'{info["r_helio"]:.5f} AU'),
            ('Δ',     f'{info["r_geo"]:.5f} AU'),
            ('Phase', f'{info["phase_angle"]:.2f}°'),
            ('RA',    f'{info["RA"]:.4f}°'),
            ('Dec',   f'{info["Dec"]:.4f}°'),
        ]:
            ax.text(0.05, y, k, color='#2a5070', fontsize=8, **kw)
            ax.text(0.95, y, v, color='#70b0e0', fontsize=8, ha='right', **kw)
            y -= dy * 0.9

        y -= dy * 0.3
        ax.text(0.05, y, 'ORBITAL ELEMENTS', color='#3a6090',
                fontsize=8, fontweight='bold', **kw)
        y -= dy * 0.6
        el = self.comet_el
        for k, v in [
            ('q',  f'{el["q"]:.5f} AU'),
            ('e',  f'{el["e"]:.6f}'),
            ('i',  f'{el["i"]:.3f}°'),
            ('Ω',  f'{el["Omega"]:.3f}°'),
            ('ω',  f'{el["omega"]:.3f}°'),
            ('T',  el.get('T', '')[:10]),
        ]:
            ax.text(0.05, y, k, color='#2a5070', fontsize=8, **kw)
            ax.text(0.95, y, v, color='#70b0e0', fontsize=8, ha='right', **kw)
            y -= dy * 0.9

        y -= dy * 0.3
        ax.text(0.05, y, f'β → GRAIN RADIUS (ρ={self.rho_d:g} g/cm³)', color='#3a6090',
                fontsize=8, fontweight='bold', **kw)
        y -= dy * 0.6
        n_s = len(self.beta_values)
        for idx, beta in enumerate(self.beta_values):
            c = self.SYNDYNE_CMAP(idx / (n_s - 1) if n_s > 1 else 0.5)
            ax.text(0.05, y, f'β={beta}', color=c, fontsize=7.5, **kw)
            ax.text(0.95, y, beta_to_size(beta, self.rho_d), color='#a0c0d8',
                    fontsize=7.5, ha='right', **kw)
            y -= dy * 0.85

        note = el.get('note', '')
        if note:
            y = max(y - dy * 0.3, 0.04)
            for line in [note[i:i + 22] for i in range(0, len(note), 22)]:
                ax.text(0.05, y, line, color='#1e3a5a', fontsize=7, **kw)
                y -= dy * 0.75

    def _add_legend(self, ax):
        handles = []
        if self.show_synd:
            handles.append(Line2D([0], [0],
                                  color=self.SYNDYNE_CMAP(0.5), lw=1.5,
                                  label='Syndynes  (const. β)'))
        if self.show_sync:
            handles.append(Line2D([0], [0],
                                  color=self.SYNC_CMAP(0.5), lw=1.2,
                                  ls='-.', label='Synchrones (const. age)'))
        if self.show_orbit:
            handles.append(Line2D([0], [0], color='#2050a0', lw=0.8,
                                  ls='--', label='Orbital path'))
        handles.append(Line2D([0], [0], color='#ffe030', lw=1.5, marker='>',
                              ms=7, label='Sun direction'))
        handles.append(Line2D([0], [0], color='#ff5078', lw=1.4, marker='>',
                              ms=6, label='Anti-velocity (−v)'))
        if handles:
            ax.legend(handles=handles, loc='lower left',
                      fontsize=7.5, framealpha=0.2, facecolor='#060b14',
                      edgecolor='#1a2a40', labelcolor='white',
                      prop={'family': 'monospace'})

    # ── Image overlay ────────────────────────────────────────────────────

    def _plot_overlay(self, save_path=None):
        if self.img_arr is None:
            print("No image loaded; switching to standalone mode.")
            self._plot_standalone(save_path)
            return

        h_img, w_img = self.img_arr.shape[:2]
        asp = h_img / w_img
        fig, ax = plt.subplots(figsize=(12, 12 * asp), facecolor='black')
        fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.01)
        fig.canvas.manager.set_window_title(f'Overlay — {self.title}')
        ax.set_facecolor('black')
        ax.imshow(self.img_arr, origin='upper', aspect='equal', zorder=0)
        ax.set_xlim(0, w_img)
        ax.set_ylim(h_img, 0)
        ax.axis('off')

        m = self.model
        self._draw_curves(ax, m, in_pixels=True)

        # Anti-solar arrow
        sun_xi, sun_eta = m['sun_dir']
        slen = np.hypot(sun_xi, sun_eta)
        if slen > 1e-10:
            arrow_au = 60 * self.au_per_px
            tip = sky_to_pixel(sun_xi / slen * arrow_au,
                               sun_eta / slen * arrow_au,
                               self.nuc_x, self.nuc_y,
                               self.au_per_px, self.north_pa)
            ax.annotate('☀', xy=tip, xytext=(self.nuc_x, self.nuc_y),
                        arrowprops=dict(arrowstyle='->',
                                        color='#ffe030', lw=1.8),
                        color='#ffe030', fontsize=12,
                        ha='center', va='center', zorder=6)

        # Anti-velocity arrow (v3.0)
        avx, avy = m.get('antivel_dir', (0.0, 0.0))
        avlen = np.hypot(avx, avy)
        if avlen > 1e-10:
            arrow_au = 60 * self.au_per_px
            tip = sky_to_pixel(avx / avlen * arrow_au,
                               avy / avlen * arrow_au,
                               self.nuc_x, self.nuc_y,
                               self.au_per_px, self.north_pa)
            ax.annotate('−v', xy=tip, xytext=(self.nuc_x, self.nuc_y),
                        arrowprops=dict(arrowstyle='->',
                                        color='#ff5078', lw=1.6),
                        color='#ff5078', fontsize=9, fontfamily='monospace',
                        ha='center', va='center', zorder=6)

        # Nucleus crosshair
        cs = 12
        ax.plot([self.nuc_x - cs, self.nuc_x + cs],
                [self.nuc_y, self.nuc_y], 'w-', lw=1.5, zorder=7)
        ax.plot([self.nuc_x, self.nuc_x],
                [self.nuc_y - cs, self.nuc_y + cs], 'w-', lw=1.5, zorder=7)
        ax.plot(self.nuc_x, self.nuc_y, 'wo', ms=4, zorder=8)

        # Orientation indicator
        aL = 40
        th = np.radians(self.north_pa)
        ox, oy = w_img - 80, h_img - 80
        ax.annotate('N', xy=(ox + np.sin(th) * aL, oy - np.cos(th) * aL),
                    xytext=(ox, oy),
                    arrowprops=dict(arrowstyle='->', color='#60c8ff', lw=2),
                    color='#60c8ff', fontsize=10,
                    fontfamily='monospace', zorder=9)

        self._add_legend(ax)

        info = m['info']
        ax.set_title(f'{self.title}  |  {info["obs_str"]}  |  '
                     f'r☉={info["r_helio"]:.4f}AU  Δ={info["r_geo"]:.4f}AU  '
                     f'Phase={info["phase_angle"]:.1f}°',
                     color='#7ab8ff', fontsize=9,
                     fontfamily='monospace', pad=4)

        if save_path:
            fig.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='black')
            print(f"  Saved → {save_path}")
        else:
            plt.show()

    # ── Export ───────────────────────────────────────────────────────────

    def save_csv(self, path: str):
        """Export syndyne and synchrone points to CSV."""
        if not self.model:
            self._compute()
        import csv
        with open(path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['type', 'parameter', 'xi_AU', 'eta_AU'])
            for synd in self.model['syndynes']:
                for xi, eta in zip(synd['xi'], synd['eta']):
                    if np.isfinite(xi) and np.isfinite(eta):
                        w.writerow(['syndyne', synd['beta'], xi, eta])
            for sync in self.model['synchrones']:
                for xi, eta in zip(sync['xi'], sync['eta']):
                    if np.isfinite(xi) and np.isfinite(eta):
                        w.writerow(['synchrone', sync['age'], xi, eta])
        print(f"  CSV saved → {path}")

    def save_fits_overlay(self, fits_input: str, fits_output: str):
        """Add WCS-aligned overlay PNG using astropy."""
        from astropy.io import fits
        from astropy.wcs import WCS
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        if not self.model:
            self._compute()
        with fits.open(fits_input) as hdul:
            wcs = WCS(hdul[0].header)
        info = self.model['info']
        comet_sc = SkyCoord(ra=info['RA'] * u.deg, dec=info['Dec'] * u.deg,
                            frame='icrs')
        nuc_px = wcs.world_to_pixel(comet_sc)
        self.nuc_x = float(nuc_px[0])
        self.nuc_y = float(nuc_px[1])
        ps = wcs.proj_plane_pixel_scales()
        self.au_per_px = float(ps[0].to(u.arcsec).value) * info['r_geo'] / 206265.0
        print(f"  WCS: nucleus at ({self.nuc_x:.1f}, {self.nuc_y:.1f}) px, "
              f"scale={self.au_per_px:.6f} AU/px")
        out = fits_output.replace('.fits', '_overlay.png')
        self._plot_overlay(save_path=out)


# ─────────────────────────────────────────────────────────────────────────────
#  COMMAND-LINE INTERFACE
# ─────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description='Comet Dust Tail Analyzer — Finson–Probstein Model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python comet_tail_analyzer.py
  python comet_tail_analyzer.py --comet "C/2020 F3 (NEOWISE)" --date 2020-07-23
  python comet_tail_analyzer.py --fetch "C/2023 A3" --date 2024-10-10
  python comet_tail_analyzer.py --list
  python comet_tail_analyzer.py --comet "Hale-Bopp" --image photo.jpg --au-per-px 0.0005 --north-pa 15
  python comet_tail_analyzer.py --comet "NEOWISE" --save output.png
  python comet_tail_analyzer.py --comet "Encke" --csv data.csv
        """)
    p.add_argument('--comet',     type=str, help='Comet name/designation — fetched from JPL Horizons')
    p.add_argument('--fetch',     type=str, help='Fetch from JPL Horizons by designation')
    p.add_argument('--date',      type=str, help='Observation date YYYY-MM-DD [default: today]')
    p.add_argument('--beta',      type=str,
                   default=','.join(map(str, CometTailVisualizer.DEFAULT_BETAS)),
                   help='Beta values comma-separated')
    p.add_argument('--ages',      type=str,
                   default=','.join(map(str, CometTailVisualizer.DEFAULT_AGES)),
                   help='Synchrone ages (days) comma-separated')
    p.add_argument('--max-age',   type=int, default=200, help='Max emission age (days)')
    p.add_argument('--npts',      type=int, default=80,  help='Curve resolution')
    p.add_argument('--image',     type=str, help='Overlay on image file (JPEG/PNG/FITS)')
    p.add_argument('--au-per-px', type=float, default=0.001, help='Image scale AU/pixel')
    p.add_argument('--north-pa',  type=float, default=0.0,   help='North PA (deg CW from up)')
    p.add_argument('--nuc-x',     type=float, help='Nucleus X pixel in image')
    p.add_argument('--nuc-y',     type=float, help='Nucleus Y pixel in image')
    p.add_argument('--save',      type=str,   help='Save plot to PNG')
    p.add_argument('--csv',       type=str,   help='Export model points to CSV')
    p.add_argument('--list',      action='store_true', help='List all comets in catalogue')
    return p.parse_args()


def interactive_menu():
    """Simple text-based comet selection if no args given."""
    print("\n" + "═" * 60)
    print("  COMET TAIL ANALYZER — Finson–Probstein Model")
    print("═" * 60)
    names = list(COMET_DB.keys())
    print("\n  SELECT A COMET (orbital elements fetched live from Horizons):\n")
    for i, name in enumerate(names):
        note = COMET_DB[name].get('note', '')[:45]
        pref = COMET_DB[name].get('obs', '')
        print(f"  [{i+1:2d}] {name:<38} {pref}  {note}")
    print("\n  [ 0] Enter any designation manually\n")
    choice = input("  Enter number (or 0): ").strip()

    if choice == '0':
        desig = input("  Designation (e.g. C/2023 A3): ").strip()
        date  = input("  Obs date YYYY-MM-DD [Enter=today]: ").strip() or None
    else:
        try:
            idx = int(choice) - 1
        except ValueError:
            print("  Invalid selection.")
            sys.exit(1)
        if not (0 <= idx < len(names)):
            print("  Invalid selection.")
            sys.exit(1)
        key   = names[idx]
        desig = key
        pref  = COMET_DB[key].get('obs', '')
        date  = input(f"  Obs date [Enter={pref}]: ").strip() or pref or None

    print(f"\n  Fetching elements for '{desig}' from JPL Horizons…")
    el = fetch_comet(desig, date=date)

    obs_jd = date_to_jd(date) if date else today_jd()
    print(f"\n  Using: {el.get('name','')}")
    print(f"  Source: {el.get('source','')}")
    print(f"  q={el['q']:.5f} AU   e={el['e']:.7f}   i={el['i']:.3f}°")
    print(f"  Ω={el['Omega']:.4f}°   ω={el['omega']:.4f}°")
    print(f"  Date:  {jd_to_str(obs_jd)}")
    return el, obs_jd


def main():
    args = parse_args()

    # ── List database ───────────────────────────────────────────────────
    if args.list:
        print(f"\n  {'Designation':<45} {'Preferred date':<14}  Note")
        print("  " + "─" * 90)
        for name, meta in COMET_DB.items():
            print(f"  {name:<45} {meta.get('obs',''):<14}  "
                  f"{meta.get('note','')[:40]}")
        print(f"\n  ({len(COMET_DB)} comets — orbital elements fetched live from JPL Horizons)\n")
        return

    # ── Select / fetch comet ────────────────────────────────────────────
    obs_jd_im: float | None = None
    if args.fetch:
        comet_el = fetch_comet(args.fetch, date=args.date)
    elif args.comet:
        comet_el = fetch_comet(args.comet, date=args.date)
    else:
        comet_el, obs_jd_im = interactive_menu()

    obs_jd = date_to_jd(args.date) if args.date else (obs_jd_im or today_jd())

    # ── Model parameters ────────────────────────────────────────────────
    betas = [float(x) for x in args.beta.split(',') if x.strip()]
    ages  = [int(x)   for x in args.ages.split(',') if x.strip()]

    vis = CometTailVisualizer(
        comet_el    = comet_el,
        obs_jd      = obs_jd,
        beta_values = betas,
        sync_ages   = ages,
        max_age     = args.max_age,
        n_pts       = args.npts,
        image_path  = args.image,
        au_per_px   = args.au_per_px,
        north_pa    = args.north_pa,
        title       = comet_el.get('name', 'Comet'),
    )

    # Override nucleus position if provided
    if args.nuc_x is not None and vis.img_arr is not None:
        vis.nuc_x = args.nuc_x
    if args.nuc_y is not None and vis.img_arr is not None:
        vis.nuc_y = args.nuc_y

    if args.csv:
        vis.save_csv(args.csv)

    vis.plot(save_path=args.save)


if __name__ == '__main__':
    main()
