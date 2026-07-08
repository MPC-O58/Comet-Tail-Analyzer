#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
  comet_tail_analyzer.py  —  Comet Tail Analyzer (CTA) Physics Engine
  Finson–Probstein + Monte Carlo Dust Tail Model
  Version 3.1.1 ·   Teerasak Thaluang (MPC O51/O58)

  SPDX-License-Identifier: MIT
  © 2024–2026 Teerasak Thaluang. See LICENSE for full terms.

  Attribution:
    The following components are ported from py_COMTAILS
    (Moreno 2025, A&A 695, A263; F. Moreno, R. Morales & N. Robles,
    IAA-CSIC; github.com/FernandoMorenoDanvila/py_COMTAILS, MIT License):
      - SCHLEICHER_COEFFS + schleicher_phase_correction()
          Source: py_COMTAILS/constants.py
      - _sample_sunward_direction()
          Source: py_COMTAILS iejec_mode == 2 branch
      - _active_area_direction()
          Source: py_COMTAILS DustTail._anisot_dir2() + iejec_mode == 3

    All other components (Finson–Probstein engine, orbital mechanics,
    COBS integration, compute_morphology_mc() core loop, contour
    extraction, Gaia/Afρ pipelines) are original to CTA.

  Cite as:
    Thaluang, T. (2026). RNAAS, doi:10.3847/2515-5172/ae6f90
=============================================================================
  Changelog:
    v3.1.1 • GUI hotfix release: no change to F-P or Monte Carlo physics.
    v3.1  • GUI workflow synchronization: F-P maximum dust age is derived
            from the largest listed synchrone; physics equations unchanged.
    v3.1  • GUI workflow synchronization: Q(t) preview is now embedded in
            Simulation > Dust production over time. No numerical physics
            or Monte Carlo sampling equations changed in this revision.
    v3.1  • Released together with CometTailGUI.py v3.1.
          • SAFETY: assess_qt_coverage() validates recent continuous COBS
            support before Q(t) is used to infer an MC release window.
            Historical points across long gaps/apparitions no longer create
            spurious multi-year window recommendations.
          • BUG FIX: sunward-hemisphere ejection now evaluates the
            comet's position at EACH PARTICLE'S OWN emission time
            (t_emits[k]) rather than at the observation time — the
            Sun-comet line rotates over the release window, so a fixed
            reference frame was misdirecting "sunward" for particles
            not emitted "now". Matches active_area's existing per-
            particle geometry pattern.
          • NEW: compute_morphology_mc() — Monte Carlo dust morphology
            simulation. Grain size power-law distribution n(r) ∝ r^κ;
            ejection velocity V = V₀·β^γ·r_H^(−κ) with optional
            sunward hemisphere or active-area (rotating nucleus)
            directionality (cos z term). Returns density grid + info dict
            including r_helio_au for downstream real-speed display.
          • NEW: extract_morphology_contours() — isophote-style contour
            extraction from MC density grid. Percentile floor now computed
            from PRE-Gaussian density (BUG FIX: smoothing previously
            lowered the threshold by populating empty bins, shifting
            contours outward — now smooth only affects line shape).
          • NEW: real_ejection_speed_ms() — compute actual grain speed
            V₀·β^γ·r_H^(−κ) for display in overlay and report.
          • NEW: _sample_sunward_direction() — port of py_COMTAILS
            iejec_mode == 2 (sunward hemisphere ejection sampling).
          • NEW: _active_area_direction() — port of py_COMTAILS
            DustTail._anisot_dir2() / iejec_mode == 3 (rotating nucleus
            active-area ejection direction).
          • NEW: SCHLEICHER_COEFFS + schleicher_phase_correction() —
            ported from py_COMTAILS constants.py (Schleicher, Millis &
            Birch 1998-style polynomial).
            (All three COMTAILS ports under MIT License — see Attribution
            block above.)
          • NEW: non-zero dust ejection velocity. dust_position() and
            compute_model() gained optional v_R0/v_T0/v_N0 (m/s, in the
            dust grain's local Radial/Transverse/Normal frame at the
            moment of release) plus γ/m_exp exponents for a power-law
            velocity scaling:
                v_R = v_R0 · β^γ · r_H(t_emit)^(−m_exp)   (and likewise v_T, v_N)
            R̂ = r̂ (sunward radial), N̂ = ĥ/|ĥ|  with ĥ = r×v (orbit
            normal — self-consistently right-handed for retrograde
            orbits too, since it's derived from the comet's own r,v
            rather than assumed), T̂ = N̂ × R̂ completes the triad.
            γ=m_exp=0 reduces to a CONSTANT ejection velocity (v_R0/
            v_T0/v_N0 used directly, no β or r_H dependence). All five
            parameters default to 0.0, which reproduces the v3.0
            zero-ejection-velocity model bit-for-bit — this is an
            additive, opt-in extension; it does NOT alter any
            previously published default-mode output.
            FIXED (same dev cycle, before release): r_H exponent's sign
            was originally r_H^(+m_exp) — wrong. Verified against two
            primary sources actually fetched and read in full: Hui et
            al. (2020, AJ 159, 77) Eq. 5, |v_ej| = V0·(β/r_H)^(1/2)
            (citing Whipple 1950, ApJ 111, 375; Ishiguro 2008, Icarus
            193, 96); and Ishiguro et al. (2007/2008; e.g. as used in
            the 67P dust-trail paper, Ishiguro 2008) Eq. 2, v_ej =
            V0·(a/a0)^(−u1)·(r_h/AU)^(−u2), the SAME two-independent-
            exponent form used here (u1,u2 typically both 0.5, but
            free in general — this is what justified keeping γ and
            m_exp independent rather than fixing m_exp=−γ). Both put
            r_H in the denominator. An earlier draft of this changelog
            cited "Probstein (1969)" and "Crifo (1995)" for this exact
            power law without having actually read either — Probstein
            (1969) is a real source (the underlying gas-dynamics, more
            complex than a simple power law) but was the wrong specific
            citation for this formula's exact form, and "Crifo (1995)"
            could not be verified at all. Replaced with the two sources
            above, which were read in full before citing.
            Design note: a single shared γ/m_exp scales all three R/T/N
            components — i.e. the model assumes one physical ejection
            mechanism with a fixed directional split, not three
            independently-fit power laws (Ishiguro et al.'s own u1=u2
            convention does the same kind of simplification, just on
            the β/r_H pair rather than R/T/N). γ=m_exp=0.5 recovers the
            literature's standard Whipple/Ishiguro/Hui-et-al. form
            exactly — reasonable starting point, but NOT verified
            against any specific comet's data — same status as β
            itself, expected to be re-fit per-object from imaging.
            CAUTION (still applies): realistic terminal velocities for
            mm-class grains at several AU are of order cm/s, not tens
            of m/s — if a visual fit needs v_R0/v_T0/v_N0 far above
            that scale, treat that as a sign of a misidentified feature
            (background structure, gas, or processing artifact) rather
            than confirmation of this model, not a license to keep
            raising the parameter.
          • NEW: check_for_update() — queries the GitHub Releases API
            (repo MPC-O58/Comet-Tail-Analyzer) and compares the latest
            tag against __version__ via tuple (not string) comparison.
            Returns None on ANY failure path (no network, GitHub's
            60 req/hr unauthenticated rate limit — observed firsthand
            during development, so budget for it — no releases, stale/
            unparseable tag, already current) so a failed check can
            never surface as an error in the GUI. __version__ was
            stale at "2.5" (out of sync with the docstring header for
            several releases); bumped to "3.0.2" mid-development to
            match the header at the time, then to "3.1" together with
            the rest of this changelog at actual release — this
            function now needs ONE authoritative version string —
            bump THIS line at every release from now on, not just the
            header text above. _HTTP_HEADERS' User-Agent now also
            tracks __version__ via f-string instead of a separately
            hardcoded (and equally stale) "2.5".
    v3.0.2 • _set_3d_equal_aspect() and draw_orbit_diagram() gained
            optional fixed_half_range/fixed_center and zoom_inner_au
            parameters respectively — a Sun-centered fixed-size view
            override for the GUI's new "Zoom to inner solar system"
            checkbox (OrbitWindow), instead of always auto-fitting to
            the whole comet orbit. For comets well beyond ~1 AU, that
            auto-fit view collapsed Earth's orbit (and the Earth/Sun
            separation) down to a barely-distinguishable point near the
            center; this trades "see the whole orbit" for "see the
            inner solar system clearly," with the comet itself simply
            off-screen if it's farther out than the chosen AU value.
    v3.0.1 • Fixed the orbit diagram's direction-of-motion arrow (added in
            v3.0) using a fixed length based on r_max_plot, a value
            computed once from the orbit's overall extent — independent
            of whatever the user has actually zoomed the 3D view to via
            OrbitWindow's plot toolbar, so the arrow could end up looking
            disproportionately large after zooming in (same underlying
            bug pattern as the Sun-direction/anti-velocity arrows in
            CometTailGUI.py's draw_model(), fixed there in v3.0.5).
            _set_3d_equal_aspect() now returns the half_range it actually
            set the view to, and the arrow is sized as a fraction of that
            fresh value instead — drawn after the aspect call now,
            having been before it. This keeps the arrow correctly
            proportioned every time draw_orbit_diagram() is called (e.g.
            on a date change) — it does NOT dynamically rescale during
            pure interactive zoom via the toolbar with no redraw in
            between, since matplotlib's 3D mouse/toolbar zoom doesn't
            call back into this function at all (unlike the 2D case,
            there's no straightforward redraw hook to attach to here).
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

__version__ = "3.1.1"
# NOTE (v3.1): was stale at "2.5" — left out of sync with the docstring
# header above for at least 3 releases. Bumped to match the header here
# because check_for_update() (below) now needs ONE authoritative version
# string to compare against GitHub Releases; from this point on, bump
# THIS line (not just the header text) at every release.

import argparse
import logging
import math
import os
import re
import sys
import time
import warnings
from datetime import datetime, timezone

import numpy as np


def _configure_frozen_https_certificates() -> None:
    """Use Certifi's CA bundle when CTA is packaged as a Windows executable."""
    try:
        import certifi
        ca_file = certifi.where()
        if ca_file and os.path.isfile(ca_file):
            os.environ.setdefault("SSL_CERT_FILE", ca_file)
            os.environ.setdefault("REQUESTS_CA_BUNDLE", ca_file)
    except Exception:
        pass


_configure_frozen_https_certificates()

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

# Unit conversion for ejection-velocity inputs (v3.1): the propagator works
# entirely in AU / day (MU above is AU³·day⁻²), but ejection speeds are far
# more naturally specified — and reported in the literature — in m/s. IAU
# (2012) exact definition: 1 AU = 1.495978707e11 m, 1 day = 86400 s.
AU_M         = 1.495978707e11           # m, IAU (2012) exact definition
MPS_TO_AUDAY = 86400.0 / AU_M           # multiply m/s by this to get AU/day

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
    "User-Agent": f"CometTailAnalyzer/{__version__}",
    "Accept":     "application/json",
}

# ─────────────────────────────────────────────────────────────────────────────
#  GITHUB UPDATE CHECK  (v3.1)
#  Queries the GitHub Releases API for the latest tagged version and
#  compares it against __version__. Designed to be called from a
#  background thread (see CometTailGUI.py's UpdateCheckWorker) — never
#  raises, ALWAYS returns either an update-info dict or None, no matter
#  what goes wrong (no network, GitHub rate-limited, repo has no
#  releases yet, version already current). This is a convenience
#  feature, not a core one: a failure here must never surface as an
#  error to the user.
# ─────────────────────────────────────────────────────────────────────────────
GITHUB_REPO = "MPC-O58/Comet-Tail-Analyzer"


def _version_tuple(v: str):
    """Parse a CTA semantic version into a normalized ``(major, minor, patch)`` tuple.

    Accepted examples include ``v3.1``, ``3.1.1`` and ``CTA-v3.2.0``.  The
    parser deliberately ignores unrelated numbers elsewhere in a tag, so a
    date or build number cannot accidentally outrank the real application
    version.  Missing patch numbers are normalized to zero.
    """
    if not v:
        return None
    match = re.search(r'(?<!\d)v?(\d+)\.(\d+)(?:\.(\d+))?', str(v), re.I)
    if not match:
        return None
    major, minor, patch = match.groups()
    return int(major), int(minor), int(patch or 0)


def get_update_status(current_version: str | None = None,
                      repo: str = GITHUB_REPO,
                      timeout: float = 7.0) -> dict:
    """Return a diagnostic GitHub Releases update-check result.

    The result always contains ``status`` with one of:

    ``"update"``
        A newer stable GitHub release exists.
    ``"current"``
        The installed version is equal to or newer than the latest release.
    ``"error"``
        The check could not be completed.  This is non-fatal and includes a
        short user-facing ``reason`` suitable for the manual update dialog.

    ``releases/latest`` excludes drafts and prereleases, which is the intended
    behaviour for normal CTA users.  This function never raises.
    """
    current_text = current_version or __version__
    current = _version_tuple(current_text)
    if current is None:
        return dict(status='error', reason=f"Invalid installed version: {current_text}")

    try:
        import requests
        headers = dict(_HTTP_HEADERS)
        headers.update({
            'Accept': 'application/vnd.github+json',
            'X-GitHub-Api-Version': '2022-11-28',
        })
        resp = requests.get(
            f"https://api.github.com/repos/{repo}/releases/latest",
            headers=headers, timeout=timeout)

        if resp.status_code == 404:
            return dict(status='error', reason='No published GitHub release was found.')
        if resp.status_code == 403:
            remaining = resp.headers.get('X-RateLimit-Remaining', '')
            reason = 'GitHub API access was denied.'
            if remaining == '0':
                reason = 'GitHub API rate limit reached. Please try again later.'
            return dict(status='error', reason=reason)
        if resp.status_code != 200:
            return dict(status='error', reason=f'GitHub returned HTTP {resp.status_code}.')

        try:
            data = resp.json()
        except Exception:
            return dict(status='error', reason='GitHub returned an invalid response.')

        tag = str(data.get('tag_name', '') or '')
        latest = _version_tuple(tag)
        if latest is None:
            return dict(status='error', reason=f'Unrecognized release tag: {tag or "(empty)"}')

        base = dict(
            latest='.'.join(map(str, latest)),
            tag=tag,
            url=data.get('html_url') or f'https://github.com/{repo}/releases/latest',
            notes=(data.get('body') or '').strip()[:500],
            installed='.'.join(map(str, current)),
        )
        if latest > current:
            base['status'] = 'update'
        else:
            base['status'] = 'current'
        return base

    except Exception as exc:
        name = exc.__class__.__name__
        if name == 'SSLError':
            reason = 'SSL certificate verification failed.'
        elif name in ('ProxyError',):
            reason = 'The configured network proxy could not be reached.'
        elif name in ('ConnectTimeout', 'ReadTimeout', 'Timeout'):
            reason = 'The GitHub request timed out.'
        elif name in ('ConnectionError',):
            reason = 'GitHub could not be reached from this computer.'
        else:
            reason = f'Update check failed: {exc}' if str(exc) else f'Update check failed ({name}).'
        return dict(status='error', reason=reason)


def check_for_update(current_version: str | None = None,
                      repo: str = GITHUB_REPO,
                      timeout: float = 7.0) -> dict | None:
    """Backward-compatible helper returning only genuine update information.

    Existing callers that expect ``None`` when no update is available can keep
    using this function.  New GUI code should use :func:`get_update_status` so
    a manual check can distinguish "up to date" from a network/API failure.
    """
    result = get_update_status(current_version=current_version,
                               repo=repo, timeout=timeout)
    return result if result.get('status') == 'update' else None


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


def _fetch_from_jpl_sbdb_api(designation: str, timeout: float = 20.0) -> dict:
    """Fetch comet elements from the official JPL SBDB JSON API.

    This is a packaging-safe fallback for frozen executables.  The normal
    path remains astroquery/JPL Horizons, but a PyInstaller bundle can fail
    before the network query if Astropy/Astroquery package data or metadata
    are incomplete.  SBDB uses plain ``requests`` + JSON and returns the six
    heliocentric ecliptic osculating elements needed by CTA: e, q, tp, om,
    w, and i.

    The fallback deliberately reports its own source in the GUI.  It is not
    presented as an observation-epoch Horizons solution.
    """
    import requests

    raw = designation.strip()
    clean = re.sub(r"\s*\([^)]*\)\s*$", "", raw).strip()
    candidates = []
    periodic = re.match(r"^(\d+[PD])", clean, re.I)
    if periodic:
        candidates.append(periodic.group(1).upper())
    candidates.extend([clean, raw])

    seen = set()
    errors = []
    for query in candidates:
        query = query.strip()
        if not query or query.lower() in seen:
            continue
        seen.add(query.lower())
        try:
            resp = requests.get(
                "https://ssd-api.jpl.nasa.gov/sbdb.api",
                params={
                    "sstr": query,
                    "full-prec": "true",
                    "alt-orbits": "false",
                    "phys-par": "false",
                },
                headers=_HTTP_HEADERS, timeout=timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            if data.get("message"):
                raise RuntimeError(str(data["message"]))
            orbit = data.get("orbit") or {}
            elements = orbit.get("elements") or []
            values = {
                str(item.get("name", "")).lower(): item.get("value")
                for item in elements
            }
            required = ("q", "e", "i", "om", "w", "tp")
            missing = [key for key in required
                       if values.get(key) in (None, "", "null")]
            if missing:
                raise RuntimeError(
                    "SBDB response did not contain: " + ", ".join(missing))

            Tp = float(values["tp"])
            obj = data.get("object") or {}
            display_name = (obj.get("fullname") or obj.get("des") or
                            designation)
            return dict(
                q=float(values["q"]),
                e=float(values["e"]),
                i=float(values["i"]),
                Omega=float(values["om"]),
                omega=float(values["w"]),
                T=jd_to_str(Tp)[:10],
                T_jd=Tp,
                name=display_name,
                source="JPL SBDB API (fallback)",
            )
        except Exception as ex:
            errors.append(f"{query!r}: {ex}")

    raise RuntimeError("; ".join(errors) if errors
                       else "No usable SBDB query candidate")


def fetch_from_horizons(designation: str, date: str | None = None) -> dict:
    """Fetch comet orbital elements, preferring JPL Horizons.

    Normal Python/source installations use ``astroquery.jplhorizons``.  If
    that path fails (including frozen-executable package-data failures), CTA
    automatically falls back to the official JPL SBDB JSON API so an
    installer build is not left unable to load orbital elements.
    """
    epoch_jd = date_to_jd(date) if date else today_jd()
    horizons_error: Exception | None = None

    try:
        from astroquery.jplhorizons import Horizons

        def _query(id_str, id_type):
            kw = dict(id=id_str, location="@10", epochs=epoch_jd)
            if id_type:
                kw["id_type"] = id_type
            el = Horizons(**kw).elements()
            if len(el) == 0:
                raise RuntimeError("Empty result")
            return el

        def _extract(el, name=None) -> dict:
            Tp = float(el["Tp_jd"][0])
            return dict(
                q=float(el["q"][0]),
                e=float(el["e"][0]),
                i=float(el["incl"][0]),
                Omega=float(el["Omega"][0]),
                omega=float(el["w"][0]),
                T=jd_to_str(Tp)[:10],
                T_jd=Tp,
                name=name or designation,
                source="JPL Horizons",
            )

        for id_str, id_type in _normalize_comet_desig(designation):
            try:
                return _extract(_query(id_str, id_type))
            except Exception as ex:
                horizons_error = ex
                if "ambiguous" in str(ex).lower():
                    rec_id = _parse_ambiguous_record(str(ex))
                    if rec_id:
                        try:
                            return _extract(_query(rec_id, None),
                                            name=designation)
                        except Exception as ex2:
                            horizons_error = ex2
    except Exception as ex:
        horizons_error = ex

    # Packaging/network-safe fallback.  This keeps the installed program
    # functional even when Astroquery/Astropy cannot locate bundled data.
    try:
        fallback = _fetch_from_jpl_sbdb_api(designation)
        logging.warning(
            "Horizons path failed for %s; using JPL SBDB fallback: %s",
            designation, horizons_error)
        return fallback
    except Exception as sbdb_error:
        raise RuntimeError(
            f"JPL orbital-element query failed for: {designation!r}\n"
            f"Horizons/Astroquery error: {horizons_error}\n"
            f"JPL SBDB fallback error: {sbdb_error}\n\n"
            "Check the internet connection or verify the designation at "
            "https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html"
        ) from sbdb_error


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

                # COBS Observation List API v1.5 returns the observation
                # array under the singular key ``object``. Older payloads
                # may use ``objects``. Accept both forms.
                recs = data.get("object")
                if recs is None:
                    recs = data.get("objects", [])
                if isinstance(recs, dict):
                    recs = [recs]
                if not isinstance(recs, list) or not recs:
                    break
                for rec in recs:
                    try:
                        date_str = str(rec.get("obs_date", ""))[:10]
                        mag = float(rec.get("magnitude", rec.get("mag", 99)))
                        if date_str and 1.0 <= mag <= 22.0:
                            # v3.1 — COBS's obs_list.api actually returns rich
                            # per-observation metadata (observer identity,
                            # method, instrument) that earlier versions threw
                            # away entirely, keeping only date+mag. Capturing
                            # obs_method/observer here lets Q(t) estimation
                            # detect/avoid a specific failure mode: a multi-day
                            # run of reports from ONE observer/method with a
                            # different effective calibration than the rest of
                            # the dataset can look exactly like a real multi-
                            # day production-rate dip to a smoothing filter
                            # (see estimate_qt_from_lightcurve()'s docstring),
                            # since smoothing alone can't tell "real multi-day
                            # event" apart from "multi-day calibration block".
                            method = rec.get("obs_method") or {}
                            observer = rec.get("observer") or {}
                            raw_obs.append({
                                "date": date_str, "mag": mag,
                                "obs_method": method.get("key", "?"),
                                "obs_method_name": method.get("name", "?"),
                                "observer": observer.get("icq_name", "?"),
                            })
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
                                r_helio=r_h, delta=dlt,
                                obs_method=obs.get("obs_method", "?"),
                                obs_method_name=obs.get("obs_method_name", "?"),
                                observer=obs.get("observer", "?")))
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
                        delta  =float(np.linalg.norm(r_C - r_E)),
                        obs_method=obs.get("obs_method", "?"),
                        obs_method_name=obs.get("obs_method_name", "?"),
                        observer=obs.get("observer", "?")))
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


# Schleicher composite dust phase function (Schleicher, Millis & Birch
# 1998-style polynomial fit), ported from py_COMTAILS's constants.py
# (Moreno, F. 2025, A&A 695, A263; python port by R. Morales & N. Robles,
# IAA-CSIC; https://github.com/FernandoMorenoDanvila/py_COMTAILS, MIT
# licence — cite per that repo's README if used in a paper). Used by
# schleicher_phase_correction() below.
SCHLEICHER_COEFFS = np.array([
    -7.4004978e-03, -1.6492566e-02, 1.0950353e-04,
    8.3640600e-07, 1.0157539e-09, -9.6882641e-11, 4.4184372e-13
])


def schleicher_phase_correction(phase_angle_deg: float) -> float:
    """
    v3.1 — Schleicher composite dust phase function Φ(α), evaluated at
    the given solar phase angle α [degrees] (e.g. comet_el-derived
    'Phase' shown in the EPHEMERIS panel). Returns a dimensionless
    scattering-efficiency correction relative to α=0° (so always ≤1 for
    α>0 — dust scatters less efficiently away from full backscatter).

    Same Horner's-method polynomial evaluation as py_COMTAILS's
    setup_coordinate_system() (see SCHLEICHER_COEFFS' own docstring for
    the citation): log10(Φ) = c0 + c1·α + c2·α² + ... + c6·α⁶.

    CLAMPED at α=130° — this degree-6 polynomial fit, like any fit,
    diverges outside the phase-angle range it was actually constrained
    by. Testing showed it blowing up to >1700× at α=180° — symptomatic
    of unconstrained extrapolation, not real forward-scattering physics
    (real cometary dust DOES brighten toward large phase angles via
    forward scattering, but not by checking against this specific
    polynomial at values this far outside its fitted regime). Values
    above 130° use the α=130° correction instead of extrapolating
    further — meaningfully better than a divergent number, but treat
    any phase angle this large as a case where this function's accuracy
    is genuinely uncertain, not silently trustworthy.

    IMPORTANT — this corrects a SINGLE NUMBER, not a per-particle or
    per-pixel quantity: every particle in one Monte Carlo run shares the
    same comet-Sun-Earth geometry (one observation epoch), so phase
    angle is the SAME for the whole simulated frame. Applying this
    therefore does NOT change a single run's contour SHAPE at all — it
    only rescales the whole frame's ABSOLUTE brightness uniformly. Its
    real use is comparing relative brightness ACROSS different runs/
    epochs (different phase angles), not refining one run's morphology.
    """
    c = SCHLEICHER_COEFFS
    x = min(float(phase_angle_deg), 130.0)
    # Horner's method, same order as the source (least->most caution
    # against floating-point error accumulation at large x):
    log_phi = c[0] + x*(c[1] + x*(c[2] + x*(c[3] + x*(c[4] + x*(c[5] + x*c[6])))))
    return float(10.0 ** log_phi)


def linear_exponential_phase_correction(phase_angle_deg: float,
                                        beta_alpha: float = 0.024,
                                        m_oe: float = 0.28,
                                        w_oe: float = 1.5) -> float:
    """
    Generic linear–exponential dust phase law, returned as a relative
    flux correction Φ(α) normalized to Φ(0)=1.

        Δm(α) = β_α α + m_oe [1 − exp(−α / w_oe)]
        Φ(α)  = 10^(−0.4 Δm)

    This is an optional generic phase-law choice for objects whose phase
    behaviour is not well represented by the default Schleicher composite
    dust phase function. Like Schleicher, it is a frame-wide scalar for a
    single epoch and does not change one run's contour shape.
    """
    a = max(0.0, float(phase_angle_deg))
    beta_alpha = float(beta_alpha)
    m_oe = float(m_oe)
    w_oe = max(1e-9, float(w_oe))
    dm = beta_alpha * a + m_oe * (1.0 - np.exp(-a / w_oe))
    return float(10.0 ** (-0.4 * dm))


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


def _rolling_median_smooth(t_jd: np.ndarray, values: np.ndarray,
                           window_days: float) -> np.ndarray:
    """
    v3.1 — robust LOCAL smoothing for a time series: each point is
    replaced by the median of all points within ±window_days/2 days of
    it. Deliberately a LOCAL median, not a global functional fit (e.g.
    the H0,n brightening-law fit elsewhere in this module): a single
    discrepant point (one observer's bad night, common in heterogeneous
    COBS data) gets outvoted by its neighbours and barely moves the
    local estimate, but a GENUINE multi-day event (a real outburst) is
    seen by enough nearby points that the local median still tracks it —
    a global fit would average a real short-lived event away just as
    readily as noise, since by construction it can't represent any
    deviation from its single fitted curve shape at all.
    window_days <= 0 returns values unchanged (smoothing off).
    """
    if window_days <= 0:
        return values.copy()
    half = window_days / 2.0
    smoothed = np.empty_like(values, dtype=float)
    for i in range(len(t_jd)):
        mask = np.abs(t_jd - t_jd[i]) <= half
        smoothed[i] = np.median(values[mask])
    return smoothed


def summarize_obs_methods(obs_list: list) -> dict:
    """
    v3.1 — count how many observations use each obs_method, for display
    (e.g. before deciding whether to filter to one method — see
    filter_obs_by_dominant_method()). Returns {method_key: (count, name)},
    sorted by count descending. method_key is '?' for anything fetched
    before this metadata was captured (older cached results).
    """
    from collections import Counter
    keys = [o.get("obs_method", "?") for o in obs_list]
    names = {o.get("obs_method", "?"): o.get("obs_method_name", "?") for o in obs_list}
    counts = Counter(keys)
    return {k: (c, names.get(k, "?")) for k, c in counts.most_common()}


def filter_obs_by_dominant_method(obs_list: list) -> tuple:
    """
    v3.1 — keep only observations using the SINGLE most common obs_method
    in obs_list, discarding the rest.

    Why this matters here specifically: different observation methods
    (visual, CCD, CCD-visual-equivalent, etc.) can carry different
    effective zero-points, and COBS observers tend to report in MULTI-DAY
    CLUSTERS (the same person observing several nights running, then
    stopping) — so a calibration offset between methods can show up as a
    multi-day brightness "step" in the raw light curve, which a smoothing
    filter alone cannot distinguish from a real multi-day production-rate
    change (see estimate_qt_from_lightcurve()'s and _rolling_median_
    smooth()'s docstrings — smoothing only rejects SINGLE-point outliers,
    not sustained multi-day blocks, by design). Restricting to one method
    removes this specific failure mode by construction, at the cost of
    discarding however many points used other methods.

    Returns (filtered_obs_list, method_counts) — method_counts is
    summarize_obs_methods()'s output for the ORIGINAL (unfiltered) list,
    so a caller can show what was kept vs. discarded.
    """
    counts = summarize_obs_methods(obs_list)
    if not counts:
        return obs_list, {}
    dominant_key = next(iter(counts))   # counts is already sorted by count desc
    filtered = [o for o in obs_list if o.get("obs_method", "?") == dominant_key]
    return filtered, counts


def estimate_qt_from_lightcurve(comet_el: dict, obs_list: list,
                                afrho_anchors: list,
                                p_v: float = 0.04,
                                v_dust_law: tuple = (0.1, -0.5),
                                smooth_window_days: float = 3.0) -> dict:
    """
    v3.1 — estimate a time-dependent dust production rate Q(t) [kg/s] for
    use as the Monte Carlo morphology sampler's release-time weighting
    (compute_morphology_mc()'s qt_weights=, replacing its default uniform
    Q(t)), by combining:

      - obs_list: heterogeneous COBS magnitude(t) (e.g. fetch_from_cobs()'s
        'obs_list' key — each dict(date, mag, r_helio, delta)). Gives the
        RELATIVE SHAPE of brightness over a long baseline, but COBS does
        NOT report a per-observation photometric aperture — most entries
        are visual estimates with no rigorously defined aperture at all —
        so this alone canNOT be converted into an absolute Afρ.
      - afrho_anchors: [(date_str, afrho_cm), …] — YOUR OWN Afρ
        measurement(s), each tied to one specific date, computed with a
        KNOWN aperture from your own image photometry (e.g. via Tycho
        Tracker). Used to calibrate the ABSOLUTE SCALE of the COBS-derived
        shape. At least one is required.

    Method (mirrors how Moreno/COMTAILS combines COBS's long-baseline
    total-magnitude data with a dedicated photometric network's aperture-
    defined Afρ measurements — see Moreno et al. 2025, MNRAS 539, 949,
    Section 2, who do exactly this with COBS + Cometas_Obs):
      1. afrho_proxy(t) ∝ Δ(t)²·r_h(t)²·10^(−0.4·m(t)) at every COBS point
         — this is Afρ's OWN defining relation (A'Hearn et al. 1984) with
         the aperture/constant terms left out, since they're unknown per-
         observation and cancel out in the calibration step below anyway.
      2. The proxy curve (log-interpolated over the COBS dates) is
         evaluated AT each anchor's own date, giving scale = Afρ_anchor /
         Afρ_proxy(anchor_date). With one anchor, that single scale is
         applied to the whole curve (shape trusted as given). With ≥2
         anchors, the (log) scale is linearly interpolated BETWEEN anchor
         dates and held constant outside their range — letting widely-
         spaced anchors correct gradual drift in the proxy, not just one
         overall multiplicative offset.
      3. Afρ_calibrated(t) = scale(t) · Afρ_proxy(t).
      4. Q(t) = compute_dust_production_rate(Afρ_calibrated(t), r_h(t),
         v_dust(t), p_v)['Qd_rough'], with v_dust(t) = v_dust_law[0] ·
         r_h(t)^v_dust_law[1] — same "uncited rule of thumb" default
         (0.1 km/s · r_h^−0.5) as compute_dust_production_rate()'s own
         GUI calculator, for consistency; override if you have a better
         v_dust(t) for this object (e.g. from the ejection-velocity fit).

    CAVEATS — real limitations, not just a disclaimer:
      - Assumes COBS's heterogeneous (mostly visual) magnitude estimates
        carry a roughly TIME-INVARIANT effective aperture/calibration
        relationship. Real visual apertures vary with observer, sky
        brightness, and instrument — exactly the scatter visible in any
        COBS light curve. This is a shape-only approximation, not a
        rigorous per-observation aperture correction.
      - The absolute scale is only as good as your anchor measurement(s)
        — garbage in, garbage out there. The relative SHAPE between
        anchors still rests on the time-invariant-aperture assumption
        above, even with several anchors.
      - This is the simpler, first-order route flagged in the v3.1
        changelog — not a substitute for full forward-fitting a Monte
        Carlo synthetic image/photometry against the data the way
        Moreno (2025)'s COMTAILS actually does (iteratively adjusting
        Q(t), size distribution, etc. until BOTH isophote shape AND the
        Afρ/magnitude time series match, across many epochs at once).

    smooth_window_days : applies _rolling_median_smooth() to the raw COBS
        magnitudes before anything else, to absorb point-to-point
        observer/method scatter without erasing a genuine multi-day
        event (see that function's docstring for why median + a SHORT
        local window does this and a global brightening-law fit
        wouldn't). Default 3.0 d; pass 0 to use the raw magnitudes
        unsmoothed.

    Returns dict(t_jd, r_helio, delta, mag, afrho_proxy, afrho_cal,
                 Q_kg_s, anchor_jd, anchor_afrho, anchor_log_scale).
    Raises ValueError if obs_list or afrho_anchors is empty.
    """
    if not obs_list:
        raise ValueError("obs_list is empty — fetch a COBS light curve first.")
    if not afrho_anchors:
        raise ValueError("Need at least one (date, afrho_cm) anchor measurement.")

    pts = sorted(obs_list, key=lambda o: o["date"])
    t_jd  = np.array([date_to_jd(o["date"]) for o in pts])
    r_h   = np.array([o["r_helio"] for o in pts])
    delta = np.array([o["delta"]   for o in pts])
    mag_raw = np.array([o["mag"]   for o in pts])
    mag     = _rolling_median_smooth(t_jd, mag_raw, smooth_window_days)

    afrho_proxy = delta ** 2 * r_h ** 2 * 10.0 ** (-0.4 * mag)

    order = np.argsort(t_jd)
    t_sorted    = t_jd[order]
    log_proxy_sorted = np.log10(np.clip(afrho_proxy[order], 1e-300, None))

    anchor_jd_raw    = np.array([date_to_jd(d) for d, _ in afrho_anchors])
    anchor_afrho_raw = np.array([float(a) for _, a in afrho_anchors])
    proxy_at_anchor  = 10.0 ** np.interp(anchor_jd_raw, t_sorted, log_proxy_sorted)
    log_scale_raw    = np.log10(np.clip(anchor_afrho_raw, 1e-300, None) / proxy_at_anchor)

    a_order = np.argsort(anchor_jd_raw)
    anchor_jd       = anchor_jd_raw[a_order]
    anchor_afrho    = anchor_afrho_raw[a_order]
    anchor_log_scale = log_scale_raw[a_order]

    if len(anchor_jd) == 1:
        log_scale_t = np.full_like(t_jd, anchor_log_scale[0])
    else:
        log_scale_t = np.interp(t_jd, anchor_jd, anchor_log_scale)

    afrho_cal = afrho_proxy * 10.0 ** log_scale_t

    v0_law, gamma_law = v_dust_law
    v_dust_kms = v0_law * r_h ** gamma_law

    Q_kg_s = np.array([
        compute_dust_production_rate(float(a), float(rh), float(vd), p_v)["Qd_rough"]
        for a, rh, vd in zip(afrho_cal, r_h, v_dust_kms)
    ])

    return dict(t_jd=t_jd, r_helio=r_h, delta=delta, mag=mag, mag_raw=mag_raw,
               afrho_proxy=afrho_proxy, afrho_cal=afrho_cal, Q_kg_s=Q_kg_s,
               anchor_jd=anchor_jd, anchor_afrho=anchor_afrho,
               anchor_log_scale=anchor_log_scale)



def assess_qt_coverage(t_jd, q_values, obs_jd: float,
                       max_age_days: float,
                       dominant_age_days: float = 0.0,
                       smooth_window_days: float = 3.0,
                       min_points_reliable: int = 8,
                       min_points_partial: int = 5) -> dict:
    """Assess whether a Q(t) curve can support an MC-window inference.

    This function is intentionally conservative.  It is used for *inference*
    (the suggested MC release window), not merely for keeping the Monte Carlo
    sampler numerically supplied with weights.  In particular it:

      1. clips the data to the F-P maximum-dust-age interval,
      2. collapses duplicate observing dates by taking the median Q,
      3. keeps only the most recent continuous observing segment,
      4. rejects interpolation across a long gap or across another apparition,
      5. grades the segment as RELIABLE / PARTIAL / INSUFFICIENT.

    The current MC sampler can flat-extrapolate edge weights.  That behaviour
    is useful as a numerical fallback, but it must not create an apparently
    precise release-window recommendation.  Callers should therefore use the
    returned ``segment_t_jd`` / ``segment_q`` for recommendation statistics,
    and enable an Apply/Use button only when ``recommendation_allowed`` is True.

    Parameters
    ----------
    t_jd, q_values
        Sampled Q(t) curve.  Duplicate dates are allowed.
    obs_jd
        Observation epoch.
    max_age_days
        F-P maximum dust age, normally the largest listed synchrone age.
    dominant_age_days
        Optional representative age of the dominant F-P dust structure.
    smooth_window_days
        Used only to set a sensible minimum gap/recency tolerance.

    Returns
    -------
    dict
        Includes quality, reason, unique_points, gap_limit_days,
        latest_age_days, reliable_lookback_days, coverage_fraction,
        largest_gap_days, segment_t_jd, segment_q, and sampling arrays.
    """
    result = dict(
        quality="INSUFFICIENT",
        recommendation_allowed=False,
        reason="",
        unique_points=0,
        gap_limit_days=0.0,
        recency_limit_days=0.0,
        latest_age_days=float("inf"),
        reliable_lookback_days=0.0,
        coverage_fraction=0.0,
        largest_gap_days=0.0,
        dominant_age_covered=False,
        segment_t_jd=np.array([], dtype=float),
        segment_q=np.array([], dtype=float),
        sampling_t_jd=np.array([], dtype=float),
        sampling_q=np.array([], dtype=float),
    )

    try:
        max_age = float(max_age_days)
        dominant_age = max(0.0, float(dominant_age_days or 0.0))
        obs = float(obs_jd)
    except Exception:
        result["reason"] = "Invalid F-P age guidance or observation epoch."
        return result

    if not np.isfinite(max_age) or max_age <= 0:
        result["reason"] = (
            "F-P maximum dust age is unavailable. Define valid synchrone ages "
            "before requesting a COBS-based MC-window recommendation.")
        return result

    t = np.asarray(t_jd, dtype=float)
    q = np.asarray(q_values, dtype=float)
    m = (np.isfinite(t) & np.isfinite(q) & (q > 0) &
         (t <= obs) & (t >= obs - max_age))
    t = t[m]
    q = q[m]
    if t.size < 2:
        result["reason"] = (
            f"Fewer than two valid COBS/Q(t) points fall within the last "
            f"{max_age:.0f} days.")
        return result

    # Collapse duplicate observing epochs.  Several observers can report on
    # the same date; treating every report as an independent time interval
    # would overweight that date in coverage diagnostics.
    order = np.argsort(t)
    t = t[order]
    q = q[order]
    unique_t, inverse = np.unique(t, return_inverse=True)
    unique_q = np.empty_like(unique_t, dtype=float)
    for i in range(unique_t.size):
        unique_q[i] = float(np.median(q[inverse == i]))
    t = unique_t
    q = unique_q
    result["unique_points"] = int(t.size)
    if t.size < 2:
        result["reason"] = "Only one unique observing date is available in the F-P age interval."
        return result

    gaps = np.diff(t)
    positive_gaps = gaps[gaps > 0]
    # Use local cadence, not a multi-month/apparition gap, to establish the
    # threshold.  The threshold is bounded so daily data are not fragmented by
    # weather gaps, while separate apparitions are never bridged.
    local_gaps = positive_gaps[positive_gaps <= 60.0]
    if local_gaps.size:
        cadence = float(np.median(local_gaps))
    elif positive_gaps.size:
        cadence = float(np.median(positive_gaps))
    else:
        cadence = 1.0
    gap_limit = float(np.clip(max(14.0,
                                  4.0 * cadence,
                                  3.0 * max(0.0, float(smooth_window_days))),
                              14.0, 45.0))
    result["gap_limit_days"] = gap_limit

    # Starting with the point nearest the observation, walk backwards until a
    # gap exceeds the adaptive limit.  Older disconnected data are retained in
    # the full Q(t) plot if desired, but excluded from inference.
    start = t.size - 1
    for i in range(t.size - 2, -1, -1):
        if (t[i + 1] - t[i]) > gap_limit:
            break
        start = i
    seg_t = t[start:]
    seg_q = q[start:]

    latest_age = max(0.0, obs - float(seg_t[-1]))
    reliable_lookback = max(0.0, obs - float(seg_t[0]))
    segment_span = max(0.0, float(seg_t[-1] - seg_t[0]))
    seg_gaps = np.diff(seg_t)
    largest_gap = float(np.max(seg_gaps)) if seg_gaps.size else 0.0
    coverage_fraction = min(1.0, reliable_lookback / max_age)
    recency_limit = float(min(30.0, max(10.0,
                                       3.0 * max(0.0, float(smooth_window_days)),
                                       0.15 * max_age)))
    recent_ok = latest_age <= recency_limit
    dominant_covered = (dominant_age <= 0.0 or
                        reliable_lookback + 1e-9 >= dominant_age)

    result.update(
        latest_age_days=latest_age,
        reliable_lookback_days=reliable_lookback,
        coverage_fraction=coverage_fraction,
        largest_gap_days=largest_gap,
        recency_limit_days=recency_limit,
        dominant_age_covered=bool(dominant_covered),
        segment_t_jd=seg_t,
        segment_q=seg_q,
    )

    n = int(seg_t.size)
    reliable = (recent_ok and dominant_covered and
                n >= int(min_points_reliable) and
                coverage_fraction >= 0.70)
    partial = (recent_ok and dominant_covered and
               n >= int(min_points_partial) and
               coverage_fraction >= 0.35)

    if reliable:
        quality = "RELIABLE"
        allowed = True
        reason = (
            f"Recent continuous COBS segment covers {coverage_fraction*100:.0f}% "
            f"of the {max_age:.0f}-day F-P interval.")
    elif partial:
        quality = "PARTIAL"
        allowed = False
        reason = (
            f"COBS constrains only the most recent {reliable_lookback:.0f} days "
            f"({coverage_fraction*100:.0f}% of the {max_age:.0f}-day F-P interval). "
            "Statistics are provisional; automatic Apply is disabled.")
    else:
        quality = "INSUFFICIENT"
        allowed = False
        problems = []
        if not recent_ok:
            problems.append(
                f"latest valid COBS point is {latest_age:.0f} days before the observation "
                f"(limit {recency_limit:.0f} days)")
        if n < int(min_points_partial):
            problems.append(f"only {n} unique dates are in the recent segment")
        if coverage_fraction < 0.35:
            problems.append(
                f"recent segment covers only {coverage_fraction*100:.0f}% of the F-P interval")
        if not dominant_covered:
            problems.append(
                f"recent segment does not reach the dominant F-P age ({dominant_age:.0f} d)")
        reason = "Insufficient continuous COBS coverage"
        if problems:
            reason += ": " + "; ".join(problems)
        reason += "."

    result["quality"] = quality
    result["recommendation_allowed"] = allowed
    result["reason"] = reason

    # Sampling arrays use only the accepted recent segment.  A small final
    # edge from the latest observation to the exact image epoch is permitted
    # only when the point passed the recency test; this avoids a several-year
    # extrapolation while keeping the sampler defined at age=0.
    sampling_t = seg_t.copy()
    sampling_q = seg_q.copy()
    if recent_ok and sampling_t.size and sampling_t[-1] < obs:
        sampling_t = np.append(sampling_t, obs)
        sampling_q = np.append(sampling_q, sampling_q[-1])
    result["sampling_t_jd"] = sampling_t
    result["sampling_q"] = sampling_q
    result["segment_span_days"] = segment_span
    return result

def find_qt_dips(t_jd: np.ndarray, Q_kg_s: np.ndarray,
                 baseline_window_days: float = 15.0,
                 dip_threshold: float = 0.5,
                 min_duration_days: float = 1.5) -> list:
    """
    v3.1 — locate calendar date ranges where Q(t) drops well below its
    own LOCAL surroundings, so a dip showing up as a disconnected loop in
    the Monte Carlo contour (see compute_morphology_mc()'s qt_weights=)
    can be checked directly against a date — e.g. against CBET/comet
    mailing-list reports — instead of reading a date back off the
    rendered contour image by eye.

    Method: a WIDE rolling median (baseline_window_days, deliberately
    much wider than estimate_qt_from_lightcurve()'s own smoothing window)
    gives a "what Q(t) would be without short-lived dips" baseline at
    each point. Points where Q(t)/baseline < dip_threshold are flagged;
    consecutive flagged points (allowing gaps of up to 3x the typical
    sample spacing, so sparse sampling doesn't split one real dip into
    several) are grouped into events, and events shorter than
    min_duration_days are dropped (too brief to be more than a single
    sparse low point even after grouping).

    Returns a list of dicts, each: date_start, date_end, duration_days,
    min_ratio (the dip's deepest Q/baseline fraction), sorted by
    duration_days descending (most significant first).
    """
    if len(t_jd) < 5:
        return []
    order = np.argsort(t_jd)
    t = t_jd[order]; Q = Q_kg_s[order]

    baseline = _rolling_median_smooth(t, Q, baseline_window_days)
    baseline = np.clip(baseline, 1e-30, None)
    ratio = Q / baseline
    flagged = ratio < dip_threshold

    if not np.any(flagged):
        return []

    # Typical sample spacing, to set how big a gap can be bridged when
    # grouping flagged points into one event (sparse sampling shouldn't
    # fragment one real dip into many tiny "events").
    spacing = np.median(np.diff(t)) if len(t) > 1 else 1.0
    max_gap = max(3.0 * spacing, 2.0)

    events = []
    idx = np.where(flagged)[0]
    start = idx[0]
    prev = idx[0]
    for i in idx[1:]:
        if t[i] - t[prev] > max_gap:
            events.append((start, prev))
            start = i
        prev = i
    events.append((start, prev))

    results = []
    for s, e in events:
        duration = t[e] - t[s]
        if duration < min_duration_days:
            continue
        results.append(dict(
            date_start=jd_to_str(t[s])[:10],
            date_end=jd_to_str(t[e])[:10],
            duration_days=round(float(duration), 1),
            min_ratio=round(float(ratio[s:e+1].min()), 2),
        ))
    results.sort(key=lambda d: d["duration_days"], reverse=True)
    return results


def _sample_weighted_ages(t_jd_grid, weight_grid, obs_jd: float,
                          max_age: float, n: int, rng) -> np.ndarray:
    """
    v3.1 — inverse-transform-sample n release "ages" (days before obs_jd,
    matching compute_morphology_mc()'s existing convention) from an
    arbitrary empirical weight(t) function — e.g. the Q_kg_s curve from
    estimate_qt_from_lightcurve() — instead of compute_morphology_mc()'s
    default uniform release-time window.

    weight_grid is linearly interpolated (np.interp, which already
    flat-extrapolates) onto a fine internal grid spanning exactly
    [obs_jd-max_age, obs_jd] — so a Q(t) estimate that doesn't fully cover
    the requested max_age window extends flatly from its nearest edge
    value rather than dropping to zero, the same spirit as estimate_qt_
    from_lightcurve()'s own anchor extrapolation.
    Falls back to plain uniform sampling (silently) if the weight is
    degenerate (all zero/negative) over the window — never raises.
    """
    t_min, t_max = obs_jd - max_age, obs_jd
    fine_t = np.linspace(t_min, t_max, 2000)
    fine_w = np.clip(np.interp(fine_t, t_jd_grid, weight_grid), 0.0, None)
    cum = np.cumsum(fine_w)
    if cum[-1] <= 0:
        return rng.uniform(0.0, max_age, n)
    cum = cum / cum[-1]
    u = rng.uniform(0.0, 1.0, n)
    sampled_t = np.interp(u, cum, fine_t)
    return obs_jd - sampled_t   # convert "absolute jd released" -> "age"


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


def _orbital_plane_basis(el: dict):
    """
    v3.1 — the perifocal basis vectors (P, Q, W) of the comet's orbit,
    expressed in heliocentric ecliptic XYZ. Same P, Q construction
    elem_to_state() already computes internally — duplicated here in
    standalone form since elem_to_state() only returns the resulting
    state vector, not these basis vectors themselves. P points toward
    periapsis (rotated into the ecliptic), Q is 90° ahead of P within
    the orbital plane, and W = P×Q is the orbit normal (same direction
    as the orbital angular momentum vector). Used by
    _active_area_direction() to convert a body-fixed active-area
    ejection direction into heliocentric XYZ.
    """
    iR, OR, oR = np.radians(el['i']), np.radians(el['Omega']), np.radians(el['omega'])
    cO, sO = np.cos(OR), np.sin(OR)
    cI, sI = np.cos(iR), np.sin(iR)
    co, so = np.cos(oR), np.sin(oR)
    P = np.array([cO*co - sO*so*cI,    sO*co + cO*so*cI,    so*sI])
    Q = np.array([-(cO*so + sO*co*cI), -(sO*so - cO*co*cI), co*sI])
    W = vcrs(P, Q)
    return P, Q, W


def _true_anomaly_at_jd(el: dict, jd: float, mu: float = MU) -> float:
    """
    v3.1 — true anomaly f [radians] at the given JD. Same Kepler-
    equation solving as elem_to_state() (duplicated here in this
    lightweight form since elem_to_state() doesn't expose f directly,
    only the resulting state vector). Needed by the active-area
    ejection model, which requires the comet's true anomaly at each
    particle's EJECTION time to know where periapsis was relative to
    the body-fixed reference used for the active area's rotation phase
    (see _active_area_direction()).
    """
    q, e = el['q'], el['e']
    T = el.get('T_jd', date_to_jd(el['T']))
    dt = jd - T
    eps_e = 1e-5
    if abs(e - 1.0) < eps_e:
        Wb = 3.0 * np.sqrt(mu / (2 * q ** 3)) * dt
        s  = np.cbrt(Wb + np.sqrt(Wb * Wb + 1))
        tf = s - 1.0 / s
        return float(2.0 * np.arctan(tf))
    elif e < 1.0:
        a = q / (1 - e)
        n = np.sqrt(mu / a ** 3)
        E = kepler_elliptic(n * dt, e)
        cf = (np.cos(E) - e) / (1 - e * np.cos(E))
        sf = np.sqrt(1 - e * e) * np.sin(E) / (1 - e * np.cos(E))
        return float(np.arctan2(sf, cf))
    else:
        ah = q / (e - 1)
        n = np.sqrt(mu / ah ** 3)
        H = kepler_hyperbolic(n * dt, e)
        cf = (e - np.cosh(H)) / (e * np.cosh(H) - 1)
        sf = np.sqrt(e * e - 1) * np.sinh(H) / (e * np.cosh(H) - 1)
        return float(np.arctan2(sf, cf))


def _active_area_direction(el: dict, t_emit: float, per_jd: float,
                           nuc_inc: float, nuc_phi: float, period: float,
                           lat_min: float, lat_max: float,
                           lon_min: float, lon_max: float,
                           isun: bool, rng, max_tries: int = 200):
    """
    v3.1 — port of py_COMTAILS's DustTail._anisot_dir2() plus the
    iejec_mode==3 ("active areas") branch of _get_ejection_velocity()
    (Moreno, F. 2025, A&A, 695, A263; python port by R. Morales &
    N. Robles, IAA-CSIC; https://github.com/FernandoMorenoDanvila/
    py_COMTAILS, MIT licence — cite per that repo's README if this is
    used in a paper). Rewritten against CTA's own orbital-element/
    state-vector conventions rather than ported line-for-line: py_
    COMTAILS builds its own helio_matrix/hpo_to_he for the perifocal-
    to-ecliptic rotation, where CTA's elem_to_state() already performs
    the equivalent rotation via the P, Q perifocal basis vectors (see
    _orbital_plane_basis()) — that part is reused, not reimplemented.

    Physical picture: the nucleus is a rotating sphere with obliquity
    nuc_inc and a known subsolar-meridian phase nuc_phi at perihelion.
    A rectangular "active area" is fixed in body coordinates between
    [lat_min,lat_max] x [lon_min,lon_max] (radians). At the particle's
    ejection time, this samples a point UNIFORMLY BY SURFACE AREA
    within that rectangle (the arccos trick on cos(colatitude) below —
    uniform-by-angle would over-sample near the poles), and — when isun
    is True — keeps re-sampling until the point lands on the
    ILLUMINATED part of the area ("particle ejection occurs only from
    the illuminated portion of the selected active area at each
    integration time step", Moreno et al. 2012, ApJ 752, 136, §4, using
    the same isun-gated-rejection approach).

    Returns (direction_xyz, cosz): direction_xyz is a heliocentric-
    ecliptic XYZ UNIT vector (not yet scaled by velocity magnitude —
    pass to dust_position_isotropic() as `direction`); cosz is the
    local solar-zenith-angle cosine at the sampled point, intended for
    the caller to apply as an additional v_ej *= cosz**expocos
    weighting (same as the source's cosz**config.expocos — illumination
    grazing the surface ejects slower than a head-on subsolar point).
    Returns (None, None) if no illuminated point turned up within
    max_tries — a real possibility (not a bug) if this particular
    active area never faces the Sun at this ejection time/rotation
    phase, not just an edge case to paper over.
    """
    f = _true_anomaly_at_jd(el, t_emit)
    cosfinu, sinfinu = np.cos(nuc_phi + f), np.sin(nuc_phi + f)
    cosi, sini = np.cos(nuc_inc), np.sin(nuc_inc)

    a11, a12, a13 = cosfinu, cosi * sinfinu, sini * sinfinu
    a21, a22, a23 = -sinfinu, cosi * cosfinu, sini * cosfinu
    a31, a32, a33 = 0.0, sini, -cosi

    theta0   = np.arctan(np.tan(f + nuc_phi) * cosi)
    theta0_t = np.arctan(np.tan(nuc_phi) * cosi)
    theta_base = (2.0 * np.pi / period) * (t_emit - per_jd) + theta0_t - theta0

    P, Q, W = _orbital_plane_basis(el)
    cf, sf = np.cos(f), np.sin(f)
    cl1, cl2 = np.cos(np.pi / 2 + lat_max), np.cos(np.pi / 2 + lat_min)

    for _ in range(max_tries):
        lon = lon_min + (lon_max - lon_min) * rng.random()
        lat = np.arccos((cl2 - cl1) * rng.random() + cl1) - np.pi / 2

        theta = theta_base + lon
        v1 = np.cos(lat) * np.cos(theta + theta0)
        v2 = np.cos(lat) * np.sin(theta + theta0)
        v3 = np.sin(lat)

        ur     = -(a11 * v1 + a12 * v2 + a13 * v3)
        utheta = -(a21 * v1 + a22 * v2 + a23 * v3)
        uz     = -(a31 * v1 + a32 * v2 + a33 * v3)
        cosz = -ur

        if cosz >= 0.0 or not isun:
            xop = ur * cf - utheta * sf
            yop = ur * sf + utheta * cf
            d_xyz = P * xop + Q * yop + W * uz
            return d_xyz, max(cosz, 0.0)
    return None, None


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
def dust_position(beta: float, t_emit: float, t_obs: float, comet_el: dict,
                   v_R0: float = 0.0, v_T0: float = 0.0, v_N0: float = 0.0,
                   gamma: float = 0.0, m_exp: float = 0.0):
    """
    Compute the position of a dust grain at t_obs, emitted at t_emit
    from the nucleus. Returns 3-vector [AU] in heliocentric ecliptic,
    or None on failure.

    Ejection velocity (v3.1, all optional, default 0.0 ⇒ exact v3.0
    zero-ejection-velocity behaviour):
        v_dust,0 = v_comet(t_emit) + v_ej
        v_ej     = v_R·R̂ + v_T·T̂ + v_N·N̂
        v_R = v_R0 · β^γ · r_H(t_emit)^(−m_exp)     (likewise v_T, v_N)
    where, evaluated at the comet's state at t_emit:
        R̂ = r̂                       (radial, sunward)
        N̂ = (r×v)/|r×v|              (orbit normal — right-handed
                                       automatically for retrograde orbits,
                                       since it's derived from this comet's
                                       own r,v rather than assumed)
        T̂ = N̂ × R̂                   (completes the right-handed triad)
    γ=m_exp=0 ⇒ constant ejection velocity (v_R0/v_T0/v_N0 used directly).
    γ=m_exp=0.5 ⇒ the standard literature form |v_ej| = V0·(β/r_H)^0.5
    (Whipple 1950, ApJ 111, 375; Ishiguro 2008, Icarus 193, 96; written
    out explicitly in Hui et al. 2020, AJ 159, 77, Eq. 5 — verified by
    reading that paper directly). Independent γ/m_exp (rather than
    forcing m_exp=γ) follows Ishiguro et al.'s own more general 2-exponent
    form, e.g. v_ej = V0·(a/a0)^(−u1)·(r_h/AU)^(−u2) as used in the 67P
    dust-trail paper, where u1, u2 are stated as independently free
    (conventionally both set to 0.5, but not required to be equal).
    See the v3.1 changelog entry above for the full literature basis and
    the caveat on physically-plausible velocity magnitudes.
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

    # ── Non-zero ejection velocity kick (v3.1) ───────────────────────────
    has_ej       = (v_R0 != 0.0 or v_T0 != 0.0 or v_N0 != 0.0)
    applied_kick = False
    if has_ej:
        if beta > 0:
            beta_term = beta ** gamma
        elif gamma == 0.0:
            beta_term = 1.0   # constant-velocity mode is valid even at β=0
        else:
            beta_term = 0.0   # power law: v_ej → 0 as β→0 for γ>0
                               # (ungasdraged/large grains aren't entrained)
        r_H_emit = vmag(r0)
        scale    = beta_term * (r_H_emit ** (-m_exp))
        if scale != 0.0:
            R_hat = vhat(r0)
            h_vec = vcrs(r0, v0)
            N_hat = vhat(h_vec)
            T_hat = vcrs(N_hat, R_hat)
            v_ej  = (v_R0 * R_hat + v_T0 * T_hat + v_N0 * N_hat) \
                    * scale * MPS_TO_AUDAY
            v0    = v0 + v_ej
            applied_kick = True

    # β = 0 (and no net ejection kick applied): particle stays on nucleus orbit
    if abs(beta) < 1e-9 and not applied_kick:
        try:
            r, _ = elem_to_state(comet_el, t_obs)
            return r
        except Exception:
            return None

    return _propagate_dust(r0, v0, beta, t_emit, t_obs)


def _propagate_dust(r0, v0, beta: float, t_emit: float, t_obs: float):
    """
    Core orbit-determination/propagation step, factored out of
    dust_position() (v3.1 refactor, no behaviour change) so the
    deterministic F-P path (dust_position(), fixed R/T/N ejection
    direction) and the Monte Carlo path (dust_position_isotropic(),
    random ejection direction — see below) share EXACTLY the same
    physics from this point on. The two differ only in how v0 was
    perturbed before this call is made; they can never diverge in how
    the resulting state is propagated, by construction.
    r0, v0 must already include any ejection-velocity kick. Returns a
    3-vector position [AU] or None on failure.
    """
    dt = t_obs - t_emit
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


def dust_position_isotropic(beta: float, t_emit: float, t_obs: float, comet_el: dict,
                            v0_coeff: float, gamma: float, m_exp: float,
                            direction):
    """
    Same trajectory physics as dust_position() (they share
    _propagate_dust() — see its docstring), but the ejection-velocity
    VECTOR's direction is supplied externally instead of being fixed by
    R̂/T̂/N̂ coefficients. Built for the Monte Carlo morphology sampler
    (v3.1, Phase 1; see compute_morphology_mc()), where each simulated
    grain gets an independently-sampled (typically isotropic-random)
    direction — dust_position()'s deterministic single direction is the
    right model for "one syndyne/synchrone test curve", but the wrong
    one for "a whole population scattered in every direction".

    Magnitude follows the same power law as dust_position():
        |v_ej| = v0_coeff · β^γ · r_H(t_emit)^(−m_exp)
    direction : any 3-vector (will be normalized here) — e.g. one row
        from _sample_isotropic_directions(). NOT required to be R/T/N-
        relative; isotropic-on-a-sphere sampling is rotation-invariant,
        so a direction sampled directly in heliocentric ecliptic XYZ is
        statistically identical to one sampled per-particle in R/T/N.
    v0_coeff=0.0 ⇒ zero ejection velocity, identical to the β-only F-P
        trajectory (direction has no effect).
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

    applied_kick = False
    if v0_coeff != 0.0:
        if beta > 0:
            beta_term = beta ** gamma
        elif gamma == 0.0:
            beta_term = 1.0
        else:
            beta_term = 0.0
        r_H_emit = vmag(r0)
        scale    = beta_term * (r_H_emit ** (-m_exp))
        if scale != 0.0:
            d_hat = vhat(np.asarray(direction, dtype=float))
            v_ej  = v0_coeff * scale * d_hat * MPS_TO_AUDAY
            v0    = v0 + v_ej
            applied_kick = True

    if abs(beta) < 1e-9 and not applied_kick:
        try:
            r, _ = elem_to_state(comet_el, t_obs)
            return r
        except Exception:
            return None

    return _propagate_dust(r0, v0, beta, t_emit, t_obs)


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


def _rotate_sky(xi, eta, rotation_offset_deg):
    """Rotate East/North sky-plane offsets by CTA's grid-rotation convention.

    ``rotation_offset_deg`` is the value from the main GUI's
    "Grid rotation (match observed tail)" control.  Positive values rotate
    the model clockwise on the sky, which is a negative angle in the usual
    mathematical (counter-clockwise-positive) x/y rotation convention.

    Parameters may be scalars or NumPy-compatible arrays.  Returned values
    preserve the broadcast shape of the inputs.
    """
    xi = np.asarray(xi)
    eta = np.asarray(eta)
    theta = np.radians(-float(rotation_offset_deg))
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    return (xi * cos_t - eta * sin_t,
            xi * sin_t + eta * cos_t)


def _orient_to_screen(xi, eta, north_pa_deg):
    """Convert sky offsets (East, North) to CTA image-screen offsets.

    ``north_pa_deg`` is the clockwise angle from image-up to celestial
    North.  The returned ``dx, dy`` follow the image convention used by
    :func:`sky_to_pixel`: x increases to the right and y increases downward.
    No pixel scale or nucleus offset is applied here, so the function can be
    reused by the Monte Carlo density-grid code in AU units.

    Parameters may be scalars or NumPy-compatible arrays.  Returned values
    preserve the broadcast shape of the inputs.
    """
    xi = np.asarray(xi)
    eta = np.asarray(eta)
    th = np.radians(float(north_pa_deg))
    cos_t, sin_t = np.cos(th), np.sin(th)
    dx = xi * (-cos_t) + eta * sin_t
    dy = xi * (-sin_t) + eta * (-cos_t)
    return dx, dy


def sky_to_pixel(xi, eta, nuc_x, nuc_y, au_per_px, north_pa_deg):
    """
    Convert sky offsets (AU East, AU North) → image pixel coordinates.
    north_pa_deg: clockwise angle from image-up to celestial North.
    """
    dx, dy = _orient_to_screen(xi, eta, north_pa_deg)
    return nuc_x + dx / au_per_px, nuc_y + dy / au_per_px


# ─────────────────────────────────────────────────────────────────────────────
#  MODEL COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────
def compute_model(comet_el: dict, obs_jd: float,
                  beta_values, sync_ages,
                  max_age: float = 200, n_pts: int = 200,
                  ejection: dict | None = None) -> dict:
    """
    Compute syndyne and synchrone curves for the given comet at obs_jd.
    Returns a dict with keys: syndynes, synchrones, sun_dir, orbit, info.

    ejection (v3.1, optional): dict with any of
        v_R0, v_T0, v_N0  — ejection-velocity coefficients [m/s]
        gamma, m_exp      — power-law exponents on β and r_H(t_emit)
    Missing keys default to 0.0. ejection=None (the default) is exactly
    the v3.0 zero-ejection-velocity model — see dust_position() for the
    full formula and literature basis.
    """
    ej    = ejection or {}
    v_R0  = ej.get('v_R0',  0.0)
    v_T0  = ej.get('v_T0',  0.0)
    v_N0  = ej.get('v_N0',  0.0)
    gamma = ej.get('gamma', 0.0)
    m_exp = ej.get('m_exp', 0.0)

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
            r_D = dust_position(beta, obs_jd - age, obs_jd, comet_el,
                                v_R0, v_T0, v_N0, gamma, m_exp)
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
            r_D = dust_position(beta, obs_jd - age, obs_jd, comet_el,
                                v_R0, v_T0, v_N0, gamma, m_exp)
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

    ej_active = bool(v_R0 or v_T0 or v_N0)
    info = dict(r_helio=r_helio, r_geo=r_geo, phase_angle=phase_ang,
                RA=RA, Dec=Dec, obs_str=jd_to_str(obs_jd), obs_jd=obs_jd,
                ejection=dict(v_R0=v_R0, v_T0=v_T0, v_N0=v_N0,
                              gamma=gamma, m_exp=m_exp, active=ej_active))

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


def _set_3d_equal_aspect(ax, points_list, fixed_half_range=None, fixed_center=None):
    """Force equal X/Y/Z scaling on a 3D Axes (mplot3d doesn't do this by
    default). Returns half_range (the half-width of the cube it set the
    view to) so callers can scale other elements (e.g. an arrow) to match
    the same fresh extent instead of a separately-computed approximation
    that could drift out of sync with what's actually being shown.

    fixed_half_range/fixed_center (v3.0.7): if given, use this Sun-
    centered fixed-size cube instead of auto-fitting to points_list at
    all — for "zoom to inner solar system" mode, where a comet many AU
    away would otherwise force the WHOLE view (including Earth's 1 AU
    orbit) down to a single indistinguishable point near the Sun."""
    if fixed_half_range is not None:
        half_range = fixed_half_range
        centers = np.array(fixed_center if fixed_center is not None else (0.0, 0.0, 0.0))
    else:
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
    return half_range


def draw_orbit_diagram(ax, diagram: dict, dark: bool = True,
                       zoom_inner_au: float | None = None):
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

    `zoom_inner_au` (v3.0.7): if given, forces a fixed Sun-centered view
    of ±zoom_inner_au AU instead of auto-fitting to the whole comet
    orbit. For comets well beyond ~1 AU, the default auto-fit view makes
    Earth's 1 AU orbit (and the Sun/Earth separation) collapse down to a
    barely-distinguishable point near the center — this trades "see the
    whole orbit" for "see the inner solar system clearly", with the
    comet itself simply off-screen if it's currently farther out than
    zoom_inner_au.
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

    # Set the view extent BEFORE sizing the motion arrow below, so the
    # arrow scales off the same fresh half_range actually being displayed
    # rather than a separately-precomputed r_max_plot that could drift
    # out of sync with it.
    half_range = _set_3d_equal_aspect(
        ax, [op, eo, np.array([[sx, sy, sz]]),
            np.array([[px, py, pz]]), np.array([[ex, ey, ez]]),
            np.array([[cx, cy, cz]])],
        fixed_half_range=zoom_inner_au,
        fixed_center=(0.0, 0.0, 0.0) if zoom_inner_au is not None else None)

    # Direction-of-motion arrow (v3.0, fixed in v3.0.5) — heliocentric
    # velocity unit vector from compute_orbit_diagram(), scaled to a
    # fraction of half_range (the actual current view's half-width, just
    # computed above) rather than the diagram's separately-precomputed
    # r_max_plot. NOTE: this keeps the arrow correctly proportioned every
    # time this function redraws (e.g. on a date change) — it does NOT
    # dynamically rescale during pure interactive zoom via the plot
    # toolbar without a redraw, since matplotlib's 3D mouse/toolbar zoom
    # doesn't call back into this function at all.
    vhx, vhy, vhz = diagram.get('velocity_dir_xyz', (0.0, 0.0, 0.0))
    if vhx or vhy or vhz:
        arrow_len = 0.36 * half_range
        ax.quiver(cx, cy, cz, vhx*arrow_len, vhy*arrow_len, vhz*arrow_len,
                  color='#00e0a0', linewidth=2.0, arrow_length_ratio=0.3,
                  zorder=7, label=f"{diagram['name']} direction of motion")

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


def _beta_to_radius_um(beta: float, rho_g_cm3: float = 0.5, Qpr: float = 1.0) -> float:
    """
    Numeric companion to beta_to_size() (which returns a formatted
    *string* for display, e.g. "57.4 µm") — this returns a plain float
    in µm, for use in calculations like the Monte Carlo size sampler
    below. Shares the same BETA_TO_RADIUS_UM constant as beta_to_size()
    so the two can never give a different radius for the same β.
    """
    if beta <= 0:
        return np.inf
    return BETA_TO_RADIUS_UM * Qpr / (beta * rho_g_cm3)


def _radius_um_to_beta(a_um, rho_g_cm3: float = 0.5, Qpr: float = 1.0):
    """Inverse of _beta_to_radius_um(). Accepts scalar or array input."""
    a_um = np.asarray(a_um, dtype=float)
    with np.errstate(divide='ignore'):
        beta = BETA_TO_RADIUS_UM * Qpr / (a_um * rho_g_cm3)
    return beta


def real_ejection_speed_ms(v0_coeff: float, beta, r_h_au: float,
                            gamma: float, m_exp: float):
    """
    The ACTUAL ejection speed (m/s) for a grain of the given β at
    heliocentric distance r_h_au — i.e. v0_coeff itself is only the
    speed of a (mostly hypothetical) β=1 grain; this is what the
    compute_morphology_mc() particles are really moving at:

        |v_ej| = v0_coeff · β^gamma · r_h_au^(−m_exp)

    Same formula as the one applied per-particle inside
    compute_morphology_mc() (see its v0_coeff/gamma/m_exp docstring) —
    kept as one shared function so the GUI's report/overlay display
    can never silently drift from what was actually simulated. Accepts
    scalar or array beta.
    """
    beta = np.asarray(beta, dtype=float)
    return v0_coeff * beta**gamma * r_h_au**(-m_exp)


# ─────────────────────────────────────────────────────────────────────────────
#  MONTE CARLO DUST MORPHOLOGY  (v3.1, Phase 1)
#
#  Where compute_model() traces ONE syndyne/synchrone curve per (β, age)
#  choice — useful for visual diagnosis, but silent on how much dust is
#  actually at each point — this samples a whole POPULATION of grains
#  from a size distribution and a release-time window, propagates each
#  one with the exact same trajectory physics (dust_position_isotropic()
#  shares _propagate_dust() with dust_position() — see those functions;
#  the two models cannot diverge in their orbit/radiation-pressure
#  physics, only in how the ejection-velocity vector is chosen), and
#  bins the results into a 2D density map: a genuine (relative)
#  brightness prediction rather than a line.
#
#  Phase 1 scope — intentionally simplified, see compute_morphology_mc()
#  docstring for the full list. No PSF convolution, no isophote-contour
#  matching against a real image yet (that's Phase 2); production rate
#  is uniform-in-time; ejection direction is fully isotropic (no day/
#  night-side restriction, no jets); 2-body Kepler only, consistent with
#  the rest of CTA (no N-body planetary perturbation as in e.g. Hui et
#  al. 2020's MERCURY6 integration — immaterial at amateur-imaging dust
#  ages of weeks to a few months).
# ─────────────────────────────────────────────────────────────────────────────

def _sample_power_law(low, high, index, n: int, rng) -> np.ndarray:
    """
    Inverse-transform sample n values from dN/dx ∝ x^index over [low, high]
    (0 < low < high). index is typically negative for a grain-size
    distribution (more small grains than large, e.g. −3.6 per Fulle 2004/
    Guzik et al. 2019 as cited in Hui et al. 2020) — but this helper is
    generic and doesn't assume a sign.

    low, high, index may each be a scalar OR a length-n array (e.g. per-
    particle values interpolated from a time-varying table — see
    compute_morphology_mc()'s size_dist_table parameter) — the formula
    below is elementwise throughout, so array inputs broadcast exactly
    like the scalar case did, no special-casing needed at the call site.

    Standard inverse-CDF result for a power law:
        x = [low^(p+1) + u·(high^(p+1) − low^(p+1))]^(1/(p+1)),  p = index
    handled separately at p = −1 (log-uniform), where the general
    formula divides by zero — done here with np.where (elementwise),
    NOT a scalar `if`, since index can be an array where only SOME
    entries happen to sit near −1.
    """
    u  = rng.uniform(0.0, 1.0, n)
    low, high, index = np.asarray(low, dtype=float), np.asarray(high, dtype=float), np.asarray(index, dtype=float)
    p1 = index + 1.0
    near_minus_one = np.abs(p1) < 1e-9
    p1_safe = np.where(near_minus_one, 1.0, p1)   # dummy value where unused, just to avoid 0**x warnings
    general = (low ** p1_safe + u * (high ** p1_safe - low ** p1_safe)) ** (1.0 / p1_safe)
    log_uniform = low * (high / low) ** u
    return np.where(near_minus_one, log_uniform, general)


def _sample_isotropic_directions(n: int, rng) -> np.ndarray:
    """
    n random unit vectors, isotropic over the full 4π sphere.
    Standard method: cos(θ) uniform in [−1,1] (NOT θ itself uniform —
    that wrongly clusters samples at the poles), φ uniform in [0,2π).

    Expressed directly in heliocentric ecliptic XYZ — isotropic-on-a-
    sphere is rotation-invariant, so sampling in this single fixed frame
    gives a statistically identical result to sampling in each
    particle's own local R/T/N frame, without needing to compute a
    per-particle frame at all.

    Returns an (n,3) array.
    """
    cos_t = rng.uniform(-1.0, 1.0, n)
    phi   = rng.uniform(0.0, 2.0 * np.pi, n)
    sin_t = np.sqrt(1.0 - cos_t ** 2)
    return np.column_stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])


def _orthonormal_basis_from_axis(axis: np.ndarray):
    """Return two unit vectors perpendicular to `axis`.

    Used by the optional CTA-only sunward-cone sampler. The choice of
    basis is arbitrary but deterministic for a given axis; the random
    azimuth below removes any preferred direction around the cone.
    """
    w = vhat(axis)
    if abs(w[0]) < 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    else:
        tmp = np.array([0.0, 1.0, 0.0])
    u = vhat(np.cross(w, tmp))
    v = np.cross(w, u)
    return u, v


def _sample_sunward_direction(r_C: np.ndarray, rng, max_tries: int = 200,
                              cone_half_angle_deg: float = 90.0,
                              *, r_comet_obs: np.ndarray | None = None,
                              r_earth_obs: np.ndarray | None = None,
                              require_projected_sunward: bool = False):
    """
    Sample a sunward-restricted ejection direction.

    r_C is the heliocentric comet position at the reference epoch:
        r_C = Sun -> comet
    The physical comet -> Sun direction is therefore -r_C.

    Default, cone_half_angle_deg=90 and require_projected_sunward=False,
    is the standard physical 3-D sunward hemisphere. This is the closest
    COMTAILS-like option and the recommended default.

    CTA diagnostics:
      • cone_half_angle_deg < 90 narrows the cone around the subsolar axis.
      • require_projected_sunward=True also requires the initial ejection
        direction to project toward the apparent Sun direction in the
        current image. This is an image-plane morphology diagnostic, not a
        new physical force and not the default physical model.
    """
    sun_dir = vhat(-np.asarray(r_C, dtype=float))
    if vmag(sun_dir) < 1e-15:
        return None, None

    cone_half_angle_deg = float(np.clip(float(cone_half_angle_deg), 0.0, 90.0))
    cos_cone_min = float(np.cos(np.radians(cone_half_angle_deg)))

    sun_axis_2d = None
    if require_projected_sunward and r_comet_obs is not None and r_earth_obs is not None:
        sun_proj = project_sky(np.zeros(3), r_comet_obs, r_earth_obs)
        if sun_proj is not None:
            sun_axis_2d = np.array([sun_proj[0], sun_proj[1]], dtype=float)
            m = float(np.hypot(sun_axis_2d[0], sun_axis_2d[1]))
            if m > 1e-15:
                sun_axis_2d /= m
            else:
                sun_axis_2d = None

    for _ in range(max_tries):
        cos_t = rng.uniform(-1.0, 1.0)
        phi   = rng.uniform(0.0, 2.0 * np.pi)
        sin_t = np.sqrt(max(0.0, 1.0 - cos_t ** 2))
        d = np.array([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t], dtype=float)

        cosz = float(np.dot(d, sun_dir))
        if cosz < cos_cone_min:
            continue

        if sun_axis_2d is not None:
            d_proj = project_vector_sky(d, r_comet_obs, r_earth_obs)
            if d_proj is None:
                continue
            d_axis_2d = np.array([d_proj[0], d_proj[1]], dtype=float)
            dm = float(np.hypot(d_axis_2d[0], d_axis_2d[1]))
            if dm <= 1e-15:
                continue
            d_axis_2d /= dm
            if float(np.dot(d_axis_2d, sun_axis_2d)) < 0.0:
                continue

        return d, cosz

    return None, None


def compute_morphology_mc(comet_el: dict, obs_jd: float,
                          beta_range: tuple, gamma_size: float,
                          max_age: float, n_particles: int,
                          v0_coeff: float = 0.0, gamma: float = 0.0, m_exp: float = 0.0,
                          grid_half_width_au: float | None = None, grid_npix: int = 200,
                          rho_g_cm3: float = 0.5, Qpr: float = 1.0, p_v: float = 0.04,
                          north_pa_deg: float = 0.0, rotation_offset_deg: float = 0.0,
                          qt_weights: dict | None = None,
                          active_area: dict | None = None,
                          sunward: bool = False, sunward_expocos: float = 1.0,
                          sunward_reference: str = "emission",
                          sunward_cone_half_angle_deg: float = 90.0,
                          require_projected_sunward: bool = False,
                          phase_law: str = "schleicher",
                          phase_linear_beta: float = 0.024,
                          phase_linear_m_oe: float = 0.28,
                          phase_linear_w_oe: float = 1.5,
                          size_dist_table: dict | None = None,
                          seed: int | None = None,
                          progress_callback=None,
                          progress_update_fraction: float = 0.005) -> dict:
    """
    Phase 1 (v3.1) Monte Carlo dust morphology model.

    IMPORTANT — beta_range is NOT the same kind of input as compute_
    model()'s beta_values. beta_values there is a handful of discrete
    points YOU choose to look at individually ("show me what β=0.01
    looks like"). beta_range here is a CONTINUOUS interval that, together
    with gamma_size, defines a whole population's size distribution — a
    genuinely different kind of input. A reasonable starting point is
    (min(beta_values), max(beta_values)) from whatever syndyne curves
    are already configured, but that's a convenience default, not an
    equivalence.

    Phase 1 simplifications (deliberate; see the module-level comment
    above and the v3.1 changelog for more):
      - Production rate Q(t) is UNIFORM over [obs_jd−max_age, obs_jd] BY
        DEFAULT. Real comets brighten toward perihelion — this default is
        a first-cut choice, not a measured production curve. Pass
        qt_weights (see below) to replace it with an estimate derived
        from a real light curve + Afρ anchor — see
        estimate_qt_from_lightcurve().
      - Ejection direction is fully ISOTROPIC (4π sr) by default — no
        sunward-hemisphere-only restriction, no jets/anisotropy. Pass
        active_area (see below) to switch to a rotating-nucleus active-
        area model instead.
      - gamma_size is a single fixed exponent (no broken power law).
      - 2-body Kepler propagation only (see module comment above).
      - No PSF convolution, no isophote-contour extraction (Phase 2).

    Parameters
    ----------
    beta_range : (beta_min, beta_max) — population size range. See the
        IMPORTANT note above before reusing compute_model()'s beta_values.
    gamma_size : differential size-distribution power-law index,
        dN/da ∝ a^gamma_size (a = grain radius). HEURISTIC DEFAULT for
        any specific comet — e.g. −3.6 sits in the middle of Fulle
        (2004)'s "Motion of cometary dust" (Comets II) review range of
        k=−3.0 to −4.1 measured across ORDINARY (non-interstellar)
        solar-system comets via the inverse-tail method — it is NOT a
        2I/Borisov-specific number (2I's own value, −3.7±1.8 per Guzik
        et al. 2020, sits comfortably inside that same ordinary-comet
        range, which is itself why Hui et al. 2020 cited both Fulle
        2004 and Guzik et al. for it). Still a heuristic default for any
        ONE specific object, though — same status as β/ejection-velocity
        γ elsewhere in CTA, expected to be refined per-object.
    max_age : outer bound (days before obs_jd) of the uniform release-
        time window — same meaning as compute_model()'s max_age. Keep
        this consistent with when THIS comet was actually active; the
        uniform-Q(t) assumption is least defensible over a window long
        enough to span a big change in heliocentric distance/activity
        (e.g. a near-parabolic comet's inbound leg).
    n_particles : Monte Carlo sample size. Larger → smoother density map,
        longer runtime (this loops in pure Python over n_particles calls
        to dust_position_isotropic(), same style as compute_model()'s
        per-point loops — not vectorized; a few thousand particles is a
        reasonable interactive-use starting point, see the v3.1
        changelog for the runtime note).
    v0_coeff, gamma, m_exp : ejection-speed power law,
        |v_ej| = v0_coeff · β^gamma · r_H(t_emit)^(−m_exp) — same
        functional form as dust_position(), but sets only the MAGNITUDE
        here; direction is isotropic-random per particle (see above),
        unlike dust_position()'s fixed v_R0/v_T0/v_N0 direction.
        v0_coeff=0.0 (default) ⇒ zero ejection velocity; every particle
        then follows the bare two-body/radiation-pressure trajectory for
        its sampled β, same as a zero-ejection-velocity F-P point.
    grid_half_width_au : half-width of the square output grid [AU].
        None (default) ⇒ auto-sized to 1.2× the 95th-percentile sky-
        plane offset of the particles that produced a valid position.
    grid_npix : pixels per side of the square output grid.
    rho_g_cm3, Qpr : passed to the β↔radius conversion (shared with
        beta_to_size()) used to sample in radius-space then convert to β.
        ALSO sets the per-particle BRIGHTNESS weight (see below).
    p_v : grain geometric albedo (py_COMTAILS's "Particle geometric
        albedo"; same default 0.04 as compute_dust_production_rate()'s
        own p_v, and as py_COMTAILS's own example TAIL_INPUTS.dat). Used
        with the phase angle at this run's obs_jd (computed
        automatically — not a parameter) to populate info[
        'relative_flux_scale'] = p_v · Schleicher_phase_correction(phase
        angle) — see schleicher_phase_correction()'s docstring for what
        this number is for (cross-run/cross-epoch brightness comparison)
        and, importantly, what it is NOT for (it cannot change a single
        run's own contour shape, since every particle in one run shares
        the same observation-epoch phase angle).
    north_pa_deg, rotation_offset_deg : optional (v3.1) — orient the
        output grid to match the main GUI canvas's image-overlay screen
        appearance instead of this module's own default East-right/
        North-up convention. rotation_offset_deg is the SAME "grid
        rotation (match observed tail)" angle as the main panel's
        slider; north_pa_deg is the image's North position angle (same
        meaning as the overlay's north_pa). Both default to 0.0, which
        keeps the original East-right/North-up output unchanged — this
        is an opt-in convenience for side-by-side comparison with a
        real-image overlay, not something compute_model() itself does
        (that one stays in sky-plane AU regardless; only the GUI's
        canvas applies sky_to_pixel() at DRAW time). When either is
        nonzero, info['oriented']=True and plot_morphology_mc() adjusts
        its axes (inverted Y, compass) to match.
    qt_weights : optional dict(t_jd=array, Q_kg_s=array) — e.g. directly
        the return value of estimate_qt_from_lightcurve(). When given,
        release times are sampled weighted by this Q(t) (inverse-
        transform sampling, see _sample_weighted_ages()) instead of
        uniformly over [obs_jd−max_age, obs_jd]. None (default) keeps
        the v3.1 Phase-1 uniform behaviour unchanged.
    active_area : optional dict — switches ejection direction from
        isotropic to a rotating-nucleus active-area model (port of py_
        COMTAILS's iejec_mode==3; see _active_area_direction()'s
        docstring for the source/citation). None (default) keeps the
        isotropic behaviour. Keys (all required if given):
          nuc_inc_deg    — obliquity: spin-axis tilt from the orbit
                           normal [deg].
          nuc_phi_deg    — argument of the subsolar meridian AT
                           PERIHELION [deg] — i.e. where, in body
                           longitude, local noon fell at perihelion.
          period_d       — nucleus rotation period [days].
          per_jd         — perihelion JD (T_jd, or date_to_jd(T) if not
                           directly available — same epoch nuc_phi is
                           referenced to).
          lat_min_deg, lat_max_deg, lon_min_deg, lon_max_deg —
                           active-area bounds in BODY-FIXED latitude/
                           longitude [deg].
          isun           — bool, default True. True = reject any
                           sampled point that isn't currently
                           illuminated (the literature's usual active-
                           area behaviour); False = use the area's
                           geometric direction regardless of
                           illumination (rarely what you want, but
                           available for comparison/debugging).
          expocos        — exponent on cos(solar zenith angle) applied
                           to v_ej for that particle (grazing incidence
                           ejects slower than a head-on subsolar point;
                           1.0 is a common literature default — see
                           Moreno et al. 2012, ApJ 752, 136, Table 2).
        A particle whose ejection time never has an illuminated point in
        the active area (isun=True can fail to find one — see
        _active_area_direction()'s own docstring) is DROPPED, same
        NaN-handling spirit as an invalid position elsewhere in this
        function — info['n_used'] reflects this.
    sunward : bool, default False — port of py_COMTAILS's iejec_mode==2
        ("sunward"; see _sample_sunward_direction()'s docstring for the
        source/citation). When True (and active_area is None — these
        two are mutually exclusive ejection-DIRECTION choices, same as
        COMTAILS's single iejec_mode selector), restricts isotropic
        kicks to the sunward hemisphere only (cosz>=0) instead of the
        full 4π sphere, weighted by cosz**sunward_expocos. Ignored if
        active_area is given (active_area's own isun handles the
        analogous restriction there). Together with active_area=None/
        sunward=False (full 4π isotropic, the default) this gives THREE
        selectable ejection-direction modes total, matching COMTAILS's
        iejec_mode 1/2/3 — see the v3.1 changelog/module notes for why
        isotropic-4π is a genuine alternative modelling choice here, not
        something sunward-restriction "corrects".
    sunward_expocos : exponent on cos(solar zenith angle) applied to
        v_ej for each particle in sunward mode (same role as active_area
        dict's own 'expocos' key for that mode). Only used when
        sunward=True.
    sunward_reference : {'emission','observation'}, default 'emission'.
        'emission' is the COMTAILS-compatible choice: the sunward
        hemisphere is evaluated from the comet position at each
        particle's own t_emit. 'observation' preserves the older CTA
        shortcut that uses the observation-epoch Sun direction for all
        particles; useful only for A/B diagnostics and reproducing older
        .mcin runs.
    sunward_cone_half_angle_deg : default 90.0. 90° = full COMTAILS-style
        sunward hemisphere. Smaller values restrict ejection to a cone
        around the subsolar direction; this is a CTA diagnostic fitting
        aid, not part of the original COMTAILS hemisphere model.
    size_dist_table : optional dict — port of py_COMTAILS's time-varying
        dM/dt/power/r_min/r_max table (its "dmdt_vel_power_rmin_rmax.dat";
        see _active_area_direction()'s docstring for source/citation).
        When given, OVERRIDES beta_range/gamma_size entirely — each
        particle's grain size is sampled from a power-law whose r_min,
        r_max, and slope are interpolated from this table at that
        particle's OWN release time, instead of one fixed distribution
        for the whole population. Keys:
          days_to_per — array of days relative to perihelion (need not
                        be sorted; sorted internally).
          power       — array, same length, size-distribution slope at
                        each days_to_per row (py_COMTAILS's "power").
          r_min_um, r_max_um — arrays, same length, grain-size bounds
                        [µm] at each row.
          per_jd      — perihelion JD this table's days_to_per is
                        referenced to (same role as active_area's own
                        per_jd key).
        Outside the table's days_to_per range, values are CLAMPED to the
        nearest row (no extrapolation) — same behaviour as np.interp's
        default and consistent with COMTAILS's own out-of-bounds
        handling. None (default): use beta_range/gamma_size as a single
        fixed distribution for the whole run, unchanged from pre-v3.1.1
        behaviour.
    seed : optional RNG seed, for reproducible sampling.

    Returns
    -------
    dict(density, xi_edges, eta_edges, info):
        density   — (grid_npix, grid_npix) array, BRIGHTNESS-weighted
                     (each particle contributes ∝ its cross-sectional
                     area a², not just a unit count — see the v3.1
                     changelog: a raw particle-count histogram is NOT
                     the same as a scattered-light brightness map, since
                     a steep size distribution means small grains vastly
                     outnumber large ones by count while contributing
                     far less light per grain). density[row=η,col=ξ]
                     in the default (unoriented) case (matplotlib imshow
                     convention with origin='lower'); if oriented (see
                     above), row/col instead follow the image-overlay's
                     screen convention (origin='upper' semantics).
        xi_edges, eta_edges — bin edges [AU], length grid_npix+1 each
                     (named for the unoriented axes regardless of
                     whether orientation was applied — they're still AU
                     extents either way, just rotated/mirrored).
        info      — sampled-parameter summary + n_used (particles with a
                     finite in-range position; some are dropped, same
                     NaN-handling spirit as compute_model()'s curves) +
                     oriented/north_pa_deg/rotation_offset_deg.
    """
    rng = np.random.default_rng(seed)
    beta_min, beta_max = beta_range
    if not (0 < beta_min <= beta_max):
        raise ValueError("beta_range must be (beta_min, beta_max) with 0 < beta_min <= beta_max")
    if n_particles < 1:
        raise ValueError("n_particles must be >= 1")
    sunward_reference = str(sunward_reference or "emission").strip().lower()
    if sunward_reference not in ("emission", "observation"):
        raise ValueError("sunward_reference must be 'emission' or 'observation'")
    sunward_cone_half_angle_deg = float(sunward_cone_half_angle_deg)
    if not np.isfinite(sunward_cone_half_angle_deg):
        raise ValueError("sunward_cone_half_angle_deg must be finite")
    sunward_cone_half_angle_deg = float(np.clip(sunward_cone_half_angle_deg, 0.0, 90.0))
    require_projected_sunward = bool(require_projected_sunward)

    phase_law = str(phase_law or "schleicher").strip().lower().replace("-", "_")
    if phase_law not in ("schleicher", "linear_exponential", "none"):
        raise ValueError("phase_law must be 'schleicher', 'linear_exponential', or 'none'")

    # Optional progress callback for GUI/CLI front ends. The callback is
    # deliberately generic and Qt-free:
    #     progress_callback(fraction_0_to_1, message)
    # It is called only every ~progress_update_fraction of the particle loop
    # to avoid slowing the Monte Carlo run down just to update the UI.
    import time as _time
    _progress_t0 = _time.perf_counter()
    _progress_last_fraction = -1.0
    try:
        progress_update_fraction = float(progress_update_fraction)
    except Exception:
        progress_update_fraction = 0.005
    if not np.isfinite(progress_update_fraction) or progress_update_fraction <= 0:
        progress_update_fraction = 0.005

    def _fmt_duration(seconds: float) -> str:
        if not np.isfinite(seconds) or seconds < 0:
            return "--:--"
        seconds = int(round(seconds))
        h, rem = divmod(seconds, 3600)
        m, sec = divmod(rem, 60)
        if h:
            return f"{h:d}:{m:02d}:{sec:02d}"
        return f"{m:02d}:{sec:02d}"

    def _progress(done: int, total: int, stage: str) -> None:
        nonlocal _progress_last_fraction
        if progress_callback is None:
            return
        total = max(int(total), 1)
        done = max(0, min(int(done), total))
        frac = done / total
        # Always report 0 and 1; throttle intermediate updates.
        if 0.0 < frac < 1.0 and (frac - _progress_last_fraction) < progress_update_fraction:
            return
        _progress_last_fraction = frac
        elapsed = _time.perf_counter() - _progress_t0
        eta = elapsed * (1.0 - frac) / frac if frac > 0 else float("nan")
        pct = 100.0 * frac
        msg = (f"{stage}: {pct:5.1f}%  "
               f"({done:,}/{total:,})  ·  "
               f"elapsed {_fmt_duration(elapsed)}  ·  "
               f"ETA {_fmt_duration(eta)}")
        try:
            progress_callback(frac, msg)
        except Exception:
            # Progress reporting must never break the physical calculation.
            pass

    _progress(0, n_particles, "Preparing particle population")

    if size_dist_table is not None:
        # v3.1 — port of py_COMTAILS's time-varying dM/dt, power, r_min,
        # r_max table (its "dmdt_vel_power_rmin_rmax.dat" / _interp5();
        # see _active_area_direction()'s docstring for the source/
        # citation). Each particle's size-distribution SHAPE (not just
        # the production RATE, which qt_weights already covers) can
        # change over the release window — e.g. grains get smaller/finer
        # near perihelion in some published models. Linear interpolation
        # between table rows, same as COMTAILS's own _interp5; clamped to
        # the table's own first/last row outside its range (np.interp's
        # default), not extrapolated.
        #
        # Release times must be known FIRST here (size now depends on
        # t_emit), unlike the single/distribution branches below where
        # age sampling order doesn't matter — hence this whole block is
        # ordered ages-before-sizes only in this branch.
        if qt_weights is not None:
            ages = _sample_weighted_ages(
                qt_weights["t_jd"], qt_weights["Q_kg_s"],
                obs_jd, max_age, n_particles, rng)
        else:
            ages = rng.uniform(0.0, max_age, n_particles)
        t_emits = obs_jd - ages

        per_jd = size_dist_table["per_jd"]
        days_to_per = np.asarray(size_dist_table["days_to_per"], dtype=float)
        order = np.argsort(days_to_per)
        days_to_per = days_to_per[order]
        power_tbl  = np.asarray(size_dist_table["power"], dtype=float)[order]
        rmin_tbl   = np.asarray(size_dist_table["r_min_um"], dtype=float)[order]
        rmax_tbl   = np.asarray(size_dist_table["r_max_um"], dtype=float)[order]

        particle_days = t_emits - per_jd
        power_p = np.interp(particle_days, days_to_per, power_tbl)
        rmin_p  = np.interp(particle_days, days_to_per, rmin_tbl)
        rmax_p  = np.interp(particle_days, days_to_per, rmax_tbl)

        # _sample_power_law()'s inverse-CDF formula is elementwise, so
        # per-particle low/high/index arrays work exactly like the
        # scalar case did — no per-particle Python loop needed.
        a_um  = _sample_power_law(rmin_p, rmax_p, power_p, n_particles, rng)
        betas = _radius_um_to_beta(a_um, rho_g_cm3, Qpr)

    elif beta_min == beta_max:
        # v3.1 — degenerate single-grain-size case (e.g. the GUI's "Single
        # grain size" mode: one β value matched by eye to the visible
        # tail, the same number you'd read off an F-P syndyne label,
        # rather than a population range). No distribution to sample —
        # every particle gets exactly this β; the only remaining
        # randomness is release time (age) and, if v0_coeff != 0, the
        # isotropic kick direction. gamma_size is unused in this case.
        a_um  = np.full(n_particles, _beta_to_radius_um(beta_min, rho_g_cm3, Qpr))
        betas = np.full(n_particles, beta_min)
        ages = (_sample_weighted_ages(qt_weights["t_jd"], qt_weights["Q_kg_s"],
                                      obs_jd, max_age, n_particles, rng)
               if qt_weights is not None else rng.uniform(0.0, max_age, n_particles))
        t_emits = obs_jd - ages
    else:
        # Sample in RADIUS space (dN/da ∝ a^gamma_size), then convert to β —
        # sampling directly in β-space would silently impose a DIFFERENT
        # (and undocumented) distribution, since β ∝ 1/a is a nonlinear
        # transform of a. Smaller β ↔ larger grain, so beta_min -> a_max.
        a_max_um = _beta_to_radius_um(beta_min, rho_g_cm3, Qpr)
        a_min_um = _beta_to_radius_um(beta_max, rho_g_cm3, Qpr)
        a_um     = _sample_power_law(a_min_um, a_max_um, gamma_size, n_particles, rng)
        betas    = _radius_um_to_beta(a_um, rho_g_cm3, Qpr)
        ages = (_sample_weighted_ages(qt_weights["t_jd"], qt_weights["Q_kg_s"],
                                      obs_jd, max_age, n_particles, rng)
               if qt_weights is not None else rng.uniform(0.0, max_age, n_particles))
        t_emits = obs_jd - ages

    r_C, _ = elem_to_state(comet_el, obs_jd)
    r_E    = earth_pos(obs_jd)

    # v3.1 — phase angle + Schleicher correction at THIS run's single
    # observation epoch. Same one number for the whole frame (see
    # schleicher_phase_correction()'s docstring for why) — stored in
    # info for cross-run/cross-epoch brightness comparison, not applied
    # per-particle (it can't change this run's own contour shape at all).
    r_helio_mc = float(vmag(r_C))
    r_geo_mc   = float(vmag(r_C - r_E))
    r_ES_mc    = float(vmag(r_E))
    cos_phase_mc = (r_helio_mc**2 + r_geo_mc**2 - r_ES_mc**2) / (2 * r_helio_mc * r_geo_mc)
    phase_angle_mc = float(np.degrees(np.arccos(np.clip(cos_phase_mc, -1, 1))))
    schleicher_corr_mc = schleicher_phase_correction(phase_angle_mc)
    if phase_law == "linear_exponential":
        phase_corr_mc = linear_exponential_phase_correction(
            phase_angle_mc, phase_linear_beta, phase_linear_m_oe, phase_linear_w_oe)
        phase_law_label = "linear-exponential"
    elif phase_law == "none":
        phase_corr_mc = 1.0
        phase_law_label = "none"
    else:
        phase_corr_mc = schleicher_corr_mc
        phase_law_label = "Schleicher composite"

    _progress(0, n_particles, "Propagating particles")

    xi  = np.full(n_particles, np.nan)
    eta = np.full(n_particles, np.nan)

    if active_area is None and not sunward:
        # Isotropic ejection directions (no effect at all if v0_coeff == 0).
        directions = _sample_isotropic_directions(n_particles, rng)
        for k in range(n_particles):
            r_D = dust_position_isotropic(betas[k], t_emits[k], obs_jd, comet_el,
                                          v0_coeff, gamma, m_exp, directions[k])
            if r_D is None:
                continue
            proj = project_sky(r_D, r_C, r_E)
            if proj is None:
                continue
            xi[k], eta[k] = proj[0], proj[1]
            _progress(k + 1, n_particles, "Propagating isotropic particles")
    elif active_area is None and sunward:
        # Sunward-restricted ejection (v3.1) — see
        # _sample_sunward_direction()'s docstring for the source/
        # citation.
        #
        # BUG FIX (v3.1.1): "sunward" must be evaluated relative to the
        # comet's position AT EACH PARTICLE'S OWN EMISSION TIME (t_emits[k]),
        # not at the observation time (obs_jd). Over a release window of
        # several tens of days, the comet moves along its orbit and the
        # Sun-comet line rotates with it — a particle emitted 40 days ago
        # was pushed "sunward" relative to where the Sun was 40 days ago,
        # not relative to today's Sun direction. Using r_C (obs_jd) for
        # every particle applies a single, slightly wrong sunward
        # reference frame to the whole population, which is exactly the
        # kind of systematic offset that would show up as a modeled
        # sunward feature pointing in a subtly different direction than
        # what COMTAILS (which correctly re-evaluates geometry at each
        # emission epoch) produces. Matches the same per-particle
        # elem_to_state() pattern the active_area branch below already
        # uses for the same reason — this was the one branch still using
        # the observation-time shortcut.
        for k in range(n_particles):
            try:
                if sunward_reference == "emission":
                    r_ref, _ = elem_to_state(comet_el, t_emits[k])
                else:   # legacy CTA / diagnostic mode
                    r_ref = r_C
            except Exception:
                continue
            d_xyz, cosz = _sample_sunward_direction(
                r_ref, rng, cone_half_angle_deg=sunward_cone_half_angle_deg,
                r_comet_obs=r_C, r_earth_obs=r_E,
                require_projected_sunward=require_projected_sunward)
            if d_xyz is None:
                continue
            v0_eff = v0_coeff * (cosz ** sunward_expocos) if sunward_expocos != 0 else v0_coeff
            r_D = dust_position_isotropic(betas[k], t_emits[k], obs_jd, comet_el,
                                          v0_eff, gamma, m_exp, d_xyz)
            if r_D is None:
                continue
            proj = project_sky(r_D, r_C, r_E)
            if proj is None:
                continue
            xi[k], eta[k] = proj[0], proj[1]
            _progress(k + 1, n_particles, "Propagating sunward particles")
    else:
        # Rotating-nucleus active-area ejection (v3.1) — see
        # _active_area_direction()'s docstring for the source/citation
        # and exactly what each active_area key means. Direction (and
        # therefore the illumination-based velocity weighting) depends
        # on t_emits[k] — unlike the isotropic case, it can't be pre-
        # sampled once for the whole population up front.
        nuc_inc = np.radians(active_area['nuc_inc_deg'])
        nuc_phi = np.radians(active_area['nuc_phi_deg'])
        period  = active_area['period_d']
        per_jd  = active_area['per_jd']
        lat_min = np.radians(active_area['lat_min_deg'])
        lat_max = np.radians(active_area['lat_max_deg'])
        lon_min = np.radians(active_area['lon_min_deg'])
        lon_max = np.radians(active_area['lon_max_deg'])
        isun    = active_area.get('isun', True)
        expocos = active_area.get('expocos', 1.0)

        for k in range(n_particles):
            d_xyz, cosz = _active_area_direction(
                comet_el, t_emits[k], per_jd, nuc_inc, nuc_phi, period,
                lat_min, lat_max, lon_min, lon_max, isun, rng)
            if d_xyz is None:
                continue   # never found an illuminated point for this particle
            v0_eff = v0_coeff * (cosz ** expocos) if expocos != 0 else v0_coeff
            r_D = dust_position_isotropic(betas[k], t_emits[k], obs_jd, comet_el,
                                          v0_eff, gamma, m_exp, d_xyz)
            if r_D is None:
                continue
            proj = project_sky(r_D, r_C, r_E)
            if proj is None:
                continue
            xi[k], eta[k] = proj[0], proj[1]
            _progress(k + 1, n_particles, "Propagating active-area particles")

    _progress(n_particles, n_particles, "Building density map")
    finite = np.isfinite(xi) & np.isfinite(eta)
    n_used = int(finite.sum())
    if n_used == 0:
        raise RuntimeError(
            "No particle produced a valid sky-plane position — "
            "check comet_el/obs_jd/beta_range for this comet/epoch.")

    xi_f, eta_f = xi[finite], eta[finite]

    # Brightness weighting (v3.1): each particle scatters light ∝ its
    # cross-sectional area a² — NOT one unit per particle. Without this,
    # a physically-standard steep NEGATIVE gamma_size (many small grains
    # by count) makes the rare-but-much-brighter large grains vanish
    # from the map entirely, which is wrong, and tempts "fixing" it by
    # making gamma_size unphysically shallow/positive instead — fixing
    # the actual bug (here) is the right move, not retuning the slope.
    # The normalization constant in a² doesn't matter for a RELATIVE
    # density map, so a_um² is used directly.
    weights = a_um[finite] ** 2

    # Orientation (v3.1, opt-in — see docstring): match the GUI's real-
    # image overlay screen convention instead of this module's own
    # default East-right/North-up. Order matters and mirrors MainWindow
    # exactly: grid-rotation first (sky-plane), THEN north_pa (screen).
    oriented = (north_pa_deg != 0.0) or (rotation_offset_deg != 0.0)
    if oriented:
        xi_f, eta_f = _rotate_sky(xi_f, eta_f, rotation_offset_deg)
        xi_f, eta_f = _orient_to_screen(xi_f, eta_f, north_pa_deg)

    if grid_half_width_au is None:
        r_sky = np.hypot(xi_f, eta_f)
        grid_half_width_au = float(1.2 * np.percentile(r_sky, 95))
        if not np.isfinite(grid_half_width_au) or grid_half_width_au <= 0:
            grid_half_width_au = 0.5

    edges = np.linspace(-grid_half_width_au, grid_half_width_au, grid_npix + 1)
    density2d, xi_edges, eta_edges = np.histogram2d(
        xi_f, eta_f, bins=[edges, edges], weights=weights)
    # histogram2d returns density2d[xi_bin, eta_bin]; transpose so
    # density[row, col] = density[η_bin, ξ_bin], matching imshow's
    # (row, col) = (y, x) convention with origin='lower'. (When oriented,
    # xi_f/eta_f are screen dx/dy in sky_to_pixel()'s sense — the row/col
    # meaning is then "down the image"/"across the image" instead, which
    # plot_morphology_mc() accounts for via invert_yaxis().)
    density = density2d.T

    info = dict(
        beta_range=(beta_min, beta_max), gamma_size=gamma_size,
        max_age=max_age, n_particles=n_particles, n_used=n_used,
        v0_coeff=v0_coeff, gamma=gamma, m_exp=m_exp,
        grid_half_width_au=grid_half_width_au, grid_npix=grid_npix,
        oriented=oriented, north_pa_deg=north_pa_deg,
        rotation_offset_deg=rotation_offset_deg,
        qt_weighted=(qt_weights is not None),
        active_area_used=(active_area is not None),
        size_dist_table_used=(size_dist_table is not None),
        sunward_used=(active_area is None and sunward),
        sunward_reference=sunward_reference,
        sunward_cone_half_angle_deg=sunward_cone_half_angle_deg,
        require_projected_sunward=require_projected_sunward,
        ejection_mode=("active_area" if active_area is not None
                       else "sunward" if sunward else "isotropic"),
        rho_g_cm3=rho_g_cm3, p_v=p_v,
        r_helio_au=r_helio_mc,
        r_geo_au=r_geo_mc,
        phase_angle_deg=phase_angle_mc,
        schleicher_correction=schleicher_corr_mc,
        phase_law=phase_law,
        phase_law_label=phase_law_label,
        phase_correction=phase_corr_mc,
        phase_linear_beta=phase_linear_beta,
        phase_linear_m_oe=phase_linear_m_oe,
        phase_linear_w_oe=phase_linear_w_oe,
        relative_flux_scale=p_v * phase_corr_mc,
        runtime_seconds=float(_time.perf_counter() - _progress_t0),
    )
    _progress(n_particles, n_particles, "Done")
    return dict(density=density, xi_edges=xi_edges, eta_edges=eta_edges, info=info)


def plot_morphology_mc(mc_result: dict, ax=None, cmap: str = 'magma',
                       overlay_model: dict | None = None, log_scale: bool = True):
    """
    Render a compute_morphology_mc() result as a 2D brightness-weighted
    density image.

    Two display modes, both driven by mc_result['info']['oriented']
    (set by compute_morphology_mc()'s north_pa_deg/rotation_offset_deg —
    see that function's docstring):
      - Default (oriented=False): this module's own East-right/North-up
        convention, axis labels say so directly.
      - Oriented (oriented=True): matches the GUI's real-image overlay
        screen convention (origin='upper' semantics — handled here via
        ax.invert_yaxis()) so this plot can be visually compared
        side-by-side against the main canvas. A small N/E compass is
        drawn instead of literal axis labels, since the axes no longer
        align with pure East/North once rotated.

    overlay_model : optional compute_model() result dict — if given,
        its syndyne/synchrone curves AND its Sun-direction/anti-velocity
        arrows are drawn on top of the density map (white curves; yellow
        ☀ / magenta −v arrows, same colors and fixed-length-normalized
        convention as the main GUI canvas) — letting you directly
        compare "where the F-P test curves/reference vectors point"
        against "where the simulated mass actually piles up". Comparing
        the two single arrows' angles between this window and the main
        canvas is a much more precise sanity check than eyeballing two
        diffuse fans/bands against each other. This is the same kind of
        comparison Hui et al. (2020) Fig. 8 makes (modelled contour over
        real image) but here it's the simplified F-P curves over OUR OWN
        simulated density, before any real-image/PSF/isophote step
        (Phase 2).
        IMPORTANT if oriented=True: overlay_model's ξ/η are assumed to
        already include any grid-rotation offset (true for MainWindow's
        self._model, which _on_model_ready() rotates in place before
        storing) — only the north_pa part is (re-)applied here, to avoid
        double-rotating it. If you pass a RAW compute_model() result
        that was never grid-rotated, set mc_result's rotation_offset_deg
        to 0 and apply that rotation to the overlay yourself first.
    log_scale : if True (default), plot log10(1+density) — even with
        brightness weighting (see compute_morphology_mc()), the central
        coma is usually orders of magnitude brighter than the faint tail,
        which would wash out on a linear scale.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    density = mc_result['density']
    xi_edges, eta_edges = mc_result['xi_edges'], mc_result['eta_edges']
    info = mc_result['info']
    oriented = info.get('oriented', False)
    img = np.log10(1.0 + density) if log_scale else density

    ax.imshow(img, origin='lower', cmap=cmap, aspect='equal',
             extent=[xi_edges[0], xi_edges[-1], eta_edges[0], eta_edges[-1]])
    # Pin the view to the density grid BEFORE adding overlay curves — those
    # curves often extend much further (a low-β syndyne can run for AU),
    # and matplotlib's autoscale would otherwise stretch the whole view to
    # fit them, shrinking the actual density map to a speck in one corner.
    ax.set_xlim(xi_edges[0], xi_edges[-1])
    ax.set_ylim(eta_edges[0], eta_edges[-1])

    if overlay_model is not None:
        npa = info.get('north_pa_deg', 0.0) if oriented else 0.0
        for s in overlay_model.get('syndynes', []):
            ox, oy = (_orient_to_screen(s['xi'], s['eta'], npa) if oriented
                      else (s['xi'], s['eta']))
            ax.plot(ox, oy, color='white', lw=0.8, alpha=0.7)
        for s in overlay_model.get('synchrones', []):
            ox, oy = (_orient_to_screen(s['xi'], s['eta'], npa) if oriented
                      else (s['xi'], s['eta']))
            ax.plot(ox, oy, color='white', lw=0.8, alpha=0.7, ls='-.')

        # Sun-direction / anti-velocity arrows (v3.1) — same fixed-length-
        # normalized convention as the main canvas (CometTailGUI.py's
        # PlotCanvas: unit vector × 13% of view width), same colors, so a
        # precise angle comparison between the two windows is a matter of
        # comparing two single arrows rather than eyeballing a diffuse
        # density band against a diffuse syndyne fan — much less error-
        # prone. Both vectors already carry the grid-rotation offset
        # baked in by MainWindow._on_model_ready() before this model was
        # stored as self._model; only the north_pa part is applied here.
        half_w = 0.5 * (xi_edges[-1] - xi_edges[0])
        arrow_len = 0.26 * half_w   # ~13% of the FULL width either side
        sun_xi, sun_eta = overlay_model.get('sun_dir', (0.0, 0.0))
        slen = np.hypot(sun_xi, sun_eta)
        if slen > 1e-10:
            sx, sy = (_orient_to_screen(np.array([sun_xi/slen]), np.array([sun_eta/slen]), npa)
                      if oriented else (sun_xi/slen, sun_eta/slen))
            sx, sy = np.ravel(sx)[0], np.ravel(sy)[0]
            ax.annotate("", xy=(sx*arrow_len, sy*arrow_len), xytext=(0, 0),
                       arrowprops=dict(arrowstyle="->", color="#ffe030", lw=1.8), zorder=6)
            ax.text(sx*arrow_len*1.15, sy*arrow_len*1.15, "☀",
                   color="#ffe030", fontsize=11, ha='center', zorder=7)
        avx, avy = overlay_model.get('antivel_dir', (0.0, 0.0))
        avlen = np.hypot(avx, avy)
        if avlen > 1e-10:
            avx_n, avy_n = (_orient_to_screen(np.array([avx/avlen]), np.array([avy/avlen]), npa)
                           if oriented else (avx/avlen, avy/avlen))
            avx_n, avy_n = np.ravel(avx_n)[0], np.ravel(avy_n)[0]
            ax.annotate("", xy=(avx_n*arrow_len, avy_n*arrow_len), xytext=(0, 0),
                       arrowprops=dict(arrowstyle="->", color="#ff3399", lw=1.8), zorder=6)
            ax.text(avx_n*arrow_len*1.15, avy_n*arrow_len*1.15, "−v",
                   color="#ff3399", fontsize=9, ha='center', zorder=7)

        ax.set_xlim(xi_edges[0], xi_edges[-1])
        ax.set_ylim(eta_edges[0], eta_edges[-1])

    if oriented:
        # origin='upper' semantics (matches the GUI's image canvas) —
        # increasing y here means DOWN the image, so flip the displayed
        # axis rather than the data (imshow/plot calls above are
        # unaffected; only the on-screen direction flips).
        ax.invert_yaxis()
        # Small N/E compass, replicated from the main canvas's own
        # (CometTailGUI.py's PlotCanvas) so the two look the same way —
        # both ultimately derive from north_pa via the same rotation.
        npa_r  = np.radians(info.get('north_pa_deg', 0.0))
        cx, cy = xi_edges[0] + 0.08 * (xi_edges[-1]-xi_edges[0]), \
                 eta_edges[-1] - 0.08 * (eta_edges[-1]-eta_edges[0])
        aL = 0.05 * (xi_edges[-1] - xi_edges[0])
        n_dx, n_dy =  np.sin(npa_r) * aL, -np.cos(npa_r) * aL
        e_dx, e_dy = -np.cos(npa_r) * aL, -np.sin(npa_r) * aL
        ax.annotate("", xy=(cx+n_dx, cy+n_dy), xytext=(cx, cy),
                   arrowprops=dict(arrowstyle="->", color="#60c8ff", lw=1.6))
        ax.text(cx+n_dx*1.3, cy+n_dy*1.3, "N", color="#60c8ff",
               fontsize=9, ha='center', va='center', fontweight='bold')
        ax.annotate("", xy=(cx+e_dx, cy+e_dy), xytext=(cx, cy),
                   arrowprops=dict(arrowstyle="->", color="#60c8ff", lw=1.3))
        ax.text(cx+e_dx*1.3, cy+e_dy*1.3, "E", color="#60c8ff",
               fontsize=9, ha='center', va='center', fontweight='bold')
        ax.set_xlabel('Δx (AU, image-aligned)')
        ax.set_ylabel('Δy (AU, image-aligned)')
    else:
        ax.set_xlabel('Δ East (AU)  →')
        ax.set_ylabel('↑  Δ North (AU)')

    ax.set_title(
        f"Monte Carlo morphology  ·  N={info['n_used']}/{info['n_particles']}"
        f"  ·  β∈[{info['beta_range'][0]:.4g}, {info['beta_range'][1]:.4g}]"
        f"  ·  γ_size={info['gamma_size']:.2g}", fontsize=9)
    return ax


def parse_tycho_photometry_file(filepath: str) -> dict:
    """
    Parse a Tycho Tracker ICQ-format photometry export (.txt) and extract
    everything needed to auto-fill the MC window's Calibrated mag/arcsec²
    mode — the Afρ summary line, and the radius-vs-cumulative-magnitude
    profile table used to derive isophote-level surface brightness.

    Tycho Tracker's export has three tables after a free-text ICQ header:
      1. A fixed-aperture magnitude table (10x10..60x60 px boxes) — not
         used here, since it doesn't carry a plate scale per column.
      2. A one-line Afρ SUMMARY table:
           COMET  UTC  DELTA  r  AP"  MAG  RSR  CM  +/-  LOG-AFRHO  OBS
         giving Δ [AU], r_h [AU], and Afρ [cm] at ONE reference aperture.
      3. A per-radius profile table:
           Radius(px) Radius(arcmin) Radius(KM) Afrho(CM) +/- Afrho(10k)
           +/- Mag
         — this is the one used to build the surface-brightness profile,
         since it gives CUMULATIVE aperture magnitude at many radii.

    Parsing is done by locating each table by its header line (rather
    than fixed line numbers), so it tolerates the free-text ICQ block at
    the top varying in length between exports.

    Returns
    -------
    dict:
        delta_au        : geocentric distance [AU]
        r_h_au          : heliocentric distance [AU]
        afrho_cm        : Afρ at the summary-table's reference aperture [cm]
        afrho_err_cm    : formal uncertainty [cm]
        aperture_arcsec : that reference aperture's radius ["]
        mag_total       : integrated magnitude at that aperture
        radii_arcsec    : np.ndarray, profile table radii [arcsec]
        cumulative_mag  : np.ndarray, cumulative aperture magnitude at
                          each radius (same length as radii_arcsec)
        obs_date        : UTC date string as written in the file, or None
        mpc_code        : observatory MPC code, or None

    Raises
    ------
    ValueError if the per-radius profile table (the one needed for
    surface-brightness calibration) cannot be found in the file.
    """
    with open(filepath, encoding='utf-8', errors='replace') as f:
        text = f.read()
    lines = text.splitlines()

    result = dict(delta_au=None, r_h_au=None, afrho_cm=None,
                  afrho_err_cm=None, aperture_arcsec=None, mag_total=None,
                  radii_arcsec=None, cumulative_mag=None,
                  obs_date=None, mpc_code=None)

    # ── MPC code (from "COD XXX" line) ──────────────────────────────────
    for ln in lines:
        m = re.match(r'^\s*COD\s+(\S+)', ln)
        if m:
            result['mpc_code'] = m.group(1)
            break

    # ── Afρ summary line — find the header row, then the next data row ──
    for i, ln in enumerate(lines):
        if 'DELTA' in ln and ('AFRHO' in ln.upper() or 'AP "' in ln or 'AP"' in ln):
            # Data row(s) follow after a "---- ---- ----" separator line
            for j in range(i + 1, min(i + 6, len(lines))):
                dl = lines[j]
                if not dl.strip() or set(dl.strip()) <= set('- '):
                    continue
                # Columns: COMET  UTC  DELTA  r  AP"  MAG  RSR  CM +/- LOGAFRHO  OBS
                parts = dl.split()
                nums = [p for p in parts if re.match(r'^-?\d+\.?\d*$', p)]
                # date/time tokens (DD/MM/YYYY and HH:MM:SS) are numeric-
                # looking but not part of the numeric measurement block —
                # strip tokens containing '/' or ':' first.
                parts_clean = [p for p in parts if '/' not in p and ':' not in p]
                nums = [p for p in parts_clean if re.match(r'^-?\d+\.?\d*$', p)]
                if len(nums) >= 6:
                    try:
                        result['delta_au']        = float(nums[0])
                        result['r_h_au']           = float(nums[1])
                        result['aperture_arcsec']  = float(nums[2])
                        result['mag_total']        = float(nums[3])
                        # nums[4] = RSR (not needed)
                        result['afrho_cm']         = float(nums[5])
                        if len(nums) >= 7:
                            result['afrho_err_cm'] = float(nums[6])
                    except (ValueError, IndexError):
                        pass
                    # obs date: first DD/MM/YYYY-looking token
                    for p in parts:
                        if re.match(r'^\d{2}/\d{2}/\d{4}$', p):
                            result['obs_date'] = p
                            break
                    break
            break

    # ── Per-radius profile table ─────────────────────────────────────────
    header_idx = None
    for i, ln in enumerate(lines):
        if 'Radius(arcmin)' in ln or ('Radius(px)' in ln and 'Mag' in ln):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(
            "Could not find the per-radius profile table (looking for a "
            "header containing 'Radius(arcmin)' and 'Mag'). This file may "
            "not be a Tycho Tracker Afρ-profile export, or its format has "
            "changed.")

    radii, mags = [], []
    for ln in lines[header_idx + 1:]:
        s = ln.strip()
        if not s:
            if radii:
                break   # blank line after data = end of table
            continue
        if set(s) <= set('-+ '):
            continue   # separator line
        parts = s.split()
        try:
            nums = [float(p) for p in parts]
        except ValueError:
            if radii:
                break
            continue
        if len(nums) < 7:
            continue
        # Columns: Radius(px) Radius(arcmin) Radius(KM) Afrho(CM) +/-
        #          Afrho(10k) +/- Mag
        radii.append(nums[1] * 60.0)   # arcmin -> arcsec
        mags.append(nums[-1])          # last column = cumulative Mag

    if len(radii) < 3:
        raise ValueError(
            f"Found the profile table header but only {len(radii)} usable "
            f"data row(s) — need at least 3 to compute a surface-brightness "
            f"profile.")

    result['radii_arcsec']   = np.array(radii)
    result['cumulative_mag'] = np.array(mags)
    return result


def surface_brightness_profile(radii_arcsec, cumulative_mag) -> list:
    """
    Convert a cumulative (aperture) magnitude profile — magnitude of ALL
    light within radius r, increasing r — into differential per-annulus
    surface brightness in mag/arcsec², directly comparable to isophote
    levels (what COMTAILS calls B_max, and what an image's own isophote
    contours are traced at).

    Method (standard photometric identities — not an approximation specific
    to comets):
        1. Cumulative magnitude -> cumulative flux (Pogson's law; flux is
           additive, magnitude is not):
               F(r) = 10^(-0.4 * m(r))
        2. Flux within the annulus between consecutive radii:
               F_ann = F(r2) - F(r1)
        3. Back to a magnitude for just that annulus's light:
               m_ann = -2.5 * log10(F_ann)
        4. Surface brightness = magnitude per unit area:
               SB = m_ann + 2.5 * log10(Area_ann)   [Area in arcsec²]

    Parameters
    ----------
    radii_arcsec   : array-like, aperture radii [arcsec], ascending order
    cumulative_mag : array-like, same length, cumulative magnitude at
                     each radius

    Returns
    -------
    list of (r1_arcsec, r2_arcsec, sb_mag_arcsec2) tuples, one per annulus
    between consecutive input radii. Annuli where the flux difference is
    non-positive (noise-dominated at large r, or non-monotonic input) are
    skipped rather than raising — a profile can still be usable even if a
    few outer points are noisy.
    """
    radii_arcsec   = np.asarray(radii_arcsec, dtype=float)
    cumulative_mag = np.asarray(cumulative_mag, dtype=float)
    flux = 10.0 ** (-0.4 * cumulative_mag)

    profile = []
    for k in range(1, len(radii_arcsec)):
        r1, r2 = radii_arcsec[k - 1], radii_arcsec[k]
        F_ann = flux[k] - flux[k - 1]
        if F_ann <= 0 or r2 <= r1:
            continue
        m_ann = -2.5 * math.log10(F_ann)
        area_ann = math.pi * (r2**2 - r1**2)
        sb = m_ann + 2.5 * math.log10(area_ann)
        profile.append((float(r1), float(r2), float(sb)))
    return profile


def suggest_sb_contour_levels(radii_arcsec, cumulative_mag,
                               margin_mag: float = 0.3) -> dict:
    """
    One-call convenience wrapper: build the surface-brightness profile
    (see surface_brightness_profile()) and return ready-to-use
    brightest/faintest contour levels for MCWindow's Calibrated
    mag/arcsec² mode, with a small safety margin so the suggested
    innermost/outermost rings sit just inside the actually-measured
    range rather than exactly at its noisy edges.

    Returns
    -------
    dict:
        profile          : the full list of (r1, r2, sb) tuples
        sb_brightest_raw : brightest (numerically smallest) SB measured
        sb_faintest_raw  : faintest (numerically largest) SB measured
        suggested_bright : sb_brightest_raw + margin_mag  (slightly fainter,
                            safely inside the measured range)
        suggested_faint  : sb_faintest_raw - margin_mag   (slightly brighter,
                            safely inside the measured range)

    Raises ValueError if surface_brightness_profile() returns no usable
    annuli (e.g. a strictly non-increasing cumulative-magnitude column).
    """
    profile = surface_brightness_profile(radii_arcsec, cumulative_mag)
    if not profile:
        raise ValueError(
            "No usable annuli in this profile — cumulative magnitude must "
            "be non-increasing (flux increasing) with radius for at least "
            "one pair of consecutive rows.")
    sbs = [p[2] for p in profile]
    sb_bright = min(sbs)
    sb_faint  = max(sbs)
    return dict(
        profile=profile,
        sb_brightest_raw=sb_bright,
        sb_faintest_raw=sb_faint,
        suggested_bright=sb_bright + margin_mag,
        suggested_faint=sb_faint - margin_mag,
    )


def calibrate_mc_to_surface_brightness(mc_result: dict,
                                        afrho_cm: float,
                                        delta_au: float,
                                        r_h_au: float,
                                        au_per_px: float,
                                        m_sun: float = -26.74) -> dict:
    """
    Convert the MC density (a²-weighted histogram) to a calibrated surface
    brightness map in mag/arcsec², using an observed Afρ measurement as the
    absolute scale anchor — the same approach COMTAILS uses to produce
    synthetic isophotes directly comparable to observed image isophotes.

    Method — fractional flux distribution:
        The total flux from the model within the grid aperture is set by
        the observed Afρ (A'Hearn et al. 1984 definition):

            F_total / F_sun = Afρ [cm] × ρ_cm / (4 × r_h_AU² × Δ_cm²)

        Each grid bin's flux is then:

            F_bin / F_sun = (density_bin / W_total) × F_total / F_sun

        where density_bin is the a²-weighted particle count (already the
        brightness-correct weight — larger grains scatter more) and W_total
        is the sum over all bins. Finally, surface brightness per bin:

            SB_bin = m_sun − 2.5 × log10(F_bin / F_sun / Ω_bin_arcsec²)

        This approach avoids explicit particle→cross-section unit conversion
        by letting Afρ anchor the absolute scale, exactly as observers do.

    Parameters
    ----------
    mc_result   : dict from compute_morphology_mc()
    afrho_cm    : observed Afρ in cm (e.g. from CTA Afρ calculator)
    delta_au    : geocentric distance at observation [AU]
    r_h_au      : heliocentric distance at observation [AU]
    au_per_px   : image plate scale [AU/pixel] — used for pixel size note only
    m_sun       : apparent V magnitude of the Sun (default −26.74)

    Returns
    -------
    dict:
        sb_map    — 2D array [mag/arcsec²], NaN where density == 0
        xi_edges  — same as mc_result['xi_edges']
        eta_edges — same as mc_result['eta_edges']
        sb_peak   — surface brightness at peak density bin [mag/arcsec²]
        sb_median — median SB of non-zero bins [mag/arcsec²]
    """
    density   = mc_result['density'].astype(float)
    xi_edges  = mc_result['xi_edges']
    eta_edges = mc_result['eta_edges']

    W_total = density.sum()
    if W_total <= 0 or (density > 0).sum() < 5:
        return dict(sb_map=np.full_like(density, np.nan),
                    xi_edges=xi_edges, eta_edges=eta_edges,
                    sb_peak=np.nan, sb_median=np.nan)

    pixel_size_au = xi_edges[1] - xi_edges[0]       # AU per grid bin
    rho_au        = abs(xi_edges[-1])                # grid half-width [AU]
    rho_cm        = rho_au * AU_M * 100.0            # [cm]
    delta_cm      = delta_au * AU_M * 100.0          # [cm]

    # Total comet flux / solar flux within grid aperture (Afρ definition):
    #   F_total / F_sun = Afρ[cm] × ρ[cm] / (4 × r_h[AU]² × Δ[cm]²)
    F_total_ratio = afrho_cm * rho_cm / (4.0 * r_h_au**2 * delta_cm**2)

    # Solid angle of one grid bin in arcsec²
    pix_angle_arcsec = pixel_size_au / delta_au * 206265.0
    pix_solid_arcsec2 = pix_angle_arcsec ** 2

    # Per-bin surface brightness
    sb_map   = np.full_like(density, np.nan)
    nonzero  = density > 0

    F_per_arcsec2 = (density[nonzero] / W_total) * F_total_ratio / pix_solid_arcsec2
    F_per_arcsec2 = np.maximum(F_per_arcsec2, 1e-30)
    sb_map[nonzero] = m_sun - 2.5 * np.log10(F_per_arcsec2)

    finite = sb_map[np.isfinite(sb_map)]
    sb_peak   = float(finite.min()) if finite.size else np.nan   # brightest (smallest mag)
    sb_median = float(np.median(finite)) if finite.size else np.nan

    return dict(sb_map=sb_map, xi_edges=xi_edges, eta_edges=eta_edges,
                sb_peak=sb_peak, sb_median=sb_median)


def extract_contours_at_magnitude_levels(sb_result: dict,
                                          mag_levels: list,
                                          min_pts: int = 15) -> list:
    """
    Extract contour paths at specific surface brightness levels [mag/arcsec²]
    from a calibrated surface brightness map (output of
    calibrate_mc_to_surface_brightness()).

    In mag units brighter = numerically SMALLER, so mag_levels should be
    ordered from faint (large number) to bright (small number), e.g.
    [23.0, 22.0, 21.0, 20.0].  The contours trace where the synthetic
    surface brightness equals each level, exactly like COMTAILS isophotes.

    Parameters
    ----------
    sb_result  : dict from calibrate_mc_to_surface_brightness()
    mag_levels : list of surface brightness levels [mag/arcsec²], faint→bright
    min_pts    : minimum contour segment length (filters noise fragments)

    Returns
    -------
    list of (N,2) arrays of [xi, eta] vertices in AU, same format as
    extract_morphology_contours() — compatible with the existing to_px()
    drawing pipeline.
    """
    sb_map   = sb_result['sb_map']
    xi_edges = sb_result['xi_edges']
    eta_edges= sb_result['eta_edges']

    xi_centers  = 0.5 * (xi_edges[:-1]  + xi_edges[1:])
    eta_centers = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    # In mag space: brighter bins have SMALLER numbers.
    # matplotlib contour expects Z[j,i] with X=xi_centers, Y=eta_centers.
    # We contour on the mag map; levels from faint→bright numerically
    # means levels are decreasing (24→20→16).  Keep the same order
    # matplotlib expects (ascending), but flip: contour at mag=20 means
    # find where sb_map = 20; fine since contour doesn't care about direction.
    levels = sorted(mag_levels)   # ascending (faint→bright numerically)

    import matplotlib.pyplot as plt
    fig_tmp, ax_tmp = plt.subplots()
    try:
        # Replace NaN with a very faint value (won't be above any level)
        sb_clean = np.where(np.isfinite(sb_map), sb_map, 99.0)
        cs = ax_tmp.contour(xi_centers, eta_centers, sb_clean, levels=levels)
        paths = [seg for level_segs in cs.allsegs
                 for seg in level_segs if len(seg) >= min_pts]
    finally:
        plt.close(fig_tmp)
    return paths


def extract_morphology_contours(mc_result: dict, n_levels: int = 3,
                                smooth_sigma: float = 1.5,
                                percentile_floor: float = 80.0,
                                smooth_sigma_au: float | None = None) -> list:
    """
    Extract isophote-style contour paths from a compute_morphology_mc()
    density grid, in the SAME pure sky-plane (East-right/North-up,
    UNROTATED) AU convention as compute_model()'s syndyne/synchrone
    output (mc_result should come from a call that left north_pa_deg/
    rotation_offset_deg at their 0.0 defaults — those parameters exist
    for the standalone-popup display case, not this one).

    v3.1, Phase 2 design choice: rather than re-deriving a second
    sky→screen orientation pipeline (the source of the Phase-1-popup
    orientation-matching headaches), these paths are meant to be fed
    straight into the GUI's EXISTING, already-trusted sky_to_pixel()/
    to_px() conversion — the SAME one syndyne/synchrone curves already
    use when drawn over a real image. One rendering path, not two, so
    there is no longer a second place for an orientation mismatch to
    hide. See CometTailGUI.py's MainWindow.set_mc_contours()/PlotCanvas
    drawing code for where these paths get converted and drawn.

    Mirrors CometTailGUI.py's _compute_isophote_levels() (the analogous
    function for a real loaded image) — Gaussian-smooth first (kills
    bin-to-bin Monte Carlo shot noise while preserving the genuine
    larger-scale tail shape), log-spaced levels (intensity isophotes ≈
    even magnitude spacing, since magnitude ∝ −2.5·log10(intensity)),
    percentile floor so contours don't trace empty/near-empty bins.

    smooth_sigma_au : v3.1 bug fix — smooth_sigma alone is in GRID
        PIXELS, but grid_npix (compute_morphology_mc()'s, default 200) is
        FIXED regardless of how big grid_half_width_au ends up being. A
        wide r_min/r_max range (or small r_min specifically) can make
        grid_half_width_au — and therefore each pixel's physical AU size
        — much larger than a narrower range would. The SAME "1.5 pixel"
        smoothing then blurs over a correspondingly larger ABSOLUTE area,
        which inflates the apparent near-nucleus "coma" size even though
        nothing about the actual near-nucleus particle distribution
        changed — purely a side effect of the grid needing to be bigger
        to capture a longer tail. Pass smooth_sigma_au to smooth by a
        FIXED ABSOLUTE AU radius instead (converted to the equivalent
        pixel sigma using this grid's own pixel scale) — this stays
        consistent across different r_min/r_max choices. Takes priority
        over smooth_sigma when given; leave at None to use the old
        pixel-based behaviour unchanged (e.g. for existing callers).

    Returns a list of (N,2) arrays, each one contour path's [ξ,η]
    vertices in AU (ready for to_px()). Empty list if there isn't
    enough signal to contour meaningfully (e.g. n_particles too low).
    """
    density = mc_result['density']
    xi_edges, eta_edges = mc_result['xi_edges'], mc_result['eta_edges']
    xi_centers  = 0.5 * (xi_edges[:-1] + xi_edges[1:])
    eta_centers = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    if smooth_sigma_au is not None:
        pixel_size_au = xi_edges[1] - xi_edges[0]   # uniform bin width
        smooth_sigma = smooth_sigma_au / max(pixel_size_au, 1e-300)

    # ── BUG FIX (v3.1) — percentile floor from PRE-smooth density ──────
    # Gaussian smoothing leaks density into neighbouring empty bins (sky
    # background). Computing percentile_floor from the POST-smooth map
    # would include those artificially-populated bins, lowering the
    # effective threshold and pushing contours outward from the coma —
    # the more smoothing applied, the worse this drift. Fixing this by
    # computing the floor from the UNSMOOTHED log-density means that
    # smoothing only changes the SHAPE of the contours (removes bin-to-
    # bin noise) without shifting WHERE the threshold sits. This is why
    # small smooth_sigma_au values appeared to bring contours closer to
    # the observed coma — they happened to avoid the threshold drift by
    # limiting the bleed, not because the physical grid scale was right.
    log_img_raw = np.log10(1.0 + density.astype(float))
    positive_raw = log_img_raw[log_img_raw > 0]
    if positive_raw.size < 10:
        return []   # not enough signal to contour meaningfully
    floor = float(np.percentile(positive_raw, percentile_floor))
    vmax_raw = float(log_img_raw.max())
    if vmax_raw <= floor:
        return []

    # NOW apply smoothing — only for improving contour LINE SHAPE,
    # not for setting thresholds.  We compute the levels entirely from
    # raw data (floor and vmax both from log_img_raw) and then trace
    # those levels on the smoothed map so the lines are less jagged.
    #
    # v3.1 BUG FIX — earlier the levels used vmax from the SMOOTHED
    # density, which is always lower than vmax_raw (Gaussian blurring
    # spreads the peak, reducing its height).  With more smoothing,
    # vmax_smooth drops further, compressing all levels downward and
    # causing the outer contour to trace increasingly sparse/faint
    # regions — the more smoothing applied, the more the outer contour
    # expanded, which is the opposite of the intended behaviour.
    # Using vmax_raw keeps levels anchored to actual particle-count
    # thresholds regardless of how much smoothing is applied.
    img = density.astype(float)
    if smooth_sigma > 0:
        try:
            from scipy.ndimage import gaussian_filter
            img = gaussian_filter(img, sigma=smooth_sigma)
        except ImportError:
            pass   # no scipy -> contour the unsmoothed grid rather than fail

    log_img = np.log10(1.0 + img)

    # Levels from RAW floor and RAW vmax — smoothing does not move them.
    levels = np.linspace(floor, vmax_raw, n_levels + 2)[1:-1]
    if len(levels) == 0:
        return []

    import matplotlib.pyplot as plt
    fig_tmp, ax_tmp = plt.subplots()
    try:
        cs = ax_tmp.contour(xi_centers, eta_centers, log_img, levels=levels)
        # Discard segments shorter than min_pts points — short segments
        # (2-4 points) from sparse outer bins appear as isolated dots
        # rather than lines. min_pts=6 keeps real contour lines while
        # eliminating noise fragments without any smoothing.
        min_pts = 15
        paths = [seg for level_segs in cs.allsegs
                 for seg in level_segs if len(seg) >= min_pts]
    finally:
        plt.close(fig_tmp)
    return paths


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
                 rho_d: float = 0.5,
                 ejection: dict | None = None):

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
        # Non-zero ejection velocity (v3.1) — see compute_model()/
        # dust_position() for the formula. None/all-zero ⇒ v3.0 behaviour.
        self.ejection    = ejection

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
        if self.ejection and any(self.ejection.get(k, 0.0) for k in ('v_R0', 'v_T0', 'v_N0')):
            print(f"  Ejection velocity (v3.1): v_R0={self.ejection.get('v_R0',0.0)} "
                  f"v_T0={self.ejection.get('v_T0',0.0)} v_N0={self.ejection.get('v_N0',0.0)} m/s, "
                  f"γ={self.ejection.get('gamma',0.0)}, m={self.ejection.get('m_exp',0.0)}")
        t0 = time.time()
        self.model = compute_model(
            self.comet_el, self.obs_jd,
            self.beta_values, self.sync_ages,
            self.max_age, self.n_pts,
            ejection=self.ejection)
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

        # ── Ejection velocity (v3.1) — only shown when active, so a
        # default (v3.0-equivalent) run never grows an extra block. ────────
        ej = info.get('ejection')
        if ej and ej.get('active'):
            y -= dy * 0.4
            ax.text(0.05, y, 'EJECTION VELOCITY (v3.1)', color='#e08030',
                    fontsize=8, fontweight='bold', **kw)
            y -= dy * 0.6
            for k, v in [
                ('v_R0', f"{ej['v_R0']:.2g} m/s"),
                ('v_T0', f"{ej['v_T0']:.2g} m/s"),
                ('v_N0', f"{ej['v_N0']:.2g} m/s"),
                ('γ, m', f"{ej['gamma']:.2g}, {ej['m_exp']:.2g}"),
            ]:
                ax.text(0.05, y, k, color='#a05020', fontsize=7.5, **kw)
                ax.text(0.95, y, v, color='#e0a060', fontsize=7.5,
                        ha='right', **kw)
                y -= dy * 0.85

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
    p.add_argument('--vR0',  type=float, default=0.0, help='Ejection velocity, radial component [m/s] (v3.1)')
    p.add_argument('--vT0',  type=float, default=0.0, help='Ejection velocity, transverse component [m/s] (v3.1)')
    p.add_argument('--vN0',  type=float, default=0.0, help='Ejection velocity, normal component [m/s] (v3.1)')
    p.add_argument('--ej-gamma', type=float, default=0.0, help='Ejection velocity β-exponent γ (v3.1; 0 = constant velocity)')
    p.add_argument('--ej-m',     type=float, default=0.0, help="Ejection velocity r_H-exponent m, used as r_H^(-m) (v3.1; 0 = no heliocentric scaling; 0.5 = literature default per Whipple 1950/Ishiguro 2008)")
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
        ejection    = dict(v_R0=args.vR0, v_T0=args.vT0, v_N0=args.vN0,
                           gamma=args.ej_gamma, m_exp=args.ej_m),
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
