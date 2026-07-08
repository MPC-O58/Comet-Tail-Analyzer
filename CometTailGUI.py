#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
  CometTailGUI.py  —  Comet Tail Analyzer (CTA)
  Finson–Probstein + Monte Carlo Dust Tail Model
  Version 3.1.1 ·   Teerasak Thaluang (MPC O51/O58)
  Native Desktop Application (PyQt6 + Matplotlib)

  SPDX-License-Identifier: MIT
  © 2024–2026 Teerasak Thaluang. See LICENSE for full terms.

  Attribution:
    Portions of the Monte Carlo module (comet_tail_analyzer.py) are ported
    from py_COMTAILS (Moreno 2025, A&A 695, A263; F. Moreno, R. Morales &
    N. Robles, IAA-CSIC; github.com/FernandoMorenoDanvila/py_COMTAILS,
    MIT License). Specifically:
      - Schleicher dust phase function coefficients (constants.py)
      - Sunward hemisphere ejection direction sampling (iejec_mode == 2)
      - Active-area rotating-nucleus ejection direction (_anisot_dir2 /
        iejec_mode == 3)
    All other components (F-P engine, orbital mechanics, COBS integration,
    GUI, MC core loop, contour extraction) are original to CTA.

  Cite as:
    Thaluang, T. (2026). RNAAS, doi:10.3847/2515-5172/ae6f90
=============================================================================
  Changelog:
    v3.1.1 • BUG FIX: CLEAR ALL MODELS now resets the hidden MC-coupled
             F-P ejection state to the classical zero-velocity baseline,
             restores F-P curve visibility, returns presentation to Analysis
             Overlay, and invalidates the current MC result. This prevents a
             later F-P run from silently inheriting the previous MC velocity
             or contour-only display state.
           • PUBLICATION AXES: Projected-distance axes are now nucleus-centred
             in kilometres, use right-positive X and up-positive Y, guarantee
             a visible 0 tick at the nucleus whenever it lies inside the
             displayed frame, and use readable scientific tick labels.
    v3.1  • NEW: MC Display now provides Analysis Overlay, Contour
            Comparison, and Publication Figure presentation modes.
            Contour-only modes suppress the image background and unrelated
            F-P annotations, render observed image isophotes in black and
            Monte Carlo contours in magenta on white, and support projected-
            distance axes plus direct PNG/TIFF/PDF/SVG figure export.
    v3.1  • UI: Main window now opens maximized and applies an initial
            screen-responsive splitter layout so the control panel, plot
            canvas, and information panel fit the user's available display.
          • NEW: File > Run F-P Model (F5) runs the same validated
            Finson–Probstein workflow as the COMPUTE MODEL button, avoiding
            repeated scrolling in the left control panel.
    v3.1  • UX FIX: MC report files now default to MC_Report_... instead
            of MC_inputs_....
          • UX FIX: Image Setup remains above the CTA main window only;
            it no longer stays on top of unrelated application windows.
    v3.1  • UX: Maximum dust age is now derived automatically from the
            largest Synchrone age; the main F-P panel now accepts one
            user-interpreted Dominant dust age. MC Q(t) guidance applies the
            dominant age as a lower bound and maximum age as an upper cap.
    v3.1  • UI refinement: Q(t) preview and MC-window guidance moved into
            Simulation > Dust production over time; MC setup now uses
            three workflow tabs and wider bounded numeric inputs.
    v3.1  • NEW: Monte Carlo dust morphology module — MCWindow (View >
            Monte Carlo morphology…) simulates grain ensembles and
            extracts isophote-style contours, auto-displayed on the
            main canvas alongside F-P syndynes/synchrones. Grain size
            power-law distribution, ejection velocity V = V₀·β^γ·r_H^(−κ),
            optional sunward or active-area (rotating nucleus) direction
            mode with cos(z) term. Density ρ, albedo p_v, phase function.
          • NEW: Save/Load .mcin input files — machine-readable JSON;
            state memory across sessions; comet-mismatch warning on load.
          • NEW: MC report — preview then save plain-text summary;
            shows actual grain ejection speed, not raw V₀.
          • NEW: Contour min. scale slider (AU) — dynamic floor =
            max(MC grid pixel, image AU/px); auto-updates per run
            and image scale; user value preserved across re-runs.
          • NEW: Sun direction / Anti-velocity (−v) / Nucleus crosshair
            toggle checkboxes in DISPLAY panel.
          • NEW: Sun/anti-velocity arrows shortened 50% (6.5% view width).
          • BUG FIX: percentile floor for contour extraction now computed
            from pre-Gaussian density — smoothing no longer shifts
            contour positions (before: more smooth → lower floor →
            contours expanded outward from coma).
          • BUG FIX: loaded smooth_sigma_au no longer overwritten by
            auto-fill on first run after loading .mcin.
          • BUG FIX: MC overlay info box now shows real grain speed
            (V₀·β^γ·r_H^(−κ)) instead of raw V₀.
          • MC contour visibility now automatic (no checkbox required);
            MC appearance controls moved to MC window Display tab.
          • V₀ spinbox range extended to 2000 m/s for interstellar
            comet models (was 500 m/s).
    v3.0.7 • NEW: Orbit Position window gained a "🔍 Zoom to inner solar
            system" checkbox + AU spinbox (default 2.0). For comets well
            beyond ~1 AU, the default view auto-fits to the whole orbit,
            which collapses Earth and the Sun down to a barely-
            distinguishable point near the center. Enabling this instead
            forces a fixed Sun-centered view of the chosen AU half-width
            — Earth's orbit (and its separation from the Sun) is then
            always clearly visible, at the cost of the comet itself
            being off-screen if it's currently farther out than that
            value. Reuses the already-computed orbit diagram data —
            toggling only changes the view, not a recomputation. Pairs
            with comet_tail_analyzer.py v3.0.2.
    v3.0.6 • BUG FIX: pre-filling the observation date from a FITS
            file's DATE-OBS header keyword truncated it to date_obs[:10]
            (date only), silently discarding the time-of-day. date_to_jd()
            already accepts the full ISO-8601 'YYYY-MM-DDTHH:MM:SS[.fff]'
            format DATE-OBS normally uses, so the whole string is now
            passed through. Also handles the older convention of a
            separate TIME-OBS keyword for FITS files that split date and
            time into two keywords instead of combining them in DATE-OBS.
          • NEW: Image Setup ▸ "VIEW FITS HEADER" button opens a read-only
            viewer (FitsHeaderDialog) showing every raw header card of the
            currently loaded FITS file, with a Copy button. No-op with an
            explanatory message for JPEG/PNG images, which have no header.
    v3.0.5 • BUG FIX: the Sun-direction and anti-velocity (−v) arrows in
            overlay mode used a FIXED length (60 × au_per_px, i.e. always
            exactly 60 pixels of the ORIGINAL unzoomed image) — once the
            new free-form Zoom feature (v3.0.2-v3.0.4) let the view get
            much smaller than the full image, that fixed 60px length
            started spanning a huge fraction of the now-tiny visible
            frame. Changed to scale with 13% of the CURRENT (possibly
            zoomed) view width instead — same percentage-of-view-width
            approach the standalone (no-image) branch already used, just
            applied to the overlay branch too. Uses a view_w_px value
            captured right when the view is set, not a fresh
            ax.get_xlim() read at arrow-draw time, so it can't be thrown
            off by matplotlib transiently auto-expanding the axes while
            plotting curves in between (the same effect the existing
            end-of-function safety reset guards against).
    v3.0.4 • BUG FIX: v3.0.3's reconnect-after-clear() fix didn't actually
            solve it either — Compute Model still reset to the full image
            regardless of toolbar pan/zoom. Replaced the whole callback-
            based approach with something architecturally simpler and
            more robust: draw_model()/draw_image_preview() now read the
            axes' CURRENT xlim/ylim directly (via ax.get_xlim()/
            get_ylim()) at the top of the function, before ax.clear()
            wipes it, and reuse that view when redrawing the overlay
            instead of resetting to the full image. No separate lock
            variables, no callbacks, nothing that can silently go stale —
            the axes' own state IS the source of truth. RESET TO FULL
            IMAGE now just sets the axes directly, which the next redraw
            picks up the same way. Removed: PlotCanvas._locked_xlim/
            _locked_ylim/_suppress_lock_capture/_connect_view_callbacks/
            _on_axes_view_changed (all dead weight from the abandoned
            v3.0.2/v3.0.3 approach).
    v3.0.3 • BUG FIX: the v3.0.2 view-lock redesign never actually worked
            — Compute Model kept resetting to the full image regardless
            of any toolbar pan/zoom. Cause: self.ax.clear() silently
            resets matplotlib's internal callback registry to empty, and
            draw_image_preview()/draw_model() both call ax.clear() on
            EVERY redraw — so the xlim_changed/ylim_changed connections
            made once in PlotCanvas.__init__() went dead the moment the
            very first redraw happened (the initial image preview),
            before the user ever got a chance to pan/zoom. No error was
            raised; the callback just silently stopped firing. Fixed by
            extracting the connection into _connect_view_callbacks() and
            calling it again after every single ax.clear() (3 call
            sites: _draw_empty, draw_image_preview, draw_model) instead
            of relying on the one-time __init__ connection.
    v3.0.2 • Redesigned the Zoom feature from v3.0/v3.0.1's discrete
            1x/2x/3x/4x presets (always a square box centered on the
            nucleus) to free-form use of the plot toolbar's own zoom
            (rubber-band) and pan (hand) tools directly on the overlay
            image. Whatever x/y range that leaves the axes at is now
            captured automatically (PlotCanvas._on_axes_view_changed,
            gated on an image actually being loaded rather than a
            zoom-level counter) and used by draw_model() for both the
            initial overlay setup and its end-of-function safety reset —
            so any arbitrary, independently-sized x and y range is
            respected, not just a square box. The 4 preset radio buttons
            are gone, replaced by a single "RESET TO FULL IMAGE" button
            that clears the lock.
    v3.0.1 • Extended Zoom from 1x/2x/3x to 1x/2x/3x/4x (same centered-on-
            nucleus, locked-across-Compute-Model behavior as v3.0 — see
            below — just one more preset level).
    v3.0  • NEW: Zoom 1x/2x/3x for the loaded overlay image (JPEG/PNG/
            FITS alike), in the Image Setup dialog next to Nucleus
            Position. Selecting 2x/3x centers a box of that size around
            the current nucleus position; drag-pan afterward with the
            plot toolbar's hand tool to fine-tune the comet's position
            within it. The resulting view is LOCKED — Compute Model no
            longer resets the overlay back to the full image every time
            it re-runs (the previous unconditional reset, including an
            end-of-draw "force back to image bounds" safety net, was
            wiping out any zoom/pan immediately). Back to 1x clears the
            lock. The lock also resets automatically whenever the
            underlying image itself changes (new file loaded, cleared,
            or rotated N-up), since a previous pixel-coordinate lock
            wouldn't correspond to a different image's dimensions.
            Implementation: PlotCanvas tracks zoom_level/_locked_xlim/
            _locked_ylim and listens for matplotlib's xlim_changed/
            ylim_changed callbacks to capture manual toolbar pans into
            the lock, guarded by a _suppress_lock_capture flag so the
            programmatic resets draw_model() itself makes don't
            re-trigger capturing themselves.
    v3.0  • REMOVED the "Activity: VERY HIGH (outburst / hyperactive
            candidate)"-style categorical line from the Dust production
            rate… calculator's output. Afρ defaults to a placeholder
            value when the dialog opens — showing a confident-sounding
            classification derived from whatever happens to be in that
            box (real measurement or untouched placeholder, no way to
            tell apart from the display) risked being read as a
            diagnosis rather than what it actually is. Q_d itself is
            still shown — it's just deterministic math given the
            current inputs, for the user to interpret themselves.
    v3.0  • Dust particle radius… calculator's β range mode used to show
            "β = 0.05 – 1 → radius = 1.1 – 23.0 µm" — two ranges the
            reader had to mentally cross-match in REVERSE order (smaller
            β gives the LARGER radius), easy to misread as a direct
            left-to-right pairing. Now shows each β explicitly paired
            with its own radius on its own line instead.
    v3.0  • Default tuning: Auto Stretch now sets Contrast=-150 (was
            200), Intensity=-300 (unchanged); Isophote levels default
            kept at 6 levels and 2.0 px smoothing by default.
    v3.0  • Contrast/Intensity now restretch the image live as the
            sliders move — REMOVED the separate "Apply Stretch" button
            entirely. Auto Stretch still exists (sets both sliders to a
            reasonable starting point in one click).
    v3.0  • Simplified the FITS image stretch from 4 technical sliders
            (Shadow %, Highlight %, Softness, Gamma) to 2 simple ones
            (Contrast, Intensity, both -1000..1000), matching the
            simpler two-knob stretch UI found in some comet-tracking
            viewers. This is a reparameterization of the same underlying
            asinh stretch — Intensity shifts the black point relative to
            the image's robust background level (negative reveals more
            faint coma/tail, positive shows only bright cores); Contrast
            sets how wide the stretched range above black is (higher =
            steeper/more contrasty). NOT a verified match to any other
            tool's exact internal formula — defaults (Auto Stretch:
            Contrast=200, Intensity=-300) are starting points to tune by
            eye, not a derived optimum.
    v3.0  • Animator's Fixed FOV now auto-syncs to the loaded overlay
            image's own angular width (image pixel width × au_per_px ×
            deg/AU at the current working date) every time Compute frames
            is clicked, if an image is loaded — so animation frames are
            shown at the same scale the overlay was built at, rather than
            whatever Fixed FOV size happened to be left over. No-op with
            no image loaded. A status-bar message reports the synced
            value (and the inputs used) each time, since it overrides
            whatever was in the Full width field.
          • DISPLAY panel: added a "Syndynes/Synchrones:" label above the
            general Width/Opacity sliders, matching the "Orbital path:"
            label added for its own sliders just below — previously only
            the orbital path ones were labelled, leaving the general
            pair's target ambiguous by comparison.
    v3.0  • Released as v3.0 — this version number now covers everything
            accumulated under the "v2.5 (unreleased, not yet pushed to
            GitHub)" working label across the prior development session:
            the embedded Animator, the Dust particle radius…/Dust
            production rate… calculators, the Animator↔Orbit↔EPHEMERIS
            date linking, the equal-aspect-ratio and N/E-compass/starfield
            fixes, the β TABLE/LIGHT CURVE/ANALYSIS tab removals, and the
            bug fixes below. Pairs with comet_tail_analyzer.py v3.0.
          • Orbital path now has its own dedicated Line width/Opacity
            sliders (DISPLAY panel), separate from the syndyne/synchrone
            ones — it was previously hardcoded at lw=0.9/alpha=0.5 with
            no way to adjust it, which against a busy starfield/isophote
            background made it nearly invisible (default raised to
            lw=1.0/alpha=0.7, both still freely adjustable down or up).
    v3.0  • EPHEMERIS (r☉, Δ, Phase, RA, Dec, Date) now also tracks the
            currently-shown Animator frame, alongside obs_date — fully
            consistent with the Orbit-view link added just above (orbital
            elements themselves don't change with date, so only that half
            of update_info() actually changes per frame).
          • BUG FIX: every frame computed by the Animator showed the
            literal title "Comet" instead of the actual comet's name.
            Cause: compute_model() itself never sets info['name'] — only
            _on_model_ready() does, for the regular Compute Model flow;
            the Animator's worker called compute_model() directly and
            never added it, so the title fell back to draw_model()'s
            info.get("name","Comet") default. Fixed by setting it in the
            worker the same way _on_model_ready() does.
          • REMOVED the duplicate in-plot "PsAng = …° (anti-solar / tail
            dir)" text label — PsAng is already shown in the plot's top
            title line, no need to repeat it as a separate annotation
            inside the plot area.
    v3.0  • The embedded Animator now re-links automatically to whichever
            comet is currently selected (new ControlPanel.comet_ready
            signal, emitted after preset/FETCH JPL/manual entry all set
            _comet_el): Start date = that comet's obs date, End date =
            Start+360d. Switching to a different comet always disables
            Play/the scrub slider and clears any previously computed
            frames (they belonged to the old comet) — Compute frames must
            be run again for the new one rather than silently keep
            showing stale results.
          • Animator ↔ Orbit linking: scrubbing/playing the Animator now
            keeps the obs_date field in lockstep with whichever frame is
            showing, and _current_obs_jd() (read by View > Orbit position
            diagram… and the Dust production rate… calculator) now
            prefers that field over the comet's static fetch-time obs_jd.
            Net effect: scrub to an interesting date in the animation,
            open Orbit, and it shows that exact date — no extra step.
          • Orbit position diagram… gained a direction-of-motion arrow at
            the comet's current position (heliocentric velocity unit
            vector from compute_orbit_diagram(), previously discarded —
            see comet_tail_analyzer.py CHANGELOG).
    v3.0  • Compute Model's standalone (no image) view now defaults to a
            600 arcmin Fixed FOV instead of the old auto-fit percentile
            box, so it's the same starting frame size the embedded
            Animator below would continue from — no jump in apparent
            zoom between looking at the static model and playing the
            animation. The Animator's own frame-size field defaults to
            600 arcmin for the same reason (both still freely editable).
            Overlay mode (image loaded) is completely unaffected — still
            always fits the plot area to the image's own native pixel
            extent exactly as in the earliest versions; this default only
            applies when there's no image to fit to.
          • Animator's default End date changed from Start+30 to
            Start+360 days.
    v3.0  • Animator moved from a separate popup window into the main
            window itself: its controls (date range, step, Fixed
            FOV/distance toggle, Auto-suggest, Compute frames, play/
            scrub) now live in a new "ANIMATOR" group box in the INFO
            tab, right below ORBITAL ELEMENTS. Frames render directly
            into the MAIN canvas (same PlotCanvas the regular Compute
            Model view uses) instead of a dedicated dialog canvas — no
            more separate window to manage. Orchestration (worker
            thread, frame cache, play timer) moved from the AnimatorWindow
            class (removed) into MainWindow's new _anim_* methods, which
            gather comet/β/age inputs fresh from the model panel on every
            Compute/Auto-suggest click rather than once at dialog-open
            time, since this panel now stays open across comet changes.
          • REMOVED the ANALYSIS tab and the "Generate analysis" menu
            item entirely. Its Afρ-based dust-production report
            duplicated the standalone Calculation > Dust production
            rate… calculator (same compute_dust_production_rate()
            formula) added earlier in v3.0, just shown a second way.
            generate_dust_analysis() is left in comet_tail_analyzer.py as
            a library function in case a text-report version is wanted
            again later; only the GUI tab/menu item/caller were removed.
            Also removed as a result: _ensure_model_computed() and the
            self._model_for/_pending_model_callback bookkeeping it used
            (dead code — Generate analysis was their only caller; the
            embedded Animator above manages its own independent state
            and never reuses the main Compute Model pipeline).
          • The right panel now has only one tab (INFO) after the β
            TABLE/LIGHT CURVE/ANALYSIS removals across v3.0 — its tab
            bar auto-hides (QTabWidget.setTabBarAutoHide) instead of
            showing a single redundant "INFO" tab button, and the
            panel's header strip was relabelled "INFO" (was "ANALYSIS",
            stale since that tab is gone).
    v3.0  • BUG FIX: the N/E compass (and to a lesser extent the now-
            removed starfield, see below) ended up positioned near the
            middle of the standalone plot instead of its intended
            bottom-right corner, as a direct side effect of the
            adjustable='datalim' fix above — set_aspect(...,
            adjustable='datalim') only computes the final, box-filling
            axis limits LAZILY at actual draw time, so the compass-
            position code's self.ax.get_xlim()/get_ylim() calls (which
            run immediately after set_aspect(), not at draw time) were
            reading the stale PRE-adjustment limits. Fixed by calling
            self.ax.apply_aspect() right after set_aspect() to force that
            computation to happen immediately — the compass formula
            itself was already correctly targeting the bottom-right
            corner, it just had the wrong numbers to work with.
          • REMOVED the random starfield (200 scattered white dots) from
            the standalone plot — purely decorative, no astrometric
            meaning, and was also reading the same stale limits as the
            compass bug above.
    v3.0  • BUG FIX: the main window's standalone (no-image) plot and the
            Animator's plot showed visibly different tail shapes for the
            EXACT same model/comet/date — main window looked stretched
            and steep, Animator looked flatter and closer to how the
            tail actually looks in real images. Cause: draw_model()'s
            standalone branch never called ax.set_aspect('equal'), so
            matplotlib used its default "auto" aspect — stretching the
            x/y degree axes independently to fill whatever pixel box the
            host canvas widget happened to have. The main window's canvas
            is much wider than the Animator's embedded one, so the same
            angular data got distorted by two different amounts. Fixed
            by adding ax.set_aspect('equal', ...) right after the
            xlim/ylim are set, with the adjustable mode split by use case:
              - Animator (fixed_xlim/fixed_ylim given): adjustable='box'
                — "Fixed FOV"/"Fixed distance" exist to show an EXACT
                requested angular size, so any aspect mismatch between
                that size and the dialog's canvas is taken up as blank
                padding, never by silently widening the requested FOV.
              - Main window auto-fit (no fixed_xlim/fixed_ylim):
                adjustable='datalim' — no precision promise to keep here,
                so instead of padding with empty space (which, for
                elongated tails far wider than tall, first showed up as
                the plot shrinking to a thin sliver with the legend now
                looming huge over it — far more confusing than the
                original distortion), matplotlib extends whichever axis
                range is needed to fill the canvas completely.
            1° in RA-direction now always renders as the same screen
            length as 1° in Dec-direction in both windows, matching true
            angular proportions as they'd appear in a real image, AND the
            main window's plot fills its canvas the way it always did.
    v3.0  • BUG FIX: Light curve…/Orbit position diagram…/Generate
            analysis/Dust production rate… kept showing the PREVIOUS
            comet after switching to a new one and clicking "Use this
            comet", unless Compute Model was pressed again first. Root
            cause: _active_comet_el() preferred self._comet_el (this
            window's own copy, only updated lazily after a compute
            finishes) over self.ctrl._comet_el (ControlPanel's copy,
            updated the instant a comet is selected/fetched — always the
            freshest). Swapped the priority; every menu action above
            goes through this one helper so the fix covers all of them.
          • BUG FIX: Orbit position diagram… stayed pinned in front of
            MainWindow at all times, including after MainWindow regained
            focus. Cause: OrbitWindow was constructed with parent=self
            (MainWindow) — a non-modal QDialog with a parent is kept
            "transient for" it by the window manager, which pins it
            above that specific parent permanently. LCWindow already
            avoided this with parent=None; applied the same fix here
            (and now stores the reference as self._orbit_win instead of
            a throwaway local, matching LCWindow's GC-safety pattern).
          • REMOVED the β TABLE and LIGHT CURVE tabs from the right-hand
            panel entirely (now just INFO + ANALYSIS). Both were fully
            redundant with dedicated menu-driven access added earlier in
            v3.0: grain radius lives in Calculation > Dust particle
            radius… (self-contained, no model needed), and the light
            curve (plot + H₀/n) is shown by the View > Light curve…
            popup, which already fetches COBS automatically. Removed
            along with them: get_selected_betas() and
            refresh_grain_radius_column() (operated on the now-gone
            table), and the now-pointless H₀/n/status mirroring in
            _fetch_cobs() (the status bar and the LCWindow popup already
            cover it). ANALYSIS is tab index 1 now, not 3.
    v3.0  • NEW: Calculation > Dust particle radius… — standalone β→radius
            calculator (single value or range), independent of any comet
            or computed model since the formula needs only β/ρ/Qpr.
            Always delegates to beta_to_size(), now extended with an
            optional Qpr parameter (default 1.0, backward compatible) so
            it can stay the single source of truth here too. The grain-
            radius section that used to print in the ANALYSIS tab (and
            the "Selected β → grain radii" use of the β TABLE's
            checkboxes) was removed — this calculator replaces it; the β
            TABLE's own radius column is unaffected.
          • NEW: Calculation > Dust production rate… — standalone Afρ→Q_d
            calculator. r_h is propagated from the selected comet's
            orbital elements for a user-entered date (read-only, like
            genuinely Horizons-sourced values) while Afρ/v_dust/p_v stay
            editable, each with its own literature-referenced default.
            Delegates to comet_tail_analyzer.py's new
            compute_dust_production_rate(), shared with the unchanged
            inline Afρ section in the ANALYSIS tab so the two formulas
            can never drift apart.
          • There is no separate combined "Physical parameters" dialog —
            ρ_d lives in Dust particle radius…, v_dust/p_v live in Dust
            production rate…, each with its own default and override,
            since each calculator is self-contained rather than reading
            from shared global state.
          • NEW: View > Animator… — steps obs date across a user-set
            start/end/step range and plays back the resulting F-P model
            frames, for spotting when a comet's tail is widest or its
            morphology changes fastest (e.g. to plan imaging sessions).
            Frame size is a toggle: "Fixed FOV" (arcmin, locks the real
            angular field of view including the effect of changing
            Earth-comet distance Δ — for imaging planning) or "Fixed
            distance" (AU, removes the Δ effect — for comparing the
            tail's actual physical growth). An Auto-suggest button samples
            the date range and proposes a size that fits the widest tail
            with ~20% margin. Reuses PlotCanvas.draw_model() for every
            frame via new optional fixed_xlim/fixed_ylim parameters
            (backward compatible — default None preserves the existing
            auto-fit behavior for the main view exactly), so animation
            frames render identically to the main view with zero
            duplicated rendering code.
          • View > Orbit position diagram… no longer requires Compute
            Model first — same auto-chain pattern as Light curve…, only
            a comet needs to be selected/fetched.
          • "Generate analysis" (Ctrl+G) replaces the old GENERATE
            ANALYSIS button.
          • View menu gained "Light curve…" (Ctrl+L) and "Animator…",
            alongside the existing Orbit diagram — consolidates every
            "open a view" action into one menu instead of splitting them
            between the menu bar and sidebar buttons.
          • REMOVED sidebar buttons: FETCH COBS, PLOT LC, GENERATE
            ANALYSIS. Their underlying actions are now reached only via
            the menu items above, which auto-chain their prerequisites:
            Light curve… fetches COBS automatically if not already
            fetched for the selected comet; Generate analysis computes
            the F-P model first if needed, then fetches COBS, then runs
            the analysis — no separate manual button presses required at
            each stage, only a comet needs to be selected/fetched first.
          • Added per-comet caching for both the computed model and COBS
            data (invalidated automatically when the selected comet
            changes) so the auto-chained actions above never silently
            reuse a previous comet's stale results.
          • Pairs with comet_tail_analyzer.py v3.0.
          • CRITICAL FIX: β→grain-radius formula was using an uncorrected
            two-parameter form (C_pr·Qpr/(ρr), C_pr=1.19e-3 kg/m²) that
            omitted a factor of 2 present in the convention that constant
            is normally paired with. Verified against Burns, Lamy & Soter
            (1979) Icarus 40, 1, Eq. 19 directly — all reported grain radii
            in v2.4 were too LARGE by a factor of 2.07. See
            comet_tail_analyzer.py CHANGELOG for the full derivation.
          • β-table column header "Size" → "Radius"; all "grain size"
            labels in the β panel and analysis report corrected to "grain
            radius" (the value was always a radius, but was unlabelled
            ambiguously — risk of being misread as diameter, doubling the
            error on top of the formula fix above).
          • Added anti-velocity arrow (−v): the negative of the comet's
            heliocentric velocity, projected onto the sky, drawn alongside
            the existing Sun-direction arrow on every plot/overlay. Shows
            how far a real (non-zero-ejection-velocity) tail is expected to
            lean from the pure antisolar line — same convention used by
            Moreno (2025, A&A 695 A263) and Mariblanca-Escalona et al.
            (2026, MNRAS) in their published tail images.
          • Added isophote overlay: optional surface-brightness contours
            traced directly from the loaded image, toggleable alongside
            the syndyne/synchrone/orbit checkboxes, for direct visual
            comparison between observed tail morphology and the model
            curves (overlay mode only). Gaussian-smooths the luminance
            and uses an 85th-percentile noise floor before contouring —
            without this, real (noisy) exposures produce solid green
            speckle across the whole frame instead of clean isophote
            curves, since contour() draws a segment at every pixel-noise
            crossing of a level. Levels-count and smoothing-radius are
            both adjustable sliders in the DISPLAY panel.
          • Added ORBIT VIEW window: a new 3D diagram of the comet's
            position on its actual orbital ellipse (Sun as a yellow
            circle at the focus, perihelion marked, Earth's orbit for
            scale, drop-lines to the ecliptic plane, click-drag to
            rotate) at the observation epoch — complements the existing
            orbital-path-on-sky overlay, which only shows the projected
            tail-axis direction, not where the comet physically sits in
            its orbit or how far out of the ecliptic plane it currently
            is.

    v2.4  • Pairs with comet_tail_analyzer.py v2.4.
          • PRESET tab: obs_date field was missing from layout (invisible) —
            fixed; _on_comet_selected(0) now called on init to seed first
            entry; "USE THIS COMET" button now triggers a live Horizons
            fetch (identical pipeline to FETCH JPL tab).
          • FETCH tab renamed "FETCH JPL"; MPC source dropdown removed
            (MPC API returns HTTP 403; JPL Horizons is the sole source).
          • FetchWorker: prefer parameter removed; always uses Horizons.
          • fetch_status / preset_status: PRESET tab now has its own status
            label; _fetch_source flag routes success/error to the correct
            label and re-enables the correct button.
          • date_to_jd: accepts decimal-day (2018-01-14.583), HH:MM,
            HH:MM:SS in all date fields.
          • LCWindow: observation-date vertical dashed line (blue) and
            perihelion vertical dashed line (orange) added to date-axis
            plot; "Now" vertical line and scatter star removed — current
            magnitude and distance shown as a text annotation instead.
          • Window title updated to v2.4; "local DB" source fallback
            replaced with "JPL Horizons"; bare print() in LC fallback
            replaced with logging.warning().

    v2.3  • Pairs with comet_tail_analyzer.py v2.3.
          • Fixed _effective_plate_scale_arcsec() to correctly delegate
            to ControlPanel attributes (self.ctrl).
          • Fixed btn_compute QSS override.
          • Fixed LCWindow._draw(): added missing timedelta import.
          • Removed dead imports and dead InfoPanel method.

    v2.2  • Pairs with comet_tail_analyzer.py v2.2.

    v2.1  • Fixed WCS auto-override; reset _wcs_ps_deg on new image load.
          • Fixed CD matrix plate-scale calculation.
          • Added auto-sync ps_spin ↔ au_px_spin.

    v2.0  • Replaced deprecated datetime.utcnow().
          • STYLE built dynamically by _build_style().
=============================================================================
  Run:
      python CometTailGUI.py

  Requires:
      pip install PyQt6 matplotlib numpy scipy astropy astroquery Pillow
=============================================================================
"""

__version__ = "3.1.1"
# Runtime patch: partial COBS coverage blocks auto-recommendation, not user-selected MC runs.

import sys, os, warnings, csv, webbrowser, logging
warnings.filterwarnings("ignore")

import numpy as np

# ── PyQt6 ────────────────────────────────────────────────────────────────────
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QLineEdit, QComboBox, QCheckBox,
    QSlider, QSpinBox, QDoubleSpinBox, QTabWidget, QGroupBox, QScrollArea,
    QSplitter, QFileDialog, QStatusBar, QProgressBar, QMessageBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QFrame, QTextEdit,
    QDialog, QDialogButtonBox, QFormLayout, QSizePolicy, QMenuBar,
    QMenu, QRadioButton, QButtonGroup, QToolBar,
)
from PyQt6.QtCore import (
    Qt, QThread, pyqtSignal, QTimer, QSize, QSettings,
)
from PyQt6.QtGui import (
    QFont, QIcon, QColor, QPalette, QAction, QPixmap, QCursor,
)

# ── Matplotlib ────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import matplotlib.patheffects as patheffects
from matplotlib.ticker import FuncFormatter

# ── Physics engine ────────────────────────────────────────────────────────────
import importlib.util

def _load_physics():
    candidates = [
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "comet_tail_analyzer.py"),
        os.path.join(os.getcwd(), "comet_tail_analyzer.py"),
    ]
    for path in candidates:
        if os.path.exists(path):
            spec = importlib.util.spec_from_file_location("cta", path)
            mod  = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    raise FileNotFoundError(
        "comet_tail_analyzer.py not found.\n"
        "Place it in the same folder as this file."
    )

try:
    cta = _load_physics()
except FileNotFoundError as e:
    app = QApplication(sys.argv)
    QMessageBox.critical(None, "Missing File", str(e))
    sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────────
#  THEME  SYSTEM  — SSTAC style: Dark + Light modes
#  (STYLE is built dynamically by _build_style() defined below)
# ─────────────────────────────────────────────────────────────────────────────

THEME_DARK = dict(
    # ── backgrounds (SSTAC-style deep dark) ──────────────────────────────
    win_bg      = "#0d1117",
    panel_bg    = "#161b22",
    input_bg    = "#1c2128",
    alt_bg      = "#1a1f27",
    group_bg    = "#161b22",
    # ── text ─────────────────────────────────────────────────────────────
    text        = "#e6edf3",
    text_muted  = "#8b949e",
    text_dim    = "#484f58",
    text_value  = "#58a6ff",
    text_sec    = "#79c0ff",
    # ── borders ──────────────────────────────────────────────────────────
    border      = "#30363d",
    border2     = "#3d4450",
    border3     = "#58a6ff",
    # ── accents (SSTAC green/blue/orange) ────────────────────────────────
    accent      = "#238636",   # SSTAC primary green
    accent2     = "#2ea043",   # lighter green
    btn_bg      = "#21262d",
    btn_hover   = "#30363d",
    highlight   = "#388bfd26",
    # ── button class colours (SSTAC multi-colour) ────────────────────────
    btn_primary   = "#1f6feb",  # blue — fetch / compute
    btn_success   = "#238636",  # green — generate / confirm
    btn_warning   = "#bb8009",  # amber — archive / save
    btn_info      = "#1158c7",  # deep blue — sky map / tools
    btn_danger    = "#b91c1c",  # red — clear / reset
    # ── matplotlib ───────────────────────────────────────────────────────
    mpl_bg      = "#0d1117",
    mpl_grid    = "#21262d",
    mpl_zero    = "#30363d",
    mpl_star    = "#c9d1d9",
    mpl_star_a  = 0.18,
    mpl_title   = "#58a6ff",
    mpl_label   = "#8b949e",
    mpl_tick    = "#6e7681",
    canvas_outer= "#0d1117",
)

THEME_LIGHT = dict(
    # ── backgrounds ──────────────────────────────────────────────────────
    win_bg      = "#f6f8fa",
    panel_bg    = "#ffffff",
    input_bg    = "#f6f8fa",
    alt_bg      = "#eaeef2",
    group_bg    = "#f6f8fa",
    # ── text ─────────────────────────────────────────────────────────────
    text        = "#24292f",
    text_muted  = "#57606a",
    text_dim    = "#8c959f",
    text_value  = "#0550ae",
    text_sec    = "#0969da",
    # ── borders ──────────────────────────────────────────────────────────
    border      = "#d0d7de",
    border2     = "#c6cdd4",
    border3     = "#0969da",
    # ── accents ──────────────────────────────────────────────────────────
    accent      = "#2da44e",
    accent2     = "#1a7f37",
    btn_bg      = "#f6f8fa",
    btn_hover   = "#eaeef2",
    highlight   = "#ddf4ff",
    # ── button class colours ─────────────────────────────────────────────
    btn_primary   = "#0969da",
    btn_success   = "#2da44e",
    btn_warning   = "#bf8700",
    btn_info      = "#1b4bab",
    btn_danger    = "#cf222e",
    # ── matplotlib ───────────────────────────────────────────────────────
    mpl_bg      = "#ffffff",
    mpl_grid    = "#eaeef2",
    mpl_zero    = "#c6cdd4",
    mpl_star    = "#424a53",
    mpl_star_a  = 0.10,
    mpl_title   = "#0969da",
    mpl_label   = "#57606a",
    mpl_tick    = "#6e7781",
    canvas_outer= "#f6f8fa",
)

# Active theme (start dark)
CURRENT_THEME = "dark"
T = THEME_DARK.copy()

def _build_style(t: dict) -> str:
    """Build QSS from theme dict — SSTAC visual style."""
    # Derive button colours (with fallback if key absent)
    bp  = t.get("btn_primary",   t["accent"])
    bs  = t.get("btn_success",   t["accent2"])
    bw  = t.get("btn_warning",   "#b08000")
    bi  = t.get("btn_info",      t["accent"])
    bd  = t.get("btn_danger",    "#c02020")
    return f"""
/* ── Global — SSTAC style ────────────────────────────────────────────── */
* {{
    font-family: "Segoe UI", "SF Pro Display", "Helvetica Neue", Arial, sans-serif;
    font-size: 12px;
    color: {t['text']};
}}
QMainWindow, QDialog {{ background: {t['win_bg']}; }}
QWidget {{ background: {t['win_bg']}; color: {t['text']}; }}

/* ── Panels ──────────────────────────────────────────────────────────── */
QGroupBox {{
    background: {t['group_bg']};
    border: 1px solid {t['border']};
    border-radius: 6px;
    margin-top: 20px;
    padding: 10px 8px 8px 8px;
    font-size: 10px;
    color: {t['text_muted']};
    letter-spacing: 1.5px;
    font-weight: 600;
    text-transform: uppercase;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 2px 8px;
    background: {t['group_bg']};
    color: {t['text_muted']};
    left: 10px; top: 3px;
    border-radius: 3px;
}}

/* ── Scroll ──────────────────────────────────────────────────────────── */
QScrollArea {{ border: none; background: transparent; }}
QScrollBar:vertical {{
    background: {t['panel_bg']}; width: 6px; border-radius: 3px; margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {t['border2']}; border-radius: 3px; min-height: 24px;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
QScrollBar:horizontal {{
    background: {t['panel_bg']}; height: 6px; border-radius: 3px;
}}
QScrollBar::handle:horizontal {{ background: {t['border2']}; border-radius: 3px; }}

/* ── Inputs ──────────────────────────────────────────────────────────── */
QLineEdit, QSpinBox, QDoubleSpinBox, QTextEdit {{
    background: {t['input_bg']};
    border: 1px solid {t['border']};
    border-radius: 6px;
    padding: 5px 9px;
    color: {t['text']};
    selection-background-color: {t['highlight']};
    min-height: 28px;
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {t['border3']};
    outline: none;
}}
QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{
    background: {t['panel_bg']}; border: 1px solid {t['border']}; width: 18px;
    border-radius: 0 6px 6px 0;
}}

/* ── Combo ───────────────────────────────────────────────────────────── */
QComboBox {{
    background: {t['input_bg']};
    border: 1px solid {t['border']};
    border-radius: 6px;
    padding: 5px 9px;
    color: {t['text']};
    min-height: 28px;
}}
QComboBox:hover {{ border-color: {t['border2']}; }}
QComboBox::drop-down {{ border: none; width: 26px; border-radius: 0 6px 6px 0; }}
QComboBox QAbstractItemView {{
    background: {t['panel_bg']};
    border: 1px solid {t['border2']};
    color: {t['text']};
    selection-background-color: {t['highlight']};
    border-radius: 6px;
    padding: 2px;
}}

/* ── Buttons — base (neutral outline) ───────────────────────────────── */
QPushButton {{
    background: {t['btn_bg']};
    border: 1px solid {t['border']};
    border-radius: 6px;
    color: {t['text']};
    padding: 5px 14px;
    min-height: 28px;
    font-size: 12px;
    font-weight: 500;
}}
QPushButton:hover {{
    background: {t['btn_hover']};
    border-color: {t['border2']};
}}
QPushButton:pressed {{ border-color: {t['border3']}; }}
QPushButton:disabled {{ color: {t['text_dim']}; border-color: {t['border']}; opacity: 0.5; }}

/* ── Button colour classes (SSTAC multi-colour) ──────────────────────── */
QPushButton[class="primary"] {{
    background: {bp};
    border: none;
    border-radius: 6px;
    color: #ffffff;
    font-size: 13px;
    font-weight: 600;
    min-height: 34px;
    letter-spacing: 0.5px;
}}
QPushButton[class="primary"]:hover {{ background: {t['accent2']}; }}

QPushButton[class="success"] {{
    background: {bs};
    border: none; border-radius: 6px; color: #ffffff;
    font-weight: 600; min-height: 28px;
}}
QPushButton[class="success"]:hover {{ filter: brightness(1.15); }}

QPushButton[class="warning"] {{
    background: {bw};
    border: none; border-radius: 6px; color: #ffffff;
    font-weight: 600; min-height: 28px;
}}

QPushButton[class="info"] {{
    background: {bi};
    border: none; border-radius: 6px; color: #ffffff;
    font-weight: 600; min-height: 28px;
}}

QPushButton[class="danger"] {{
    background: {bd};
    border: none; border-radius: 6px; color: #ffffff;
    font-weight: 600;
}}
QPushButton[class="danger"]:hover {{ background: #e02020; }}

/* ── COMPUTE button override (always green like SSTAC Generate Plan) ─── */
QPushButton#btn_compute {{
    background: {bs};
    border: none; border-radius: 6px;
    color: #ffffff; font-size: 13px; font-weight: 700;
    min-height: 36px; letter-spacing: 0.5px;
}}
QPushButton#btn_compute:hover {{ background: {t.get('accent2', bs)}; }}

/* ── Tabs ────────────────────────────────────────────────────────────── */
QTabWidget::pane {{
    border: 1px solid {t['border']};
    border-radius: 6px;
    background: {t['group_bg']};
    top: -1px;
}}
QTabBar::tab {{
    background: transparent;
    color: {t['text_muted']};
    border: none;
    border-bottom: 2px solid transparent;
    padding: 6px 14px;
    margin-right: 2px;
    font-size: 11px;
    font-weight: 500;
    min-width: 55px;
}}
QTabBar::tab:selected {{
    color: {t['border3']};
    border-bottom: 2px solid {t['border3']};
    background: transparent;
}}
QTabBar::tab:hover:!selected {{ color: {t['text']}; }}

/* ── Checkboxes ──────────────────────────────────────────────────────── */
QCheckBox {{ spacing: 7px; color: {t['text']}; font-size: 12px; }}
QCheckBox::indicator {{
    width: 15px; height: 15px;
    border: 1px solid {t['border2']};
    border-radius: 3px;
    background: {t['input_bg']};
}}
QCheckBox::indicator:checked {{
    background: {t.get('btn_success', t['accent'])};
    border-color: {t.get('btn_success', t['accent'])};
}}

/* ── Sliders ─────────────────────────────────────────────────────────── */
QSlider::groove:horizontal {{
    height: 4px; background: {t['border']}; border-radius: 2px;
}}
QSlider::handle:horizontal {{
    width: 14px; height: 14px; margin: -5px 0;
    background: {t.get('btn_primary', t['accent'])};
    border-radius: 7px;
}}
QSlider::sub-page:horizontal {{
    background: {t.get('btn_primary', t['accent'])}; border-radius: 2px;
}}

/* ── Labels ──────────────────────────────────────────────────────────── */
QLabel[class="section"] {{
    color: {t['text_muted']};
    font-size: 10px;
    letter-spacing: 1.5px;
    font-weight: 600;
    padding: 2px 0;
    text-transform: uppercase;
}}
QLabel[class="muted"]  {{ color: {t['text_muted']}; font-size: 11px; }}
QLabel[class="value"]  {{ color: {t['text_value']}; font-weight: 600; font-size: 13px; }}
QLabel[class="title"]  {{
    color: {t['text_value']};
    font-size: 15px;
    font-weight: 700;
    letter-spacing: 0.5px;
}}

/* ── Table ───────────────────────────────────────────────────────────── */
QTableWidget {{
    background: {t['input_bg']};
    border: 1px solid {t['border']};
    gridline-color: {t['border']};
    color: {t['text']};
    alternate-background-color: {t['alt_bg']};
    border-radius: 6px;
}}
QTableWidget::item:selected {{ background: {t['highlight']}; color: {t['text']}; }}
QHeaderView::section {{
    background: {t['panel_bg']};
    color: {t['text_muted']};
    border: none;
    border-bottom: 1px solid {t['border']};
    border-right: 1px solid {t['border']};
    padding: 5px 8px;
    font-weight: 600;
    font-size: 11px;
    letter-spacing: 0.5px;
}}

/* ── Status bar ──────────────────────────────────────────────────────── */
QStatusBar {{
    background: {t['panel_bg']};
    border-top: 1px solid {t['border']};
    color: {t['text_muted']};
    font-size: 11px;
}}
QStatusBar::item {{ border: none; }}

/* ── Menu ────────────────────────────────────────────────────────────── */
QMenuBar {{
    background: {t['panel_bg']};
    border-bottom: 1px solid {t['border']};
    color: {t['text']};
    padding: 2px;
}}
QMenuBar::item {{ padding: 4px 10px; border-radius: 4px; }}
QMenuBar::item:selected {{ background: {t['btn_hover']}; color: {t['text']}; }}
QMenu {{
    background: {t['panel_bg']};
    border: 1px solid {t['border2']};
    color: {t['text']};
    border-radius: 6px;
    padding: 4px;
}}
QMenu::item {{ padding: 5px 20px; border-radius: 4px; }}
QMenu::item:selected {{ background: {t['highlight']}; }}
QMenu::separator {{ background: {t['border']}; height: 1px; margin: 4px 10px; }}

/* ── Toolbar ─────────────────────────────────────────────────────────── */
QToolBar {{
    background: {t['panel_bg']};
    border: none;
    border-bottom: 1px solid {t['border']};
    spacing: 4px; padding: 3px 6px;
}}
QToolBar::separator {{
    background: {t['border']}; width: 1px; margin: 4px 2px;
}}

/* ── Progress bar ────────────────────────────────────────────────────── */
QProgressBar {{
    background: {t['input_bg']};
    border: 1px solid {t['border']};
    border-radius: 4px;
    text-align: center;
    color: {t['text_muted']};
    height: 6px;
    max-height: 6px;
}}
QProgressBar::chunk {{
    background: {t.get('btn_success', t['accent'])};
    border-radius: 4px;
}}

/* ── Splitter ────────────────────────────────────────────────────────── */
QSplitter::handle {{ background: {t['border']}; width: 1px; height: 1px; }}
QSplitter::handle:hover {{ background: {t['border3']}; }}

/* ── Dialog ──────────────────────────────────────────────────────────── */
QDialog {{ background: {t['win_bg']}; }}
"""


STYLE = _build_style(THEME_DARK)

# ─────────────────────────────────────────────────────────────────────────────
#  COLORS  (matplotlib — also theme-aware, updated at switch)
# ─────────────────────────────────────────────────────────────────────────────
SYNDYNE_COLORS = ["#ff4444","#ff6b35","#ff8c42","#ffaa33","#ffcc22","#ffdd60","#ffe89a"]
SYNC_COLORS    = ["#ffd700","#f0b800","#e09a00","#cc8000","#b86800","#a05000","#883a00"]
ISOPHOTE_COLOR = "#60e0a0"   # cool green — distinct from all syndyne/sync/sun/orbit colors
BG      = T["mpl_bg"]
PANEL_BG= T["panel_bg"]
MF      = "DejaVu Sans"   # clean sans-serif for plot labels


def _compute_isophote_levels(img_arr, n_levels: int = 6, smooth_sigma: float = 2.0):
    """
    Build surface-brightness contour levels directly from a loaded image,
    without requiring photometric calibration (v3.0 isophote overlay).

    Since we have no zeropoint for an arbitrary user image, we approximate
    the conventional "isophotes spaced evenly in magnitude" look (e.g.
    Mariblanca-Escalona et al. 2026, Fig. 1: 0.5-mag spacing) by using
    LOG-spaced levels in raw pixel intensity — magnitude is itself
    -2.5·log10(intensity), so even log-intensity spacing is a reasonable,
    calibration-free stand-in for even magnitude spacing.

    v3.0: real astronomical exposures have substantial pixel-to-pixel
    shot/read noise. Contouring the raw frame makes matplotlib draw a tiny
    line segment at every noise fluctuation that crosses a level — over an
    entire noisy frame this looks like solid speckle rather than clean
    isophote curves (this is exactly what happens without smoothing; no
    isophote-fitting tool contours a raw frame for this reason). Two
    changes fix it:
      1. Gaussian-smooth the luminance first (kills sub-PSF pixel noise
         while preserving the coma/tail's real, much larger-scale shape).
      2. Raise the lower percentile floor well above the noise floor
         (85th, not 55th) so contours don't trace the dark sky background
         at all — only the genuinely bright coma/tail structure.

    Returns (luminance_2d, levels_array), or (None, None) if the image
    doesn't have enough usable signal to contour meaningfully.
    """
    if img_arr is None:
        return None, None
    arr = img_arr
    if arr.ndim == 3:
        # Rec.601 luminance from RGB(A); alpha channel (if any) is ignored
        weights = np.array([0.299, 0.587, 0.114])
        lum = arr[..., :3].astype(float) @ weights
    else:
        lum = arr.astype(float)

    if smooth_sigma > 0:
        try:
            from scipy.ndimage import gaussian_filter
            lum = gaussian_filter(lum, sigma=smooth_sigma)
        except Exception:
            pass   # scipy missing/failed — fall back to unsmoothed (still works, just noisier)

    finite = lum[np.isfinite(lum)]
    if finite.size < 100:
        return None, None
    lo, hi = np.percentile(finite, [85.0, 99.7])
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        return lum, None
    levels = np.geomspace(max(lo, 1e-6), hi, n_levels)
    return lum, levels


def apply_theme(theme_name: str, app: "QApplication"):
    """Switch application theme at runtime."""
    global CURRENT_THEME, T, BG, PANEL_BG, STYLE
    CURRENT_THEME = theme_name
    T = THEME_DARK.copy() if theme_name == "dark" else THEME_LIGHT.copy()
    BG       = T["mpl_bg"]
    PANEL_BG = T["panel_bg"]
    STYLE    = _build_style(T)
    app.setStyleSheet(STYLE)
    # Update palette
    from PyQt6.QtGui import QPalette, QColor
    pal = app.palette()
    pal.setColor(QPalette.ColorRole.Window,          QColor(T["win_bg"]))
    pal.setColor(QPalette.ColorRole.WindowText,      QColor(T["text"]))
    pal.setColor(QPalette.ColorRole.Base,            QColor(T["input_bg"]))
    pal.setColor(QPalette.ColorRole.AlternateBase,   QColor(T["alt_bg"]))
    pal.setColor(QPalette.ColorRole.Text,            QColor(T["text"]))
    pal.setColor(QPalette.ColorRole.Button,          QColor(T["panel_bg"]))
    pal.setColor(QPalette.ColorRole.ButtonText,      QColor(T["accent2"]))
    pal.setColor(QPalette.ColorRole.Highlight,       QColor(T["highlight"]))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor(T["text"]))
    app.setPalette(pal)


# ─────────────────────────────────────────────────────────────────────────────
#  WORKER THREAD  (keeps UI responsive during computation)
# ─────────────────────────────────────────────────────────────────────────────
class ComputeWorker(QThread):
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)

    def __init__(self, comet_el, obs_jd, betas, ages, max_age, n_pts, ejection=None):
        super().__init__()
        self.comet_el = comet_el
        self.obs_jd   = obs_jd
        self.betas    = betas
        self.ages     = ages
        self.max_age  = max_age
        self.n_pts    = n_pts
        self.ejection = ejection

    def run(self):
        try:
            self.progress.emit(10, "Computing dust trajectories…")
            model = cta.compute_model(
                self.comet_el, self.obs_jd,
                self.betas, self.ages,
                self.max_age, self.n_pts,
                ejection=self.ejection)
            self.progress.emit(100, "Done")
            self.finished.emit(model)
        except Exception as ex:
            self.error.emit(str(ex))


class FetchWorker(QThread):
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)

    def __init__(self, desig, date):
        super().__init__()
        self.desig = desig
        self.date  = date

    def run(self):
        try:
            el = cta.fetch_comet(self.desig, date=self.date or None)
            self.finished.emit(el)
        except Exception as ex:
            self.error.emit(str(ex))


class UpdateCheckWorker(QThread):
    """Run the GitHub Releases update check without blocking the GUI.

    The worker returns a diagnostic dictionary with status ``update``,
    ``current`` or ``error``. Silent startup checks ignore non-update results,
    while an explicit manual check can now distinguish an up-to-date install
    from a network, SSL, API, or rate-limit failure.
    """
    finished = pyqtSignal(object)

    def run(self):
        try:
            if hasattr(cta, "get_update_status"):
                result = cta.get_update_status()
            else:
                legacy = cta.check_for_update()
                result = legacy or {"status": "current", "installed": cta.__version__}
        except Exception as exc:
            result = {"status": "error", "reason": str(exc) or "Unknown update-check error."}
        self.finished.emit(result)


class MCWorker(QThread):
    """v3.1 — runs cta.compute_morphology_mc() off the UI thread. Unlike
    ComputeWorker's per-curve loop (a handful of syndynes/synchrones,
    fast), this loops over n_particles (potentially thousands) of
    individual dust_position_isotropic() calls — slow enough on the UI
    thread to visibly freeze the window, hence its own worker rather
    than reusing ComputeWorker."""
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)

    def __init__(self, comet_el, obs_jd, beta_range, gamma_size, max_age,
                n_particles, v0_coeff, gamma, m_exp, seed=None,
                north_pa_deg=0.0, rotation_offset_deg=0.0, qt_weights=None,
                active_area=None, rho_g_cm3=0.5, sunward=False, sunward_expocos=1.0,
                sunward_reference="emission", sunward_cone_half_angle_deg=90.0,
                require_projected_sunward=False,
                phase_law="schleicher", phase_linear_beta=0.024,
                phase_linear_m_oe=0.28, phase_linear_w_oe=1.5,
                size_dist_table=None, p_v=0.04, grid_npix=300):
        super().__init__()
        self.comet_el, self.obs_jd = comet_el, obs_jd
        self.beta_range, self.gamma_size = beta_range, gamma_size
        self.max_age, self.n_particles = max_age, n_particles
        self.v0_coeff, self.gamma, self.m_exp = v0_coeff, gamma, m_exp
        self.seed = seed
        self.north_pa_deg, self.rotation_offset_deg = north_pa_deg, rotation_offset_deg
        self.qt_weights = qt_weights
        self.active_area = active_area
        self.rho_g_cm3 = rho_g_cm3
        self.sunward = sunward
        self.sunward_expocos = sunward_expocos
        self.sunward_reference = sunward_reference
        self.sunward_cone_half_angle_deg = sunward_cone_half_angle_deg
        self.require_projected_sunward = require_projected_sunward
        self.phase_law = phase_law
        self.phase_linear_beta = phase_linear_beta
        self.phase_linear_m_oe = phase_linear_m_oe
        self.phase_linear_w_oe = phase_linear_w_oe
        self.size_dist_table = size_dist_table
        self.p_v = p_v
        self.grid_npix = grid_npix

    def _progress_callback(self, fraction, message):
        try:
            pct = int(round(float(fraction) * 100.0))
        except Exception:
            pct = 0
        pct = max(0, min(100, pct))
        self.progress.emit(pct, str(message))

    def run(self):
        try:
            result = cta.compute_morphology_mc(
                self.comet_el, self.obs_jd, self.beta_range, self.gamma_size,
                self.max_age, self.n_particles,
                v0_coeff=self.v0_coeff, gamma=self.gamma, m_exp=self.m_exp,
                seed=self.seed, north_pa_deg=self.north_pa_deg,
                rotation_offset_deg=self.rotation_offset_deg,
                qt_weights=self.qt_weights, active_area=self.active_area,
                rho_g_cm3=self.rho_g_cm3,
                sunward=self.sunward, sunward_expocos=self.sunward_expocos,
                sunward_reference=self.sunward_reference,
                sunward_cone_half_angle_deg=self.sunward_cone_half_angle_deg,
                require_projected_sunward=self.require_projected_sunward,
                phase_law=self.phase_law,
                phase_linear_beta=self.phase_linear_beta,
                phase_linear_m_oe=self.phase_linear_m_oe,
                phase_linear_w_oe=self.phase_linear_w_oe,
                size_dist_table=self.size_dist_table, p_v=self.p_v,
                grid_npix=self.grid_npix,
                progress_callback=self._progress_callback)
            self.finished.emit(result)
        except Exception as ex:
            self.error.emit(str(ex))

# ─────────────────────────────────────────────────────────────────────────────
#  MATPLOTLIB CANVAS WIDGET
# ─────────────────────────────────────────────────────────────────────────────
class PlotCanvas(QWidget):
    nucleus_clicked = pyqtSignal(float, float)   # px, py on image

    def __init__(self, parent=None):
        super().__init__(parent)
        self.fig  = Figure(facecolor=BG)
        self.ax   = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setIconSize(QSize(18, 18))

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)
        vbox.addWidget(self.toolbar)
        vbox.addWidget(self.canvas)

        self._picking_nucleus = False
        self._model  = None
        self._imgArr = None
        self.nuc_x   = 0.0
        self.nuc_y   = 0.0
        self.au_per_px = 0.001
        self.north_pa  = 0.0
        self._vis = dict(synd=True, sync=True, orbit=True, isophote=False,
                         isophote_levels=6, isophote_smooth=2.0, lw=1.5, alpha=0.85,
                         orbit_lw=1.5, orbit_alpha=1.0, mc_contour=False)
        # Phase 2 (v3.1) — paths from extract_morphology_contours(), pushed
        # in by MainWindow.set_mc_contours() (called from the MC popup, via
        # an explicit reference — NOT Qt parent/child, since that window
        # deliberately uses parent=None; see MCWindow's own comment on why).
        # Each path is an (N,2) AU [ξ,η] array, drawn through the SAME
        # to_px()/prep() pipeline as syndyne/synchrone curves below — one
        # rendering path, so it can't develop its own separate orientation
        # bug the way the standalone popup's plot did.
        self.mc_contours = []
        self.mc_info = None   # v3.1 — info dict for the optional parameter box

        self._draw_empty()

    # ── Nucleus picking ───────────────────────────────────────────────────
    def set_picking_nucleus(self, val: bool):
        self._picking_nucleus = val
        self.canvas.setCursor(
            QCursor(Qt.CursorShape.CrossCursor) if val else
            QCursor(Qt.CursorShape.ArrowCursor))

    def _on_click(self, event):
        if self._picking_nucleus and event.inaxes and event.button == 1:
            self.nuc_x = event.xdata
            self.nuc_y = event.ydata
            self._picking_nucleus = False
            self.canvas.setCursor(QCursor(Qt.CursorShape.ArrowCursor))
            self.nucleus_clicked.emit(self.nuc_x, self.nuc_y)
            if self._model:
                self.draw_model(self._model, self._imgArr)

    # ── Empty state ───────────────────────────────────────────────────────
    def _draw_empty(self):
        self.ax.clear()
        self.ax.set_facecolor(BG)
        self.fig.patch.set_facecolor(BG)
        for sp in self.ax.spines.values():
            sp.set_edgecolor("#1a2540")
        self.ax.tick_params(colors="#1a2540")
        self.ax.set_xlabel("Δ East (AU)", color="#1a2540", fontfamily=MF)
        self.ax.set_ylabel("Δ North (AU)", color="#1a2540", fontfamily=MF)
        self.ax.text(0.5, 0.5, "☄", transform=self.ax.transAxes,
                     ha="center", va="center", fontsize=48, color="#0e1c38", alpha=0.6)
        self.ax.text(0.5, 0.38, "SELECT A COMET AND CLICK  COMPUTE",
                     transform=self.ax.transAxes, ha="center", va="center",
                     fontsize=10, color="#0e1c38", fontfamily=MF, alpha=0.8)
        self.canvas.draw_idle()

    # ── Image preview (before model computed) ────────────────────────────
    def draw_image_preview(self, img_arr: np.ndarray, filename: str = "",
                           wcs_info: dict = None):
        """Show image immediately on open, with compass and WCS info overlay."""
        prev_xlim = prev_ylim = None
        if self._imgArr is not None:
            prev_xlim = self.ax.get_xlim()
            prev_ylim = self.ax.get_ylim()

        self._imgArr = img_arr
        self.ax.clear()
        self.fig.patch.set_facecolor("black")
        self.ax.set_facecolor("black")
        for sp in self.ax.spines.values():
            sp.set_edgecolor("#1a2540")

        h, w = img_arr.shape[:2]
        self.ax.imshow(img_arr, origin="upper", aspect="equal",
                       extent=[0, w, h, 0], zorder=0)

        # Preserve observational isophotes when model overlays are cleared.
        # Isophotes are derived from the image, not from either the F-P or
        # Monte Carlo model, so CLEAR ALL MODELS should leave them visible.
        if self._vis.get("isophote", False):
            lum, lvl = _compute_isophote_levels(
                img_arr,
                n_levels=self._vis.get("isophote_levels", 6),
                smooth_sigma=self._vis.get("isophote_smooth", 2.0))
            if lum is not None and lvl is not None and len(lvl) >= 2:
                try:
                    self.ax.contour(
                        lum, levels=lvl, origin="upper",
                        extent=[0, w, h, 0], colors=ISOPHOTE_COLOR,
                        linewidths=0.6, alpha=0.75, zorder=1)
                except Exception:
                    pass
        if prev_xlim is not None and prev_ylim is not None:
            self.ax.set_xlim(*prev_xlim)
            self.ax.set_ylim(*prev_ylim)
        else:
            self.ax.set_xlim(0, w)
            self.ax.set_ylim(h, 0)
        self.ax.tick_params(labelcolor="#3a5070", labelsize=8, color="#1a2540")
        self.ax.set_xlabel("X (pixels)", color="#3a6080", fontsize=9, fontfamily=MF)
        self.ax.set_ylabel("Y (pixels)", color="#3a6080", fontsize=9, fontfamily=MF)

        # ── Nucleus crosshair (if set) ────────────────────────────────────
        if self.nuc_x > 0 or self.nuc_y > 0:
            cs = max(h, w) * 0.012
            self.ax.plot([self.nuc_x - cs, self.nuc_x + cs],
                         [self.nuc_y, self.nuc_y], "w-", lw=1.5, zorder=6)
            self.ax.plot([self.nuc_x, self.nuc_x],
                         [self.nuc_y - cs, self.nuc_y + cs], "w-", lw=1.5, zorder=6)
            self.ax.plot(self.nuc_x, self.nuc_y, "wo", ms=4, zorder=7)

        # ── N/E Compass ───────────────────────────────────────────────────
        npa = np.radians(self.north_pa)
        cx, cy = w * 0.07, h * 0.07
        aL = min(w, h) * 0.055
        # North arrow
        n_dx =  np.sin(npa) * aL;  n_dy = -np.cos(npa) * aL
        self.ax.annotate("", xy=(cx + n_dx, cy + n_dy), xytext=(cx, cy),
                         arrowprops=dict(arrowstyle="->", color="#60c8ff", lw=2.2),
                         zorder=8)
        self.ax.text(cx + n_dx * 1.32, cy + n_dy * 1.32, "N",
                     color="#60c8ff", fontsize=10, ha="center", va="center",
                     fontfamily=MF, fontweight="bold", zorder=9)
        # East arrow (90° CW from North on screen)
        e_dx = -np.cos(npa) * aL;  e_dy = -np.sin(npa) * aL
        self.ax.annotate("", xy=(cx + e_dx, cy + e_dy), xytext=(cx, cy),
                         arrowprops=dict(arrowstyle="->", color="#60c8ff", lw=1.8),
                         zorder=8)
        self.ax.text(cx + e_dx * 1.32, cy + e_dy * 1.32, "E",
                     color="#60c8ff", fontsize=10, ha="center", va="center",
                     fontfamily=MF, fontweight="bold", zorder=9)

        # ── WCS / file info (top-right) ───────────────────────────────────
        lines = [filename]
        if wcs_info:
            lines += [
                f"PS = {wcs_info.get('ps_arcsec',0):.3f}\"/px",
                f"NPA = {wcs_info.get('npa',0):.2f}°",
                f"Nucleus = ({self.nuc_x:.0f}, {self.nuc_y:.0f}) px",
                f"Size = {w}×{h} px",
            ]
            if wcs_info.get("object"):
                lines.insert(1, wcs_info["object"])
            if wcs_info.get("date_obs"):
                lines.insert(2, wcs_info["date_obs"][:19])
        for j, line in enumerate(lines):
            if not line: continue
            self.ax.text(w - 8, 10 + j * 18, line,
                         color="white", fontsize=8, fontfamily=MF,
                         ha="right", va="top",
                         bbox=dict(boxstyle="round,pad=0.2", facecolor="black",
                                   alpha=0.55, edgecolor="none"),
                         zorder=9)

        # ── Stretch hint ──────────────────────────────────────────────────
        self.ax.text(0.5, 0.01,
                     "IMAGE LOADED  ·  Adjust stretch in left panel  ·  "
                     "Then COMPUTE MODEL",
                     transform=self.ax.transAxes, ha="center", va="bottom",
                     fontsize=9, color="#3a7080", fontfamily=MF, alpha=0.9)

        self.canvas.draw_idle()

    # ── Publication projected-distance axes ─────────────────────────────
    @staticmethod
    def _nice_projected_tick_values(vmin: float, vmax: float,
                                    target_intervals: int = 5) -> list[float]:
        """Return readable kilometre ticks and include zero when visible.

        The plot data stay in image-pixel coordinates so the image, WCS, and
        MC/F-P overlay geometry are not modified.  This helper instead builds
        a fixed set of physically labelled ticks around the nucleus.  Using a
        nucleus-centred locator is important because Matplotlib's default
        pixel locator does not guarantee that the nucleus pixel will be a
        tick, even though its projected coordinate must be exactly 0 km.
        """
        lo, hi = sorted((float(vmin), float(vmax)))
        span = hi - lo
        if not np.isfinite(span) or span <= 0:
            return [0.0]

        raw = span / max(int(target_intervals), 1)
        power = 10.0 ** np.floor(np.log10(max(raw, 1e-12)))
        fraction = raw / power
        if fraction <= 1.0:
            factor = 1.0
        elif fraction <= 2.0:
            factor = 2.0
        elif fraction <= 2.5:
            factor = 2.5
        elif fraction <= 5.0:
            factor = 5.0
        else:
            factor = 10.0
        step = factor * power

        start = np.ceil((lo - 1e-12 * step) / step) * step
        stop  = np.floor((hi + 1e-12 * step) / step) * step
        if stop >= start:
            count = int(round((stop - start) / step)) + 1
            ticks = [float(start + i * step) for i in range(max(count, 0))]
        else:
            ticks = []

        # The nucleus coordinate is the physical origin.  Always include it
        # when it falls in the visible crop, even if the automatic nice-step
        # sequence would otherwise skip it.
        if lo <= 0.0 <= hi and not any(abs(v) <= step * 1e-8 for v in ticks):
            ticks.append(0.0)
        ticks = sorted(v for v in ticks if lo - step*1e-8 <= v <= hi + step*1e-8)
        return ticks or ([0.0] if lo <= 0.0 <= hi else [lo, hi])

    @staticmethod
    def _format_projected_km_tick(value: float) -> str:
        """Format a projected-distance tick using compact ×10 superscripts."""
        value = float(value)
        if not np.isfinite(value) or abs(value) < 0.5:
            return "0"

        av = abs(value)
        if av >= 1.0e4:
            exp = int(np.floor(np.log10(av)))
            coeff = value / (10.0 ** exp)
            if abs(coeff - round(coeff)) < 1e-9:
                coeff_text = str(int(round(coeff)))
            else:
                coeff_text = f"{coeff:.1f}".rstrip("0").rstrip(".")
            coeff_text = coeff_text.replace("-", "−")
            superscript = str(exp).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
            return f"{coeff_text}×10{superscript}"

        if av >= 100.0:
            text = f"{value:.0f}"
        elif av >= 10.0:
            text = f"{value:.1f}".rstrip("0").rstrip(".")
        else:
            text = f"{value:.2f}".rstrip("0").rstrip(".")
        return text.replace("-", "−")

    def _apply_projected_km_axes(self, white_canvas: bool) -> None:
        """Apply nucleus-centred projected-kilometre ticks to an image plot.

        Sign convention follows the requested publication layout:
          X < 0 left of the nucleus; X > 0 right of the nucleus.
          Y > 0 above the nucleus; Y < 0 below the nucleus.
        """
        km_per_px = float(self.au_per_px) * 149597870.7
        if not np.isfinite(km_per_px) or km_per_px <= 0:
            return

        x0_px, x1_px = self.ax.get_xlim()
        y0_px, y1_px = self.ax.get_ylim()
        x0_km = (x0_px - self.nuc_x) * km_per_px
        x1_km = (x1_px - self.nuc_x) * km_per_px
        y0_km = (self.nuc_y - y0_px) * km_per_px
        y1_km = (self.nuc_y - y1_px) * km_per_px

        x_ticks_km = self._nice_projected_tick_values(x0_km, x1_km)
        y_ticks_km = self._nice_projected_tick_values(y0_km, y1_km)

        x_pairs = sorted(
            ((self.nuc_x + value / km_per_px, value) for value in x_ticks_km),
            key=lambda item: item[0])
        y_pairs = sorted(
            ((self.nuc_y - value / km_per_px, value) for value in y_ticks_km),
            key=lambda item: item[0])

        self.ax.set_xticks([p for p, _ in x_pairs])
        self.ax.set_xticklabels(
            [self._format_projected_km_tick(v) for _, v in x_pairs],
            fontfamily=MF, fontsize=8)
        self.ax.set_yticks([p for p, _ in y_pairs])
        self.ax.set_yticklabels(
            [self._format_projected_km_tick(v) for _, v in y_pairs],
            fontfamily=MF, fontsize=8)

        axis_color = "black" if white_canvas else "#3a6080"
        self.ax.set_xlabel("Distance to nucleus projected (km)",
                           color=axis_color, fontfamily=MF)
        self.ax.set_ylabel("Distance to nucleus projected (km)",
                           color=axis_color, fontfamily=MF)

    # ── Main draw ─────────────────────────────────────────────────────────
    def draw_model(self, model, img_arr=None, fixed_xlim=None, fixed_ylim=None):
        # Capture the CURRENT view — whatever a previous draw_model() call
        # set, or whatever the user's toolbar pan/zoom has since changed
        # it to — directly from the axes, BEFORE ax.clear() wipes it. This
        # is what makes the overlay view persist across Compute Model
        # re-runs. (An earlier attempt used matplotlib's xlim_changed/
        # ylim_changed callbacks instead; ax.clear() silently resets that
        # registry on every redraw, so the callback kept going dead with
        # no visible error. Reading the axes' own current state directly
        # has no such failure mode — it's always accurate by definition.)
        prev_xlim = prev_ylim = None
        if self._imgArr is not None:
            prev_xlim = self.ax.get_xlim()
            prev_ylim = self.ax.get_ylim()

        self._model  = model
        self._imgArr = img_arr
        self.ax.clear()
        overlay = img_arr is not None
        view_style = self._vis.get("mc_view_style", "analysis")
        # Contour-only presentation is meaningful only when an MC result is
        # actually present.  After CLEAR ALL MODELS, a new F-P run must return
        # to a normal diagnostic display rather than remain silently hidden by
        # the last MC publication preset.
        contour_only = (view_style in ("contour", "publication") and
                        bool(self.mc_contours) and
                        bool(self._vis.get("mc_contour", False)))
        background_mode = self._vis.get(
            "mc_background", "white" if contour_only else "image")
        white_canvas = background_mode == "white"
        grayscale_safe = bool(self._vis.get("grayscale_safe", False))
        observed_color = self._vis.get(
            "observed_color", "black" if contour_only else ISOPHOTE_COLOR)
        observed_ls = self._vis.get("observed_ls", "-")
        model_color = self._vis.get("model_color", "#ff3399")
        model_ls = self._vis.get("model_ls", "--")
        if grayscale_safe:
            observed_color, observed_ls = "black", "-"
            model_color, model_ls = "black", "--"

        if white_canvas:
            self.ax.set_facecolor("white")
            self.fig.patch.set_facecolor("white")
            for sp in self.ax.spines.values():
                sp.set_edgecolor("black")
            self.ax.tick_params(labelcolor="black", labelsize=8, color="black",
                                direction="in", top=True, right=True)
            self.ax.grid(False)
        else:
            self.ax.set_facecolor(BG)
            self.fig.patch.set_facecolor(BG)
            for sp in self.ax.spines.values():
                sp.set_edgecolor("#1a2540")
            self.ax.tick_params(labelcolor="#2a4060", labelsize=8,
                                color="#1a2540")
            self.ax.grid(color="#0d1a2e", lw=0.4, zorder=0)

        # ── AU → degree conversion factor (standalone only) ──────────────
        r_geo = model["info"].get("r_geo", 1.0)
        K = (180.0 / np.pi) / r_geo   # deg per AU  (small-angle, accurate for <5°)

        if overlay:
            # Draw the raster image only in analysis mode.  Contour-only
            # modes still use img_arr as the source of the observed
            # isophotes and retain its pixel/WCS geometry, but deliberately
            # suppress the background so the comparison can be exported as
            # a clean paper figure.
            h_img, w_img = img_arr.shape[:2]
            if background_mode == "image":
                self.ax.imshow(img_arr, origin="upper", aspect="equal",
                               extent=[0, w_img, h_img, 0], zorder=0)
            if prev_xlim is not None and prev_ylim is not None:
                self.ax.set_xlim(*prev_xlim)
                self.ax.set_ylim(*prev_ylim)
                view_w_px = abs(prev_xlim[1] - prev_xlim[0])
            else:
                self.ax.set_xlim(0, w_img)
                self.ax.set_ylim(h_img, 0)
                view_w_px = w_img
            axis_color = "black" if white_canvas else "#3a6080"
            self.ax.set_xlabel("X (pixels)", color=axis_color, fontfamily=MF)
            self.ax.set_ylabel("Y (pixels)", color=axis_color, fontfamily=MF)

            # ── Isophote overlay (v3.0) — traced from the loaded image ────
            show_observed = bool(self._vis.get("observed_show", True)) and (
                self._vis.get("isophote", False) or contour_only)
            if show_observed:
                lum, lvl = _compute_isophote_levels(
                    img_arr,
                    n_levels=self._vis.get("isophote_levels", 6),
                    smooth_sigma=self._vis.get("isophote_smooth", 2.0))
                if lum is not None and lvl is not None and len(lvl) >= 2:
                    try:
                        self.ax.contour(lum, levels=lvl, origin="upper",
                                        extent=[0, w_img, h_img, 0],
                                        colors=observed_color,
                                        linewidths=self._vis.get("observed_lw", 0.6),
                                        alpha=self._vis.get("observed_alpha", 0.75),
                                        linestyles=observed_ls, zorder=2)
                    except Exception:
                        pass   # malformed/flat image — skip silently, never block the plot

            def to_px(xi, eta):
                return cta.sky_to_pixel(xi, eta, self.nuc_x, self.nuc_y,
                                        self.au_per_px, self.north_pa)

            def prep(xi_arr, eta_arr):
                """
                Convert (xi, eta) AU → pixel, with limited NaN interpolation.
                Only interpolate short NaN gaps (≤10 consecutive points) to avoid
                unrealistic straight lines across large missing data regions.
                """
                xo, yo = [], []
                for xi, eta in zip(xi_arr, eta_arr):
                    if np.isfinite(xi) and np.isfinite(eta):
                        px, py = to_px(xi, eta)
                        xo.append(px); yo.append(py)
                    else:
                        xo.append(np.nan); yo.append(np.nan)
                
                xo, yo = np.array(xo, dtype=float), np.array(yo, dtype=float)
                
                # Interpolate only SHORT NaN gaps to preserve curve shape
                mask = np.isfinite(xo) & np.isfinite(yo)
                if not mask.any():
                    return xo, yo  # All NaN, nothing to do
                
                # Find NaN runs (consecutive NaN segments)
                diff = np.diff(np.concatenate(([True], mask, [True])).astype(int))
                gap_starts = np.where(diff == -1)[0]
                gap_ends = np.where(diff == 1)[0]
                
                # Interpolate gaps shorter than max_gap_size
                max_gap_size = 10  # Don't interpolate gaps > 10 points
                indices = np.arange(len(xo))
                for start, end in zip(gap_starts, gap_ends):
                    gap_len = end - start
                    if 0 < gap_len <= max_gap_size:
                        # Short gap: interpolate
                        gap_idx = indices[start:end]
                        valid_before = indices[:start][mask[:start]]
                        valid_after = indices[end:][mask[end:]]
                        if len(valid_before) > 0 and len(valid_after) > 0:
                            # Have valid points on both sides: interpolate
                            valid_idx = np.concatenate([valid_before, valid_after])
                            valid_x = np.concatenate([xo[valid_before], xo[valid_after]])
                            valid_y = np.concatenate([yo[valid_before], yo[valid_after]])
                            xo[gap_idx] = np.interp(gap_idx, valid_idx, valid_x)
                            yo[gap_idx] = np.interp(gap_idx, valid_idx, valid_y)
                    # else: Long gap → leave as NaN → matplotlib will break the line
                
                return xo, yo

            def clip_to_image(xv, yv, w, h, margin=20):
                """
                Replace points outside image bounds with NaN so matplotlib
                doesn't draw those segments. Adds a break (NaN) at exit points.
                """
                out_x = np.array(xv, dtype=float)
                out_y = np.array(yv, dtype=float)
                outside = ((out_x < -margin) | (out_x > w + margin) |
                           (out_y < -margin) | (out_y > h + margin))
                out_x[outside] = np.nan
                out_y[outside] = np.nan
                return out_x, out_y
        else:
            # ── Convert all AU offsets to angular degrees ─────────────────
            all_xi  = [0.0]; all_eta = [0.0]
            for s in model["syndynes"]:
                all_xi  += list(s["xi"] [np.isfinite(s["xi"])]  * K)
                all_eta += list(s["eta"][np.isfinite(s["eta"])] * K)
            for s in model["synchrones"]:
                all_xi  += list(s["xi"] [np.isfinite(s["xi"])]  * K)
                all_eta += list(s["eta"][np.isfinite(s["eta"])] * K)
            if fixed_xlim is not None and fixed_ylim is not None:
                # v3.0 Animator: lock the frame to a caller-specified size
                # (Fixed FOV in degrees, or Fixed distance in AU passed
                # through the same K-scaled xi/eta units) instead of
                # auto-fitting per frame — auto-fit would zoom out as the
                # tail grows, hiding the very change the animation exists
                # to show.
                self.ax.set_xlim(*fixed_xlim)
                self.ax.set_ylim(*fixed_ylim)
            elif len(all_xi) > 1:
                xr  = np.percentile(all_xi,  [2, 98])
                yr  = np.percentile(all_eta, [2, 98])
                pad = 0.12
                dx  = max((xr[1]-xr[0])*pad, 0.01)
                dy  = max((yr[1]-yr[0])*pad, 0.01)
                # East LEFT = inverted x-axis (E is positive xi, appears left)
                self.ax.set_xlim(xr[1]+dx, xr[0]-dx)
                self.ax.set_ylim(yr[0]-dy, yr[1]+dy)
            # Equal aspect so 1° in RA-direction always renders as the same
            # screen length as 1° in Dec-direction (see CHANGELOG v3.0 for
            # the bug this fixes: without it, the main window and Animator
            # rendered the same model at different distortions because
            # their canvas widgets have different pixel aspect ratios).
            #
            # adjustable mode differs by use case:
            #   - Animator (fixed_xlim/fixed_ylim given): 'box' — the whole
            #     point of "Fixed FOV"/"Fixed distance" is an EXACT
            #     requested angular size, so any aspect mismatch is taken
            #     up as blank padding, never by silently widening the FOV.
            #   - Main window auto-fit (no fixed_xlim/fixed_ylim): 'datalim'
            #     — there's no precision promise to keep here, so instead
            #     of padding with empty space, extend whichever axis range
            #     is needed to fill the canvas completely. Far less
            #     confusing than a tiny plot adrift in a sea of black with
            #     the legend now looming huge over it.
            adj = 'box' if (fixed_xlim is not None and fixed_ylim is not None) else 'datalim'
            self.ax.set_aspect('equal', adjustable=adj)
            # set_aspect() with adjustable='datalim' only computes the
            # final, box-filling data limits LAZILY at actual draw time —
            # so anything below that calls self.ax.get_xlim()/get_ylim()
            # right after (the N/E compass position, the starfield, label
            # placement) would otherwise read the stale PRE-adjustment
            # limits, not the final ones. apply_aspect() forces that
            # computation to happen now instead of later.
            self.ax.apply_aspect()
            self.ax.axhline(0, color="#1a2a40", lw=0.7, zorder=1)
            self.ax.axvline(0, color="#1a2a40", lw=0.7, zorder=1)
            self.ax.set_xlabel("← East  ·  Δ (°)  ·  West →",
                               color="#3a6080", fontsize=9, fontfamily=MF)
            self.ax.set_ylabel("↑  Δ North (°)",
                               color="#3a6080", fontsize=9, fontfamily=MF)

            def prep(xi_arr, eta_arr):
                # Convert AU → degrees for standalone display
                return xi_arr * K, eta_arr * K

        lw  = self._vis["lw"]
        alp = self._vis["alpha"]

        # ── Helper: find best label position along a curve ────────────────
        def _best_label_pos(xv, yv, xl_lim, yl_lim, prefer_end=True):
            """
            Find a visible point on the curve for label placement.
            Strategy: find the LAST visible point (near where curve exits frame).
            Always returns a position — uses boundary clamp if nothing is ideal.
            """
            mask = np.isfinite(xv) & np.isfinite(yv)
            if not mask.any():
                return 0, 0, False
            xs, ys = xv[mask], yv[mask]
            n = len(xs)
            if n == 0:
                return 0, 0, False

            x0, x1 = min(xl_lim), max(xl_lim)
            y0, y1 = min(yl_lim), max(yl_lim)
            # Small padding so label doesn't touch exact edge
            px_pad = (x1 - x0) * 0.02
            py_pad = abs(y1 - y0) * 0.02

            def _inside(px, py):
                return (x0 + px_pad < px < x1 - px_pad and
                        y0 + py_pad < py < y1 - py_pad)

            # Collect all inside indices
            inside_idx = [i for i in range(n) if _inside(xs[i], ys[i])]

            if not inside_idx:
                # Curve entirely outside — place at closest point to plot centre
                cx, cy = (x0+x1)/2, (y0+y1)/2
                dist = [(xs[i]-cx)**2 + (ys[i]-cy)**2 for i in range(n)]
                best = int(np.argmin(dist))
                return xs[best], ys[best], False

            # prefer_end=True → label at far/exit end  (syndynes)
            # prefer_end=False → label at near end     (synchrones)
            if prefer_end:
                best_i = inside_idx[-1]          # last visible
            else:
                # Pick ~55% along the inside segment
                mid = inside_idx[len(inside_idx) * 55 // 100]
                best_i = mid

            return xs[best_i], ys[best_i], True

        # Orbital path
        if (not contour_only) and self._vis["orbit"] and model["orbit"]:
            op = np.array(model["orbit"])
            if overlay:
                xs, ys = zip(*[to_px(p[0], p[1]) for p in op])
            else:
                xs = op[:,0] * K
                ys = op[:,1] * K
            self.ax.plot(xs, ys, "--", color="#2050a0",
                         lw=self._vis.get("orbit_lw", 1.0),
                         alpha=self._vis.get("orbit_alpha", 0.7),
                         zorder=2, label="Orbital path")

        def _truncate_at_turnaround(xv, yv, nuc_x, nuc_y):
            """
            Keep only the portion of the curve that moves AWAY from nucleus.
            Truncate at the first point where cumulative distance starts
            decreasing (curve loops back). Preserves the physically meaningful
            part and removes the spurious 'return' arc.
            """
            n = len(xv)
            if n < 3:
                return xv, yv
            dists = np.sqrt((xv - nuc_x)**2 + (yv - nuc_y)**2)
            # Find the global maximum distance index
            peak_i = int(np.argmax(dists))
            # Only keep up to the peak
            return xv[:peak_i+1], yv[:peak_i+1]

        def _truncate_syndyne(xv, yv, nuc_x, nuc_y):
            """
            For syndynes: keep only the outward portion (0 → peak distance).
            Removes loop-backs where curve returns toward nucleus.
            """
            n = len(xv)
            if n < 3: return xv, yv
            mask = np.isfinite(xv) & np.isfinite(yv)
            if not mask.any(): return xv, yv
            dists = np.where(mask, np.sqrt((xv-nuc_x)**2+(yv-nuc_y)**2), 0.0)
            peak_i = int(np.argmax(dists))
            return xv[:peak_i+1], yv[:peak_i+1]

        def _truncate_synchrone(xv, yv, nuc_x, nuc_y):
            """
            For synchrones: the curve spans from β_small (old comet pos, sunward)
            to β_large (anti-solar, tail direction).
            Keep only the anti-solar half: from the point CLOSEST to nucleus
            outward toward the tail.
            Physics: β_small → dust stayed where comet WAS (sunward of now)
                     β_large → dust pushed anti-solar (tail, we want this)
            """
            n = len(xv)
            if n < 3: return xv, yv
            mask = np.isfinite(xv) & np.isfinite(yv)
            if not mask.any(): return xv, yv
            dists = np.where(mask, np.sqrt((xv-nuc_x)**2+(yv-nuc_y)**2), 1e9)
            min_i = int(np.argmin(dists))
            # Keep from the minimum-distance point onward (toward tail)
            # The tail end has larger β = earlier indices in the array
            # (β is log-spaced from small→large in compute_model)
            # min_i divides sunward portion (after min_i) from tail portion (before)
            if min_i <= n // 2:
                # min is near start → tail is after min_i
                return xv[min_i:], yv[min_i:]
            else:
                # min is near end → tail is before min_i
                return xv[:min_i+1], yv[:min_i+1]

        xl_lim = self.ax.get_xlim()
        yl_lim = self.ax.get_ylim()

        # Nucleus position for label-distance calculation
        if overlay:
            nuc_x = float(self.nuc_x) if self.nuc_x else 0.0
            nuc_y = float(self.nuc_y) if self.nuc_y else 0.0
        else:
            nuc_x, nuc_y = 0.0, 0.0

        _lbl_kw = dict(fontfamily=MF, fontsize=7.5, va="center",
                       fontweight="bold", zorder=6, clip_on=True)

        def _label_far_from_nucleus(xv, yv, xl_lim, yl_lim, nuc_x, nuc_y,
                                    min_dist_frac=0.12):
            """
            Find best label position: inside frame AND far from nucleus.
            Returns (x, y, ok). ok=False means curve too short to label.
            """
            mask = np.isfinite(xv) & np.isfinite(yv)
            if not mask.any(): return 0, 0, False
            xs, ys = xv[mask], yv[mask]
            n = len(xs)
            if n < 2: return 0, 0, False

            x0, x1 = min(xl_lim), max(xl_lim)
            y0, y1 = min(yl_lim), max(yl_lim)
            frame_w = x1 - x0
            frame_h = abs(y1 - y0)
            min_dist = min_dist_frac * max(frame_w, frame_h)
            px_pad = frame_w * 0.02
            py_pad = frame_h * 0.02

            def _inside(px, py):
                return (x0+px_pad < px < x1-px_pad and
                        min(y0,y1)+py_pad < py < max(y0,y1)-py_pad)

            def _dist(px, py):
                return np.sqrt((px-nuc_x)**2 + (py-nuc_y)**2)

            # All inside points with distance from nucleus
            inside_pts = [(i, _dist(xs[i], ys[i]))
                          for i in range(n) if _inside(xs[i], ys[i])]
            if not inside_pts:
                return 0, 0, False

            # Far enough from nucleus
            far_pts = [(i, d) for i, d in inside_pts if d >= min_dist]

            if far_pts:
                best_i = max(far_pts, key=lambda x: x[1])[0]
                return xs[best_i], ys[best_i], True

            # Relax: farthest inside point, but require at least half min_dist
            best_i, best_d = max(inside_pts, key=lambda x: x[1])
            return xs[best_i], ys[best_i], best_d >= min_dist * 0.5

        # Syndynes — draw with proper clipping after NaN interpolation
        # prep() interpolates NaN gaps, then matplotlib clips to axes bounds
        if (not contour_only) and self._vis["synd"]:
            for idx, synd in enumerate(model["syndynes"]):
                col  = SYNDYNE_COLORS[idx % len(SYNDYNE_COLORS)]
                xv, yv = prep(synd["xi"], synd["eta"])
                self.ax.plot(xv, yv, color=col, lw=lw, alpha=alp,
                             zorder=3, clip_on=True)  # Enable clip after interpolation
                lx, ly, ok = _label_far_from_nucleus(
                    xv, yv, xl_lim, yl_lim, nuc_x, nuc_y, min_dist_frac=0.08)
                if ok:
                    self.ax.text(lx, ly, f"β={synd['beta']}",
                                 color=col, **_lbl_kw)

        # Synchrones — draw with clipping after interpolation
        if (not contour_only) and self._vis["sync"]:
            for idx, sync in enumerate(model["synchrones"]):
                col  = SYNC_COLORS[idx % len(SYNC_COLORS)]
                xv, yv = prep(sync["xi"], sync["eta"])
                # Draw the curve with proper clipping
                self.ax.plot(xv, yv, color=col, lw=lw*0.85,
                             alpha=alp, ls="-.", zorder=3, clip_on=True)  # Enable clip
                # Add label only when a suitable far point exists
                lx, ly, ok = _label_far_from_nucleus(
                    xv, yv, xl_lim, yl_lim, nuc_x, nuc_y, min_dist_frac=0.12)
                if ok:
                    self.ax.text(lx, ly, f"t−{sync['age']}d",
                                 color=col, **_lbl_kw)

        # MC morphology contour (v3.1, Phase 2) — same to_px()/prep()
        # pipeline as syndynes/synchrones above, so it shares their
        # already-correct orientation handling rather than needing its own.
        # Own width/opacity sliders (not tied to the syndyne ones). Plain
        # dashed line — an earlier version added a black+white double
        # outline (path_effects) for visibility on any background, but with
        # several nested contour levels overlapping it read as a blurry
        # band rather than a clean line, so that's been dropped in favour
        # of this simpler look.
        if self._vis.get("mc_contour") and self.mc_contours:
            mc_lw    = self._vis.get("mc_lw", 1.2)
            mc_alpha = self._vis.get("mc_alpha", 0.95)
            for path in self.mc_contours:
                xv, yv = prep(path[:, 0], path[:, 1])
                self.ax.plot(xv, yv, color=model_color, lw=mc_lw,
                            alpha=mc_alpha, ls=model_ls, zorder=3, clip_on=True)

            # Optional parameter summary box (v3.1) — off by default (see
            # the DISPLAY checkbox's tooltip for why); bottom-right corner,
            # same visual style as the WCS info box (top-right, see
            # draw_image_preview()) but on the opposite corner so the two
            # never overlap.
            if self._vis.get("mc_info_box") and self.mc_info:
                info  = self.mc_info
                b_lo, b_hi = info['beta_range']
                rho_used   = info.get('rho_g_cm3', 0.5)
                single_size = (b_lo == b_hi)

                def _fmt_um(a):
                    return f"{a:.0f}" if a >= 10 else f"{a:.2g}"

                if single_size:
                    a_um = cta._beta_to_radius_um(b_lo, rho_used)
                    grain_str = f"{_fmt_um(a_um)} µm"
                else:
                    a_max = cta._beta_to_radius_um(b_lo, rho_used)
                    a_min = cta._beta_to_radius_um(b_hi, rho_used)
                    grain_str = f"{_fmt_um(a_min)}–{_fmt_um(a_max)} µm"

                speed_str = ""
                if info['v0_coeff'] != 0:
                    r_h = info.get('r_helio_au')
                    if r_h:
                        v_hi_raw = cta.real_ejection_speed_ms(
                            info['v0_coeff'], b_hi, r_h, info['gamma'], info['m_exp'])
                        v_lo_raw = cta.real_ejection_speed_ms(
                            info['v0_coeff'], b_lo, r_h, info['gamma'], info['m_exp'])
                        v_lo_s, v_hi_s = sorted([float(v_lo_raw), float(v_hi_raw)])
                        speed_str = (f"{v_lo_s:.1f} m/s" if single_size
                                     else f"{v_lo_s:.1f}–{v_hi_s:.1f} m/s")

                # Ejection mode label
                if info.get('active_area_used'):
                    mode_str = "Active area"
                elif info.get('sunward_used'):
                    cone = float(info.get('sunward_cone_half_angle_deg', 90.0))
                    ref  = info.get('sunward_reference', 'emission')
                    mode_str = ("Sunward hemisphere" if cone >= 89.999
                                else f"Sunward cone {cone:.0f}°")
                    if ref == 'observation':
                        mode_str += " (legacy ref)"
                else:
                    mode_str = "Isotropic"

                lines = ["MONTE CARLO MODEL"]
                lines.append(f"Grain size: {grain_str}")
                if not single_size:
                    lines.append(f"Size index: {info['gamma_size']:.2g}")
                lines.append(f"Particles: {info['n_used']:,}")
                lines.append(f"Dust release window: \u2264 {info['max_age']:.0f} d")
                if speed_str:
                    lines.append(f"Ejection speed: {speed_str}")
                lines.append(f"Model mode: {mode_str}")

                # Build block with title first - single ax.text call, no overlap
                block = "\n".join(lines)
                x_pos = 0.98
                y_pos = 0.02 if self._imgArr is not None else 0.98
                va    = "bottom" if self._imgArr is not None else "top"

                self.ax.text(x_pos, y_pos, block,
                             transform=self.ax.transAxes,
                             color="#ff99cc", fontsize=8, fontfamily=MF,
                             ha="right", va=va,
                             linespacing=1.7,
                             bbox=dict(boxstyle="round,pad=0.6",
                                       facecolor="#0d0d0d",
                                       alpha=0.82,
                                       edgecolor="#ff3399",
                                       linewidth=0.8),
                             zorder=9)

        # ── Sun direction arrow ───────────────────────────────────────────
        # sun_dir = (xi, eta) pointing FROM COMET TOWARD SUN
        sun_xi, sun_eta = model["sun_dir"]
        slen = np.sqrt(sun_xi**2 + sun_eta**2)
        if (not contour_only) and slen > 1e-10 and self._vis.get("show_sun", True):

            if overlay:
                sc_px = view_w_px * 0.065               # 6.5% of view width (50% of old 13%)
                sc = sc_px * self.au_per_px            # → AU, for to_px() below
                ax_end, ay_end = to_px(sun_xi/slen*sc, sun_eta/slen*sc)
                ax_start, ay_start = self.nuc_x, self.nuc_y
            else:
                xl2 = self.ax.get_xlim()
                sc  = abs(xl2[1] - xl2[0]) * 0.065  # 50% of old 0.13
                ax_end   = (sun_xi/slen) * sc
                ay_end   = (sun_eta/slen) * sc
                ax_start = 0; ay_start = 0

            self.ax.annotate("", xy=(ax_end, ay_end),
                             xytext=(ax_start, ay_start),
                             arrowprops=dict(arrowstyle="->",
                                             color="#ffe030", lw=2.0), zorder=6)
            self.ax.text(ax_end + (ax_end - ax_start)*0.12,
                         ay_end + (ay_end - ay_start)*0.12,
                         "☀", color="#ffe030", fontsize=13,
                         ha="center", zorder=7)

        # ── Anti-velocity arrow (v3.0) ──────────────────────────────────────
        avx, avy = model.get("antivel_dir", (0.0, 0.0))
        avlen = np.sqrt(avx**2 + avy**2)
        if (not contour_only) and avlen > 1e-10 and self._vis.get("show_antivel", True):
            if overlay:
                sc_px = view_w_px * 0.065               # 50% of old 0.13
                sc = sc_px * self.au_per_px
                avx_end, avy_end = to_px(avx/avlen*sc, avy/avlen*sc)
                avx_start, avy_start = self.nuc_x, self.nuc_y
            else:
                xl2 = self.ax.get_xlim()
                sc  = abs(xl2[1] - xl2[0]) * 0.065
                avx_end   = (avx/avlen) * sc
                avy_end   = (avy/avlen) * sc
                avx_start = 0; avy_start = 0

            self.ax.annotate("", xy=(avx_end, avy_end),
                             xytext=(avx_start, avy_start),
                             arrowprops=dict(arrowstyle="->",
                                             color="#ff5078", lw=1.8), zorder=6)
            self.ax.text(avx_end + (avx_end - avx_start)*0.14,
                         avy_end + (avy_end - avy_start)*0.14,
                         "−v", color="#ff5078", fontsize=10,
                         fontfamily=MF, ha="center", zorder=7)

        # Nucleus crosshair
        nx = self.nuc_x if overlay else 0.0
        ny = self.nuc_y if overlay else 0.0
        if (not contour_only) and self._vis.get("show_crosshair", True):
            self.ax.plot(nx, ny, "+", color="white", ms=14, mew=1.8, zorder=10)
            self.ax.plot(nx, ny, "o", color="white", ms=4,  zorder=11)

        # ── N/E Compass ──────────────────────────────────────────────────
        # Position: lower-right corner in SCREEN axes (fraction 0–1)
        # We use annotate with xycoords='axes fraction' for the base point,
        # and xycoords='data' for the arrow tip (converted below).
        #
        # North PA: angle CW from image-up to celestial North
        #   standalone always uses north_pa=0 (N=up, E=left after inversion)
        #   overlay uses self.north_pa (0 after N-up rotation)
        #
        # In DATA coordinates for the compass:
        #   standalone: x inverted, y normal → N=(0, +aL_c), E=(−aL_c, 0) in data
        #   overlay:    normal x,y image coords → N=(0, −aL_c), E=(−aL_c, 0) in px

        npa_r = np.radians(self.north_pa)

        if not overlay:
            # Standalone: data coords, x-axis inverted (E=left = positive xi = right in data
            # displayed left). N=up = positive eta direction.
            #   North direction in data: (sin(NPA), +cos(NPA))
            #     NPA=0: (0, +1) → up ✓
            #   East direction in data (90° CW from North IN SKY = LEFT on screen):
            #     E always left in astronomical convention.
            #     In inverted-x standalone: data x increases rightward, but display is reversed.
            #     E should appear LEFT on screen = LARGER data x (since axis is inverted).
            #     E data direction: (+cos(NPA), −sin(NPA)) → for NPA=0: (+1, 0) = right in data
            #     which appears LEFT on screen after inversion ✓
            xl  = self.ax.get_xlim()   # (max, min) due to inversion
            yl  = self.ax.get_ylim()
            # compass center: 90% toward the low-xi (West=right-data) side + 10% up
            aL_c = abs(xl[0] - xl[1]) * 0.055   # positive length
            cx  = xl[1] + (xl[0] - xl[1]) * 0.10  # near min-data (= right on screen = West)
            cy  = yl[0] + (yl[1] - yl[0]) * 0.08
            # North arrow: +aL_c in y (up), NPA rotation included
            n_dx =  np.sin(npa_r) * aL_c
            n_dy =  np.cos(npa_r) * aL_c   # positive = upward in data ✓
            # East arrow (90° CCW from N on sky = left on screen):
            # In inverted x-axis, LEFT on screen = POSITIVE data x direction
            e_dx =  np.cos(npa_r) * aL_c   # positive data x → LEFT on screen after inversion ✓
            e_dy = -np.sin(npa_r) * aL_c
        else:
            # Overlay: image pixel coords, y increases downward (origin=
            # "upper" convention, so get_ylim() returns (bottom, top) with
            # the TOP of the current view at the smaller value).
            # BUG FIX (v3.1): previously positioned at a fixed fraction of
            # the FULL image's pixel dimensions (w_img*0.07, h_img*0.07) —
            # correct only when the view happened to show the whole image.
            # Zooming into a sub-region (the normal way to actually look
            # at the comet/tail) moved that fixed point outside the
            # visible area entirely, making the compass silently vanish
            # with nothing to explain why. Use the CURRENT axis limits
            # instead — the same fix the standalone branch above already
            # had, for the same reason.
            xl = self.ax.get_xlim()
            yl = self.ax.get_ylim()
            aL_c = min(abs(xl[1] - xl[0]), abs(yl[1] - yl[0])) * 0.045
            cx = xl[0] + (xl[1] - xl[0]) * 0.10
            cy = yl[1] + (yl[0] - yl[1]) * 0.10   # 10% down from the TOP of the current view
            # North: screen-up = negative py
            n_dx =  np.sin(npa_r) * aL_c
            n_dy = -np.cos(npa_r) * aL_c   # negative py = up in image ✓
            # East: screen-left = negative px
            e_dx = -np.cos(npa_r) * aL_c   # negative px = left ✓
            e_dy = -np.sin(npa_r) * aL_c

        if self._vis.get("show_compass", True):
            # Draw N arrow
            self.ax.annotate("", xy=(cx + n_dx, cy + n_dy), xytext=(cx, cy),
                             arrowprops=dict(arrowstyle="->", color="#60c8ff", lw=2.0),
                             zorder=12)
            self.ax.text(cx + n_dx * 1.35, cy + n_dy * 1.35, "N",
                         color="#60c8ff", fontsize=9, ha="center", va="center",
                         fontfamily=MF, fontweight="bold", zorder=12)
            # Draw E arrow
            self.ax.annotate("", xy=(cx + e_dx, cy + e_dy), xytext=(cx, cy),
                             arrowprops=dict(arrowstyle="->", color="#60c8ff", lw=1.6),
                             zorder=12)
            self.ax.text(cx + e_dx * 1.38, cy + e_dy * 1.38, "E",
                         color="#60c8ff", fontsize=9, ha="center", va="center",
                         fontfamily=MF, fontweight="bold", zorder=12)

        # Legend
        handles = []
        if (not contour_only) and self._vis["synd"]:
            handles.append(Line2D([0],[0], color=SYNDYNE_COLORS[0], lw=1.5,
                                  label="Syndynes  (β = const)"))
        if (not contour_only) and self._vis["sync"]:
            handles.append(Line2D([0],[0], color=SYNC_COLORS[0], lw=1.2,
                                  ls="-.", label="Synchrones (age = const)"))
        if (not contour_only) and self._vis["orbit"]:
            handles.append(Line2D([0],[0], color="#2050a0", lw=0.9,
                                  ls="--", label="Orbital path"))
        show_observed_legend = overlay and bool(self._vis.get("observed_show", True)) and (
            self._vis.get("isophote", False) or contour_only)
        if show_observed_legend:
            handles.append(Line2D([0],[0], color=observed_color,
                                  lw=self._vis.get("observed_lw", 0.6),
                                  ls=observed_ls,
                                  label="Observed isophotes"))
        if self._vis.get("mc_contour", False) and self.mc_contours:
            handles.append(Line2D([0],[0], color=model_color,
                                  lw=self._vis.get("mc_lw", 1.2), ls=model_ls,
                                  label="Monte Carlo model"))
        if not contour_only:
            handles.append(Line2D([0],[0], color="#ffe030", lw=1.5,
                                  marker=">", ms=7, label="Sun direction"))
            handles.append(Line2D([0],[0], color="#ff5078", lw=1.4,
                                  marker=">", ms=6, label="Anti-velocity (−v)"))
        if handles and self._vis.get("show_legend", True):
            if white_canvas:
                self.ax.legend(handles=handles, loc="lower left", fontsize=8,
                               framealpha=0.9, facecolor="white",
                               edgecolor="black", labelcolor="black",
                               prop={"family": MF})
            else:
                self.ax.legend(handles=handles, loc="lower left", fontsize=8,
                               framealpha=0.2, facecolor="#060b14",
                               edgecolor="#1a2540", labelcolor="white",
                               prop={"family": MF})

        info = model["info"]
        r_g  = info.get("r_geo", 1.0)
        # PsAng from model for title
        _sxi, _set = model["sun_dir"]
        _sl = np.sqrt(_sxi**2 + _set**2)
        _PsAng = ((np.degrees(np.arctan2(_sxi/_sl, _set/_sl)) + 180.0) % 360) if _sl > 1e-10 else 0.0
        ovr_tag = "  ★obs" if info.get("obs_override") else ""
        if self._vis.get("show_title", True):
            if contour_only:
                self.ax.set_title(
                    f'{info.get("name","Comet")}  ·  {info["obs_str"]}',
                    color="black" if white_canvas else "#7ab8ff",
                    fontsize=9, fontfamily=MF, pad=8)
            else:
                self.ax.set_title(
                    f'{info.get("name","Comet")}   ·   {info["obs_str"]}   ·   Finson–Probstein\n'
                    f'r☉={info["r_helio"]:.4f} AU  ·  Δ={r_g:.4f} AU  ·  '
                    f'Phase={info["phase_angle"]:.1f}°  ·  PsAng={_PsAng:.1f}°{ovr_tag}',
                    color="#7ab8ff", fontsize=9, fontfamily=MF, pad=8)
        else:
            self.ax.set_title("")

        # Force axes back to the intended view after drawing (matplotlib
        # may auto-expand axes when curves extend beyond) — same view
        # chosen at the top of this function (the image's previous pan/
        # zoom state if there was one, else the full image bounds).
        if overlay:
            if prev_xlim is not None and prev_ylim is not None:
                self.ax.set_xlim(*prev_xlim)
                self.ax.set_ylim(*prev_ylim)
            else:
                self.ax.set_xlim(0, w_img)
                self.ax.set_ylim(h_img, 0)

            # Publication coordinate display.  Data remain in image pixels,
            # but the fixed tick locations are generated in physical units
            # around the nucleus.  This guarantees a 0-km tick at the nucleus
            # whenever it lies inside the visible crop.
            coord_mode = self._vis.get("mc_coordinates", "pixels")
            if coord_mode == "projected_km":
                self._apply_projected_km_axes(white_canvas)
            elif coord_mode == "hide":
                self.ax.set_axis_off()

        self.canvas.draw_idle()

# ─────────────────────────────────────────────────────────────────────────────
#  LEFT CONTROL PANEL
# ─────────────────────────────────────────────────────────────────────────────
class ControlPanel(QScrollArea):
    compute_requested = pyqtSignal(dict, float, list, list, int, int, dict)
    fetch_requested   = pyqtSignal(str, str)
    image_loaded      = pyqtSignal(object)   # object allows None or ndarray
    comet_ready       = pyqtSignal(dict)     # emitted whenever _comet_el is freshly set

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWidgetResizable(True)
        self.setMinimumWidth(320)
        self.setMaximumWidth(400)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        # v3.1 — ejection-velocity state, set EXCLUSIVELY via
        # set_ejection_params() from the Monte Carlo window's Run button.
        # See get_ejection_params()/set_ejection_params() below.
        self._current_ejection = dict(v_R0=0.0, v_T0=0.0, v_N0=0.0,
                                      gamma=0.0, m_exp=0.0)

        # v3.1 — replaces the old "MC morphology contour" checkbox: there's
        # no manual show/hide step anymore, the contour just appears as
        # soon as a Monte Carlo run is sent to the main canvas, and is
        # cleared automatically when switching to a different comet. See
        # MainWindow._on_mc_done() / _on_comet_changed().
        self._mc_contour_visible = False
        # v3.1 — MC presentation state, controlled from MCWindow's Display
        # tab.  The original analysis overlay remains the default; the two
        # contour-only modes use the loaded image only as the SOURCE of the
        # observed isophotes while suppressing the raster background.
        self._mc_style = dict(
            lw=1.2, alpha=0.95, info_box=True,
            view_style="analysis",
            background="image",
            coordinates="pixels",
            observed_show=True,
            observed_color=ISOPHOTE_COLOR,
            observed_lw=0.6,
            observed_alpha=0.75,
            observed_ls="-",
            model_color="#ff3399",
            model_ls="--",
            show_legend=True,
            show_title=True,
            show_compass=True,
            grayscale_safe=False,
        )

        inner = QWidget()
        self.setWidget(inner)
        vbox  = QVBoxLayout(inner)
        vbox.setContentsMargins(8, 8, 8, 8)
        vbox.setSpacing(8)

        # ── Header ───────────────────────────────────────────────────────
        hdr_row = QHBoxLayout()
        title = QLabel("☄  COMET TAIL ANALYZER")
        title.setProperty("class", "title")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        hdr_row.addWidget(title, 1)
        self.btn_theme = QPushButton("☀")
        self.btn_theme.setToolTip("Switch to Light mode")
        self.btn_theme.setFixedSize(28, 28)
        self.btn_theme.setStyleSheet("font-size:14px; padding:0; min-height:0; border-radius:14px;")
        self.btn_theme.clicked.connect(self._toggle_theme)
        hdr_row.addWidget(self.btn_theme)
        vbox.addLayout(hdr_row)
        sub = QLabel("FINSON–PROBSTEIN MODEL  ·  1968")
        sub.setProperty("class", "muted")
        sub.setAlignment(Qt.AlignmentFlag.AlignCenter)
        sub.setStyleSheet("font-size:9px; letter-spacing:2px;")
        vbox.addWidget(sub)

        sep = QFrame(); sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color:#1a2540;")
        vbox.addWidget(sep)

        # ── Comet tabs ────────────────────────────────────────────────────
        self.comet_tabs = QTabWidget()
        vbox.addWidget(self.comet_tabs)

        # Preset
        preset_w = QWidget()
        pv = QVBoxLayout(preset_w); pv.setSpacing(6); pv.setContentsMargins(6,8,6,6)
        pv.addWidget(self._lbl("SELECT COMET"))
        self.combo_comet = QComboBox()
        self.combo_comet.addItems(list(cta.COMET_DB.keys()))
        self.combo_comet.currentIndexChanged.connect(self._on_comet_selected)
        pv.addWidget(self.combo_comet)
        self.note_label = QLabel()
        self.note_label.setWordWrap(True)
        self.note_label.setStyleSheet("color:#1e4060; font-size:10px;")
        pv.addWidget(self.note_label)
        pv.addWidget(self._lbl("OBSERVATION DATE (UT)"))
        self.obs_date = QLineEdit()
        self.obs_date.setPlaceholderText("YYYY-MM-DD  or  YYYY-MM-DD.fff")
        self.obs_date.setToolTip(
            "Observation date in UT.\n"
            "Accepted formats:\n"
            "  2018-01-14          (midnight)\n"
            "  2018-01-14.583      (decimal day — astronomical standard)\n"
            "  2018-01-14 13:55    (HH:MM)\n"
            "  2018-01-14 13:55:12 (HH:MM:SS)"
        )
        pv.addWidget(self.obs_date)          # was missing — field was invisible
        self.btn_use = QPushButton("USE THIS COMET")
        self.btn_use.clicked.connect(self._use_preset)
        pv.addWidget(self.btn_use)
        self.preset_status = QLabel("")
        self.preset_status.setWordWrap(True)
        self.preset_status.setStyleSheet("color:#2a5060; font-size:10px;")
        pv.addWidget(self.preset_status)
        pv.addStretch()
        self.comet_tabs.addTab(preset_w, "PRESET")
        self._on_comet_selected(0)           # seed note + preferred date for first item

        # Manual
        manual_w = QWidget()
        mv = QFormLayout(manual_w); mv.setSpacing(5); mv.setContentsMargins(6,8,6,6)
        self.m_q  = self._dspin(0.2946, 7, 0.00001)
        self.m_e  = self._dspin(0.9992, 7, 0.000001)
        self.m_i  = self._dspin(128.938, 4, 0.001)
        self.m_Om = self._dspin(61.021, 4, 0.001)
        self.m_om = self._dspin(37.281, 4, 0.001)
        self.m_T  = QLineEdit("2020-07-03")
        self.m_obs= QLineEdit("2020-07-23")
        self.m_name= QLineEdit("Custom Comet")
        for label, w in [("Name",self.m_name),("q (AU)",self.m_q),("e",self.m_e),
                          ("i (°)",self.m_i),("Ω (°)",self.m_Om),("ω (°)",self.m_om),
                          ("T perihelion",self.m_T),("Obs date",self.m_obs)]:
            mv.addRow(QLabel(label), w)
        btn_man = QPushButton("USE THESE ELEMENTS")
        btn_man.clicked.connect(self._use_manual)
        mv.addRow(btn_man)
        self.comet_tabs.addTab(manual_w, "MANUAL")

        # Fetch
        fetch_w = QWidget()
        fv = QVBoxLayout(fetch_w); fv.setSpacing(6); fv.setContentsMargins(6,8,6,6)
        fv.addWidget(self._lbl("DESIGNATION"))
        self.fetch_desig = QLineEdit()
        self.fetch_desig.setPlaceholderText("e.g.  C/2023 A3")
        fv.addWidget(self.fetch_desig)
        fv.addWidget(self._lbl("OBS DATE (optional)"))
        self.fetch_date = QLineEdit()
        self.fetch_date.setPlaceholderText("YYYY-MM-DD  or  YYYY-MM-DD.fff")
        fv.addWidget(self.fetch_date)
        self.btn_fetch = QPushButton("FETCH ORBITAL ELEMENTS")
        self.btn_fetch.clicked.connect(self._do_fetch)
        fv.addWidget(self.btn_fetch)
        self.fetch_status = QLabel("")
        self.fetch_status.setWordWrap(True)
        self.fetch_status.setStyleSheet("color:#2a5060; font-size:10px;")
        fv.addWidget(self.fetch_status)
        fv.addStretch()
        self.comet_tabs.addTab(fetch_w, "FETCH JPL")

        # ── Model params ──────────────────────────────────────────────────
        grp_model = QGroupBox("⚙  MODEL")
        gm = QVBoxLayout(grp_model); gm.setSpacing(5)

        gm.addWidget(self._lbl("β VALUES  (comma-separated)"))
        self.beta_str = QLineEdit("0.001,0.01,0.05,0.1,0.3,0.6,1.0")
        gm.addWidget(self.beta_str)
        self.beta_hint = QLabel()
        self.beta_hint.setStyleSheet("color:#1a3a50; font-size:9px;")
        self.beta_hint.setWordWrap(True)
        gm.addWidget(self.beta_hint)
        self.beta_str.textChanged.connect(self._update_beta_hint)
        self._update_beta_hint()

        gm.addWidget(self._lbl("SYNCHRONE AGES  (days, comma-sep.)"))
        self.age_str = QLineEdit("10,30,60,90,120,180")
        self.age_str.setToolTip(
            "Enter the F-P synchrone ages to draw. The largest valid age is used\n"
            "automatically as the maximum dust age and the upper limit for MC guidance.")
        gm.addWidget(self.age_str)

        # Maximum dust age is derived from the largest listed synchrone age.
        # Keep a hidden numeric widget because the existing compute/animation
        # pipeline already reads self.max_age.value().
        self.max_age = QSpinBox()
        self.max_age.setRange(1, 3650)
        self.max_age.setValue(180)
        self.max_age.setVisible(False)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Maximum dust age"))
        self.lbl_max_age_auto = QLabel("180 d  — automatic")
        self.lbl_max_age_auto.setStyleSheet("color:#90d8ff; font-weight:600;")
        self.lbl_max_age_auto.setToolTip(
            "Automatically set from the largest valid Synchrone age.\n"
            "The same value is used as the upper F-P age limit for the\n"
            "Monte Carlo release-window recommendation.")
        row1.addWidget(self.lbl_max_age_auto, 1)
        gm.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Dominant dust age (d)"))
        self.dominant_age = QDoubleSpinBox()
        self.dominant_age.setRange(0.0, 3650.0)
        self.dominant_age.setDecimals(1)
        self.dominant_age.setSingleStep(5.0)
        self.dominant_age.setSuffix(" d")
        self.dominant_age.setSpecialValueText("Not specified")
        self.dominant_age.setValue(0.0)
        self.dominant_age.setMinimumWidth(135)
        self.dominant_age.setToolTip(
            "Approximate age of the single F-P synchrone that most closely follows\n"
            "the main visible dust structure. This is a user-interpreted morphology\n"
            "estimate, not an automatic or formal best-fit result.")
        row2.addWidget(self.dominant_age)
        gm.addLayout(row2)

        # Resolution is a technical sampling parameter, not a morphology input.
        # Keep it internal at a robust routine default.
        self.n_pts = QSpinBox()
        self.n_pts.setRange(20, 300)
        self.n_pts.setValue(100)
        self.n_pts.setVisible(False)

        self.age_str.textChanged.connect(self._sync_max_age_from_synchrone_text)
        self._sync_max_age_from_synchrone_text()

        # ── PA Rotation Offset ────────────────────────────────────────────
        gm.addWidget(self._lbl("GRID ROTATION  (match observed tail)"))
        rot_hint = QLabel(
            "Rotate entire model grid around nucleus to align\n"
            "with observed tail. Corrects for uncertain ω / Ω\n"
            "in orbital elements (= orbital trail direction).")
        rot_hint.setStyleSheet("color:#3a6070; font-size:9px;")
        rot_hint.setWordWrap(True)
        gm.addWidget(rot_hint)

        row_rot = QHBoxLayout()
        self.rot_slider = QSlider(Qt.Orientation.Horizontal)
        self.rot_slider.setRange(-300, 300)   # ± 30.0°
        self.rot_slider.setValue(0)
        self.rot_slider.setTickInterval(50)
        self.rot_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.rot_label = QLabel("0.0°")
        self.rot_label.setMinimumWidth(42)
        self.rot_label.setStyleSheet("color:#90d8ff; font-size:11px;")
        self.rot_slider.valueChanged.connect(
            lambda v: self.rot_label.setText(f"{v/10:+.1f}°"))
        row_rot.addWidget(self.rot_slider, 1)
        row_rot.addWidget(self.rot_label)
        gm.addLayout(row_rot)

        row_rot2 = QHBoxLayout()
        btn_rot_m = QPushButton("−1°"); btn_rot_m.setMaximumWidth(42)
        btn_rot_0 = QPushButton("0°");  btn_rot_0.setMaximumWidth(42)
        btn_rot_p = QPushButton("+1°"); btn_rot_p.setMaximumWidth(42)
        btn_rot_m.clicked.connect(lambda: self.rot_slider.setValue(self.rot_slider.value()-10))
        btn_rot_0.clicked.connect(lambda: self.rot_slider.setValue(0))
        btn_rot_p.clicked.connect(lambda: self.rot_slider.setValue(self.rot_slider.value()+10))
        row_rot2.addStretch(); row_rot2.addWidget(btn_rot_m)
        row_rot2.addWidget(btn_rot_0); row_rot2.addWidget(btn_rot_p)
        gm.addLayout(row_rot2)
        self.rot_info_lbl = QLabel("Orbital trail PA: — °  |  PsAng: — °")
        self.rot_info_lbl.setStyleSheet("color:#2a5060; font-size:9px;")
        self.rot_info_lbl.setWordWrap(True)
        gm.addWidget(self.rot_info_lbl)

        vbox.addWidget(grp_model)

        # Ejection velocity (v3.1) is now set EXCLUSIVELY from the Monte
        # Carlo window (View > Monte Carlo morphology…) — see ControlPanel.
        # get_ejection_params()/set_ejection_params() below. There used to
        # be a second, separate set of v_R0/v_T0/v_N0/γ/m fields here, but
        # that meant keying the same kind of parameter into two different
        # windows (and the two used different ejection models — fixed R/T/N
        # direction here vs. isotropic-random in the MC window — which
        # invited them to silently disagree). One input location now.

        # ── Display options ───────────────────────────────────────────────
        grp_disp = QGroupBox("DISPLAY")
        gd = QVBoxLayout(grp_disp); gd.setSpacing(4)
        self.chk_synd  = QCheckBox("Syndynes (β = const)")
        self.chk_sync  = QCheckBox("Synchrones (age = const)")
        self.chk_orbit = QCheckBox("Orbital path")
        self.chk_isophote = QCheckBox("Isophotes (from image)")
        self.chk_sun   = QCheckBox("Sun direction ☀")
        self.chk_antivel = QCheckBox("Anti-velocity (−v)")
        self.chk_crosshair = QCheckBox("Nucleus crosshair (+)")
        self.chk_synd.setChecked(True)
        self.chk_sync.setChecked(True)
        self.chk_orbit.setChecked(True)
        self.chk_isophote.setChecked(False)
        self.chk_sun.setChecked(True)
        self.chk_antivel.setChecked(True)
        self.chk_crosshair.setChecked(True)
        self.chk_synd.setStyleSheet("color:#ff7070;")
        self.chk_sync.setStyleSheet("color:#ffd060;")
        self.chk_orbit.setStyleSheet("color:#4a8adf;")
        self.chk_isophote.setStyleSheet("color:#60e0a0;")
        self.chk_sun.setStyleSheet("color:#ffe030;")
        self.chk_antivel.setStyleSheet("color:#ff5078;")
        self.chk_isophote.setToolTip(
            "Trace surface-brightness contours directly from the loaded\n"
            "image (overlay mode only) for direct comparison between the\n"
            "observed tail morphology and the model curves. (v3.0)")
        self.chk_sun.setToolTip("Show/hide the yellow Sun-direction arrow and ☀ label.")
        self.chk_antivel.setToolTip("Show/hide the red –v (anti-velocity) arrow and label.")
        self.chk_crosshair.setToolTip("Show/hide the white + crosshair marking the nucleus position.")
        # MC morphology contour no longer has its own checkbox — it is
        # drawn automatically on the main canvas as soon as a Monte Carlo
        # run is sent over (see ControlPanel._mc_contour_visible below),
        # the same way Run already updates the canvas with no extra step.
        for w in [self.chk_synd, self.chk_sync, self.chk_orbit, self.chk_isophote,
                  self.chk_sun, self.chk_antivel, self.chk_crosshair]:
            gd.addWidget(w)

        gd.addWidget(QLabel("Syndynes/Synchrones:"))
        row_lw = QHBoxLayout()
        row_lw.addWidget(QLabel("  Width"))
        self.lw_slider = QSlider(Qt.Orientation.Horizontal)
        self.lw_slider.setRange(5, 50); self.lw_slider.setValue(15)
        self.lw_label  = QLabel("1.5")
        self.lw_slider.valueChanged.connect(lambda v: self.lw_label.setText(f"{v/10:.1f}"))
        row_lw.addWidget(self.lw_slider); row_lw.addWidget(self.lw_label)
        gd.addLayout(row_lw)

        row_op = QHBoxLayout()
        row_op.addWidget(QLabel("  Opacity"))
        self.op_slider = QSlider(Qt.Orientation.Horizontal)
        self.op_slider.setRange(10, 100); self.op_slider.setValue(50)
        self.op_label  = QLabel("0.85")
        self.op_slider.valueChanged.connect(lambda v: self.op_label.setText(f"{v/100:.2f}"))
        row_op.addWidget(self.op_slider); row_op.addWidget(self.op_label)
        gd.addLayout(row_op)

        # ── Orbital-path-specific width/opacity (v3.0) ────────────────────
        # Previously hardcoded at lw=0.9/alpha=0.5 regardless of the
        # general sliders above (which only affect syndynes/synchrones) —
        # against a busy starfield/isophote background that combination
        # is nearly invisible, with no way to turn it up.
        gd.addWidget(QLabel("Orbital path:"))
        row_orbit_lw = QHBoxLayout()
        row_orbit_lw.addWidget(QLabel("  Width"))
        self.orbit_lw_slider = QSlider(Qt.Orientation.Horizontal)
        self.orbit_lw_slider.setRange(5, 50); self.orbit_lw_slider.setValue(15)
        self.orbit_lw_label  = QLabel("1.5")
        self.orbit_lw_slider.valueChanged.connect(
            lambda v: self.orbit_lw_label.setText(f"{v/10:.1f}"))
        row_orbit_lw.addWidget(self.orbit_lw_slider); row_orbit_lw.addWidget(self.orbit_lw_label)
        gd.addLayout(row_orbit_lw)

        row_orbit_op = QHBoxLayout()
        row_orbit_op.addWidget(QLabel("  Opacity"))
        self.orbit_op_slider = QSlider(Qt.Orientation.Horizontal)
        self.orbit_op_slider.setRange(10, 100); self.orbit_op_slider.setValue(100)
        self.orbit_op_label  = QLabel("1.00")
        self.orbit_op_slider.valueChanged.connect(
            lambda v: self.orbit_op_label.setText(f"{v/100:.2f}"))
        row_orbit_op.addWidget(self.orbit_op_slider); row_orbit_op.addWidget(self.orbit_op_label)
        gd.addLayout(row_orbit_op)

        # ── Isophote-specific controls (v3.0) ─────────────────────────────
        row_isolvl = QHBoxLayout()
        row_isolvl.addWidget(QLabel("Isophote levels"))
        self.isolvl_slider = QSlider(Qt.Orientation.Horizontal)
        self.isolvl_slider.setRange(3, 15); self.isolvl_slider.setValue(6)
        self.isolvl_label  = QLabel("6")
        self.isolvl_slider.valueChanged.connect(
            lambda v: self.isolvl_label.setText(str(v)))
        row_isolvl.addWidget(self.isolvl_slider); row_isolvl.addWidget(self.isolvl_label)
        gd.addLayout(row_isolvl)

        row_isosm = QHBoxLayout()
        row_isosm.addWidget(QLabel("Isophote smoothing (px)"))
        self.isosm_slider = QSlider(Qt.Orientation.Horizontal)
        self.isosm_slider.setRange(0, 60); self.isosm_slider.setValue(20)
        self.isosm_label  = QLabel("2.0")
        self.isosm_slider.valueChanged.connect(
            lambda v: self.isosm_label.setText(f"{v/10:.1f}"))
        row_isosm.addWidget(self.isosm_slider); row_isosm.addWidget(self.isosm_label)
        gd.addLayout(row_isosm)
        isosm_hint = QLabel(
            "Real exposures need smoothing >0 — contouring raw pixel noise\n"
            "looks like speckle, not clean isophote curves. Increase if\n"
            "contours still look noisy; decrease if the coma looks blobby.")
        isosm_hint.setStyleSheet("color:#2a5060; font-size:9px;")
        isosm_hint.setWordWrap(True)
        gd.addWidget(isosm_hint)

        # v3.1 — MC contour width/opacity/info-box controls used to live
        # here, but with the visibility checkbox gone (contour now shows
        # automatically) these were the last MC-specific bits left in the
        # main DISPLAY panel — moved out to keep that panel free of
        # anything you only ever touch from inside the Monte Carlo
        # window. Width/opacity/info-box are now set on the MC window's
        # own "3. Display" tab instead (see MCWindow), next to the other
        # re-thresholding controls that already lived there.

        vbox.addWidget(grp_disp)

        # ── Image overlay (compact – full config in dialog) ──────────────
        grp_img = QGroupBox("🖼  IMAGE OVERLAY")
        gi = QVBoxLayout(grp_img); gi.setSpacing(5)

        row_btns = QHBoxLayout()
        self.btn_img_open = QPushButton("⚙  OPEN / SETUP IMAGE…")
        self.btn_img_open.clicked.connect(self._open_image_dialog)
        row_btns.addWidget(self.btn_img_open)
        gi.addLayout(row_btns)

        self.img_label = QLabel("No image loaded")
        self.img_label.setStyleSheet("color:#5a90b0; font-size:10px;")
        self.img_label.setWordWrap(True)
        gi.addWidget(self.img_label)

        # Hidden stub widgets — parented to self so they stay hidden and
        # are destroyed with the panel (never become floating windows)
        _stub = QWidget(self)   # invisible container with proper parent
        _stub.setVisible(False)
        _stub_lay = QVBoxLayout(_stub)

        self.nuc_x_spin = self._dspin(0, 1, 1.0)
        self.nuc_y_spin = self._dspin(0, 1, 1.0)
        self.au_px_spin = self._dspin(0.001, 9, 0.000001)
        self.npa_spin   = self._dspin(0.0, 2, 0.5)
        self.npa_hint   = QLabel("")
        self.sl_contrast  = QSlider(Qt.Orientation.Horizontal)
        self.sl_intensity = QSlider(Qt.Orientation.Horizontal)
        self.sl_contrast.setRange(-1000, 1000);  self.sl_contrast.setValue(0)
        self.sl_intensity.setRange(-1000, 1000); self.sl_intensity.setValue(0)
        self.btn_pick_nuc  = QPushButton("⊕  CLICK NUCLEUS ON PLOT")
        self.btn_rot_nup   = QPushButton("↻  ROTATE N-UP / E-LEFT")
        self.rot_status    = QLabel("")
        # View persistence (v3.0.4) — replaces the discrete 1x-4x preset
        # buttons from v3.0/v3.0.1, and the callback-based capture from
        # v3.0.2/v3.0.3 (which turned out unreliable — see CHANGELOG). Use
        # the plot toolbar's own zoom (rubber-band) and pan (hand) tools
        # directly to frame the overlay however you want — draw_model()
        # reads whatever x/y range that leaves the axes at directly off
        # the axes itself each time it redraws, so it persists across
        # Compute Model re-runs automatically, with no separate "lock"
        # state to go stale. This button is the only explicit control
        # needed: it snaps the view back to showing the full image.
        self.btn_reset_view = QPushButton("⟲  RESET TO FULL IMAGE")
        self.btn_reset_view.clicked.connect(self._reset_overlay_view)
        self.btn_auto_stretch = QPushButton("✨  AUTO STRETCH")
        self.btn_auto_stretch.setToolTip("Automatically compute a reasonable starting Contrast/Intensity from the image histogram")
        self.btn_auto_stretch.clicked.connect(self._auto_stretch)
        # v3.0: live restretch — no more separate Apply button. Re-stretches
        # and redraws on every slider tick (not just on release), which is
        # cheap enough for typical comet-image sizes; if it ever feels
        # sluggish on very large images, debounce with a short QTimer here.
        self.sl_contrast.valueChanged.connect(self._restretch)
        self.sl_intensity.valueChanged.connect(self._restretch)
        self.btn_rot_nup.clicked.connect(self._rotate_north_up)
        self.npa_spin.valueChanged.connect(self._update_npa_hint)

        # Give all stub widgets the hidden container as parent
        for w in [self.nuc_x_spin, self.nuc_y_spin, self.au_px_spin,
                  self.npa_spin, self.npa_hint, self.sl_contrast, self.sl_intensity,
                  self.btn_pick_nuc, self.btn_reset_view,
                  self.btn_rot_nup, self.rot_status]:
            _stub_lay.addWidget(w)

        # btn_clear_img — must also have a parent, not be orphaned
        self.btn_clear_img = QPushButton("CLEAR", _stub)
        self.btn_clear_img.setProperty("class","danger")
        self.btn_clear_img.clicked.connect(self._clear_image)

        # overlay_params stub
        self.overlay_params = QWidget(_stub)
        self.btn_img = self.btn_img_open   # alias

        vbox.addWidget(grp_img)

        # ── Monte Carlo overlay control ───────────────────────────────────
        # The MC window may be closed while its contour remains on the main
        # canvas, so provide a persistent, obvious way to remove that result.
        grp_mc_overlay = QGroupBox("🧹  MODEL RESULTS")
        gmc = QVBoxLayout(grp_mc_overlay)
        gmc.setContentsMargins(8, 7, 8, 8)
        self.btn_clear_mc = QPushButton("✕  CLEAR ALL MODELS")
        self.btn_clear_mc.setProperty("class", "danger")
        self.btn_clear_mc.setToolTip(
            "Clear all computed model results from the main canvas: "
            "Finson–Probstein curves/vectors and Monte Carlo contours. "
            "The loaded image and editable MC inputs are preserved; hidden "
            "MC-coupled F-P velocity is reset to the zero-velocity baseline.")
        self.btn_clear_mc.setEnabled(False)
        self.btn_clear_mc.clicked.connect(
            lambda: getattr(self.window(), "_confirm_clear_all_models")())
        gmc.addWidget(self.btn_clear_mc)
        vbox.addWidget(grp_mc_overlay)
        # v3.1.1: model actions now live in the dedicated Model menu and
        # compact toolbar. Keep the button object for backward-compatible
        # signal/state handling, but remove the duplicate left-panel section.
        grp_mc_overlay.setVisible(False)

        # ── Observed Ephemeris Override (also in dialog) ──────────────────
        self.chk_use_obs = QCheckBox("★ Ephemeris override active")
        self.chk_use_obs.setStyleSheet("color:#5a8060; font-size:10px;")
        self.chk_use_obs.setEnabled(False)    # enabled via dialog
        vbox.addWidget(self.chk_use_obs)
        self.obs_rhelio  = self._dspin(0.576, 5, 0.001)
        self.obs_rgeo    = self._dspin(1.007, 5, 0.001)
        self.obs_phase   = self._dspin(72.8,  2, 0.1)
        self.obs_psang   = self._dspin(288.3, 2, 0.1)
        self.psang_hint  = QLabel("")
        self.obs_psang.valueChanged.connect(self._update_psang_hint)

        # ── Compute button ────────────────────────────────────────────────
        self.btn_compute = QPushButton("▶   COMPUTE MODEL")
        self.btn_compute.setObjectName("btn_compute")
        self.btn_compute.setProperty("class", "primary")
        self.btn_compute.setMinimumHeight(40)
        self.btn_compute.clicked.connect(self._emit_compute)
        vbox.addWidget(self.btn_compute)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumHeight(6)
        vbox.addWidget(self.progress_bar)

        vbox.addStretch()

        # ── Credit footer ─────────────────────────────────────────────────
        credit = QLabel(
            "☄ Teerasak Thaluang\n"
            "MPC O51 / O58\n"
            "Finson–Probstein Model · 2026")
        credit.setAlignment(Qt.AlignmentFlag.AlignCenter)
        credit.setStyleSheet(
            "color:#1a3050; font-size:9px; font-family:'Segoe UI',Arial,sans-serif;"
            "border-top: 1px solid #1a2540; padding:8px 0 4px 0;")
        vbox.addWidget(credit)

        # State
        self._comet_el   = None
        self._img_arr    = None
        self._fits_raw   = None
        self._fits_header= None
        self._wcs_ps_deg = None
        self._wcs_npa    = None
        self._on_comet_selected(0)

    # ── Helpers ───────────────────────────────────────────────────────────
    def _lbl(self, text):
        l = QLabel(text)
        l.setProperty("class", "section")
        return l

    def _dspin(self, val, decimals, step):
        w = QDoubleSpinBox()
        w.setDecimals(decimals); w.setRange(-1e9, 1e9)
        w.setValue(val); w.setSingleStep(step)
        return w

    def _open_image_dialog(self):
        """Open the Image Setup floating dialog."""
        if not hasattr(self, '_img_dialog') or self._img_dialog is None:
            self._img_dialog = ImageSetupDialog(self)
        self._img_dialog.show()
        # NOTE: raise_()/activateWindow() called synchronously right after
        # show() can run before the OS has finished realizing the native
        # window (HWND on Windows), leaving the dialog visible but not yet
        # registered with the window manager as a normal draggable/
        # focusable top-level window — opening a modal QFileDialog later
        # would force the event loop to flush and "fix" it retroactively,
        # which is why the symptom only appeared before an image was
        # selected. Deferring via singleShot(0, ...) lets the event loop
        # finish realizing the window first.
        QTimer.singleShot(0, self._img_dialog.raise_)
        QTimer.singleShot(0, self._img_dialog.activateWindow)

    def _toggle_theme(self):
        """Toggle dark/light theme — called from the ☀/☾ button in the header."""
        mw = self.window()
        if hasattr(mw, '_set_theme'):
            new = "light" if CURRENT_THEME == "dark" else "dark"
            mw._set_theme(new)

    def _update_beta_hint(self):
        betas = self._parse_floats(self.beta_str.text())
        rho_d = getattr(self.window(), 'phys_params', {}).get('rho_d', 0.5)
        hints = [f"β{b}→{cta.beta_to_size(b, rho_d)}" for b in betas[:4]]
        self.beta_hint.setText("  ".join(hints))

    def _update_psang_hint(self):
        psang = self.obs_psang.value()
        anti  = (psang + 180.0) % 360.0
        self.psang_hint.setText(
            f"Anti-solar (tail) dir: PA = {anti:.1f}°\n"
            f"(PsAng {psang:.1f}° − 180°)  from N toward E")

    def get_obs_override(self):
        """Return observed ephemeris override values if checkbox is checked, else None."""
        if not self.chk_use_obs.isChecked():
            return None
        return dict(
            r_helio    = self.obs_rhelio.value(),
            r_geo      = self.obs_rgeo.value(),
            phase_angle= self.obs_phase.value(),
            psang      = self.obs_psang.value(),   # JPL Horizons PsAng (deg)
        )

    def _parse_floats(self, s):
        try:
            return [float(x) for x in s.split(",") if x.strip() and float(x.strip()) > 0]
        except: return []

    def _parse_ints(self, s):
        try:
            return [int(x) for x in s.split(",") if x.strip()]
        except: return []

    def _sync_max_age_from_synchrone_text(self, *_):
        """Set Maximum dust age = max(valid Synchrone ages)."""
        ages = self._parse_floats(self.age_str.text())
        if not ages:
            if hasattr(self, 'lbl_max_age_auto'):
                self.lbl_max_age_auto.setText("—  no valid synchrone age")
                self.lbl_max_age_auto.setStyleSheet(
                    "color:#ffb020; font-weight:600;")
            return
        max_age = max(1, int(round(max(ages))))
        self.max_age.setValue(max_age)
        if hasattr(self, 'dominant_age'):
            self.dominant_age.setMaximum(float(max_age))
        if hasattr(self, 'lbl_max_age_auto'):
            self.lbl_max_age_auto.setText(f"{max_age:g} d  — automatic")
            self.lbl_max_age_auto.setStyleSheet(
                "color:#90d8ff; font-weight:600;")


    # ── Preset comet ──────────────────────────────────────────────────────
    def _on_comet_selected(self, idx):
        key = list(cta.COMET_DB.keys())[idx]
        meta = cta.COMET_DB[key]
        self.note_label.setText(meta.get("note", ""))
        self.obs_date.setText(meta.get("obs", ""))

    def _use_preset(self):
        """Fetch live orbital elements from Horizons for the selected catalogue comet."""
        key  = self.combo_comet.currentText()
        date = self.obs_date.text().strip()
        if not key:
            return
        self._pending_obs_date = date        # consumed by on_fetch_done
        self._fetch_source     = 'preset'    # so callbacks know which status to update
        self.btn_use.setEnabled(False)
        self.preset_status.setText("Fetching from JPL Horizons…")
        self.fetch_requested.emit(key, date)

    def _use_manual(self):
        try:
            T_jd = cta.date_to_jd(self.m_T.text())
            self._comet_el = dict(
                q=self.m_q.value(), e=self.m_e.value(), i=self.m_i.value(),
                Omega=self.m_Om.value(), omega=self.m_om.value(),
                T=self.m_T.text(), T_jd=T_jd, name=self.m_name.text(),
                obs_jd=cta.date_to_jd(self.m_obs.text()))
            self.comet_ready.emit(self._comet_el)
        except Exception as ex:
            QMessageBox.warning(self, "Input Error", str(ex))

    # ── Fetch ─────────────────────────────────────────────────────────────
    def _do_fetch(self):
        desig = self.fetch_desig.text().strip()
        if not desig:
            self.fetch_status.setText("Enter a designation.")
            return
        self._fetch_source = 'fetch'
        self.btn_fetch.setEnabled(False)
        self.fetch_status.setText("Fetching from JPL Horizons…")
        self.fetch_requested.emit(desig, self.fetch_date.text())

    def on_fetch_done(self, el):
        self._comet_el = el
        pending = getattr(self, '_pending_obs_date', '')
        self._pending_obs_date = ''
        date_str = pending or self.fetch_date.text() or el.get("T", "")[:10]
        self._comet_el["obs_jd"] = cta.date_to_jd(date_str)
        msg = (f"✓  {el.get('name','')[:40]}\n"
               f"q={el['q']:.5f}  e={el['e']:.6f}  i={el['i']:.2f}°")
        if getattr(self, '_fetch_source', '') == 'preset':
            self.preset_status.setText(msg)
            self.btn_use.setEnabled(True)
        else:
            self.fetch_status.setText(msg)
            self.btn_fetch.setEnabled(True)
        self.comet_ready.emit(self._comet_el)

    def on_fetch_error(self, msg):
        if getattr(self, '_fetch_source', '') == 'preset':
            self.preset_status.setText(f"✗  {msg}")
            self.btn_use.setEnabled(True)
        else:
            self.fetch_status.setText(f"✗  {msg}")
            self.btn_fetch.setEnabled(True)

    # ── Image ──────────────────────────────────────────────────────────────
    def _load_image(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Comet Image", "",
            "All Supported (*.fits *.fit *.fts *.jpg *.jpeg *.png *.tif *.tiff);;"
            "FITS Images (*.fits *.fit *.fts);;"
            "Raster Images (*.jpg *.jpeg *.png *.tif *.tiff)")
        if not path: return
        self._last_image_path = path   # store for MainWindow
        ext = os.path.splitext(path)[1].lower()
        self._fits_header = None   # reset
        self._wcs_ps_deg  = None   # reset FITS WCS metadata to prevent stale state
        self._wcs_npa     = None

        if ext in (".fits", ".fit", ".fts"):
            self._img_arr, self._fits_header = self._load_fits(path)
        else:
            from PIL import Image as PILImage
            pil = PILImage.open(path).convert("RGB")
            self._img_arr = np.array(pil)

        if self._img_arr is None:
            return

        h, w = self._img_arr.shape[:2]

        # ── Auto-calibrate from FITS WCS ─────────────────────────────────
        if self._fits_header is not None:
            hdr = self._fits_header
            # Nucleus position (from plate-solver keyword)
            nx = float(hdr.get("OBJPOSX", w/2))
            ny = float(hdr.get("OBJPOSY", h/2))
            # Plate scale: sqrt(|det(CD)|) → deg/px (correct area element)
            # CD matrix maps pixel coords to sky coords (RA, Dec)
            cd11 = float(hdr.get("CD1_1", 0))
            cd12 = float(hdr.get("CD1_2", 0))
            cd21 = float(hdr.get("CD2_1", 0))
            cd22 = float(hdr.get("CD2_2", 1e-4))
            det = abs(cd11*cd22 - cd12*cd21)
            # Use determinant if valid, fallback to single-column norm for legacy headers
            if det > 1e-14:
                ps_deg = float(np.sqrt(det))
            else:
                ps_deg = float(np.sqrt(cd12**2 + cd22**2)) or 1e-4  # fallback
            # North PA: angle of North vector from image-up (row-decreasing dir)
            npa = float(np.degrees(np.arctan2(cd21, cd22))) # CW from top
            # Store for caller; AU/px filled when model computed
            self._wcs_ps_deg = ps_deg    # deg/px
            self._wcs_npa    = npa       # North PA

            self.nuc_x_spin.setValue(nx)
            self.nuc_y_spin.setValue(ny)
            self.npa_spin.setValue(round(npa, 2))
            # Preliminary AU/px estimate (1 AU ~ 8.3 lmin; use 1 AU default)
            self.au_px_spin.setValue(round(ps_deg * (3.14159/180) * 1.0, 8))

            # Pre-fill observation date+time from FITS header. date_to_jd()
            # already accepts the full ISO-8601 'YYYY-MM-DDTHH:MM:SS[.fff]'
            # format DATE-OBS normally uses, so the WHOLE string is passed
            # through (previously this truncated to date_obs[:10], silently
            # discarding the time-of-day entirely). Some older FITS files
            # split the time into a separate TIME-OBS keyword instead of
            # embedding it in DATE-OBS — combine them in that case.
            date_obs = hdr.get("DATE-OBS", "") or hdr.get("DATE_OBS", "") or hdr.get("DATE-AVG", "")
            if date_obs and len(date_obs) <= 10 and hdr.get("TIME-OBS"):
                date_obs = f"{date_obs}T{hdr['TIME-OBS']}"
            if date_obs:
                # BUG FIX (v3.1): there are THREE separate date fields, one
                # per comet_tabs tab (PRESET=obs_date, MANUAL=m_obs,
                # FETCH JPL=fetch_date) — only obs_date was being updated
                # here, so loading a FITS image while on the FETCH JPL tab
                # (the common case — JPL-fetched elements + your own image)
                # left ITS date field looking untouched even though the
                # actual computation (see below) was already fixed. Update
                # all three that exist, regardless of which tab is active,
                # so switching tabs afterward still shows the right date
                # everywhere rather than just in whichever was visible at
                # load time.
                for attr in ("obs_date", "fetch_date", "m_obs"):
                    w = getattr(self, attr, None)
                    if w is not None:
                        w.setText(date_obs)
                # _emit_compute() reads obs_jd as `self._comet_el.get(
                # "obs_jd", <fallback from whichever tab's date field>)`,
                # and on_fetch_done() (run once, when the comet was
                # fetched) permanently bakes a fixed "obs_jd" key into
                # self._comet_el at that time. That baked-in key takes
                # priority over ANY of the text fields above
                # UNCONDITIONALLY — so updating the field text alone was
                # never enough; also overwrite comet_el's obs_jd here, so
                # a freshly loaded image's own date wins, the way you'd
                # expect, no matter which tab is currently selected.
                try:
                    new_obs_jd = cta.date_to_jd(date_obs)
                    if self._comet_el is not None:
                        self._comet_el["obs_jd"] = new_obs_jd
                except Exception:
                    pass   # malformed DATE-OBS string — leave obs_jd alone

            # Pre-fill object name if recognised. blockSignals is
            # important here (v3.1 bug fix) — without it, setCurrentIndex()
            # fires currentIndexChanged -> _on_comet_selected(), which does
            # self.obs_date.setText(meta.get("obs","")) and overwrote the
            # FITS file's own DATE-OBS we *just* set above with the DB's
            # static reference date instead — this was the actual, more
            # common trigger for the "date doesn't get pulled in" symptom
            # (it WAS pulled in, then silently clobbered a few lines later
            # by this auto-match, which is only meant to visually highlight
            # a matching preset, not to reset other fields as a side effect).
            obj = hdr.get("OBJECT","")
            if obj:
                # Try to match against DB
                self.combo_comet.blockSignals(True)
                try:
                    for i, name in enumerate(list(cta.COMET_DB.keys())):
                        if any(tok in name for tok in obj.split()):
                            self.combo_comet.setCurrentIndex(i)
                            break
                finally:
                    self.combo_comet.blockSignals(False)

            wcs_info = (f"WCS OK · {ps_deg*3600:.3f}\"/px · "
                        f"NPA={npa:.1f}° · nuc=({nx:.0f},{ny:.0f})")
            self.img_label.setText(
                f"✓  {os.path.basename(path)}  ({w}×{h})  |  {wcs_info}")
        else:
            self.nuc_x_spin.setValue(w/2)
            self.nuc_y_spin.setValue(h/2)
            self.img_label.setText(f"✓  {os.path.basename(path)}  ({w}×{h} px)")

        self.overlay_params.setVisible(True)
        self.btn_clear_img.setVisible(True)

        # Default behaviour for real image workflows: as soon as an
        # image is loaded, show its image-derived isophote contours.
        # Users can still turn this off from DISPLAY.
        try:
            self.chk_isophote.setChecked(True)
            self.isolvl_slider.setValue(6)
            self.isosm_slider.setValue(20)   # 2.0 px
        except Exception:
            pass

        self._reset_overlay_view()
        self.image_loaded.emit(self._img_arr)
        # Sync dialog
        if hasattr(self, '_img_dialog') and self._img_dialog:
            try:
                self._img_dialog.sync_label(self.img_label.text())
                if hasattr(self, '_wcs_ps_deg') and self._wcs_ps_deg:
                    self._img_dialog.ps_spin.setValue(round(self._wcs_ps_deg * 3600.0, 4))
                self._img_dialog._update_scale_preview()
            except Exception:
                pass

    def _load_fits(self, path: str):
        """Load a FITS file → (RGB uint8 ndarray, header). Saves raw data for restretch."""
        try:
            from astropy.io import fits as astrofits
            hdul = astrofits.open(path, memmap=False)
            data, hdr = None, None
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim >= 2:
                    data = hdu.data.copy(); hdr = hdu.header.copy(); break
            hdul.close()
            if data is None:
                QMessageBox.warning(self, "FITS Error", "No 2D image data found.")
                return None, None
            while data.ndim > 2:
                data = data[0]
            self._fits_raw = data.astype(np.float32)   # store for restretch
            return self._stretch_raw(self._fits_raw), hdr
        except Exception as ex:
            QMessageBox.critical(self, "FITS Load Error", f"Could not read FITS:\n{ex}")
            return None, None

    def _view_fits_header(self):
        """Open a read-only viewer showing the raw FITS header cards of
        the currently loaded file (v3.0.6). No-op with an explanatory
        message for non-FITS images (JPEG/PNG), which have no header."""
        if self._fits_header is None:
            QMessageBox.information(self, "No FITS Header",
                "The currently loaded image isn't a FITS file (or has no "
                "header) — JPEG/PNG images don't carry a FITS header.")
            return
        name = ""
        if hasattr(self, "_last_image_path") and self._last_image_path:
            name = os.path.basename(self._last_image_path)
        dlg = FitsHeaderDialog(self, self._fits_header, name)
        dlg.exec()

    def _stretch_raw(self, data: np.ndarray) -> np.ndarray:
        """
        Apply the current Contrast/Intensity stretch to raw float32 FITS
        data → RGB uint8.

        v3.0: simplified from 4 technical sliders (Shadow%/Highlight%/
        Softness/Gamma) to 2 simple ones (Contrast, Intensity, both
        -1000..1000), inspired by the simpler two-knob stretch UI in
        comet-tracking viewers like Tycho Tracker. This is Claude's own
        reparameterization of the same underlying asinh stretch, NOT a
        reverse-engineered match of any other tool's exact internal
        formula (no access to that) — tune the formula below if the
        defaults don't feel right for your images.

        Mapping:
          Background level = median; noise scale = robust σ (1.4826×MAD).
          INTENSITY shifts the black point relative to the background:
            negative → black point pulled well below background, reveals
            faint signal/coma (image looks brighter/more "exposed");
            positive → black point pushed above background, only bright
            cores survive. ±1000 spans roughly ∓15σ.
          CONTRAST sets how wide the stretched range above black is:
            higher → narrower range → steeper, more contrasty; lower/
            negative → wider range → flatter, more gradual. ±1000 spans
            a ×0.1 to ×10 multiplier on a fixed 40σ base width.
          A fixed mild asinh softness (0.05) is kept internally — not
          user-exposed, since the 2-knob goal was simplicity.
        """
        contrast  = self.sl_contrast.value()
        intensity = self.sl_intensity.value()

        med = float(np.median(data))
        mad = float(np.median(np.abs(data - med)))
        sigma = mad * 1.4826 + 1e-6

        black = med + (intensity / 1000.0) * 15.0 * sigma
        range_factor = 10.0 ** (-contrast / 1000.0)   # +1000→×0.1, -1000→×10
        white = black + 40.0 * sigma * range_factor
        if white <= black:
            white = black + 1e-6

        d = np.clip(data, black, white)
        d = (d - black) / (white - black)
        d = np.arcsinh(d / 0.05) / np.arcsinh(1.0 / 0.05)
        d = np.clip(d, 0.0, 1.0)
        img8 = (d * 255).astype(np.uint8)
        from PIL import Image as PILImage
        return np.array(PILImage.fromarray(img8, "L").convert("RGB"))

    def _auto_stretch(self):
        """
        Apply a reasonable default Contrast/Intensity for comet images:
        push Intensity negative to reveal faint coma/tail above the
        background, with Contrast also negative for a wider/flatter
        stretch range (avoids clipping the faint signal Intensity just
        revealed). These are starting points, not a derived optimum —
        adjust by eye afterward.
        """
        if not hasattr(self, "_fits_raw") or self._fits_raw is None:
            QMessageBox.information(self, "Auto Stretch",
                "No FITS image loaded. Auto stretch only works with FITS files.")
            return

        self.sl_contrast.setValue(-150)
        self.sl_intensity.setValue(-300)

        self._restretch()

    def _restretch(self):
        """Re-apply stretch using current slider values and redraw."""
        if not hasattr(self, "_fits_raw") or self._fits_raw is None:
            return
        self._img_arr = self._stretch_raw(self._fits_raw)
        self.image_loaded.emit(self._img_arr)

    def _update_npa_hint(self):
        npa = self.npa_spin.value()
        if abs(npa) < 1:
            hint = "North ↑  East ←  (standard)"
        else:
            hint = f"North is {abs(npa):.1f}° {'CW' if npa>0 else 'CCW'} from image top"
        self.npa_hint.setText(hint)

    def _rotate_north_up(self):
        """
        Rotate image so North is up and East is left (standard astronomical orientation).

        scipy.ndimage.rotate with axes=(0,1) and POSITIVE angle θ rotates the image
        content CCW in screen view (y-down convention):
          - A pixel to the RIGHT  (North at NPA≈90°) → moves to the TOP ✓
          - A pixel at the TOP    (East when NPA≈90°) → moves to the LEFT ✓

        Therefore angle = +NPA  (NOT  -NPA).

        Nucleus coordinate forward transform under scipy (axes=(0,1), angle=+θ):
          col_out = sin(θ)·row_in + cos(θ)·col_in
          row_out = cos(θ)·row_in − sin(θ)·col_in
        """
        if self._img_arr is None:
            return
        npa = self.npa_spin.value()
        if abs(npa) < 0.5:
            self.rot_status.setText("✓ Already North-up  (NPA < 0.5°)")
            return
        try:
            from scipy.ndimage import rotate as nd_rotate

            # ── CORRECT sign: +npa for N-up E-left ──────────────────────
            angle = npa          # positive → CCW in screen → Right→Top ✓
            h, w  = self._img_arr.shape[:2]

            rotated = nd_rotate(self._img_arr, angle, axes=(0,1),
                                reshape=True, order=1, mode="constant", cval=0)

            # ── Correct nucleus position (scipy forward transform) ────────
            cx_old, cy_old = w / 2.0, h / 2.0
            col_in = self.nuc_x_spin.value() - cx_old   # rightward offset
            row_in = self.nuc_y_spin.value() - cy_old   # downward offset
            th = np.radians(angle)
            col_out = np.sin(th) * row_in + np.cos(th) * col_in
            row_out = np.cos(th) * row_in - np.sin(th) * col_in

            h2, w2 = rotated.shape[:2]
            self.nuc_x_spin.setValue(round(col_out + w2 / 2.0, 1))
            self.nuc_y_spin.setValue(round(row_out + h2 / 2.0, 1))

            # ── Rotate raw FITS data too (for restretch) ─────────────────
            if hasattr(self, "_fits_raw") and self._fits_raw is not None:
                self._fits_raw = nd_rotate(self._fits_raw, angle, axes=(0,1),
                                           reshape=True, order=1,
                                           mode="constant", cval=0)
                rotated = self._stretch_raw(self._fits_raw)

            self._img_arr = rotated

            if hasattr(self, "_wcs_npa"):
                self._wcs_npa = 0.0
            self.npa_spin.setValue(0.0)

            dir_str = "CCW" if npa > 0 else "CW"
            self.rot_status.setText(
                f"✓ Rotated {abs(npa):.2f}° {dir_str}\n"
                f"North ↑  East ←  (N-up / E-left)\n"
                f"New size: {w2}×{h2} px\n"
                f"Nucleus: ({col_out+w2/2:.0f}, {row_out+h2/2:.0f}) px")
            self.image_loaded.emit(self._img_arr)

        except Exception as ex:
            self.rot_status.setText(f"✗ Rotation failed: {ex}")
        # Rotation changes the image's pixel dimensions — any previous
        # view lock was in the OLD image's pixel coordinates and is now
        # meaningless, so reset to the full (newly rotated) image.
        self._reset_overlay_view()

    def _reset_overlay_view(self):
        """Snap the overlay back to showing the full image — call
        whenever the loaded image itself changes (new file, cleared, or
        rotated), since whatever pan/zoom was set before no longer
        corresponds to the new image's dimensions, and also wired
        directly to the RESET TO FULL IMAGE button. Nothing to "clear" as
        separate state — draw_model() reads the view directly off the
        axes each time, so setting it here is itself enough; the next
        Compute Model run will read this back as "the previous view"."""
        mw = self.window()
        canvas = getattr(mw, "canvas", None) if mw else None
        if canvas is None or self._img_arr is None:
            return
        h_img, w_img = self._img_arr.shape[:2]
        canvas.ax.set_xlim(0, w_img)
        canvas.ax.set_ylim(h_img, 0)
        canvas.canvas.draw_idle()

    def _clear_image(self):
        self._img_arr     = None
        self._fits_raw    = None
        self._fits_header = None
        self._wcs_ps_deg  = None
        self._wcs_npa     = None
        self.img_label.setText("No image loaded")
        # overlay_params and btn_clear_img are stubs in new layout — safe to call setVisible
        try: self.overlay_params.setVisible(False)
        except Exception: pass
        try: self.btn_clear_img.setVisible(False)
        except Exception: pass
        # Reset image-derived display state.
        try: self.chk_isophote.setChecked(False)
        except Exception: pass

        # Reset ephemeris override checkbox
        try: self.chk_use_obs.setChecked(False)
        except Exception: pass
        # Update dialog label if open
        if hasattr(self, '_img_dialog') and self._img_dialog is not None:
            try: self._img_dialog.sync_label("No image loaded")
            except Exception: pass
        self.image_loaded.emit(None)

    # ── Compute ───────────────────────────────────────────────────────────
    def _emit_compute(self):
        if self._comet_el is None:
            QMessageBox.warning(self, "No Comet",
                                "Select or fetch a comet first.")
            return
        betas = self._parse_floats(self.beta_str.text())
        ages  = self._parse_ints(self.age_str.text())
        self._sync_max_age_from_synchrone_text()
        if not betas or not ages:
            QMessageBox.warning(self, "Bad Input",
                                "Enter at least one valid β value and synchrone age.")
            return
        # obs_jd: prefer whatever's already baked into _comet_el (set by
        # on_fetch_done()/_use_preset()/_use_manual(), and now also kept
        # current by _load_image()'s DATE-OBS fix — see that method).
        # Fallback below only matters if that key is somehow absent; it
        # now checks EACH tab's own date field (previously only tab 0's,
        # so tabs 1/2 always silently fell back to today() even when
        # their own field had a perfectly good date in it — v3.1 fix).
        obs_jd = self._comet_el.get("obs_jd")
        if obs_jd is None:
            idx = self.comet_tabs.currentIndex()
            date_text = (self.obs_date.text()   if idx == 0 else
                        self.m_obs.text()       if idx == 1 else
                        self.fetch_date.text()  if idx == 2 else "")
            try:
                obs_jd = cta.date_to_jd(date_text) if date_text else cta.today_jd()
            except Exception:
                obs_jd = cta.today_jd()
        self.compute_requested.emit(
            self._comet_el, obs_jd, betas, ages,
            self.max_age.value(), self.n_pts.value(),
            self.get_ejection_params())

    def get_rotation_offset(self) -> float:
        """Return grid rotation offset in degrees (positive = CW = more toward West)."""
        return self.rot_slider.value() / 10.0

    def get_ejection_params(self) -> dict:
        """Return the v3.1 ejection-velocity dict for compute_model()/
        dust_position(). Set EXCLUSIVELY via set_ejection_params() below,
        called from the Monte Carlo window's "Run Monte Carlo" — there is
        no longer a separate input for this in the main panel (see the
        comment above grp_model). All-zero (the initial default)
        reproduces the v3.0 zero-ejection-velocity model exactly."""
        return dict(self._current_ejection)

    def set_ejection_params(self, d: dict):
        """Called from MCWindow when its Run Monte Carlo button is
        clicked, so the main F-P model picks up the SAME ejection
        velocity the Monte Carlo population is using, without it being
        keyed in twice in two different places."""
        self._current_ejection = dict(d)

    def reset_ejection_params(self):
        """Return the hidden F-P ejection state to classical F-P defaults.

        MC Run synchronizes its velocity law into this internal state.  Because
        those values are intentionally not duplicated as visible controls on
        the main panel, preserving them after CLEAR ALL MODELS would leave an
        invisible, stale assumption in the next F-P run.
        """
        self._current_ejection = dict(
            v_R0=0.0, v_T0=0.0, v_N0=0.0, gamma=0.0, m_exp=0.0)

    def reset_mc_presentation(self):
        """Return the main canvas to the normal analysis presentation."""
        self._mc_style.update(dict(
            lw=1.2, alpha=0.95, info_box=True,
            view_style="analysis", background="image", coordinates="pixels",
            observed_show=True, observed_color=ISOPHOTE_COLOR,
            observed_lw=0.6, observed_alpha=0.75, observed_ls="-",
            model_color="#ff3399", model_ls="--",
            show_legend=True, show_title=True, show_compass=True,
            grayscale_safe=False))

    def get_vis(self):
        mc = dict(self._mc_style)
        return dict(
            synd=self.chk_synd.isChecked(),
            sync=self.chk_sync.isChecked(),
            orbit=self.chk_orbit.isChecked(),
            isophote=self.chk_isophote.isChecked(),
            isophote_levels=self.isolvl_slider.value(),
            isophote_smooth=self.isosm_slider.value()/10.0,
            lw=self.lw_slider.value()/10,
            alpha=self.op_slider.value()/100,
            orbit_lw=self.orbit_lw_slider.value()/10,
            orbit_alpha=self.orbit_op_slider.value()/100,
            mc_contour=self._mc_contour_visible,
            mc_lw=mc.get("lw", 1.2),
            mc_alpha=mc.get("alpha", 0.95),
            mc_info_box=mc.get("info_box", True),
            mc_view_style=mc.get("view_style", "analysis"),
            mc_background=mc.get("background", "image"),
            mc_coordinates=mc.get("coordinates", "pixels"),
            observed_show=mc.get("observed_show", True),
            observed_color=mc.get("observed_color", ISOPHOTE_COLOR),
            observed_lw=mc.get("observed_lw", 0.6),
            observed_alpha=mc.get("observed_alpha", 0.75),
            observed_ls=mc.get("observed_ls", "-"),
            model_color=mc.get("model_color", "#ff3399"),
            model_ls=mc.get("model_ls", "--"),
            show_legend=mc.get("show_legend", True),
            show_title=mc.get("show_title", True),
            show_compass=mc.get("show_compass", True),
            grayscale_safe=mc.get("grayscale_safe", False),
            show_sun=self.chk_sun.isChecked(),
            show_antivel=self.chk_antivel.isChecked(),
            show_crosshair=self.chk_crosshair.isChecked(),
        )

    def get_overlay(self):
        return dict(
            img_arr=self._img_arr,
            nuc_x=self.nuc_x_spin.value(),
            nuc_y=self.nuc_y_spin.value(),
            au_per_px=self.au_px_spin.value(),
            north_pa=self.npa_spin.value(),
        )

    def set_computing(self, val: bool):
        self.btn_compute.setEnabled(not val)
        self.progress_bar.setVisible(val)
        if val: self.progress_bar.setValue(0)

# ─────────────────────────────────────────────────────────────────────────────
#  RIGHT INFO PANEL  — tabbed to prevent overflow
# ─────────────────────────────────────────────────────────────────────────────
class InfoPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(240)
        self.setMaximumWidth(300)

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        # Thin header strip (stored so update_theme() can re-style it)
        self.hdr = QLabel("  INFO")
        self.hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_muted']}; font-size:10px;"
            f"letter-spacing:2px; font-weight:bold; padding:6px 8px;"
            f"border-left:1px solid {T['border']};"
            f"border-bottom:1px solid {T['border']};")
        vbox.addWidget(self.hdr)

        self.tabs = QTabWidget()
        self.tabs.setTabBarAutoHide(True)
        self.tabs.setStyleSheet(
            f"QTabWidget::pane {{ border-left:1px solid {T['border']}; }}")
        vbox.addWidget(self.tabs)

        # ── Tab 1 : INFO (Ephemeris + Orbital Elements) ──────────────────
        sa1 = QScrollArea(); sa1.setWidgetResizable(True)
        sa1.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        w1 = QWidget(); sa1.setWidget(w1)
        v1 = QVBoxLayout(w1); v1.setContentsMargins(8,8,8,8); v1.setSpacing(6)

        grp_eph = QGroupBox("EPHEMERIS")
        ge = QGridLayout(grp_eph); ge.setSpacing(5); ge.setColumnStretch(1,1)
        self.eph_labels = {}
        self._eph_key_lbs = []
        for i,(k,tip) in enumerate([("r☉","heliocentric dist."),("Δ","geocentric dist."),
                ("Phase","solar phase angle"),("RA","right ascension J2000"),
                ("Dec","declination J2000"),("Date","observation (UT)")]):
            lk = QLabel(k); lk.setToolTip(tip)
            lv = QLabel("—"); lv.setWordWrap(True)
            self._eph_key_lbs.append(lk)
            ge.addWidget(lk,i,0); ge.addWidget(lv,i,1); self.eph_labels[k] = lv
        v1.addWidget(grp_eph)

        grp_orb = QGroupBox("ORBITAL ELEMENTS")
        go = QGridLayout(grp_orb); go.setSpacing(5); go.setColumnStretch(1,1)
        self.orb_labels = {}
        self._orb_key_lbs = []
        for i,(k,tip) in enumerate([("q","Perihelion dist. (AU)"),("e","Eccentricity"),
                ("i","Inclination (°)"),("Ω","Long. ascending node (°)"),
                ("ω","Arg. of perihelion (°)"),("T","Perihelion date"),("Source","Data source")]):
            lk = QLabel(k); lk.setToolTip(tip)
            lv = QLabel("—"); lv.setWordWrap(True)
            self._orb_key_lbs.append(lk)
            go.addWidget(lk,i,0); go.addWidget(lv,i,1); self.orb_labels[k] = lv
        v1.addWidget(grp_orb)

        # ── ANIMATOR (v3.0: embedded in main window, not a popup) ────────
        # Steps obs date across a range and plays back the resulting F-P
        # model frames directly into the MAIN canvas — for spotting when
        # the tail is widest/changing fastest (e.g. to plan imaging
        # sessions). Orchestration (worker thread, frame cache, play
        # timer, drawing) lives in MainWindow; this group box only holds
        # the input widgets.
        grp_anim = QGroupBox("ANIMATOR")
        av = QVBoxLayout(grp_anim); av.setSpacing(6)

        af = QFormLayout(); af.setSpacing(5)
        self.anim_start = QLineEdit(cta.jd_to_str(cta.today_jd())[:10])
        self.anim_end   = QLineEdit(cta.jd_to_str(cta.today_jd() + 360)[:10])
        self.anim_step  = QDoubleSpinBox()
        self.anim_step.setRange(0.1, 100.0)
        self.anim_step.setDecimals(1)
        self.anim_step.setValue(2.0)
        self.anim_step.setSuffix(" d")
        af.addRow("Start:", self.anim_start)
        af.addRow("End:", self.anim_end)
        af.addRow("Step:", self.anim_step)
        av.addWidget(QLabel("Date range (YYYY-MM-DD):"))
        av.addLayout(af)

        size_hdr = QHBoxLayout()
        self.anim_rb_fov  = QRadioButton("Fixed FOV")
        self.anim_rb_dist = QRadioButton("Fixed dist.")
        self.anim_rb_fov.setChecked(True)
        self.anim_rb_fov.setToolTip(
            "Locks the angular field of view (arcmin) — what you'd see "
            "through a fixed camera/eyepiece FOV every night, including "
            "the real effect of changing Earth-comet distance Δ.")
        self.anim_rb_dist.setToolTip(
            "Locks the physical size (AU) shown on screen, removing the "
            "effect of changing Δ — for comparing the tail's actual "
            "physical growth.")
        size_hdr.addWidget(self.anim_rb_fov)
        size_hdr.addWidget(self.anim_rb_dist)
        av.addLayout(size_hdr)

        size_row = QHBoxLayout()
        self.anim_size = QDoubleSpinBox()
        self.anim_size.setRange(0.01, 100000.0)
        self.anim_size.setDecimals(3)
        self.anim_size.setValue(600.0)
        self.anim_size.setSuffix(" arcmin")
        self.anim_size.setMaximumWidth(130)
        self.anim_btn_auto = QPushButton("Auto")
        self.anim_btn_auto.setMaximumWidth(52)
        self.anim_btn_auto.setToolTip(
            "Samples the date range and suggests a size that fits the "
            "widest tail with ~20% margin.")
        size_row.addWidget(QLabel("Full width:"))
        size_row.addWidget(self.anim_size)
        size_row.addWidget(self.anim_btn_auto)
        av.addLayout(size_row)
        self.anim_rb_fov.toggled.connect(
            lambda c: self.anim_size.setSuffix(" arcmin" if c else " AU"))

        self.anim_btn_compute = QPushButton("▶  COMPUTE FRAMES")
        self.anim_btn_compute.setProperty("class", "primary")
        av.addWidget(self.anim_btn_compute)
        self.anim_progress = QProgressBar()
        self.anim_progress.setTextVisible(False)
        self.anim_progress.setMaximumHeight(8)
        av.addWidget(self.anim_progress)

        play_row = QHBoxLayout()
        self.anim_btn_play = QPushButton("▶")
        self.anim_btn_play.setEnabled(False)
        self.anim_btn_play.setMaximumWidth(36)
        self.anim_slider = QSlider(Qt.Orientation.Horizontal)
        self.anim_slider.setEnabled(False)
        play_row.addWidget(self.anim_btn_play)
        play_row.addWidget(self.anim_slider, 1)
        av.addLayout(play_row)

        self.anim_lbl_frame = QLabel("No frames computed yet.")
        self.anim_lbl_frame.setWordWrap(True)
        self.anim_lbl_frame.setStyleSheet("font-size:10px;")
        av.addWidget(self.anim_lbl_frame)

        v1.addWidget(grp_anim)
        v1.addStretch()
        self._apply_info_label_styles()   # initial colour pass
        self.tabs.addTab(sa1, "INFO")

        # v3.0: β TABLE and LIGHT CURVE tabs removed — both are now fully
        # covered elsewhere with no loss of function: grain radius lives
        # in Calculation > Dust particle radius… (self-contained, doesn't
        # need a per-comet model), and the light curve (plot + H0/n) is
        # shown directly by the View > Light curve… popup window, which
        # already fetches COBS automatically. Keeping a thin duplicate
        # status display here risked it silently drifting out of sync
        # with whichever the user actually looked at.

        # v3.0: ANALYSIS tab removed. Its Afρ-based dust-production
        # report duplicated the standalone Calculation > Dust production
        # rate… calculator added earlier in v3.0 — same formula
        # (compute_dust_production_rate(), still in comet_tail_analyzer.py
        # if a text-report version is wanted again later), just shown a
        # second way. generate_dust_analysis() itself is left in place as
        # a library function; only this GUI tab and its caller
        # (_run_analysis/_menu_generate_analysis/the Calculation menu's
        # "Generate analysis" item) were removed.

    def update_info(self, model, comet_el):
        info = model["info"]
        self.eph_labels["r☉"].setText(f"{info['r_helio']:.5f} AU")
        self.eph_labels["Δ"].setText(f"{info['r_geo']:.5f} AU")
        self.eph_labels["Phase"].setText(f"{info['phase_angle']:.2f}°")
        self.eph_labels["RA"].setText(f"{info['RA']:.4f}°")
        self.eph_labels["Dec"].setText(f"{info['Dec']:.4f}°")
        self.eph_labels["Date"].setText(info["obs_str"][:16])

        self.orb_labels["q"].setText(f"{comet_el['q']:.6f} AU")
        self.orb_labels["e"].setText(f"{comet_el['e']:.7f}")
        self.orb_labels["i"].setText(f"{comet_el['i']:.4f}°")
        self.orb_labels["Ω"].setText(f"{comet_el['Omega']:.4f}°")
        self.orb_labels["ω"].setText(f"{comet_el['omega']:.4f}°")
        self.orb_labels["T"].setText(comet_el.get("T","")[:10])
        self.orb_labels["Source"].setText(comet_el.get("source", "JPL Horizons"))
        self.tabs.setCurrentIndex(0)

    # ── Theme helpers ──────────────────────────────────────────────────────
    def _apply_info_label_styles(self):
        """Re-apply eph / orb key-value label colours from the active theme."""
        for lk in getattr(self, '_eph_key_lbs', []):
            lk.setStyleSheet(f"color:{T['text_muted']}; font-size:11px;")
        for lv in self.eph_labels.values():
            lv.setStyleSheet(f"color:{T['text_value']}; font-weight:bold; font-size:12px;")
        for lk in getattr(self, '_orb_key_lbs', []):
            lk.setStyleSheet(f"color:{T['text_muted']}; font-size:11px;")
        for lv in self.orb_labels.values():
            lv.setStyleSheet(f"color:{T['text_value']}; font-weight:bold; font-size:11px;")

    def _apply_misc_label_styles(self):
        """Re-apply the header strip colour from the active theme."""
        # Header strip
        self.hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_muted']}; font-size:10px;"
            f"letter-spacing:2px; font-weight:bold; padding:6px 8px;"
            f"border-bottom:1px solid {T['border']};")

    def update_theme(self):
        """Re-apply ALL theme-sensitive styles in the right panel."""
        self._apply_info_label_styles()
        self._apply_misc_label_styles()



# ─────────────────────────────────────────────────────────────────────────────
#  FITS HEADER VIEWER  (v3.0.6) — read-only display of every raw header
#  card in the currently loaded FITS file, opened from the Image Setup
#  dialog's "VIEW FITS HEADER" button.
# ─────────────────────────────────────────────────────────────────────────────
class FitsHeaderDialog(QDialog):
    def __init__(self, parent, header, filename: str = ""):
        super().__init__(parent)
        title = f"FITS Header — {filename}" if filename else "FITS Header"
        self.setWindowTitle(title)
        self.setMinimumSize(560, 600)
        v = QVBoxLayout(self)

        info = QLabel(f"{len(header)} card(s)" +
                      (f"  ·  {filename}" if filename else ""))
        info.setStyleSheet("color:#3a6070; font-size:10px;")
        v.addWidget(info)

        text = QTextEdit()
        text.setReadOnly(True)
        text.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        text.setStyleSheet(
            "font-family:monospace; font-size:11px;"
            "background:#0a0f1a; color:#c0d0e0; border:1px solid #1a2540;")
        # tostring() renders the standard 80-column FITS card format, one
        # card per line — the exact raw keyword/value/comment content,
        # not a re-summarized version of it.
        try:
            body = header.tostring(sep='\n', padding=False)
        except Exception:
            body = "\n".join(str(card) for card in header.cards)
        text.setPlainText(body)
        v.addWidget(text)

        btns = QHBoxLayout()
        btn_copy = QPushButton("📋  COPY")
        btn_copy.clicked.connect(
            lambda: QApplication.clipboard().setText(text.toPlainText()))
        btns.addWidget(btn_copy)
        btns.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btns.addWidget(btn_close)
        v.addLayout(btns)


# ─────────────────────────────────────────────────────────────────────────────
#  IMAGE SETUP DIALOG — floating window for all image calibration
# ─────────────────────────────────────────────────────────────────────────────
class ImageSetupDialog(QDialog):
    """
    Floating non-modal dialog for:
      - Open / Clear image
      - Nucleus position (manual + click)
      - Scale & Orientation (AU/px, North PA)
      - Rotate N-up / E-left
      - FITS stretch (Contrast, Intensity)
      - Observed Ephemeris Override (PsAng, r☉, Δ, Phase)
    """
    def __init__(self, ctrl: "ControlPanel"):
        super().__init__(ctrl.window() if ctrl.window() else ctrl,
                         Qt.WindowType.Dialog |
                         Qt.WindowType.WindowCloseButtonHint |
                         Qt.WindowType.WindowTitleHint)
        # Keep this as a modeless child dialog of MainWindow.  A parented
        # QDialog remains above its CTA main window, but without
        # WindowStaysOnTopHint it no longer stays above unrelated
        # applications or every other desktop window.
        # Qt.WindowType.Dialog also provides a normal draggable title bar.
        self.ctrl = ctrl
        self._sync_in_progress = False  # Flag to prevent infinite loops in ps_spin ↔ au_px_spin sync
        self.setWindowTitle("Image Setup & Calibration")
        self.setMinimumWidth(420)
        self.setMaximumWidth(480)
        self.resize(440, 820)

        root = QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(8)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        inner = QWidget()
        scroll.setWidget(inner)
        vb = QVBoxLayout(inner)
        vb.setContentsMargins(6, 6, 6, 6)
        vb.setSpacing(8)
        root.addWidget(scroll, 1)

        def sec(text):
            l = QLabel(text)
            l.setStyleSheet("color:#6ab0d8;font-size:10px;letter-spacing:2px;"
                            "font-weight:bold;padding:2px 0;border-bottom:1px solid #1a2540;")
            return l

        def inp_row(label, widget, unit=""):
            row = QHBoxLayout()
            lbl = QLabel(label); lbl.setMinimumWidth(90)
            lbl.setStyleSheet("color:#7aA8c8;font-size:11px;")
            row.addWidget(lbl); row.addWidget(widget, 1)
            if unit:
                ul = QLabel(unit); ul.setStyleSheet("color:#3a6070;font-size:10px;")
                row.addWidget(ul)
            return row

        # ── File ─────────────────────────────────────────────────────────
        vb.addWidget(sec("📂  IMAGE FILE"))
        fb = QHBoxLayout()
        self.btn_open  = QPushButton("OPEN IMAGE…")
        self.btn_open.clicked.connect(self._on_open)
        self.btn_clear = QPushButton("CLEAR")
        self.btn_clear.setProperty("class","danger")
        self.btn_clear.clicked.connect(self._on_clear)
        fb.addWidget(self.btn_open, 3); fb.addWidget(self.btn_clear, 1)
        vb.addLayout(fb)
        self.btn_view_header = QPushButton("📄  VIEW FITS HEADER")
        self.btn_view_header.setToolTip(
            "Shows every raw header card of the loaded FITS file. "
            "No-op for JPEG/PNG (no FITS header to show).")
        self.btn_view_header.clicked.connect(self.ctrl._view_fits_header)
        vb.addWidget(self.btn_view_header)
        self.file_lbl = QLabel("No image loaded")
        self.file_lbl.setStyleSheet("color:#5a90b0;font-size:10px;")
        self.file_lbl.setWordWrap(True)
        vb.addWidget(self.file_lbl)

        # ── Nucleus ───────────────────────────────────────────────────────
        # ── Scale & Orientation ───────────────────────────────────────────
        vb.addWidget(sec("📐  SCALE & ORIENTATION"))

        # ── Method A: enter arcsec/px directly (most intuitive) ──────────
        scl_hint = QLabel(
            "Enter plate scale in arcsec/px  (from telescope specs or FITS header)\n"
            "AU/px is computed automatically using geocentric distance Δ.")
        scl_hint.setStyleSheet("color:#3a6070;font-size:10px;"); scl_hint.setWordWrap(True)
        vb.addWidget(scl_hint)

        self.ps_spin = QDoubleSpinBox()
        self.ps_spin.setDecimals(4); self.ps_spin.setRange(0.01, 500.0)
        self.ps_spin.setValue(1.441); self.ps_spin.setSingleStep(0.1)
        self.ps_spin.setSuffix("  \"/px")
        self.ps_spin.setToolTip(
            "Plate scale in arc-seconds per pixel.\n"
            "Typical values:\n"
            "  Phone/DSLR on wide lens : 50–200\"/px\n"
            "  DSLR on 200mm tele      : 10–30\"/px\n"
            "  DSLR on 1000mm refractor:  2–5\"/px\n"
            "  CCD on dedicated SCT/RC :  0.5–2\"/px\n"
            "  Same telescope as FITS  :  1.441\"/px")
        vb.addLayout(inp_row("Plate scale", self.ps_spin))

        # FOV readout
        self.fov_lbl = QLabel("")
        self.fov_lbl.setStyleSheet("color:#3a6070;font-size:10px;")
        self.fov_lbl.setWordWrap(True)
        vb.addWidget(self.fov_lbl)

        # AU/px readout (computed, semi-read-only but still editable)
        vb.addLayout(inp_row("AU / pixel", ctrl.au_px_spin))
        self.aupx_lbl = QLabel("")
        self.aupx_lbl.setStyleSheet("color:#3a6070;font-size:10px;")
        vb.addWidget(self.aupx_lbl)

        # Apply button
        self.btn_apply_scale = QPushButton("▶  APPLY PLATE SCALE → COMPUTE AU/PX")
        self.btn_apply_scale.setProperty("class","primary")
        self.btn_apply_scale.clicked.connect(self._apply_plate_scale)
        vb.addWidget(self.btn_apply_scale)
        
        # Scale result label
        self.scale_result_lbl = QLabel("")
        self.scale_result_lbl.setStyleSheet("color:#4a9060;font-size:10px;")
        self.scale_result_lbl.setWordWrap(True)
        vb.addWidget(self.scale_result_lbl)

        # ── North PA (moved from IMAGE ROTATION section) ──────────────────
        vb.addLayout(inp_row("North PA (°)", ctrl.npa_spin, "° CW from top"))
        ctrl.npa_hint.setStyleSheet("color:#4a8080;font-size:9px;")
        vb.addWidget(ctrl.npa_hint)
        vb.addWidget(ctrl.btn_rot_nup)
        ctrl.rot_status.setStyleSheet("color:#4a8060;font-size:9px;")
        ctrl.rot_status.setWordWrap(True)
        vb.addWidget(ctrl.rot_status)

        # Connect signals
        self.ps_spin.valueChanged.connect(self._update_scale_preview)
        self._update_scale_preview()   # initial fill
        
        # Auto-sync: ps_spin (arcsec/px) ↔ au_px_spin (AU/px)
        # When user changes ps_spin, update au_px_spin and vice versa
        self.ps_spin.valueChanged.connect(self._on_ps_changed)
        self.ctrl.au_px_spin.valueChanged.connect(self._on_au_px_changed)

        # ── Stretch ───────────────────────────────────────────────────────
        vb.addWidget(sec("🎚  IMAGE STRETCH  (FITS only)"))

        def sl_row(label, slider, lo, hi, val, scale, fmt):
            sl = slider
            sl.setRange(lo, hi); sl.setValue(val)
            vl = QLabel(fmt.format(val/scale))
            vl.setMinimumWidth(50)
            vl.setStyleSheet("color:#90d8ff;font-size:11px;")
            sl.valueChanged.connect(lambda v, lv=vl, sc=scale, f=fmt: lv.setText(f.format(v/sc)))
            r = QHBoxLayout()
            lbl = QLabel(label); lbl.setMinimumWidth(90)
            lbl.setStyleSheet("color:#7aA8c8;font-size:11px;")
            r.addWidget(lbl); r.addWidget(sl, 1); r.addWidget(vl)
            return r

        vb.addLayout(sl_row("Contrast",  ctrl.sl_contrast,  -1000, 1000, 0, 1, "{:.0f}"))
        vb.addLayout(sl_row("Intensity", ctrl.sl_intensity, -1000, 1000, 0, 1, "{:.0f}"))
        
        # Auto Stretch + Apply Stretch buttons
        vb.addWidget(ctrl.btn_auto_stretch)

        # ── Observed Ephemeris Override ───────────────────────────────────
        # ── Nucleus Position (moved to bottom) ────────────────────────────
        vb.addWidget(sec("⊕  NUCLEUS POSITION"))
        nuc_hint = QLabel("Click on the comet nucleus in the plot to set its position,\n"
                          "or enter X,Y pixel coordinates manually.")
        nuc_hint.setStyleSheet("color:#3a6070;font-size:10px;"); nuc_hint.setWordWrap(True)
        vb.addWidget(nuc_hint)
        nrow = QHBoxLayout()
        nrow.addWidget(QLabel("X (px)")); nrow.addWidget(ctrl.nuc_x_spin, 1)
        nrow.addWidget(QLabel("Y (px)")); nrow.addWidget(ctrl.nuc_y_spin, 1)
        vb.addLayout(nrow)
        vb.addWidget(ctrl.btn_pick_nuc)

        vb.addWidget(ctrl.btn_reset_view)
        zoom_hint = QLabel(
            "Use the plot toolbar's own zoom (rubber-band) and pan (hand) "
            "tools directly on the image above to frame it however you "
            "want. Whatever view that leaves you at is locked automatically "
            "and used for the overlay every time Compute Model re-runs, "
            "instead of resetting to the full image. Click RESET TO FULL "
            "IMAGE to clear the lock.")
        zoom_hint.setStyleSheet("color:#3a6070;font-size:10px;")
        zoom_hint.setWordWrap(True)
        vb.addWidget(zoom_hint)

        vb.addStretch()

        # ── Close button ──────────────────────────────────────────────────
        close_btn = QPushButton("CLOSE  (changes auto-saved)")
        close_btn.clicked.connect(self.hide)
        root.addWidget(close_btn)

    def _update_scale_preview(self):
        """Live preview: show FOV and AU/px estimate from current arcsec/px."""
        ps = self.ps_spin.value()
        # Try to get r_geo from main window model
        r_geo = self._get_r_geo()
        au_px = ps * (np.pi / 180.0 / 3600.0) * r_geo

        img  = self.ctrl._img_arr
        if img is not None:
            h, w = img.shape[:2]
            fov_x = w * ps / 3600.0
            fov_y = h * ps / 3600.0
            self.fov_lbl.setText(
                f"FOV: {fov_x*60:.1f}′ × {fov_y*60:.1f}′  "
                f"({w}×{h} px  @  Δ={r_geo:.4f} AU)\n"
                f"AU/px preview: {au_px:.3e}")
        else:
            self.fov_lbl.setText(f"AU/px preview: {au_px:.3e}  (Δ={r_geo:.4f} AU)")

        self.aupx_lbl.setText(
            f"→ AU/px = {ps:.4f}\" × π/180/3600 × {r_geo:.4f} AU = {au_px:.4e}")

    def _on_ps_changed(self, ps_arcsec):
        """
        Auto-sync handler: when ps_spin changes, update au_px_spin.
        Converts arcsec/px → AU/px using current geocentric distance.
        """
        if self._sync_in_progress:
            return  # Prevent infinite loop
        self._sync_in_progress = True
        try:
            r_geo = self._get_r_geo()
            if ps_arcsec > 0 and r_geo > 0:
                # arcsec/px → rad/px → AU/px
                au_per_px = ps_arcsec * (np.pi / 180.0 / 3600.0) * r_geo
                # Only update if different (avoid rounding loop)
                current = self.ctrl.au_px_spin.value()
                if abs(current - au_per_px) > 1e-10:
                    self.ctrl.au_px_spin.setValue(round(au_per_px, 9))
        finally:
            self._sync_in_progress = False

    def _on_au_px_changed(self, au_per_px):
        """
        Auto-sync handler: when au_px_spin changes, update ps_spin.
        Converts AU/px → arcsec/px using current geocentric distance.
        """
        if self._sync_in_progress:
            return  # Prevent infinite loop
        self._sync_in_progress = True
        try:
            r_geo = self._get_r_geo()
            if au_per_px > 0 and r_geo > 0:
                # AU/px → rad/px → arcsec/px
                ps_arcsec = au_per_px / r_geo * (180.0 * 3600.0 / np.pi)
                # Only update if different
                current = self.ps_spin.value()
                if abs(current - ps_arcsec) > 1e-6:
                    self.ps_spin.setValue(round(ps_arcsec, 4))
        finally:
            self._sync_in_progress = False

    def _get_r_geo(self) -> float:
        """Retrieve geocentric distance from main window model or obs override."""
        try:
            mw = self.ctrl.window()
            if mw and hasattr(mw, '_model') and mw._model:
                ovr = self.ctrl.get_obs_override()
                return ovr["r_geo"] if ovr else mw._model["info"]["r_geo"]
        except Exception:
            pass
        return 1.0   # fallback

    def _apply_plate_scale(self):
        """Compute AU/px from entered arcsec/px + r_geo, write to ctrl.au_px_spin."""
        ps    = self.ps_spin.value()
        r_geo = self._get_r_geo()
        au_px = ps * (np.pi / 180.0 / 3600.0) * r_geo
        self.ctrl.au_px_spin.setValue(au_px)
        self.scale_result_lbl.setText(
            f"✓ Scale applied: {ps:.4f}\"/px → AU/px = {au_px:.4e}\n"
            f"(Δ = {r_geo:.4f} AU)  Re-COMPUTE MODEL to update overlay.")

    def _on_open(self):
        self.ctrl._load_image()
        if self.ctrl._img_arr is not None:
            self.file_lbl.setText(self.ctrl.img_label.text())
            # If FITS with WCS: read plate scale back into ps_spin
            if hasattr(self.ctrl, "_wcs_ps_deg") and self.ctrl._wcs_ps_deg:
                ps_as = self.ctrl._wcs_ps_deg * 3600.0
                self.ps_spin.setValue(round(ps_as, 4))
            self._update_scale_preview()

    def _on_clear(self):
        self.ctrl._clear_image()
        self.file_lbl.setText("No image loaded")
        self.scale_result_lbl.setText("")
        self._update_scale_preview()

    def sync_label(self, text):
        self.file_lbl.setText(text)


# ─────────────────────────────────────────────────────────────────────────────
#  LC POPUP WINDOW  — full-size light curve viewer
# ─────────────────────────────────────────────────────────────────────────────
class LCWindow(QDialog):
    """Standalone light curve window — independent, safe matplotlib canvas."""

    def __init__(self, cobs_data: dict, obs_r=None, obs_delta=None,
                 obs_jd=None, T_jd=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(
            f"☄ Light Curve — {cobs_data.get('comet_name','')}")
        self.setMinimumSize(860, 580)
        self.resize(960, 660)
        self.setWindowFlags(
            Qt.WindowType.Window |
            Qt.WindowType.WindowCloseButtonHint |
            Qt.WindowType.WindowMinMaxButtonsHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(10, 10, 10, 10)
        vbox.setSpacing(6)

        # ── Header ─────────────────────────────────────────────────────
        H0 = cobs_data.get("H0"); n = cobs_data.get("n")
        hdr_txt = cobs_data.get("comet_name","Comet")
        if H0 is not None:
            hdr_txt += f"   |   H₀ = {H0:.2f}   n = {n:.2f}"
        hdr = QLabel(hdr_txt)
        hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_value']};"
            f"font-size:13px; font-weight:700; padding:8px 12px;"
            f"border-radius:6px; border:1px solid {T['border']};")
        vbox.addWidget(hdr)

        # ── Matplotlib canvas ───────────────────────────────────────────
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg as FC, NavigationToolbar2QT as NT)
        import matplotlib
        matplotlib.rcParams['axes.facecolor'] = T["mpl_bg"]
        matplotlib.rcParams['figure.facecolor'] = T["mpl_bg"]

        self._fig = Figure(figsize=(9.5, 5.5))
        self._fig.set_facecolor(T["mpl_bg"])
        self._ax = self._fig.add_subplot(111)
        self._canvas = FC(self._fig)
        toolbar = NT(self._canvas, self)
        vbox.addWidget(toolbar)
        vbox.addWidget(self._canvas)

        # ── Source + close ──────────────────────────────────────────────
        bot = QHBoxLayout()
        src_lbl = QLabel(f"Source: {cobs_data.get('source','')}")
        src_lbl.setStyleSheet(f"color:{T['text_muted']}; font-size:10px;")
        btn_close = QPushButton("Close"); btn_close.setMaximumWidth(90)
        btn_close.clicked.connect(self.close)
        bot.addWidget(src_lbl); bot.addStretch(); bot.addWidget(btn_close)
        vbox.addLayout(bot)

        # Draw — wrap in try/except
        try:
            self._draw(cobs_data, obs_r, obs_delta, obs_jd, T_jd)
        except Exception as e:
            self._ax.text(0.5, 0.5, f"Plot error:\n{e}",
                          ha='center', va='center',
                          transform=self._ax.transAxes,
                          color='#ff4444', fontsize=10)
        self._canvas.draw()

    def _draw(self, d, obs_r, obs_delta, obs_jd=None, T_jd=None):
        import numpy as np
        from datetime import datetime, timedelta, timezone

        ax = self._ax
        ax.set_facecolor(T["mpl_bg"])
        self._fig.set_facecolor(T["mpl_bg"])
        for sp in ax.spines.values(): sp.set_edgecolor(T["border"])
        ax.tick_params(colors=T["mpl_tick"], labelsize=9)
        ax.grid(True, color=T["mpl_grid"], lw=0.5, alpha=0.4)
        ax.set_ylabel("Apparent Magnitude", color=T["mpl_label"], fontsize=10)

        H0 = d.get("H0"); n = d.get("n")
        obs_list = d.get("obs_list", [])  # has r_helio + delta
        raw_obs  = d.get("raw_obs", [])   # date + mag only

        use_date_axis = bool(raw_obs)

        if use_date_axis:
            # ── X-axis: date (like COBS website) ───────────────────────
            dates, mags = [], []
            for o in raw_obs:
                try:
                    dt = datetime.strptime(o["date"], "%Y-%m-%d")
                    m  = float(o["mag"])
                    if 1 <= m <= 22:
                        dates.append(dt); mags.append(m)
                except Exception: continue

            if dates:
                ax.scatter(dates, mags, s=8, color="#58a6ff",
                           alpha=0.45, zorder=2,
                           label=f"{len(dates)} COBS observations")
                self._fig.autofmt_xdate()

            # Fit curve over time range (smooth)
            if obs_list and H0 and n:
                # Build predicted curve using r from obs_with_eph
                sorted_oph = sorted(obs_list, key=lambda o: o["date"])
                try:
                    dt_arr = [datetime.strptime(o["date"],"%Y-%m-%d") for o in sorted_oph]
                    
                    if len(dt_arr) >= 2:
                        # Smooth curve: interpolate between observations
                        from scipy.interpolate import interp1d
                        
                        # Convert dates to numbers for interpolation
                        date_min = min(dt_arr)
                        dates_num = [(d - date_min).days for d in dt_arr]
                        r_helio_vals = [o["r_helio"] for o in sorted_oph]
                        delta_vals = [o["delta"] for o in sorted_oph]
                        
                        # Create interpolation functions
                        f_r = interp1d(dates_num, r_helio_vals, kind='linear')
                        f_delta = interp1d(dates_num, delta_vals, kind='linear')
                        
                        # Generate 200 smooth points
                        date_range = max(dates_num) - min(dates_num)
                        if date_range > 0:
                            smooth_nums = np.linspace(min(dates_num), max(dates_num), 200)
                            r_smooth = f_r(smooth_nums)
                            delta_smooth = f_delta(smooth_nums)
                            dt_smooth = [date_min + timedelta(days=float(d)) for d in smooth_nums]
                            
                            # Calculate smooth magnitude
                            m_pred = [H0 + 5*np.log10(d) + 2.5*n*np.log10(r)
                                     for r, d in zip(r_smooth, delta_smooth)]
                            
                            ax.plot(dt_smooth, m_pred, color="#ffaa00", lw=2.5,
                                   label=f"H₀={H0:.2f}  n={n:.2f}", zorder=4)
                        else:
                            # Same date: just plot points
                            m_pred = [H0 + 5*np.log10(o["delta"]) + 2.5*n*np.log10(o["r_helio"])
                                     for o in sorted_oph]
                            ax.plot(dt_arr, m_pred, color="#ffaa00", lw=2.5, marker='o',
                                   label=f"H₀={H0:.2f}  n={n:.2f}", zorder=4)
                    else:
                        # Single point
                        m_pred = [H0 + 5*np.log10(o["delta"]) + 2.5*n*np.log10(o["r_helio"])
                                 for o in sorted_oph]
                        ax.plot(dt_arr, m_pred, color="#ffaa00", lw=2.5, marker='o', markersize=8,
                               label=f"H₀={H0:.2f}  n={n:.2f}", zorder=4)
                except Exception as e:
                    # Fallback: simple plot without interpolation
                    logging.warning("Light curve interpolation failed: %s — using simple plot", e)
                    dt_arr = [datetime.strptime(o["date"],"%Y-%m-%d") for o in sorted_oph]
                    m_pred = [H0 + 5*np.log10(o["delta"]) + 2.5*n*np.log10(o["r_helio"])
                             for o in sorted_oph]
                    ax.plot(dt_arr, m_pred, color="#ffaa00", lw=2.5,
                           label=f"H₀={H0:.2f}  n={n:.2f}", zorder=4)

            # ── Observation date vertical line ─────────────────────────
            if obs_jd is not None:
                obs_dt = datetime(1858, 11, 17) + timedelta(days=obs_jd - 2400000.5)
                ax.axvline(obs_dt, color="#58a6ff", lw=1.2, ls="--", alpha=0.85,
                           zorder=5, label=f"Obs date: {obs_dt.strftime('%Y-%m-%d')}")

            # ── Perihelion date vertical line ──────────────────────────
            if T_jd is not None:
                peri_dt = datetime(1858, 11, 17) + timedelta(days=T_jd - 2400000.5)
                ax.axvline(peri_dt, color="#ffa657", lw=1.2, ls="--", alpha=0.85,
                           zorder=5, label=f"Perihelion: {peri_dt.strftime('%Y-%m-%d')}")

            # ── "Now" info — text only, no vertical line ───────────────
            if obs_r and obs_delta and H0 and n:
                now = datetime.now(timezone.utc).replace(tzinfo=None)
                m_now = H0 + 5*np.log10(obs_delta) + 2.5*n*np.log10(obs_r)
                # Show as a small annotation in the lower-right corner
                ax.annotate(
                    f"Now:  m = {m_now:.1f}    r = {obs_r:.3f} AU    Δ = {obs_delta:.3f} AU",
                    xy=(1.0, 0.02), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=8,
                    color="#f85149",
                    bbox=dict(boxstyle="round,pad=0.3",
                              facecolor=T["panel_bg"], edgecolor="#f85149",
                              alpha=0.7))

            ax.set_xlabel("Date", color=T["mpl_label"], fontsize=10)

        else:
            # ── X-axis: r☉ (AU) — fallback if no raw obs ───────────────
            ax.set_xlabel("r☉ (AU)", color=T["mpl_label"], fontsize=10)
            if obs_list:
                pts = [(o["r_helio"],o["mag"]) for o in obs_list
                       if 1<=o["mag"]<=22 and o["r_helio"]>0.1]
                if pts:
                    rs=np.array([p[0] for p in pts])
                    ms=np.array([p[1] for p in pts])
                    ax.scatter(rs, ms, s=10, color="#58a6ff", alpha=0.5,
                               label=f"{len(pts)} obs")
            if H0 and n:
                r_fit = np.linspace(0.5, 6.0, 300)
                d_ref = obs_delta or 1.5
                m_fit = H0 + 5*np.log10(d_ref) + 2.5*n*np.log10(r_fit)
                ax.plot(r_fit, m_fit, color="#ffaa00", lw=2.0,
                        label=f"H₀={H0:.2f}  n={n:.2f}")
            if obs_r and obs_delta and H0 and n:
                m_now = H0 + 5*np.log10(obs_delta) + 2.5*n*np.log10(obs_r)
                ax.annotate(
                    f"Now:  m = {m_now:.1f}    r = {obs_r:.3f} AU",
                    xy=(1.0, 0.02), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=8, color="#f85149",
                    bbox=dict(boxstyle="round,pad=0.3",
                              facecolor=T["panel_bg"], edgecolor="#f85149",
                              alpha=0.7))

        ax.set_title(d.get("comet_name",""), color=T["mpl_title"],
                     fontsize=11, pad=6)
        ax.invert_yaxis()
        ax.legend(fontsize=9, labelcolor=T["text"], framealpha=0.35,
                  facecolor=T["panel_bg"], edgecolor=T["border"])
        self._fig.tight_layout(pad=1.2)


# ─────────────────────────────────────────────────────────────────────────────
#  ORBIT-POSITION WINDOW  (v3.0)
# ─────────────────────────────────────────────────────────────────────────────
class OrbitWindow(QDialog):
    """
    Standalone 3D orbit-position diagram window (v3.0).

    Distinct from the "Orbital path" overlay already drawn on the main
    canvas: that overlay shows the comet's path PROJECTED ONTO THE SKY
    (used to judge tail-axis direction). This window shows where the
    comet physically SITS in its orbit right now, in full 3D — Sun at
    the focus (yellow circle), perihelion marked, Earth's orbit for
    scale, drag-to-rotate — independent, safe matplotlib canvas, same
    pattern as LCWindow.
    """

    def __init__(self, comet_el: dict, obs_jd: float, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"☄ Orbit Position — {comet_el.get('name','Comet')}")
        self.setMinimumSize(680, 680)
        self.resize(760, 760)
        self.setWindowFlags(
            Qt.WindowType.Window |
            Qt.WindowType.WindowCloseButtonHint |
            Qt.WindowType.WindowMinMaxButtonsHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(10, 10, 10, 10)
        vbox.setSpacing(6)

        hdr = QLabel(f"{comet_el.get('name','Comet')}   |   "
                     f"q={comet_el.get('q',0):.4f} AU   e={comet_el.get('e',0):.5f}")
        hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_value']};"
            f"font-size:13px; font-weight:700; padding:8px 12px;"
            f"border-radius:6px; border:1px solid {T['border']};")
        vbox.addWidget(hdr)

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg as FC, NavigationToolbar2QT as NT)
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 — registers '3d' projection

        self._fig = Figure(figsize=(7, 7))
        self._fig.set_facecolor(T["mpl_bg"])
        self._ax = self._fig.add_subplot(111, projection='3d')
        self._canvas = FC(self._fig)
        toolbar = NT(self._canvas, self)
        vbox.addWidget(toolbar)
        vbox.addWidget(self._canvas)

        # Zoom to inner solar system (v3.0.7) — for comets well beyond
        # ~1 AU, the default view auto-fits to the WHOLE orbit, which
        # collapses Earth's 1 AU orbit (and the Earth/Sun separation) down
        # to a barely-distinguishable point near the center. This swaps
        # to a fixed Sun-centered view instead — the comet itself may end
        # up off-screen if it's currently farther out than the AU value
        # below, which is the explicit trade-off being made.
        zoom_row = QHBoxLayout()
        self.chk_zoom_inner = QCheckBox("🔍  Zoom to inner solar system")
        self.chk_zoom_inner.toggled.connect(lambda _: self._redraw())
        self.sp_zoom_au = QDoubleSpinBox()
        self.sp_zoom_au.setRange(0.1, 50.0)
        self.sp_zoom_au.setDecimals(2)
        self.sp_zoom_au.setValue(2.0)
        self.sp_zoom_au.setSuffix("  AU")
        self.sp_zoom_au.setToolTip(
            "Half-width of the fixed Sun-centered view. 2 AU comfortably "
            "shows Earth's full 1 AU orbit with some margin.")
        self.sp_zoom_au.valueChanged.connect(
            lambda _: self.chk_zoom_inner.isChecked() and self._redraw())
        zoom_row.addWidget(self.chk_zoom_inner)
        zoom_row.addWidget(self.sp_zoom_au)
        zoom_row.addStretch()
        vbox.addLayout(zoom_row)

        bot = QHBoxLayout()
        note_lbl = QLabel("3D heliocentric view  ·  click-drag to rotate  ·  "
                          "physical orbit position, not the sky-projected tail")
        note_lbl.setStyleSheet(f"color:{T['text_muted']}; font-size:10px;")
        btn_close = QPushButton("Close"); btn_close.setMaximumWidth(90)
        btn_close.clicked.connect(self.close)
        bot.addWidget(note_lbl); bot.addStretch(); bot.addWidget(btn_close)
        vbox.addLayout(bot)

        self._diagram = None
        try:
            self._diagram = cta.compute_orbit_diagram(comet_el, obs_jd)
            cta.draw_orbit_diagram(self._ax, self._diagram, dark=(CURRENT_THEME == "dark"))
            self._fig.tight_layout(pad=1.2)
        except Exception as e:
            self._ax.text2D(0.5, 0.5, f"Diagram error:\n{e}",
                            ha='center', va='center',
                            transform=self._ax.transAxes,
                            color='#ff4444', fontsize=10)
        self._canvas.draw()

    def _redraw(self):
        """Re-render with the current Zoom-to-inner-solar-system setting.
        Reuses the already-computed self._diagram — only the VIEW changes,
        not the underlying orbit data, so there's no need to recompute it."""
        if self._diagram is None:
            return
        zoom_au = self.sp_zoom_au.value() if self.chk_zoom_inner.isChecked() else None
        self._ax.clear()
        try:
            cta.draw_orbit_diagram(self._ax, self._diagram,
                                   dark=(CURRENT_THEME == "dark"),
                                   zoom_inner_au=zoom_au)
            self._fig.tight_layout(pad=1.2)
        except Exception as e:
            self._ax.text2D(0.5, 0.5, f"Diagram error:\n{e}",
                            ha='center', va='center',
                            transform=self._ax.transAxes,
                            color='#ff4444', fontsize=10)
        self._canvas.draw()


class MCWindow(QDialog):
    """
    Monte Carlo dust morphology control window (v3.1) — population-level
    extension of the F-P syndyne/synchrone model. Where the main canvas
    traces one curve per (β, age) choice, this samples a whole grain
    population (size distribution + release-time window) and bins the
    result into a density map, then contours it straight onto the MAIN
    canvas (Phase 2) — through that canvas's own already-correct to_px()
    pipeline, not a separate one here. See cta.compute_morphology_mc()'s
    docstring for the full model and the Phase 1 simplification list
    (uniform Q(t), isotropic ejection direction, fixed size-distribution
    slope, no PSF/isophote-vs-real-image fitting yet).

    v3.1 design change: this window has NO plot of its own anymore (the
    contour on the main canvas IS the result — keeping a second, separate
    plot here was redundant once Phase 2 existed, and Run now sends
    automatically rather than needing a separate Send step). It is also
    now the ONLY place ejection velocity is entered — see ControlPanel.
    get_ejection_params()/set_ejection_params(): clicking Run here pushes
    v0→v_N0 (others 0), γ, m into the main panel's ejection state AND
    triggers a main-model recompute, so the F-P curves and the Monte
    Carlo population always agree on what ejection velocity is in effect
    instead of risking two different answers kept in two different
    windows. v0→v_N0 specifically (not v_R0 or v_T0) because the normal-
    to-orbit-plane direction has been this project's actual diagnostic
    interest throughout (out-of-plane broadening of a tail feature) —
    if you need a different single-direction mapping, edit _on_done()
    below; the comment there says exactly where.

    beta_range here is deliberately NOT the same input as the main
    panel's β list — prefilled from min/max of that list as a
    convenience starting point, but it's a continuous population range,
    not a handful of discrete diagnostic test points (see the v3.1
    changelog for why these aren't interchangeable).
    """

    # v3.1 — class attribute (not instance): survives this window being
    # closed and a brand-new MCWindow being constructed next time it's
    # opened (WA_DeleteOnClose means the old instance is gone). See
    # closeEvent()/_collect_state()/_apply_state() below.
    _saved_state: dict | None = None

    def __init__(self, comet_el: dict, obs_jd: float,
                 beta_values=None, max_age: float = 200,
                 fp_sync_ages=None, fp_dominant_age: float = 0.0,
                 fp_model: dict | None = None, parent=None,
                 main_window=None):
        super().__init__(parent)
        self.comet_el  = comet_el
        self.obs_jd    = obs_jd
        self.fp_model  = fp_model
        try:
            self.fp_sync_ages = [float(a) for a in (fp_sync_ages or [])
                                 if float(a) > 0]
        except Exception:
            self.fp_sync_ages = []
        # Maximum dust age is authoritative from the largest listed
        # synchrone whenever that list is available.  The legacy max_age
        # argument is only a fallback for callers that do not supply ages.
        self.fp_max_age = (max(self.fp_sync_ages) if self.fp_sync_ages else
                           (float(max_age) if float(max_age) > 0 else 0.0))
        try:
            self.fp_dominant_age = max(0.0, float(fp_dominant_age))
        except Exception:
            self.fp_dominant_age = 0.0
        if self.fp_max_age > 0 and self.fp_dominant_age > self.fp_max_age:
            self.fp_dominant_age = self.fp_max_age
        # v3.1, Phase 2 — plain Python reference (NOT Qt parent; this window
        # uses parent=None, same reasoning as OrbitWindow/LCWindow) so Run
        # can push contours + ejection params back without going through
        # Qt's parent/child signal plumbing.
        self.main_window = main_window
        self._mc_result = None
        self._smooth_au_initialized = False
        self._cobs_obs_list = None   # v3.1 — populated by _fetch_cobs_for_qt()
        # A small persistent cache used by the editable-input workflow.
        # Some GUI paths can redraw/refresh tabs after load; the cache prevents
        # a stale/cleared _cobs_obs_list from blocking Run MC after the status
        # line already reports that COBS is ready.
        self._cobs_ready_cache = None
        self._pending_run_after_cobs = False
        self._last_qt_result = None
        self._last_qt_dips = []
        # Q(t) inference/sampling coverage state.  COBS can contain old
        # apparitions separated by long gaps; those data must never create a
        # multi-year MC-window recommendation by interpolation.
        self._qt_coverage = None
        self._qt_recommendation_quality = ""
        self._qt_provisional_window = None
        # Editable input-file workflow.  A loaded .mcin file is treated as a
        # transparent template: it restores visible GUI fields only.  Users may
        # edit any field before running MC, then save the edited state as a new
        # .mcin file.
        self._loaded_input_path = None
        self._loaded_input_profile = {}
        self._loaded_input_comet = ""
        self._worker     = None
        # Async-result generation. Clear All increments this value so a late
        # Monte Carlo worker cannot restore a result after it was cleared.
        self._run_generation = 0

        # First-use guided workflow.  The tab content is always available, but
        # new users start on F-P Guide with Back/Next navigation.  After one
        # successful MC run (or Skip guide), later sessions open directly on
        # Simulation.  This is a UI preference only and never changes physics.
        self._guide_settings = QSettings("MPC-O58", "CometTailAnalyzer")
        guide_done_raw = self._guide_settings.value("mc_guided_completed", False)
        self._guided_first_use = str(guide_done_raw).strip().lower() not in (
            "1", "true", "yes", "on")

        self.setWindowTitle(f"🎲 Monte Carlo Morphology — {comet_el.get('name','Comet')}")
        self.setMinimumSize(940, 720)
        self.resize(1180, 900)
        self.setWindowFlags(
            Qt.WindowType.Window |
            Qt.WindowType.WindowCloseButtonHint |
            Qt.WindowType.WindowMinMaxButtonsHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(10, 10, 10, 10)
        vbox.setSpacing(6)

        hdr = QLabel(f"{comet_el.get('name','Comet')}   |   "
                     f"q={comet_el.get('q',0):.4f} AU   e={comet_el.get('e',0):.5f}")
        hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_value']};"
            f"font-size:13px; font-weight:700; padding:8px 12px;"
            f"border-radius:6px; border:1px solid {T['border']};")
        vbox.addWidget(hdr)

        note = QLabel(
            "Monte Carlo dust morphology — qualitative comparison with observed image. "
            "β range is a continuous population range (not the F-P syndyne list). "
            "Ejection velocity set here also updates the main F-P model.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{T['text_muted']}; font-size:9px; padding:0 4px 4px 4px;")
        vbox.addWidget(note)

        # First-use guide banner. Users can skip it without changing any input.
        self.guide_banner = QFrame()
        self.guide_banner.setStyleSheet(
            "QFrame { background:#102338; border:1px solid #4a8adf; border-radius:6px; }"
            "QLabel { border:none; background:transparent; }")
        guide_banner_l = QHBoxLayout(self.guide_banner)
        guide_banner_l.setContentsMargins(10, 7, 10, 7)
        self.lbl_guide_banner = QLabel(
            "Guided setup: review the F-P starting range, continue to Simulation, "
            "then run the model. Display is optional and can be adjusted afterward.")
        self.lbl_guide_banner.setWordWrap(True)
        self.lbl_guide_banner.setStyleSheet("color:#eaf4ff; font-size:10px;")
        self.btn_skip_guide = QPushButton("Skip guide")
        self.btn_skip_guide.setMaximumWidth(100)
        self.btn_skip_guide.setAutoDefault(False)
        self.btn_skip_guide.setDefault(False)
        guide_banner_l.addWidget(self.lbl_guide_banner, 1)
        guide_banner_l.addWidget(self.btn_skip_guide)
        self.guide_banner.setVisible(self._guided_first_use)
        vbox.addWidget(self.guide_banner)

        def _dspin(val, lo, hi, decimals=4, step=0.001, suffix=""):
            w = QDoubleSpinBox()
            w.setRange(lo, hi); w.setDecimals(decimals); w.setSingleStep(step)
            w.setValue(val)
            if suffix:
                w.setSuffix(suffix)
            return w

        # v3.1 — three-tab layout, separating three genuinely different
        # KINDS of decision so they don't all sit in one undifferentiated
        # wall of fields:
        #   Tab 1 "F-P Guide"  — read what the main panel's F-P syndynes
        #     already suggest (β range, release window) and convert that
        #     into a starting grain-size range for the simulation below.
        #     Nothing here is itself fed into the simulation directly —
        #     it only WRITES into Tab 2's fields when you click Apply.
        #   Tab 2 "Simulate"   — every input that actually changes what
        #     gets SIMULATED (grain sizes, particle count, release
        #     window, ejection speed/direction, production-rate shape).
        #     Changing anything here means Run must be clicked again.
        #   Tab 3 "Display"    — re-thresholding an ALREADY-COMPUTED
        #     result into a contour (floor/levels/smoothing). Re-extract
        #     here is fast and needs no new simulation at all.
        # Labels throughout favour plain wording over raw Python variable
        # names, and follow py_COMTAILS's TAIL_INPUTS.dat naming where
        # there's a direct equivalent (Moreno, F. 2025, A&A 695, A263;
        # python port by R. Morales & N. Robles, IAA-CSIC; https://
        # github.com/FernandoMorenoDanvila/py_COMTAILS, MIT licence).

        self._smooth_au_initialized = False

        # sp_smooth_au — hidden state carrier for smooth sigma used in
        # extract_morphology_contours(). Always 0 (no smoothing); kept
        # as a QDoubleSpinBox so _collect_state() / save / load work
        # without special-casing it. _on_done() writes 0.0 here each run.
        self.sp_smooth_au = QDoubleSpinBox()
        self.sp_smooth_au.setRange(0.0, 1.0)
        self.sp_smooth_au.setDecimals(7)
        self.sp_smooth_au.setValue(0.0)
        self.sp_smooth_au.hide()

        tabs = QTabWidget()
        self.tabs = tabs
        self._visited_tabs = {0}
        tabs.setDocumentMode(True)
        tabs.setStyleSheet(
            "QTabBar::tab { padding: 6px 18px; font-size: 11px; font-weight: bold; }"
            "QTabBar::tab:selected { border-bottom: 3px solid #4a8adf; }")

        page_guide = QWidget(); guide_v = QVBoxLayout(page_guide)
        page_sim   = QWidget(); sim_v   = QVBoxLayout(page_sim); sim_v.setSpacing(10)
        page_disp  = QWidget(); disp_v  = QVBoxLayout(page_disp)

        # Input widths are intentionally bounded: wide enough to show the
        # complete value + suffix, but not so wide that numeric controls
        # dominate the form on large screens.
        W_SMALL = 112
        W_MED   = 145
        W_WIDE  = 172

        def _sp(w, width=W_MED):
            w.setMinimumWidth(width)
            w.setMaximumWidth(width + 18)
            return w

        # ════════════════════════════════════════════════════════════════
        # GRAIN DENSITY (shared by Guide tab and Simulate tab)
        # ════════════════════════════════════════════════════════════════
        self.sp_rho = _dspin(0.5, 0.05, 10.0, 2, 0.05, " g/cm³")
        _sp(self.sp_rho, W_MED)
        self.sp_rho.setToolTip(
            "Grain bulk density — converts between grain radius and β.\n"
            "CTA default 0.5 g/cm³; many published models use 1.0 g/cm³.")

        # ════════════════════════════════════════════════════════════════
        # TAB 2 — SIMULATE
        # ════════════════════════════════════════════════════════════════

        # ── Grain size ───────────────────────────────────────────────────
        grp_grain = QGroupBox("Grain size distribution")
        gg = QGridLayout(grp_grain); gg.setHorizontalSpacing(10); gg.setVerticalSpacing(6)

        beta_lo = min(beta_values) if beta_values else 0.001
        beta_hi = max(beta_values) if beta_values else 0.1
        r_max_default = cta._beta_to_radius_um(beta_lo)
        r_min_default = cta._beta_to_radius_um(beta_hi)

        self.sp_r_min = _sp(_dspin(r_min_default, 0.0001, 1.0e5, 4, 0.01, " µm"), W_WIDE)
        self.sp_r_max = _sp(_dspin(r_max_default, 0.0001, 1.0e5, 2, 1.0,  " µm"), W_WIDE)
        self.sp_gamma_size = _sp(_dspin(-3.6, -10.0, 5.0, 2, 0.1), W_SMALL)
        self.sp_gamma_size.setToolTip(
            "Size-distribution slope: dN/da ∝ a^κ.  More negative = more small grains.")

        self.lbl_beta_equiv = QLabel()
        self.lbl_beta_equiv.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        def _update_beta_equiv():
            try:
                rho = self.sp_rho.value()
                b_max = cta._radius_um_to_beta(self.sp_r_min.value(), rho)
                b_min = cta._radius_um_to_beta(self.sp_r_max.value(), rho)
                self.lbl_beta_equiv.setText(
                    f"≡ β [{b_min:.4g} – {b_max:.4g}]  (ρ={rho:.2g} g/cm³, Qpr=1)")
            except Exception:
                self.lbl_beta_equiv.setText("")
        self.sp_r_min.valueChanged.connect(_update_beta_equiv)
        self.sp_r_max.valueChanged.connect(_update_beta_equiv)
        self.sp_rho.valueChanged.connect(_update_beta_equiv)
        _update_beta_equiv()

        gg.addWidget(QLabel("r_min (µm)"), 0, 0); gg.addWidget(self.sp_r_min, 0, 1)
        gg.addWidget(QLabel("r_max (µm)"), 0, 2); gg.addWidget(self.sp_r_max, 0, 3)
        gg.addWidget(self.lbl_beta_equiv, 1, 0, 1, 4)
        gg.addWidget(QLabel("Slope κ"), 2, 0); gg.addWidget(self.sp_gamma_size, 2, 1)
        sim_v.addWidget(grp_grain)

        # ── Simulation parameters ────────────────────────────────────────
        grp_sim = QGroupBox("Simulation")
        sg = QGridLayout(grp_sim); sg.setHorizontalSpacing(10); sg.setVerticalSpacing(8)

        self.sp_n = QSpinBox(); self.sp_n.setRange(100, 2000000)
        self.sp_n.setValue(50000); self.sp_n.setSingleStep(50000)
        _sp(self.sp_n, W_MED)
        self.sp_n.setToolTip(
            "Number of simulated particles.\n"
            "50,000  — quick preview (~5s)\n"
            "500,000 — good quality (~1 min)\n"
            "1,000,000 — publication quality (~2 min)\n"
            "2,000,000 — maximum quality (~5 min)")

        # Grid resolution fixed at 400 — not user-adjustable.
        # Grid is a computational parameter (always higher=better); users
        # should control quality/speed via particle count instead.
        # 400×400 = 160,000 bins; with 1M particles → ~6 particles/bin
        # average, enough for smooth contours without smoothing.
        self._grid_npix = 400

        self.sp_max_age = _sp(_dspin(max_age, 1.0, 2000.0, 1, 10.0, " d"), W_MED)
        self.sp_max_age.setToolTip(
            "Dust release window: how many days before the observation\n"
            "activity is assumed to span. Not the same as the main panel\n"
            "Max age — estimate this independently from the image size.")

        self.sp_v0    = _sp(_dspin(0.0, 0.0, 2000.0, 2, 1.0, " m/s"), W_MED)
        self.sp_gamma = _sp(_dspin(0.5, 0.0, 5.0, 2, 0.1), W_SMALL)
        self.sp_mexp  = _sp(_dspin(0.5, -5.0, 5.0, 2, 0.1), W_SMALL)
        self.sp_v0.setToolTip(
            "V₀ is the reference coefficient at β=1 and r_H=1 AU — it is NOT "
            "the actual speed of every grain.\n"
            "CTA evaluates each grain with V = V₀·β^γ·r_H^−κ.")
        self.sp_gamma.setToolTip("γ: V ∝ β^γ  (size dependence)")
        self.sp_mexp.setToolTip("κ: V ∝ r_H^−κ  (heliocentric distance dependence)")

        # Live translation of the abstract V₀ coefficient into the actual
        # terminal-speed range represented by the current grain-size limits.
        # This prevents V₀ from being mistaken for the speed of every particle.
        self.lbl_actual_speed = QLabel()
        self.lbl_actual_speed.setWordWrap(True)
        self.lbl_actual_speed.setStyleSheet(
            f"color:{T['text_muted']}; font-size:9px; padding:2px 0 4px 0;")

        def _update_actual_speed_label():
            try:
                v0 = float(self.sp_v0.value())
                rho = float(self.sp_rho.value())
                beta_fast = cta._radius_um_to_beta(self.sp_r_min.value(), rho)
                beta_slow = cta._radius_um_to_beta(self.sp_r_max.value(), rho)
                r_c, _ = cta.elem_to_state(self.comet_el, self.obs_jd)
                r_h = float(cta.vmag(r_c))
                v_fast = float(cta.real_ejection_speed_ms(
                    v0, beta_fast, r_h,
                    self.sp_gamma.value(), self.sp_mexp.value()))
                v_slow = float(cta.real_ejection_speed_ms(
                    v0, beta_slow, r_h,
                    self.sp_gamma.value(), self.sp_mexp.value()))
                v_lo, v_hi = sorted((v_fast, v_slow))
                if abs(v_hi - v_lo) < 0.005:
                    speed_text = f"{v_lo:.2f} m/s"
                elif v_hi < 10.0:
                    speed_text = f"{v_lo:.2f}–{v_hi:.2f} m/s"
                else:
                    speed_text = f"{v_lo:.1f}–{v_hi:.1f} m/s"
                self.lbl_actual_speed.setText(
                    f"Actual grain-speed range: {speed_text} at r_H={r_h:.3f} AU "
                    f"for the current r_min–r_max range.  "
                    "V₀ is a β=1, r_H=1 AU reference coefficient, not the speed "
                    "assigned to every grain.")
            except Exception:
                self.lbl_actual_speed.setText(
                    "V₀ is a reference coefficient. Actual grain speed is evaluated "
                    "from V = V₀·β^γ·r_H^−κ.")

        for _w in (self.sp_v0, self.sp_gamma, self.sp_mexp,
                   self.sp_r_min, self.sp_r_max, self.sp_rho):
            _w.valueChanged.connect(_update_actual_speed_label)
        _update_actual_speed_label()

        self.chk_seed = QCheckBox("Fixed seed")
        self.sp_seed  = QSpinBox(); self.sp_seed.setRange(0, 999999)
        self.sp_seed.setValue(42); _sp(self.sp_seed, W_SMALL); self.sp_seed.setEnabled(False)
        self.chk_seed.toggled.connect(self.sp_seed.setEnabled)
        self.chk_seed.setToolTip("Fix random seed for reproducible / comparable runs.")

        # ── Row layout — each field gets ONE unbroken grid row, no shared
        # row indices between unrelated fields (previous version had V0/γ
        # silently overwriting the Release-window/grid-note row).
        r = 0
        sg.addWidget(QLabel("Particles"), r, 0)
        sg.addWidget(self.sp_n, r, 1)
        r += 1
        # Quick-select buttons on their OWN full-width row, roomy enough
        # to not clip labels (previous version squeezed 4 buttons into a
        # 2-column grid cell meant for one spinbox, cutting off text).
        n_row = QHBoxLayout(); n_row.setSpacing(6)
        n_row.addWidget(QLabel(""))  # small left indent to align under field
        for label, val in [("50k", 50000), ("500k", 500000),
                           ("1M", 1000000), ("2M", 2000000)]:
            btn = QPushButton(label)
            btn.setMinimumWidth(56)
            btn.setToolTip(f"{val:,} particles")
            btn.clicked.connect(lambda _, v=val: self.sp_n.setValue(v))
            n_row.addWidget(btn)
        n_row.addStretch()
        sg.addLayout(n_row, r, 0, 1, 4)
        r += 1

        sg.addWidget(QLabel("Release window"), r, 0)
        sg.addWidget(self.sp_max_age, r, 1)
        grid_note = QLabel("Grid: 400×400 (fixed)")
        grid_note.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        sg.addWidget(grid_note, r, 2, 1, 2)
        r += 1

        release_hint = QLabel(
            "A Q(t)-based recommendation will appear in Dust production over time below "
            "when COBS or a manual dM/dt table is selected.")
        release_hint.setWordWrap(True)
        release_hint.setStyleSheet(f"color:{T['text_muted']}; font-size:9px; padding:2px 0 4px 0;")
        sg.addWidget(release_hint, r, 0, 1, 4)
        r += 1

        sg.addWidget(QLabel("V₀ reference coefficient (m/s)"), r, 0)
        sg.addWidget(self.sp_v0, r, 1)
        sg.addWidget(QLabel("γ (speed-size)"), r, 2)
        sg.addWidget(self.sp_gamma, r, 3)
        r += 1

        sg.addWidget(self.lbl_actual_speed, r, 0, 1, 4)
        r += 1

        sg.addWidget(QLabel("κ (speed-dist)"), r, 0)
        sg.addWidget(self.sp_mexp, r, 1)
        seed_row = QHBoxLayout(); seed_row.setSpacing(6)
        seed_row.addWidget(self.chk_seed); seed_row.addWidget(self.sp_seed)
        seed_row.addStretch()
        sg.addLayout(seed_row, r, 2, 1, 2)
        sim_v.addWidget(grp_sim)

        # ── Ejection direction ───────────────────────────────────────────
        grp_dir = QGroupBox("Ejection direction")
        dir_v = QVBoxLayout(grp_dir); dir_v.setSpacing(5)

        dir_row = QHBoxLayout()
        self.rad_mode_isotropic = QRadioButton("Isotropic")
        self.rad_mode_sunward   = QRadioButton("Sunward hemisphere")
        self.rad_mode_active    = QRadioButton("Active area")
        self.rad_mode_isotropic.setChecked(True)
        self._mode_group = QButtonGroup(self)
        for rb in (self.rad_mode_isotropic, self.rad_mode_sunward, self.rad_mode_active):
            self._mode_group.addButton(rb); dir_row.addWidget(rb)
        dir_row.addStretch()
        dir_v.addLayout(dir_row)

        self.sunward_group = QWidget()
        swg = QGridLayout(self.sunward_group)
        swg.setContentsMargins(0, 0, 0, 0)
        swg.setHorizontalSpacing(10); swg.setVerticalSpacing(4)

        self.sp_sunward_expocos = _sp(_dspin(1.0, 0.0, 10.0, 2, 0.1), W_SMALL)
        self.sp_sunward_expocos.setToolTip(
            "Illumination-angle exponent ε: V_ej ∝ cos(z)^ε.\n"
            "ε=1 follows the usual cosine-weighted sunward emission.\n"
            "Higher values concentrate ejection closer to the subsolar direction.")

        # The physical Sun direction is always evaluated at each particle's
        # emission time. Keep the legacy combo as a hidden state carrier so
        # older .mcin files still load, but do not expose the observation-time
        # shortcut in the standard workflow.
        self.cmb_sunward_reference = QComboBox()
        self.cmb_sunward_reference.addItems([
            "Emission time — physical",
            "Observation time — legacy diagnostic",
        ])
        self.cmb_sunward_reference.setCurrentIndex(0)
        self.cmb_sunward_reference.hide()

        self.sp_sunward_cone = _sp(_dspin(90.0, 0.0, 90.0, 1, 5.0, "°"), W_MED)
        self.sp_sunward_cone.setToolTip(
            "Half-angle around the subsolar direction.\n"
            "90° = full physical sunward hemisphere.\n"
            "Smaller values are diagnostic/fitting aids, not the default model.")

        sunward_note = QLabel(
            "Dust is emitted relative to the Sun direction at each particle's "
            "emission time (physical default).")
        sunward_note.setWordWrap(True)
        sunward_note.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")

        self.chk_projected_sunward = QCheckBox(
            "Projection gate: retain particles projected toward the apparent Sun")
        self.chk_projected_sunward.setChecked(False)
        self.chk_projected_sunward.setToolTip(
            "Diagnostic only. This image-plane filter can create or enhance an\n"
            "apparent sunward structure. It is not a physical force and should\n"
            "not be interpreted as a measured emission constraint.")

        # Retained internally for backward compatibility with older .mcin
        # files, but intentionally hidden from the standard v3.1.1 UI. The
        # projection gate is a diagnostic image-plane filter rather than a
        # physical emission parameter and may be reconsidered in a later
        # advanced/research workflow.
        self.grp_projection_diag = QGroupBox("Advanced projection diagnostics")
        self.grp_projection_diag.setCheckable(True)
        self.grp_projection_diag.setChecked(False)
        self.grp_projection_diag.setVisible(False)
        diag_v = QVBoxLayout(self.grp_projection_diag)
        diag_v.setContentsMargins(8, 5, 8, 7)
        diag_v.addWidget(self.chk_projected_sunward)
        diag_warn = QLabel(
            "Diagnostic only — keep Off for normal physical modeling. "
            "Any use is recorded in the MC report.")
        diag_warn.setWordWrap(True)
        diag_warn.setStyleSheet("color:#ffb020; font-size:9px;")
        diag_v.addWidget(diag_warn)
        self.chk_projected_sunward.setVisible(False)
        diag_warn.setVisible(False)

        def _toggle_projection_diag(opened):
            self.chk_projected_sunward.setVisible(bool(opened))
            diag_warn.setVisible(bool(opened))
        self.grp_projection_diag.toggled.connect(_toggle_projection_diag)

        swg.addWidget(QLabel("cos(z) exponent ε"), 0, 0)
        swg.addWidget(self.sp_sunward_expocos, 0, 1)
        swg.addWidget(QLabel("Cone half-angle"), 0, 2)
        swg.addWidget(self.sp_sunward_cone, 0, 3)
        swg.addWidget(sunward_note, 1, 0, 1, 4)
        self.sunward_group.setVisible(False)
        dir_v.addWidget(self.sunward_group)

        self.active_area_group = QWidget()
        aag = QGridLayout(self.active_area_group)
        aag.setContentsMargins(0, 0, 0, 0); aag.setHorizontalSpacing(10); aag.setVerticalSpacing(4)
        self.sp_nuc_inc  = _sp(_dspin(45.0, 0.0, 90.0, 1, 1.0, "°"))
        self.sp_nuc_phi  = _sp(_dspin(0.0, -360.0, 360.0, 1, 5.0, "°"))
        self.sp_period   = _sp(_dspin(0.5, 0.001, 1000.0, 4, 0.01, " d"))
        self.sp_lat_min  = _sp(_dspin(-15.0, -90.0, 90.0, 1, 1.0, "°"))
        self.sp_lat_max  = _sp(_dspin(15.0, -90.0, 90.0, 1, 1.0, "°"))
        self.sp_lon_min  = _sp(_dspin(-20.0, -360.0, 360.0, 1, 1.0, "°"))
        self.sp_lon_max  = _sp(_dspin(20.0, -360.0, 360.0, 1, 1.0, "°"))
        self.chk_isun    = QCheckBox("Sunlit ground only")
        self.chk_isun.setChecked(True)
        self.sp_expocos  = _sp(_dspin(1.0, 0.0, 10.0, 2, 0.1))
        aag.addWidget(QLabel("  Obliquity I"), 0, 0); aag.addWidget(self.sp_nuc_inc, 0, 1)
        aag.addWidget(QLabel("Subsolar phase Φ (perihelion)"), 0, 2); aag.addWidget(self.sp_nuc_phi, 0, 3)
        aag.addWidget(QLabel("  Rotation period"), 1, 0); aag.addWidget(self.sp_period, 1, 1)
        aag.addWidget(self.chk_isun, 1, 2, 1, 2)
        aag.addWidget(QLabel("  Latitude range"), 2, 0)
        lat_row = QHBoxLayout(); lat_row.setSpacing(4)
        lat_row.addWidget(self.sp_lat_min); lat_row.addWidget(QLabel("to")); lat_row.addWidget(self.sp_lat_max)
        aag.addLayout(lat_row, 2, 1)
        aag.addWidget(QLabel("Longitude range"), 2, 2)
        lon_row = QHBoxLayout(); lon_row.setSpacing(4)
        lon_row.addWidget(self.sp_lon_min); lon_row.addWidget(QLabel("to")); lon_row.addWidget(self.sp_lon_max)
        aag.addLayout(lon_row, 2, 3)
        aag.addWidget(QLabel("  cos(z) exponent"), 3, 0); aag.addWidget(self.sp_expocos, 3, 1)
        self.active_area_group.setVisible(False)
        dir_v.addWidget(self.active_area_group)

        def _on_mode_changed():
            self.sunward_group.setVisible(self.rad_mode_sunward.isChecked())
            self.active_area_group.setVisible(self.rad_mode_active.isChecked())
        for rb in (self.rad_mode_isotropic, self.rad_mode_sunward, self.rad_mode_active):
            rb.toggled.connect(_on_mode_changed)
        sim_v.addWidget(grp_dir)

        # ── Dust production ───────────────────────────────────────────────
        grp_qt = QGroupBox("Dust production over time")
        qt_v = QVBoxLayout(grp_qt); qt_v.setSpacing(8)

        # Source selector — 3 visually distinct toggle-style buttons
        src_lbl = QLabel("Source:")
        src_lbl.setStyleSheet("font-weight: bold;")
        qt_v.addWidget(src_lbl)

        self.rad_dmdt_steady = QRadioButton("① Steady  (constant rate)")
        self.rad_dmdt_cobs   = QRadioButton("② COBS light curve  (scaled by your Afρ)")
        self.rad_dmdt_manual = QRadioButton("③ Manual table  (dM/dt vs. days to perihelion)")
        self.rad_dmdt_steady.setChecked(True)
        _INACTIVE = f"color:{T['text_muted']};"
        _ACTIVE   = "color:#ffffff; font-weight:bold;"
        def _style_src(*_):
            self.rad_dmdt_steady.setStyleSheet(_ACTIVE if self.rad_dmdt_steady.isChecked() else _INACTIVE)
            self.rad_dmdt_cobs  .setStyleSheet(_ACTIVE if self.rad_dmdt_cobs  .isChecked() else _INACTIVE)
            self.rad_dmdt_manual.setStyleSheet(_ACTIVE if self.rad_dmdt_manual.isChecked() else _INACTIVE)
        self._dmdt_mode_group = QButtonGroup(self)
        for rb in (self.rad_dmdt_steady, self.rad_dmdt_cobs, self.rad_dmdt_manual):
            self._dmdt_mode_group.addButton(rb)
            rb.toggled.connect(_style_src)
            rb.toggled.connect(lambda *_: self._update_run_button_state())
            rb.toggled.connect(lambda *_: self._update_qt_suggestion_summary_labels())
            qt_v.addWidget(rb)
        _style_src()

        self.lbl_steady_qt_note = QLabel(
            "Constant production is used across the selected release window. "
            "No Q(t) curve or light-curve-based MC-window recommendation is required.")
        self.lbl_steady_qt_note.setWordWrap(True)
        self.lbl_steady_qt_note.setStyleSheet(
            f"color:{T['text_muted']}; background:{T['input_bg']}; "
            f"border:1px solid {T['border']}; border-radius:5px; padding:8px;")
        qt_v.addWidget(self.lbl_steady_qt_note)

        # ── Steady: no extra UI needed ──────────────────────────────────

        # ── COBS group ──────────────────────────────────────────────────
        self.dmdt_cobs_group = QWidget()
        cobs_v = QVBoxLayout(self.dmdt_cobs_group)
        cobs_v.setContentsMargins(0, 4, 0, 0); cobs_v.setSpacing(4)
        qt_hint = QLabel(
            "Uses the COBS curve for the SHAPE of production over time,\n"
            "and your own Afρ measurement(s) for the absolute scale.")
        qt_hint.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        cobs_v.addWidget(qt_hint)
        fetch_row = QHBoxLayout()
        self.btn_fetch_cobs = QPushButton("📡 Fetch COBS light curve")
        self.btn_fetch_cobs.clicked.connect(self._fetch_cobs_for_qt)
        self.btn_fetch_cobs.setEnabled(self.main_window is not None)
        self.lbl_cobs_status = QLabel("No COBS data fetched yet.")
        self.lbl_cobs_status.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        fetch_row.addWidget(self.btn_fetch_cobs); fetch_row.addWidget(self.lbl_cobs_status)
        fetch_row.addStretch(); cobs_v.addLayout(fetch_row)
        # Show at least 3 editable anchor rows by default.  Afρ anchors are
        # usually 2–3 measurements, so a single visible row makes input
        # awkward even though the table can scroll.
        self.tbl_anchors = QTableWidget(3, 2)
        self.tbl_anchors.setHorizontalHeaderLabels(["Date (YYYY-MM-DD)", "Afρ [cm]"])
        self.tbl_anchors.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_anchors.verticalHeader().setVisible(False)
        self.tbl_anchors.verticalHeader().setDefaultSectionSize(30)
        self.tbl_anchors.setMinimumWidth(300)
        self.tbl_anchors.setMinimumHeight(120)
        self.tbl_anchors.setMaximumHeight(180)
        for _r in range(3):
            self.tbl_anchors.setItem(_r, 0, QTableWidgetItem(""))
            self.tbl_anchors.setItem(_r, 1, QTableWidgetItem(""))
        cobs_v.addWidget(self.tbl_anchors)
        self.tbl_anchors.itemChanged.connect(lambda *_: self._update_run_button_state())
        self.tbl_anchors.itemChanged.connect(lambda *_: self._update_qt_suggestion_summary_labels())
        anchor_btn_row = QHBoxLayout()
        btn_add_anchor = QPushButton("+ Row"); btn_add_anchor.setMaximumWidth(70)
        btn_add_anchor.clicked.connect(lambda: self.tbl_anchors.insertRow(self.tbl_anchors.rowCount()))
        btn_del_anchor = QPushButton("− Row"); btn_del_anchor.setMaximumWidth(70)
        btn_del_anchor.clicked.connect(self._remove_anchor_row)
        anchor_btn_row.addWidget(btn_add_anchor); anchor_btn_row.addWidget(btn_del_anchor)
        anchor_btn_row.addStretch(); cobs_v.addLayout(anchor_btn_row)
        smooth_row2 = QHBoxLayout()
        smooth_row2.addWidget(QLabel("Smoothing window:"))
        self.sp_qt_smooth = QDoubleSpinBox()
        self.sp_qt_smooth.setRange(0.0, 30.0); self.sp_qt_smooth.setValue(3.0)
        self.sp_qt_smooth.setSingleStep(0.5); self.sp_qt_smooth.setSuffix(" d")
        _sp(self.sp_qt_smooth, W_SMALL)
        smooth_row2.addWidget(self.sp_qt_smooth); smooth_row2.addStretch()
        cobs_v.addLayout(smooth_row2)
        qt_v.addWidget(self.dmdt_cobs_group)
        self.dmdt_cobs_group.setVisible(False)

        # ── Manual table group ──────────────────────────────────────────
        self.dmdt_manual_group = QWidget()
        man_v = QVBoxLayout(self.dmdt_manual_group)
        man_v.setContentsMargins(0, 4, 0, 0); man_v.setSpacing(4)
        man_hint = QLabel(
            "Days to perihelion (negative = pre-perihelion) vs. relative dM/dt.\n"
            "Values are interpolated and normalised to your Afρ scale.")
        man_hint.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        man_v.addWidget(man_hint)
        self.tbl_dmdt_manual = QTableWidget(3, 2)
        self.tbl_dmdt_manual.setHorizontalHeaderLabels(["Days to perihelion", "Relative dM/dt"])
        self.tbl_dmdt_manual.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_dmdt_manual.verticalHeader().setVisible(False)
        self.tbl_dmdt_manual.verticalHeader().setDefaultSectionSize(30)
        self.tbl_dmdt_manual.setMinimumWidth(300)
        self.tbl_dmdt_manual.setMinimumHeight(100)
        self.tbl_dmdt_manual.setMaximumHeight(160)
        for _r in range(3):
            self.tbl_dmdt_manual.setItem(_r, 0, QTableWidgetItem(""))
            self.tbl_dmdt_manual.setItem(_r, 1, QTableWidgetItem(""))
        man_v.addWidget(self.tbl_dmdt_manual)
        man_btn_row = QHBoxLayout()
        btn_add_man = QPushButton("+ Row"); btn_add_man.setMaximumWidth(70)
        btn_add_man.clicked.connect(lambda: self.tbl_dmdt_manual.insertRow(self.tbl_dmdt_manual.rowCount()))
        btn_del_man = QPushButton("− Row"); btn_del_man.setMaximumWidth(70)
        btn_del_man.clicked.connect(
            lambda: self.tbl_dmdt_manual.removeRow(self.tbl_dmdt_manual.currentRow())
                   if self.tbl_dmdt_manual.currentRow() >= 0 else None)
        man_btn_row.addWidget(btn_add_man); man_btn_row.addWidget(btn_del_man)
        man_btn_row.addStretch(); man_v.addLayout(man_btn_row)

        self.tbl_anchors_manual = QTableWidget(1, 2)
        self.tbl_anchors_manual.setHorizontalHeaderLabels(["Date (YYYY-MM-DD)", "Afρ [cm]"])
        self.tbl_anchors_manual.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_anchors_manual.verticalHeader().setVisible(False)
        self.tbl_anchors_manual.verticalHeader().setDefaultSectionSize(30)
        self.tbl_anchors_manual.setMinimumWidth(300)
        self.tbl_anchors_manual.setMinimumHeight(60)
        self.tbl_anchors_manual.setMaximumHeight(90)
        self.tbl_anchors_manual.setItem(0, 0, QTableWidgetItem(""))
        self.tbl_anchors_manual.setItem(0, 1, QTableWidgetItem(""))
        man_v.addWidget(QLabel("Afρ anchor (for absolute scale):"))
        man_v.addWidget(self.tbl_anchors_manual)
        qt_v.addWidget(self.dmdt_manual_group)
        self.dmdt_manual_group.setVisible(False)

        def _on_dmdt_mode():
            is_cobs = self.rad_dmdt_cobs.isChecked()
            is_manual = self.rad_dmdt_manual.isChecked()
            show_preview = is_cobs or is_manual
            self.dmdt_cobs_group.setVisible(is_cobs)
            self.dmdt_manual_group.setVisible(is_manual)
            self.lbl_steady_qt_note.setVisible(not show_preview)
            if hasattr(self, 'qt_analysis_group'):
                self.qt_analysis_group.setVisible(show_preview)
                if show_preview:
                    QTimer.singleShot(0, self._refresh_qt_plot)
        self.rad_dmdt_steady.toggled.connect(_on_dmdt_mode)
        self.rad_dmdt_cobs.toggled.connect(_on_dmdt_mode)
        self.rad_dmdt_manual.toggled.connect(_on_dmdt_mode)

        sim_v.addWidget(grp_qt)
        sim_v.addStretch()

        # ── Size-over-time optional table (keep, at bottom of sim_v) ───
        grp_size_t = QGroupBox("Variable grain size over time (optional)")
        size_t_v = QVBoxLayout(grp_size_t)
        self.chk_size_table = QCheckBox("Enable time-varying grain size")
        size_t_v.addWidget(self.chk_size_table)
        self.size_table_group = QWidget()
        size_tg_v = QVBoxLayout(self.size_table_group)
        size_tg_v.setContentsMargins(0, 0, 0, 0)
        size_hint2 = QLabel(
            "Each row overrides r_min / r_max for particles released in that\n"
            "time window.  Leave empty to use the values above throughout.")
        size_hint2.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        size_tg_v.addWidget(size_hint2)
        self.tbl_size_over_time = QTableWidget(0, 4)
        self.tbl_size_over_time.setHorizontalHeaderLabels(
            ["Start (d to peri)", "End (d to peri)", "r_min (µm)", "r_max (µm)"])
        self.tbl_size_over_time.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_size_over_time.verticalHeader().setVisible(False)
        self.tbl_size_over_time.verticalHeader().setDefaultSectionSize(30)
        self.tbl_size_over_time.setMinimumWidth(300)
        self.tbl_size_over_time.setMinimumHeight(70)
        self.tbl_size_over_time.setMaximumHeight(150)
        size_tg_v.addWidget(self.tbl_size_over_time)
        size_t_btn_row = QHBoxLayout()
        btn_add_size = QPushButton("+ Row"); btn_add_size.setMaximumWidth(70)
        btn_add_size.clicked.connect(
            lambda: self.tbl_size_over_time.insertRow(self.tbl_size_over_time.rowCount()))
        btn_del_size = QPushButton("− Row"); btn_del_size.setMaximumWidth(70)
        btn_del_size.clicked.connect(
            lambda: self.tbl_size_over_time.removeRow(self.tbl_size_over_time.currentRow())
                   if self.tbl_size_over_time.currentRow() >= 0 else None)
        size_t_btn_row.addWidget(btn_add_size); size_t_btn_row.addWidget(btn_del_size)
        size_t_btn_row.addStretch(); size_tg_v.addLayout(size_t_btn_row)
        self.size_table_group.setVisible(False)
        self.chk_size_table.toggled.connect(self.size_table_group.setVisible)
        size_t_v.addWidget(self.size_table_group)
        sim_v.addWidget(grp_size_t)

        # ════════════════════════════════════════════════════════════════
        # TAB 1 — F-P GUIDE
        # ════════════════════════════════════════════════════════════════
        guide_text = QLabel(
            "Workflow:\n"
            "  1. On the main panel, match F-P syndynes/synchrones to the tail by eye\n"
            "     → gives a rough β range and release window estimate (fast).\n"
            "  2. Convert that β range to grain sizes below (using your grain density ρ).\n"
            "  3. Click Apply → go to Simulate tab → Run to see a real particle distribution.")
        guide_text.setStyleSheet(f"color:{T['text_muted']}; font-size:10px;")
        guide_v.addWidget(guide_text)

        fp_box = QGroupBox("Current F-P model settings (main panel)")
        fp_v = QGridLayout(fp_box)
        beta_list_str = ", ".join(f"{b:g}" for b in beta_values) if beta_values else "(none set)"
        self.lbl_fp_betas  = QLabel(f"β values: {beta_list_str}")
        self.lbl_fp_maxage = QLabel(f"Maximum dust age: {self.fp_max_age:g} d")
        dom_text = (f"{self.fp_dominant_age:g} d" if self.fp_dominant_age > 0
                    else "Not specified")
        self.lbl_fp_dominant = QLabel(f"Dominant dust age: {dom_text}")
        fp_v.addWidget(self.lbl_fp_betas, 0, 0)
        fp_v.addWidget(self.lbl_fp_maxage, 0, 1)
        fp_v.addWidget(self.lbl_fp_dominant, 1, 0, 1, 2)
        guide_v.addWidget(fp_box)

        convert_box = QGroupBox("Convert F-P β range → grain-size range")
        conv_v = QGridLayout(convert_box)
        conv_v.addWidget(QLabel("Grain density ρ"), 0, 0)
        conv_v.addWidget(self.sp_rho, 0, 1)
        self.lbl_guide_suggestion = QLabel()
        self.lbl_guide_suggestion.setWordWrap(True)
        self.lbl_guide_suggestion.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        conv_v.addWidget(self.lbl_guide_suggestion, 1, 0, 1, 2)
        self.btn_apply_fp_guess = QPushButton("Apply to Simulate tab →")
        conv_v.addWidget(self.btn_apply_fp_guess, 2, 0, 1, 2)
        guide_v.addWidget(convert_box)

        assumptions_box = QGroupBox("Brightness / phase assumptions")
        bright_v = QGridLayout(assumptions_box)
        bright_v.setHorizontalSpacing(10); bright_v.setVerticalSpacing(6)

        self.sp_pv = _sp(_dspin(0.04, 0.001, 1.0, 3, 0.005), W_SMALL)
        self.sp_pv.setToolTip("Grain geometric albedo p_v (default 0.04).")

        self.cmb_phase_law = QComboBox()
        self.cmb_phase_law.addItems([
            "Schleicher composite dust phase function",
            "Linear–exponential custom",
            "None / relative morphology only",
        ])
        self.cmb_phase_law.setToolTip(
            "Phase correction is a single frame-wide brightness scalar.\n"
            "It changes absolute brightness comparison, not one run's contour shape.")

        self.sp_phase_beta = _sp(_dspin(0.024, 0.000, 0.200, 3, 0.001, " mag/deg"), W_WIDE)
        self.sp_phase_moe  = _sp(_dspin(0.28,  0.000, 2.000, 3, 0.010, " mag"), W_MED)
        self.sp_phase_woe  = _sp(_dspin(1.5,   0.050, 30.000, 2, 0.100, "°"), W_MED)
        for w in (self.sp_phase_beta, self.sp_phase_moe, self.sp_phase_woe):
            w.setToolTip("Used only when Phase law = Linear–exponential custom.")

        self.phase_custom_group = QWidget()
        pcg = QHBoxLayout(self.phase_custom_group)
        pcg.setContentsMargins(0, 0, 0, 0); pcg.setSpacing(6)
        pcg.addWidget(QLabel("β_α")); pcg.addWidget(self.sp_phase_beta)
        pcg.addWidget(QLabel("m_oe")); pcg.addWidget(self.sp_phase_moe)
        pcg.addWidget(QLabel("w_oe")); pcg.addWidget(self.sp_phase_woe)
        pcg.addStretch()
        self.phase_custom_group.setVisible(False)

        def _on_phase_law_changed(idx):
            self.phase_custom_group.setVisible(idx == 1)
        self.cmb_phase_law.currentIndexChanged.connect(_on_phase_law_changed)

        bright_v.addWidget(QLabel("Geometric albedo p_v"), 0, 0)
        bright_v.addWidget(self.sp_pv, 0, 1)
        bright_v.addWidget(QLabel("Phase law"), 1, 0)
        bright_v.addWidget(self.cmb_phase_law, 1, 1, 1, 3)
        bright_v.addWidget(self.phase_custom_group, 2, 1, 1, 3)
        guide_v.addWidget(assumptions_box)

        input_box = QGroupBox("Input-file driven configurations")
        input_v = QVBoxLayout(input_box)
        input_note = QLabel(
            "Recommended cases should be distributed as .mcin input files, not hidden GUI presets.\n"
            "Use Load inputs… to apply a literature-guided configuration, then edit any visible field.\n"
            "Save inputs… writes the same transparent JSON format for sharing or publication support.")
        input_note.setWordWrap(True)
        input_note.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        input_v.addWidget(input_note)
        fmt_note = QLabel(
            "Format stores: schema, profile metadata, target/date, interpretation notes, "
            "and a generic snapshot of all visible MC inputs. Files without schema='CTA_MC_INPUT' are rejected.")
        fmt_note.setWordWrap(True)
        fmt_note.setStyleSheet(f"color:{T['text_value']}; font-size:9px;")
        input_v.addWidget(fmt_note)
        guide_v.addWidget(input_box)

        guide_v.addStretch()

        def _refresh_guide_suggestion():
            try:
                rho = self.sp_rho.value()
                a_max_um = cta._beta_to_radius_um(beta_lo, rho)
                a_min_um = cta._beta_to_radius_um(beta_hi, rho)
                self.lbl_guide_suggestion.setText(
                    f"→ Suggested: r_min={a_min_um:.3g} µm,  r_max={a_max_um:.3g} µm "
                    f"(β=[{beta_lo:g}, {beta_hi:g}], ρ={rho:.2g} g/cm³)")
            except Exception:
                self.lbl_guide_suggestion.setText("")
        self.sp_rho.valueChanged.connect(_refresh_guide_suggestion)
        _refresh_guide_suggestion()

        def _apply_fp_guess():
            try:
                rho = self.sp_rho.value()
                self.sp_r_min.setValue(cta._beta_to_radius_um(beta_hi, rho))
                self.sp_r_max.setValue(cta._beta_to_radius_um(beta_lo, rho))
                tabs.setCurrentWidget(page_sim)
            except Exception:
                pass
        self.btn_apply_fp_guess.clicked.connect(_apply_fp_guess)

        # ════════════════════════════════════════════════════════════════
        # TAB 3 — DISPLAY
        # ════════════════════════════════════════════════════════════════
        grp_contour = QGroupBox("Contour extraction")
        cv = QGridLayout(grp_contour); cv.setHorizontalSpacing(10); cv.setVerticalSpacing(6)

        self.sp_floor = _sp(_dspin(80.0, 0.0, 99.9, 1, 1.0, " %ile"), W_MED)
        self.sp_floor.setToolTip(
            "Percentile of non-zero density bins used as the innermost\n"
            "contour threshold. Used when 'Percentile mode' is selected.\n"
            "Higher = fewer, tighter contours near the dense coma.")
        self.sp_nlevels = _sp(QSpinBox(), W_SMALL); self.sp_nlevels.setRange(1, 20)
        self.sp_nlevels.setValue(3)
        self.sp_nlevels.setToolTip(
            "Number of contour rings to trace — applies in both\n"
            "Percentile and Calibrated mag/arcsec² modes.")

        # Calibrated surface brightness mode — fixed surface-brightness isophote levels
        self.chk_sb_mode = QCheckBox("Calibrated mag/arcsec² mode")
        self.chk_sb_mode.setChecked(False)
        self.chk_sb_mode.setStyleSheet("color:#66ccff;")
        self.chk_sb_mode.setToolTip(
            "When checked: extract contours at specific surface brightness\n"
            "levels [mag/arcsec²] calibrated using your Afρ measurement.\n"
            "Use fixed surface-brightness isophotes at the same magnitude\n"
            "levels as the observed image.\n\n"
            "Requires: Afρ value from Image Overlay → Analysis panel,\n"
            "and image plate scale (AU/px) set correctly.\n\n"
            "When unchecked: percentile mode (original CTA method).")
        cv.addWidget(self.chk_sb_mode, 0, 0, 1, 2)

        self.sb_mode_group = QWidget()
        sbg = QGridLayout(self.sb_mode_group); sbg.setContentsMargins(0,0,0,0); sbg.setHorizontalSpacing(10)

        btn_load_tycho = QPushButton("📂 Load Tycho photometry…")
        btn_load_tycho.setToolTip(
            "Parse a Tycho Tracker ICQ-format export (.txt) and auto-fill\n"
            "Afρ + Faintest/Brightest levels from its radius-vs-magnitude\n"
            "profile table. If the file can't be parsed, current values\n"
            "(or the defaults below) are left untouched.")
        btn_load_tycho.clicked.connect(self._load_tycho_photometry)
        self.lbl_tycho_status = QLabel(
            "No file loaded — using defaults below (or COBS-anchor Afρ).")
        self.lbl_tycho_status.setWordWrap(True)
        self.lbl_tycho_status.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        sbg.addWidget(btn_load_tycho, 0, 0, 1, 2)
        sbg.addWidget(self.lbl_tycho_status, 1, 0, 1, 2)

        self.sp_afrho_cal = _sp(_dspin(500.0, 0.1, 1e7, 1, 10.0, " cm"), W_MED)
        self.sp_afrho_cal.setToolTip(
            "Observed Afρ [cm] at the observation date.\n"
            "Auto-filled from COBS anchor table or a loaded Tycho file\n"
            "when available. Otherwise enter manually.")
        self.lbl_afrho_source = QLabel("(enter manually)")
        self.lbl_afrho_source.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        self.sp_mag_faint = _sp(_dspin(23.0, 10.0, 35.0, 1, 0.5, " mag/\"²"), W_WIDE)
        self.sp_mag_faint.setToolTip(
            "Faintest isophote level [mag/arcsec²].\n"
            "Match the outermost visible isophote of your image.\n"
            "Auto-filled by 'Load Tycho photometry…' above, or enter\n"
            "manually by reading it off your image processing software.")
        self.sp_mag_bright = _sp(_dspin(19.0, 10.0, 35.0, 1, 0.5, " mag/\"²"), W_WIDE)
        self.sp_mag_bright.setToolTip(
            "Brightest isophote level [mag/arcsec²].\n"
            "Match the innermost isophote of your image.\n"
            "Auto-filled by 'Load Tycho photometry…' above.")
        afrho_row = QHBoxLayout()
        afrho_row.addWidget(self.sp_afrho_cal)
        afrho_row.addWidget(self.lbl_afrho_source)
        sbg.addWidget(QLabel("Observed Afρ"), 2, 0); sbg.addLayout(afrho_row, 2, 1)
        sbg.addWidget(QLabel("Faintest level"), 3, 0); sbg.addWidget(self.sp_mag_faint, 3, 1)
        sbg.addWidget(QLabel("Brightest level"), 4, 0); sbg.addWidget(self.sp_mag_bright, 4, 1)
        self.sb_mode_group.setVisible(False)
        cv.addWidget(self.sb_mode_group, 1, 0, 1, 2)

        self.percentile_group = QWidget()
        pg = QGridLayout(self.percentile_group); pg.setContentsMargins(0,0,0,0); pg.setHorizontalSpacing(10)
        pg.addWidget(QLabel("Sensitivity (percentile floor)"), 0, 0); pg.addWidget(self.sp_floor, 0, 1)
        cv.addWidget(self.percentile_group, 1, 0, 1, 2)

        def _on_sb_mode(checked):
            self.sb_mode_group.setVisible(checked)
            self.percentile_group.setVisible(not checked)
            if checked:
                self._autofill_afrho()
        self.chk_sb_mode.toggled.connect(_on_sb_mode)

        # sp_nlevels applies to BOTH modes — placed once here so it's
        # never fought over between the two mode-specific sub-layouts
        # (a Qt widget can only live in one parent layout at a time;
        # adding it to both silently strips it from whichever group
        # loses, leaving that group's "Contour rings" row with no field).
        cv.addWidget(QLabel("Number of contour rings"), 2, 0)
        cv.addWidget(self.sp_nlevels, 2, 1)

        auto_note = QLabel("Contour smoothing: auto (1× grid pixel size)")
        auto_note.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        cv.addWidget(auto_note, 3, 0, 1, 2)
        disp_v.addWidget(grp_contour)

        self.btn_reextract = QPushButton("↻ Re-extract contour  (fast, no re-sample)")
        self.btn_reextract.setEnabled(False)
        self.btn_reextract.setToolTip(
            "Re-apply the threshold/levels above to the LAST MC density grid\n"
            "without re-simulating particles — typically < 0.1s.")
        self.btn_reextract.clicked.connect(self._reextract_contours)
        disp_v.addWidget(self.btn_reextract)

        # ── Presentation preset ─────────────────────────────────────────
        grp_present = QGroupBox("PRESENTATION MODE")
        pv = QGridLayout(grp_present); pv.setHorizontalSpacing(10); pv.setVerticalSpacing(6)
        self.cmb_mc_view_style = QComboBox()
        self.cmb_mc_view_style.addItem("Analysis Overlay", "analysis")
        self.cmb_mc_view_style.addItem("Contour Comparison", "contour")
        self.cmb_mc_view_style.addItem("Publication Figure", "publication")
        self.cmb_mc_view_style.setToolTip(
            "Analysis Overlay keeps the comet image and all diagnostic overlays.\n"
            "Contour Comparison hides the raster image and compares observed\n"
            "isophotes with the MC model on white. Publication Figure applies\n"
            "paper-oriented axes, line weights, and annotation defaults.")
        self.lbl_mc_view_help = QLabel(
            "Analysis Overlay: full diagnostic view on the image background.")
        self.lbl_mc_view_help.setWordWrap(True)
        self.lbl_mc_view_help.setStyleSheet(
            f"color:{T['text_muted']}; font-size:9px;")
        pv.addWidget(QLabel("View style"), 0, 0)
        pv.addWidget(self.cmb_mc_view_style, 0, 1)
        pv.addWidget(self.lbl_mc_view_help, 1, 0, 1, 2)
        disp_v.addWidget(grp_present)

        # ── Observed + model contour appearance ────────────────────────
        grp_appear = QGroupBox("CONTOUR APPEARANCE")
        av = QGridLayout(grp_appear); av.setHorizontalSpacing(10); av.setVerticalSpacing(6)

        def _add_color_items(combo, items):
            for label, color in items:
                combo.addItem(label, color)

        def _add_line_items(combo):
            combo.addItem("Solid", "-")
            combo.addItem("Dashed", "--")
            combo.addItem("Dash-dot", "-.")
            combo.addItem("Dotted", ":")

        self.chk_observed_show = QCheckBox("Show observed isophotes")
        self.chk_observed_show.setChecked(True)
        self.chk_observed_show.setToolTip(
            "In contour-only modes this directly shows/hides isophotes derived\n"
            "from the loaded image. In Analysis Overlay the main DISPLAY\n"
            "checkbox 'Isophotes (from image)' remains the master switch.")
        self.cmb_observed_color = QComboBox()
        _add_color_items(self.cmb_observed_color, [
            ("Green (analysis)", ISOPHOTE_COLOR),
            ("Black", "black"),
            ("White", "white"),
            ("Blue", "#1f5fbf"),
        ])
        self.cmb_observed_ls = QComboBox(); _add_line_items(self.cmb_observed_ls)
        self.sp_observed_lw = _dspin(0.6, 0.2, 5.0, 1, 0.1)
        self.sp_observed_lw.setSuffix(" pt")
        self.sp_observed_lw.setMaximumWidth(W_SMALL + 18)
        self.obs_op_slider = QSlider(Qt.Orientation.Horizontal)
        self.obs_op_slider.setRange(10, 100); self.obs_op_slider.setValue(75)
        self.obs_op_label = QLabel("0.75"); self.obs_op_label.setMinimumWidth(32)

        av.addWidget(QLabel("Observed isophotes"), 0, 0, 1, 3)
        av.addWidget(self.chk_observed_show, 1, 0, 1, 3)
        av.addWidget(QLabel("Color"), 2, 0); av.addWidget(self.cmb_observed_color, 2, 1, 1, 2)
        av.addWidget(QLabel("Line style"), 3, 0); av.addWidget(self.cmb_observed_ls, 3, 1, 1, 2)
        av.addWidget(QLabel("Width"), 4, 0); av.addWidget(self.sp_observed_lw, 4, 1, 1, 2)
        av.addWidget(QLabel("Opacity"), 5, 0); av.addWidget(self.obs_op_slider, 5, 1); av.addWidget(self.obs_op_label, 5, 2)

        self.cmb_model_color = QComboBox()
        _add_color_items(self.cmb_model_color, [
            ("Magenta", "#ff3399"),
            ("Red", "#e31a1c"),
            ("Black", "black"),
            ("Blue", "#1f78b4"),
        ])
        self.cmb_model_ls = QComboBox(); _add_line_items(self.cmb_model_ls)
        self.mc_lw_slider = QSlider(Qt.Orientation.Horizontal)
        self.mc_lw_slider.setRange(5, 50); self.mc_lw_slider.setValue(12)
        self.mc_lw_label = QLabel("1.2"); self.mc_lw_label.setMinimumWidth(30)
        self.mc_op_slider = QSlider(Qt.Orientation.Horizontal)
        self.mc_op_slider.setRange(10, 100); self.mc_op_slider.setValue(95)
        self.mc_op_label = QLabel("0.95"); self.mc_op_label.setMinimumWidth(30)

        av.addWidget(QLabel("Monte Carlo model contours"), 6, 0, 1, 3)
        av.addWidget(QLabel("Color"), 7, 0); av.addWidget(self.cmb_model_color, 7, 1, 1, 2)
        av.addWidget(QLabel("Line style"), 8, 0); av.addWidget(self.cmb_model_ls, 8, 1, 1, 2)
        av.addWidget(QLabel("Width"), 9, 0); av.addWidget(self.mc_lw_slider, 9, 1); av.addWidget(self.mc_lw_label, 9, 2)
        av.addWidget(QLabel("Opacity"), 10, 0); av.addWidget(self.mc_op_slider, 10, 1); av.addWidget(self.mc_op_label, 10, 2)

        self.chk_grayscale_safe = QCheckBox("Grayscale-safe line styles")
        self.chk_grayscale_safe.setToolTip(
            "Force observed contours to black solid and model contours to\n"
            "black dashed, so the figure remains distinguishable in grayscale.")
        av.addWidget(self.chk_grayscale_safe, 11, 0, 1, 3)
        disp_v.addWidget(grp_appear)

        # ── Figure layout + export ──────────────────────────────────────
        grp_layout = QGroupBox("FIGURE LAYOUT")
        fv = QGridLayout(grp_layout); fv.setHorizontalSpacing(10); fv.setVerticalSpacing(6)
        self.cmb_mc_background = QComboBox()
        self.cmb_mc_background.addItem("Comet image", "image")
        self.cmb_mc_background.addItem("White", "white")
        self.cmb_mc_coordinates = QComboBox()
        self.cmb_mc_coordinates.addItem("Image pixels", "pixels")
        self.cmb_mc_coordinates.addItem("Projected distance (km)", "projected_km")
        self.cmb_mc_coordinates.addItem("Hide axes", "hide")
        self.chk_mc_legend = QCheckBox("Show legend"); self.chk_mc_legend.setChecked(True)
        self.chk_mc_title = QCheckBox("Show title"); self.chk_mc_title.setChecked(True)
        self.chk_mc_compass = QCheckBox("Show N/E compass"); self.chk_mc_compass.setChecked(True)
        self.chk_mc_info = QCheckBox("Show MC parameters box")
        self.chk_mc_info.setChecked(True)
        self.btn_export_mc_figure = QPushButton("Export figure…")
        self.btn_export_mc_figure.setToolTip(
            "Export the current main-canvas figure as PNG, TIFF, PDF, or SVG.\n"
            "The plot toolbar and application interface are not included.")
        self.btn_export_mc_figure.clicked.connect(self._export_current_mc_figure)

        fv.addWidget(QLabel("Background"), 0, 0); fv.addWidget(self.cmb_mc_background, 0, 1)
        fv.addWidget(QLabel("Coordinates"), 1, 0); fv.addWidget(self.cmb_mc_coordinates, 1, 1)
        fv.addWidget(self.chk_mc_legend, 2, 0)
        fv.addWidget(self.chk_mc_title, 2, 1)
        fv.addWidget(self.chk_mc_compass, 3, 0)
        fv.addWidget(self.chk_mc_info, 3, 1)
        fv.addWidget(self.btn_export_mc_figure, 4, 0, 1, 2)
        disp_v.addWidget(grp_layout)

        def _set_combo_data(combo, value):
            idx = combo.findData(value)
            if idx >= 0:
                combo.setCurrentIndex(idx)

        self._applying_presentation_preset = False

        def _apply_presentation_preset(index):
            if self._applying_presentation_preset:
                return
            self._applying_presentation_preset = True
            try:
                mode = self.cmb_mc_view_style.itemData(index)
                if mode == "analysis":
                    self.lbl_mc_view_help.setText(
                        "Analysis Overlay: full diagnostic view on the image background.")
                    _set_combo_data(self.cmb_mc_background, "image")
                    _set_combo_data(self.cmb_mc_coordinates, "pixels")
                    _set_combo_data(self.cmb_observed_color, ISOPHOTE_COLOR)
                    _set_combo_data(self.cmb_observed_ls, "-")
                    _set_combo_data(self.cmb_model_color, "#ff3399")
                    _set_combo_data(self.cmb_model_ls, "--")
                    self.sp_observed_lw.setValue(0.6)
                    self.obs_op_slider.setValue(75)
                    self.mc_lw_slider.setValue(12)
                    self.mc_op_slider.setValue(95)
                    self.chk_mc_legend.setChecked(True)
                    self.chk_mc_title.setChecked(True)
                    self.chk_mc_compass.setChecked(True)
                    self.chk_mc_info.setChecked(True)
                    self.chk_grayscale_safe.setChecked(False)
                elif mode == "contour":
                    self.lbl_mc_view_help.setText(
                        "Contour Comparison: observed isophotes and MC contours only, on white.")
                    _set_combo_data(self.cmb_mc_background, "white")
                    _set_combo_data(self.cmb_mc_coordinates, "pixels")
                    _set_combo_data(self.cmb_observed_color, "black")
                    _set_combo_data(self.cmb_observed_ls, "-")
                    _set_combo_data(self.cmb_model_color, "#ff3399")
                    _set_combo_data(self.cmb_model_ls, "-")
                    self.sp_observed_lw.setValue(1.0)
                    self.obs_op_slider.setValue(100)
                    self.mc_lw_slider.setValue(12)
                    self.mc_op_slider.setValue(100)
                    self.chk_mc_legend.setChecked(True)
                    self.chk_mc_title.setChecked(True)
                    self.chk_mc_compass.setChecked(False)
                    self.chk_mc_info.setChecked(False)
                    self.chk_grayscale_safe.setChecked(False)
                else:
                    self.lbl_mc_view_help.setText(
                        "Publication Figure: paper-oriented white background and projected-distance axes.")
                    _set_combo_data(self.cmb_mc_background, "white")
                    _set_combo_data(self.cmb_mc_coordinates, "projected_km")
                    _set_combo_data(self.cmb_observed_color, "black")
                    _set_combo_data(self.cmb_observed_ls, "-")
                    _set_combo_data(self.cmb_model_color, "#ff3399")
                    _set_combo_data(self.cmb_model_ls, "-")
                    self.sp_observed_lw.setValue(1.2)
                    self.obs_op_slider.setValue(100)
                    self.mc_lw_slider.setValue(12)
                    self.mc_op_slider.setValue(100)
                    self.chk_mc_legend.setChecked(True)
                    self.chk_mc_title.setChecked(False)
                    self.chk_mc_compass.setChecked(False)
                    self.chk_mc_info.setChecked(False)
                    self.chk_grayscale_safe.setChecked(False)
            finally:
                self._applying_presentation_preset = False
            _push_mc_style()

        def _push_mc_style(*_):
            self.mc_lw_label.setText(f"{self.mc_lw_slider.value()/10:.1f}")
            self.mc_op_label.setText(f"{self.mc_op_slider.value()/100:.2f}")
            self.obs_op_label.setText(f"{self.obs_op_slider.value()/100:.2f}")
            if self.main_window is None:
                return
            self.main_window.ctrl._mc_style = dict(
                lw=self.mc_lw_slider.value()/10,
                alpha=self.mc_op_slider.value()/100,
                info_box=self.chk_mc_info.isChecked(),
                view_style=self.cmb_mc_view_style.currentData(),
                background=self.cmb_mc_background.currentData(),
                coordinates=self.cmb_mc_coordinates.currentData(),
                observed_show=self.chk_observed_show.isChecked(),
                observed_color=self.cmb_observed_color.currentData(),
                observed_lw=self.sp_observed_lw.value(),
                observed_alpha=self.obs_op_slider.value()/100,
                observed_ls=self.cmb_observed_ls.currentData(),
                model_color=self.cmb_model_color.currentData(),
                model_ls=self.cmb_model_ls.currentData(),
                show_legend=self.chk_mc_legend.isChecked(),
                show_title=self.chk_mc_title.isChecked(),
                show_compass=self.chk_mc_compass.isChecked(),
                grayscale_safe=self.chk_grayscale_safe.isChecked())
            mw = self.main_window
            mw.canvas._vis = mw.ctrl.get_vis()
            if mw._model:
                ov = mw.ctrl.get_overlay()
                mw.canvas.draw_model(
                    mw._model, ov["img_arr"] if mw.ctrl._img_arr is not None else None)

        self.cmb_mc_view_style.currentIndexChanged.connect(_apply_presentation_preset)
        for sig in [
                self.mc_lw_slider.valueChanged,
                self.mc_op_slider.valueChanged,
                self.obs_op_slider.valueChanged,
                self.sp_observed_lw.valueChanged,
                self.chk_mc_info.toggled,
                self.chk_observed_show.toggled,
                self.chk_grayscale_safe.toggled,
                self.chk_mc_legend.toggled,
                self.chk_mc_title.toggled,
                self.chk_mc_compass.toggled,
                self.cmb_mc_background.currentIndexChanged,
                self.cmb_mc_coordinates.currentIndexChanged,
                self.cmb_observed_color.currentIndexChanged,
                self.cmb_observed_ls.currentIndexChanged,
                self.cmb_model_color.currentIndexChanged,
                self.cmb_model_ls.currentIndexChanged]:
            sig.connect(_push_mc_style)

        # Restore the current main-window style without forcing a preset.
        if self.main_window is not None:
            st = self.main_window.ctrl._mc_style
            _set_combo_data(self.cmb_mc_view_style, st.get("view_style", "analysis"))
            _set_combo_data(self.cmb_mc_background, st.get("background", "image"))
            _set_combo_data(self.cmb_mc_coordinates, st.get("coordinates", "pixels"))
            _set_combo_data(self.cmb_observed_color, st.get("observed_color", ISOPHOTE_COLOR))
            _set_combo_data(self.cmb_observed_ls, st.get("observed_ls", "-"))
            _set_combo_data(self.cmb_model_color, st.get("model_color", "#ff3399"))
            _set_combo_data(self.cmb_model_ls, st.get("model_ls", "--"))
            self.chk_observed_show.setChecked(st.get("observed_show", True))
            self.sp_observed_lw.setValue(st.get("observed_lw", 0.6))
            self.obs_op_slider.setValue(round(st.get("observed_alpha", 0.75)*100))
            self.mc_lw_slider.setValue(round(st.get("lw", 1.2)*10))
            self.mc_op_slider.setValue(round(st.get("alpha", 0.95)*100))
            self.chk_mc_info.setChecked(st.get("info_box", True))
            self.chk_mc_legend.setChecked(st.get("show_legend", True))
            self.chk_mc_title.setChecked(st.get("show_title", True))
            self.chk_mc_compass.setChecked(st.get("show_compass", True))
            self.chk_grayscale_safe.setChecked(st.get("grayscale_safe", False))

        self.mc_lw_label.setText(f"{self.mc_lw_slider.value()/10:.1f}")
        self.mc_op_label.setText(f"{self.mc_op_slider.value()/100:.2f}")
        self.obs_op_label.setText(f"{self.obs_op_slider.value()/100:.2f}")
        _push_mc_style()

        disp_v.addStretch()

        # ════════════════════════════════════════════════════════════════
        # ASSEMBLE TABS + BOTTOM BAR
        # ════════════════════════════════════════════════════════════════
        # ════════════════════════════════════════════════════════════════
        # Q(t) PREVIEW + MC-WINDOW GUIDANCE
        # Embedded inside Dust production over time.  It is shown only for
        # COBS or manual dM/dt sources; steady production needs no preview.
        # ════════════════════════════════════════════════════════════════
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

        self.qt_analysis_group = QGroupBox("Q(t) preview and MC-window guidance")
        qt_prev_v = QVBoxLayout(self.qt_analysis_group)
        qt_prev_v.setSpacing(8)

        qt_intro = QLabel(
            "The curve below is the release-time weighting used by the Monte Carlo sampler. "
            "Use the recommendation as a starting point, then validate the final window "
            "against the observed morphology.")
        qt_intro.setWordWrap(True)
        qt_intro.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        qt_prev_v.addWidget(qt_intro)

        # F-P morphology guidance comes from the main panel. The user
        # supplies one representative dominant age; maximum dust age is
        # automatically derived from the largest listed synchrone age.
        grp_fp_age = QGroupBox("F-P morphology guidance  (automatic)")
        fpag = QGridLayout(grp_fp_age)
        fpag.setHorizontalSpacing(12); fpag.setVerticalSpacing(6)

        fpag.addWidget(QLabel("Dominant dust age"), 0, 0)
        self.lbl_fp_dom_guidance = QLabel()
        self.lbl_fp_dom_guidance.setStyleSheet(
            "color:#f4f8ff; font-weight:600;")
        fpag.addWidget(self.lbl_fp_dom_guidance, 0, 1)

        fpag.addWidget(QLabel("Maximum dust age"), 1, 0)
        self.lbl_fp_max_guidance = QLabel()
        self.lbl_fp_max_guidance.setStyleSheet(
            "color:#f4f8ff; font-weight:600;")
        fpag.addWidget(self.lbl_fp_max_guidance, 1, 1)

        fp_note = QLabel(
            "Dominant age is your F-P morphology interpretation. Maximum dust age "
            "is the largest synchrone age and is applied as the upper limit for "
            "MC-window guidance. Q(t) supplies the release-time weighting.")
        fp_note.setWordWrap(True)
        fp_note.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        fpag.addWidget(fp_note, 2, 0, 1, 2)
        qt_prev_v.addWidget(grp_fp_age)
        self._update_fp_guidance_labels()

        self._qt_suggested_window = None
        self._qt_suggested_window_range = None
        self._qt_suggested_window_note = ""
        self._qt_qonly_window = None
        self._qt_fp_adjustment = ""
        self._qt_coverage = None
        self._qt_recommendation_quality = ""
        self._qt_provisional_window = None

        suggest_card = QFrame()
        suggest_card.setStyleSheet(
            "QFrame { background:#13202b; border:1px solid #4a8adf; "
            "border-radius:7px; } QLabel { border:none; background:transparent; }")
        suggest_layout = QHBoxLayout(suggest_card)
        suggest_layout.setContentsMargins(12, 10, 12, 10)
        suggest_layout.setSpacing(12)
        self.lbl_qt_suggest_summary = QLabel(
            "MC-window recommendation not calculated yet.\n"
            "Enter/fetch Q(t) data, then refresh the preview.")
        self.lbl_qt_suggest_summary.setWordWrap(True)
        self.lbl_qt_suggest_summary.setMinimumHeight(54)
        self.lbl_qt_suggest_summary.setStyleSheet(
            "color:#f4f8ff; font-size:13px; font-weight:700;")
        suggest_layout.addWidget(self.lbl_qt_suggest_summary, 1)

        self.btn_apply_qt_suggested_window = QPushButton("Use recommended window")
        self.btn_apply_qt_suggested_window.setEnabled(False)
        self.btn_apply_qt_suggested_window.setAutoDefault(False)
        self.btn_apply_qt_suggested_window.setDefault(False)
        self.btn_apply_qt_suggested_window.setMinimumWidth(190)
        self.btn_apply_qt_suggested_window.setToolTip(
            "Copy the recommended value into Simulation → Release window. "
            "The value remains editable and is not a unique physical solution.")
        self.btn_apply_qt_suggested_window.clicked.connect(
            self._apply_qt_suggested_window)
        suggest_layout.addWidget(self.btn_apply_qt_suggested_window, 0, Qt.AlignmentFlag.AlignVCenter)
        qt_prev_v.addWidget(suggest_card)

        self._qt_fig  = Figure(figsize=(5, 3.2), facecolor=T['panel_bg'])
        self._qt_canvas = FigureCanvasQTAgg(self._qt_fig)
        self._qt_canvas.setMinimumHeight(300)
        qt_prev_v.addWidget(self._qt_canvas, 1)

        self._qt_info_lbl = QLabel("")
        self._qt_info_lbl.setWordWrap(True)
        self._qt_info_lbl.setStyleSheet(f"color:{T['text_muted']}; font-size:9px;")
        qt_prev_v.addWidget(self._qt_info_lbl)

        btn_qt_refresh = QPushButton("↻  Refresh Q(t) plot")
        btn_qt_refresh.setToolTip(
            "Re-draw the Q(t) preview using current COBS data and\n"
            "anchor values. Fetch COBS data first in the COBS source section above.")
        btn_qt_refresh.clicked.connect(self._refresh_qt_plot)
        qt_prev_v.addWidget(btn_qt_refresh, 0, Qt.AlignmentFlag.AlignLeft)

        qt_v.addWidget(self.qt_analysis_group)
        self.qt_analysis_group.setVisible(False)

        # Debounced live refresh: source inputs, F-P constraints, and the
        # current release window all affect the recommendation.  Waiting a
        # short moment avoids recomputing on every keystroke while a table
        # cell is still being edited.
        self._qt_refresh_timer = QTimer(self)
        self._qt_refresh_timer.setSingleShot(True)
        self._qt_refresh_timer.setInterval(350)
        self._qt_refresh_timer.timeout.connect(self._refresh_qt_plot)

        def _schedule_qt_refresh(*_):
            self._update_qt_suggestion_summary_labels()
            if self.rad_dmdt_cobs.isChecked() or self.rad_dmdt_manual.isChecked():
                self._qt_refresh_timer.start()

        self.sp_max_age.valueChanged.connect(_schedule_qt_refresh)
        self.sp_qt_smooth.valueChanged.connect(_schedule_qt_refresh)
        self.tbl_anchors.itemChanged.connect(_schedule_qt_refresh)
        self.tbl_dmdt_manual.itemChanged.connect(_schedule_qt_refresh)

        _on_dmdt_mode()

        tabs.addTab(page_guide, "1 · F-P Guide")
        tabs.addTab(page_sim,   "2 · Simulation")
        tabs.addTab(page_disp,  "3 · Display")
        tabs.currentChanged.connect(self._on_mc_tab_changed)

        tabs_scroll = QScrollArea()
        tabs_scroll.setWidgetResizable(True)
        tabs_scroll.setWidget(tabs)
        tabs_scroll.setFrameShape(QFrame.Shape.NoFrame)
        vbox.addWidget(tabs_scroll, 1)

        # Guided Back/Next navigation. Tab headers remain directly clickable,
        # but first-time users no longer need to discover the workflow by trial.
        nav = QHBoxLayout()
        self.btn_step_back = QPushButton("← Back")
        self.btn_step_back.setAutoDefault(False)
        self.btn_step_back.setDefault(False)
        self.btn_step_back.clicked.connect(self._navigate_mc_back)
        self.lbl_step = QLabel("Step 1 of 3 — F-P Guide")
        self.lbl_step.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_step.setStyleSheet(
            f"color:{T['text_value']}; font-size:10px; font-weight:600;")
        self.btn_step_next = QPushButton("Next: Simulation →")
        self.btn_step_next.setAutoDefault(False)
        self.btn_step_next.setDefault(False)
        self.btn_step_next.clicked.connect(self._navigate_mc_next)
        nav.addWidget(self.btn_step_back)
        nav.addStretch()
        nav.addWidget(self.lbl_step)
        nav.addStretch()
        nav.addWidget(self.btn_step_next)
        vbox.addLayout(nav)

        # Run is intentionally shown only on Simulation. Display is optional
        # and no longer has to be visited before a valid model can run.
        self.btn_run = QPushButton(
            "▶  Run Monte Carlo  (Ctrl+Enter)  →  contour appears on main canvas")
        self.btn_run.setAutoDefault(False)
        self.btn_run.setDefault(False)
        self.btn_run.setMinimumHeight(38)
        self.btn_run.setStyleSheet(
            "QPushButton { background:#1a6e2e; color:white; font-weight:bold; border-radius:4px; }"
            "QPushButton:hover { background:#21883a; }"
            "QPushButton:disabled { background:#2a3a2e; color:#666; }")
        self.btn_run.clicked.connect(self._run)
        vbox.addWidget(self.btn_run)

        # Dialog-local shortcuts: Ctrl+Enter runs only when validation passes;
        # Ctrl+R re-extracts contours only when a cached MC result exists.
        self._act_run_mc_shortcut = QAction(self)
        self._act_run_mc_shortcut.setShortcuts(["Ctrl+Return", "Ctrl+Enter"])
        self._act_run_mc_shortcut.setShortcutContext(
            Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self._act_run_mc_shortcut.triggered.connect(self._run_mc_from_shortcut)
        self.addAction(self._act_run_mc_shortcut)

        self._act_reextract_shortcut = QAction(self)
        self._act_reextract_shortcut.setShortcut("Ctrl+R")
        self._act_reextract_shortcut.setShortcutContext(
            Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self._act_reextract_shortcut.triggered.connect(self._reextract_contours)
        self.addAction(self._act_reextract_shortcut)

        self.btn_skip_guide.clicked.connect(self._skip_guided_setup)

        # ── Progress bar + status (below run button) ─────────────────────
        self.progress = QProgressBar()
        self.progress.setRange(0, 100); self.progress.setValue(0)
        self.progress.setFormat("0%")
        self.progress.setMaximumHeight(18); self.progress.setTextVisible(True)
        vbox.addWidget(self.progress)

        self.lbl_status = QLabel("Set parameters, then click Run.")
        self.lbl_status.setWordWrap(True)
        self.lbl_status.setStyleSheet(f"color:{T['text_muted']}; font-size:10px;")
        vbox.addWidget(self.lbl_status)

        # ── Bottom action bar ─────────────────────────────────────────────
        bot = QHBoxLayout()
        btn_reset = QPushButton("↺ Reset")
        btn_reset.setToolTip("Reset all fields to factory defaults.")
        btn_reset.clicked.connect(self._reset_to_defaults)
        btn_save_input = QPushButton("💾 Save inputs as…")
        btn_save_input.setToolTip("Save the current edited GUI inputs as a new .mcin file.")
        btn_save_input.clicked.connect(self._save_input_file)
        btn_load_input = QPushButton("📂 Load inputs…")
        btn_load_input.setToolTip("Load a .mcin file as an editable template. You can edit before running MC.")
        btn_load_input.clicked.connect(self._load_input_file)
        btn_report = QPushButton("📄 Report…")
        btn_report.setToolTip("Preview and save a plain-text parameter summary.")
        btn_report.clicked.connect(self._preview_report)
        btn_close = QPushButton("Close"); btn_close.setMaximumWidth(80)
        btn_close.clicked.connect(self.close)
        bot.addWidget(btn_reset); bot.addWidget(btn_save_input)
        bot.addWidget(btn_load_input); bot.addWidget(btn_report)
        bot.addStretch(); bot.addWidget(btn_close)
        vbox.addLayout(bot)

        self._default_state = self._collect_state()
        self._default_state["mc_lw_slider"] = ("slider", 12)
        self._default_state["mc_op_slider"] = ("slider", 95)
        self._default_state["chk_mc_info"]  = ("check", True)
        if MCWindow._saved_state is not None:
            self._apply_state(MCWindow._saved_state)

        # v3.1.1 standardizes the sunward reference to the physical
        # emission-time direction even when an older .mcin/state stored the
        # legacy observation-time option.
        self.cmb_sunward_reference.setCurrentIndex(0)
        if self.chk_projected_sunward.isChecked():
            self.grp_projection_diag.setChecked(True)

        # Keep validation live while the user edits the physical inputs.
        for w in (self.sp_r_min, self.sp_r_max, self.sp_gamma_size,
                  self.sp_n, self.sp_max_age, self.sp_v0, self.sp_gamma,
                  self.sp_mexp, self.sp_sunward_cone, self.sp_sunward_expocos,
                  self.sp_lat_min, self.sp_lat_max, self.sp_lon_min,
                  self.sp_lon_max, self.sp_period):
            try:
                w.valueChanged.connect(self._update_run_button_state)
            except Exception:
                pass
        for rb in (self.rad_mode_isotropic, self.rad_mode_sunward,
                   self.rad_mode_active, self.rad_dmdt_steady,
                   self.rad_dmdt_cobs, self.rad_dmdt_manual):
            rb.toggled.connect(self._update_run_button_state)
        self.tbl_dmdt_manual.itemChanged.connect(self._update_run_button_state)
        self.tbl_size_over_time.itemChanged.connect(self._update_run_button_state)
        self.chk_size_table.toggled.connect(self._update_run_button_state)

        self.tabs.setCurrentIndex(0 if self._guided_first_use else 1)
        self._update_qt_suggestion_summary_labels()
        self._update_guided_navigation()
        self._update_run_button_state()

    def _update_fp_guidance_labels(self):
        """Refresh the read-only F-P guidance shown in the MC window."""
        dom = float(getattr(self, 'fp_dominant_age', 0.0) or 0.0)
        max_age = float(getattr(self, 'fp_max_age', 0.0) or 0.0)
        if hasattr(self, 'lbl_fp_dom_guidance'):
            self.lbl_fp_dom_guidance.setText(
                f"{dom:g} d" if dom > 0
                else "Not specified — Q(t) sets the lower guidance")
        if hasattr(self, 'lbl_fp_max_guidance'):
            self.lbl_fp_max_guidance.setText(
                f"{max_age:g} d  (largest synchrone age)"
                if max_age > 0 else "Unavailable")
        if hasattr(self, 'lbl_fp_maxage'):
            self.lbl_fp_maxage.setText(
                f"Maximum dust age: {max_age:g} d"
                if max_age > 0 else "Maximum dust age: unavailable")
        if hasattr(self, 'lbl_fp_dominant'):
            self.lbl_fp_dominant.setText(
                f"Dominant dust age: {dom:g} d"
                if dom > 0 else "Dominant dust age: Not specified")

    def update_fp_guidance(self, fp_sync_ages, fp_dominant_age=0.0):
        """Synchronize an already-open MC window with the main F-P panel."""
        try:
            ages = sorted({float(a) for a in (fp_sync_ages or [])
                           if float(a) > 0})
        except Exception:
            ages = []
        self.fp_sync_ages = ages
        self.fp_max_age = max(ages) if ages else 0.0
        try:
            self.fp_dominant_age = max(0.0, float(fp_dominant_age))
        except Exception:
            self.fp_dominant_age = 0.0
        if self.fp_max_age > 0 and self.fp_dominant_age > self.fp_max_age:
            self.fp_dominant_age = self.fp_max_age
        self._update_fp_guidance_labels()
        if hasattr(self, '_qt_fig'):
            QTimer.singleShot(0, self._refresh_qt_plot)

    def _update_qt_suggestion_summary_labels(self, assessment: str = ""):
        """Update the Q(t) recommendation/coverage card."""
        val = getattr(self, '_qt_suggested_window', None)
        rng = getattr(self, '_qt_suggested_window_range', None)
        note = getattr(self, '_qt_suggested_window_note', '') or ''
        cur = float(self.sp_max_age.value()) if hasattr(self, 'sp_max_age') else 0.0
        coverage = getattr(self, '_qt_coverage', None)

        if val is None:
            enabled = False
            if self.rad_dmdt_cobs.isChecked() and isinstance(coverage, dict):
                quality = coverage.get('quality', 'INSUFFICIENT')
                lookback = float(coverage.get('reliable_lookback_days', 0.0) or 0.0)
                frac = 100.0 * float(coverage.get('coverage_fraction', 0.0) or 0.0)
                msg = f"MC-window recommendation unavailable   |   COBS coverage: {quality}"
                msg += (f"\nReliable recent lookback: {lookback:.0f} d "
                        f"({frac:.0f}% of F-P maximum age)")
                provisional = getattr(self, '_qt_provisional_window', None)
                if provisional is not None:
                    msg += f"\nProvisional Q(t)-only value: {float(provisional):.0f} d — not applied automatically"
                reason = coverage.get('reason', '') or note
                if reason:
                    msg += f"\n{reason}"
                if quality in ('PARTIAL', 'INSUFFICIENT') and lookback > 0:
                    msg += ("\nRun MC may still be used with a manually selected window. "
                            "If the window exceeds the supported lookback, CTA will "
                            "ask before applying flat edge extrapolation.")
            else:
                msg = ("MC-window recommendation not calculated yet.\n"
                       "Enter/fetch Q(t) data, then refresh the preview.")
                if note:
                    msg += f"\n{note}"
        else:
            enabled = True
            msg = f"Recommended MC window: {float(val):.0f} d"
            if rng is not None:
                msg += f"   |   Suggested range: {float(rng[0]):.0f}–{float(rng[1]):.0f} d"
            msg += f"\nCurrent release window: {cur:.1f} d"
            qonly = getattr(self, '_qt_qonly_window', None)
            if qonly is not None:
                msg += f"\nQ(t)-only recommendation: {float(qonly):.0f} d"
            adjustment = getattr(self, '_qt_fp_adjustment', '') or ''
            if adjustment:
                msg += f"   |   F-P adjustment: {adjustment}"
            if isinstance(coverage, dict):
                msg += (f"\nCOBS coverage: {coverage.get('quality','')} — "
                        f"{float(coverage.get('reliable_lookback_days',0.0)):.0f} d reliable lookback")
            if assessment:
                msg += f"\n{assessment}"
            elif note:
                msg += f"\n{note}"

        if hasattr(self, 'lbl_qt_suggest_summary'):
            self.lbl_qt_suggest_summary.setText(msg)
        if hasattr(self, 'btn_apply_qt_suggested_window'):
            self.btn_apply_qt_suggested_window.setEnabled(enabled)
            if enabled:
                self.btn_apply_qt_suggested_window.setText(
                    f"Use recommended ({float(val):.0f} d)")
            else:
                self.btn_apply_qt_suggested_window.setText("Recommendation unavailable")

    def _usable_qt_sampling_arrays(self, qt_result=None):
        """Return finite, time-sorted Q(t) arrays usable by the MC sampler.

        Preference order:
          1. the recent continuous segment selected by the coverage analyser;
          2. the calibrated raw Q(t) curve itself.

        F-P age guidance and the coverage quality gate control only whether CTA
        may issue an automatic release-window recommendation.  They must not
        erase an otherwise valid COBS Q(t) curve or prevent a user-selected MC
        window from running.  The fallback is therefore deliberate: it keeps
        automatic Apply disabled, while allowing explicit simulation with the
        documented flat-edge extrapolation used by _sample_weighted_ages().
        """
        coverage = getattr(self, '_qt_coverage', None)
        if isinstance(coverage, dict):
            t = np.asarray(coverage.get('sampling_t_jd', []), dtype=float)
            q = np.asarray(coverage.get('sampling_q', []), dtype=float)
            if (t.size >= 2 and q.size == t.size and
                    np.all(np.isfinite(t)) and np.all(np.isfinite(q))):
                order = np.argsort(t)
                return t[order], q[order], 'coverage-segment'

        result = qt_result if isinstance(qt_result, dict) else getattr(self, '_last_qt_result', None)
        if not isinstance(result, dict):
            return np.array([], dtype=float), np.array([], dtype=float), 'none'

        t = np.asarray(result.get('t_jd', []), dtype=float)
        q = np.asarray(result.get('Q_kg_s', []), dtype=float)
        if t.size != q.size or t.size < 2:
            return np.array([], dtype=float), np.array([], dtype=float), 'none'

        mask = np.isfinite(t) & np.isfinite(q) & (q >= 0.0)
        # Release times are sampled only at or before the observation epoch.
        mask &= (t <= float(self.obs_jd) + 1e-9)
        t = t[mask]; q = q[mask]
        if t.size < 2:
            return np.array([], dtype=float), np.array([], dtype=float), 'none'

        order = np.argsort(t)
        t = t[order]; q = q[order]
        # np.interp expects increasing abscissae.  Keep one value per timestamp;
        # repeated dates are common in COBS and the calibrated curve is already
        # smoothed, so averaging duplicates is the least surprising choice.
        unique_t, inverse = np.unique(t, return_inverse=True)
        if unique_t.size != t.size:
            q_sum = np.zeros(unique_t.size, dtype=float)
            q_n = np.zeros(unique_t.size, dtype=float)
            np.add.at(q_sum, inverse, q)
            np.add.at(q_n, inverse, 1.0)
            q = q_sum / np.maximum(q_n, 1.0)
            t = unique_t
        if t.size < 2:
            return np.array([], dtype=float), np.array([], dtype=float), 'none'
        return t, q, 'raw-qt-fallback'

    def _qt_source_ready_for_run(self):
        """Return (ready, reason) for the selected dust-production source.

        Coverage quality controls automatic MC-window inference only.  A valid
        calibrated COBS Q(t) curve remains runnable even when F-P maximum age is
        unavailable or the continuous-coverage classifier is INSUFFICIENT.
        """
        try:
            self._qt_run_coverage_warning = ""
            self._qt_run_extrapolation_days = 0.0
            if self.rad_dmdt_cobs.isChecked():
                anchors = self._read_anchor_table()
                has_curve = bool(getattr(self, '_cobs_obs_list', None) or
                                 getattr(self, '_cobs_ready_cache', None) or
                                 getattr(self, '_last_qt_result', None))
                if not anchors:
                    return False, "COBS Q(t) selected: add at least one Afρ anchor."
                if not has_curve:
                    return False, "COBS Q(t) selected: fetch/generate COBS Q(t) first."

                t, q, source = self._usable_qt_sampling_arrays()
                if t.size < 2 or q.size != t.size:
                    coverage = getattr(self, '_qt_coverage', None)
                    reason = (coverage.get('reason', '') if isinstance(coverage, dict) else '')
                    return False, reason or "No usable calibrated COBS Q(t) curve is available."

                lookback = max(0.0, float(self.obs_jd) - float(np.min(t)))
                current_window = float(self.sp_max_age.value())
                extra = max(0.0, current_window - lookback)
                self._qt_run_extrapolation_days = extra
                coverage = getattr(self, '_qt_coverage', None)
                quality = coverage.get('quality', 'UNASSESSED') if isinstance(coverage, dict) else 'UNASSESSED'

                if extra > 1.0:
                    self._qt_run_coverage_warning = (
                        f"Run available with COBS edge extrapolation: the oldest "
                        f"{extra:.0f} d of the {current_window:.0f}-d release window "
                        f"lie outside the available {lookback:.0f}-d Q(t) support.")
                elif source == 'raw-qt-fallback':
                    self._qt_run_coverage_warning = (
                        "Run available from the calibrated COBS Q(t) curve. "
                        "F-P/coverage guidance is unavailable, so automatic window "
                        "Apply remains disabled.")
                elif quality != 'RELIABLE':
                    self._qt_run_coverage_warning = (
                        f"Run available with {quality.lower()} COBS coverage; "
                        "automatic window Apply remains disabled.")
            # Manual table validation still happens in _run(); do not over-block here.
            return True, ""
        except Exception as exc:
            return False, f"Q(t) input check failed: {exc}"

    def _validate_mc_inputs(self):
        """Return ``(errors, warnings)`` for the current physical inputs.

        Validation is based on actual values, not on whether a user happened
        to click every tab. Display choices never block a scientifically valid
        Monte Carlo run.
        """
        errors, warnings = [], []
        try:
            rmin = float(self.sp_r_min.value())
            rmax = float(self.sp_r_max.value())
            if not (rmin > 0.0 and rmax > 0.0):
                errors.append("Grain radii must be positive.")
            elif rmin >= rmax:
                errors.append("Set r_min smaller than r_max.")

            if float(self.sp_max_age.value()) <= 0.0:
                errors.append("Release window must be greater than zero.")
            if int(self.sp_n.value()) < 100:
                errors.append("Particle count must be at least 100.")

            if self.rad_mode_active.isChecked():
                if self.sp_lat_min.value() > self.sp_lat_max.value():
                    errors.append("Active-area latitude minimum exceeds maximum.")
                if self.sp_period.value() <= 0.0:
                    errors.append("Active-area rotation period must be positive.")

            if self.chk_size_table.isChecked():
                try:
                    self._read_size_over_time_table()
                except Exception as exc:
                    errors.append(f"Size-over-time table: {exc}")

            if self.rad_dmdt_manual.isChecked():
                try:
                    self._read_dmdt_manual_table()
                except Exception as exc:
                    errors.append(f"Manual dM/dt table: {exc}")

            qt_ready, reason = self._qt_source_ready_for_run()
            if not qt_ready:
                errors.append(reason)
            warning = getattr(self, '_qt_run_coverage_warning', '')
            if warning:
                warnings.append(warning)

            if (self.rad_mode_sunward.isChecked() and
                    self.grp_projection_diag.isChecked() and
                    self.chk_projected_sunward.isChecked()):
                warnings.append(
                    "Advanced apparent-Sun projection gate is enabled; report it as a diagnostic filter.")
        except Exception as exc:
            errors.append(f"Input validation failed: {exc}")
        return errors, warnings

    def _set_mc_status_style(self, level="neutral"):
        colors = {
            "ready": "#60d394",
            "warning": "#ffb020",
            "error": "#ff6b6b",
            "neutral": T['text_muted'],
        }
        if hasattr(self, 'lbl_status'):
            self.lbl_status.setStyleSheet(
                f"color:{colors.get(level, colors['neutral'])}; font-size:10px;")

    def _update_run_button_state(self, *_args):
        """Enable Run whenever the actual MC inputs are valid.

        The former v3.1 behaviour required visiting all three tabs. That was a
        discovery aid masquerading as validation and confused first-time users.
        v3.1.1 validates values directly; Display is optional.
        """
        if not hasattr(self, 'btn_run'):
            return
        try:
            if self._worker is not None and self._worker.isRunning():
                self.btn_run.setEnabled(False)
                self._set_mc_status_style("neutral")
                return
        except Exception:
            pass

        errors, warnings = self._validate_mc_inputs()
        if errors:
            self.btn_run.setEnabled(False)
            if hasattr(self, 'lbl_status'):
                shown = errors[:4]
                more = len(errors) - len(shown)
                text = "Complete the following:\n• " + "\n• ".join(shown)
                if more > 0:
                    text += f"\n• …and {more} more issue(s)"
                self.lbl_status.setText(text)
            self._set_mc_status_style("error")
        else:
            self.btn_run.setEnabled(True)
            if warnings:
                self.lbl_status.setText(
                    "Ready to run Monte Carlo.\nNote: " + " ".join(warnings))
                self._set_mc_status_style("warning")
            else:
                self.lbl_status.setText("Ready to run Monte Carlo.")
                self._set_mc_status_style("ready")
        self._update_guided_navigation()

    def _update_guided_navigation(self):
        if not hasattr(self, 'tabs') or not hasattr(self, 'btn_step_next'):
            return
        idx = int(self.tabs.currentIndex())
        labels = [
            "Step 1 of 3 — F-P Guide",
            "Step 2 of 3 — Simulation",
            "Step 3 of 3 — Display (optional)",
        ]
        self.lbl_step.setText(labels[max(0, min(2, idx))])
        self.btn_step_back.setEnabled(idx > 0)
        if idx == 0:
            self.btn_step_next.setText("Next: Simulation →")
            self.btn_step_next.setEnabled(True)
        elif idx == 1:
            self.btn_step_next.setText("Next: Display (optional) →")
            self.btn_step_next.setEnabled(True)
        else:
            self.btn_step_next.setText("Display is the final optional step")
            self.btn_step_next.setEnabled(False)
        self.btn_run.setVisible(idx == 1)

    def _navigate_mc_back(self):
        self.tabs.setCurrentIndex(max(0, self.tabs.currentIndex() - 1))

    def _navigate_mc_next(self):
        self.tabs.setCurrentIndex(min(2, self.tabs.currentIndex() + 1))

    def _skip_guided_setup(self):
        self._guided_first_use = False
        self._guide_settings.setValue("mc_guided_completed", True)
        self.guide_banner.hide()
        self.tabs.setCurrentIndex(1)
        self._update_guided_navigation()
        self._update_run_button_state()

    def _mark_guided_setup_complete(self):
        if getattr(self, '_guided_first_use', False):
            self._guided_first_use = False
            self._guide_settings.setValue("mc_guided_completed", True)
            if hasattr(self, 'guide_banner'):
                self.guide_banner.hide()

    def _run_mc_from_shortcut(self):
        # Ctrl+Enter always means "review/run Simulation", never an accidental
        # run from a hidden tab with unclear context.
        if self.tabs.currentIndex() != 1:
            self.tabs.setCurrentIndex(1)
        self._update_run_button_state()
        if self.btn_run.isEnabled():
            self._run()

    def _on_mc_tab_changed(self, index: int):
        self._visited_tabs.add(int(index))
        if index == 1 and (self.rad_dmdt_cobs.isChecked() or
                           self.rad_dmdt_manual.isChecked()):
            QTimer.singleShot(0, self._refresh_qt_plot)
        self._update_guided_navigation()
        self._update_run_button_state()

    def _refresh_qt_plot(self):
        """Safe Q(t) preview refresh wrapper.

        Qt signal handlers should not allow uncaught exceptions to escape.
        In some PyQt/Python builds an exception raised while switching tabs can
        terminate the dialog or the whole application.  Keep the Monte Carlo
        window alive and show the error inside the Q(t) panel instead.
        """
        try:
            return self._refresh_qt_plot_impl()
        except Exception as exc:
            import traceback
            tb = traceback.format_exc(limit=6)
            msg = f"Q(t) preview error: {exc}"

            # Reset the suggestion state so the Apply button cannot use stale
            # values from an earlier successful refresh.
            self._qt_suggested_window = None
            self._qt_suggested_window_range = None
            self._qt_suggested_window_note = ""
            self._qt_qonly_window = None
            self._qt_fp_adjustment = ""
            self._qt_coverage = None
            self._qt_recommendation_quality = ""
            self._qt_provisional_window = None
            if hasattr(self, 'btn_apply_qt_suggested_window'):
                self.btn_apply_qt_suggested_window.setEnabled(False)
            if hasattr(self, 'btn_step2_apply_qt_suggested_window'):
                self.btn_step2_apply_qt_suggested_window.setEnabled(False)
            if hasattr(self, '_update_qt_suggestion_summary_labels'):
                self._update_qt_suggestion_summary_labels()

            # Draw an in-panel diagnostic instead of letting the window close.
            try:
                fig = self._qt_fig
                fig.clear()
                ax = fig.add_subplot(111)
                ax.set_facecolor('#1a1a1a')
                ax.text(0.5, 0.56, 'Q(t) preview could not be drawn',
                        transform=ax.transAxes, ha='center', va='center',
                        color='#ffb020', fontsize=10, fontweight='bold')
                ax.text(0.5, 0.43, str(exc),
                        transform=ax.transAxes, ha='center', va='center',
                        color='#cccccc', fontsize=8, wrap=True)
                ax.set_xticks([]); ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_color('#444')
                self._qt_canvas.draw()
            except Exception:
                pass

            if hasattr(self, '_qt_info_lbl'):
                self._qt_info_lbl.setText(msg + "\n\n" + tb)
            if hasattr(self, 'lbl_status'):
                self.lbl_status.setText(msg)
            return None

    def _refresh_qt_plot_impl(self):
        """Re-draw the Q(t) preview using the SAME Q(t) logic as the MC run.

        Important: this plot is deliberately tied to estimate_qt_from_lightcurve()
        / manual dM/dt parsing rather than an independent magnitude-to-Afρ
        shortcut.  The preview should show the release-time weights that the
        Monte Carlo sampler will actually use.
        """
        import numpy as np
        fig = self._qt_fig
        fig.clear()
        ax = fig.add_subplot(111)
        T_colors = dict(text='#e0e0e0', text_muted='#888', grid='#333',
                        cobs='#4a8adf', anchor='#ffe030', window='#ff3399',
                        qw='#60e0a0', warn='#ffb020', fp='#c084ff',
                        fp_band='#8b5cf6')

        window_d = float(self.sp_max_age.value())
        self._qt_suggested_window = None
        self._qt_suggested_window_range = None
        self._qt_suggested_window_note = ""
        self._qt_qonly_window = None
        self._qt_fp_adjustment = ""
        self._qt_coverage = None
        self._qt_recommendation_quality = ""
        self._qt_provisional_window = None
        if hasattr(self, 'btn_apply_qt_suggested_window'):
            self.btn_apply_qt_suggested_window.setEnabled(False)
        if hasattr(self, 'btn_step2_apply_qt_suggested_window'):
            self.btn_step2_apply_qt_suggested_window.setEnabled(False)
        self._update_qt_suggestion_summary_labels()

        # Automatic F-P morphology guidance: one optional dominant age
        # plus a maximum age derived from the largest listed synchrone.
        fp_dom = float(getattr(self, 'fp_dominant_age', 0.0) or 0.0)
        fp_outer = float(getattr(self, 'fp_max_age', 0.0) or 0.0)
        use_fp_age = fp_outer > 0.0
        if use_fp_age and fp_dom > fp_outer:
            fp_dom = fp_outer

        info_parts = []
        has_data = False

        def _weighted_age_stats(t_jd, q_w, obs_jd, t_min=None, t_max=None):
            """Stats for a Q(t) curve, integrated on a uniform time grid.

            The COBS sampling cadence is highly non-uniform.  Using raw point
            sums would overweight nights with many reports.  This function
            mirrors the MC sampler's behaviour: interpolate Q(t) onto a uniform
            grid, then compute age statistics from that continuous curve.
            """
            t_jd = np.asarray(t_jd, dtype=float)
            q_w  = np.asarray(q_w, dtype=float)
            m = np.isfinite(t_jd) & np.isfinite(q_w) & (q_w > 0)
            if m.sum() < 2:
                return None
            t = t_jd[m]
            q = q_w[m]
            order = np.argsort(t)
            t = t[order]
            q = q[order]
            if t_min is None:
                t_min = max(float(t[0]), obs_jd - window_d)
            if t_max is None:
                t_max = min(float(t[-1]), obs_jd)
            t_min = float(t_min); t_max = float(t_max)
            if not (t_min < t_max):
                return None
            fine_t = np.linspace(t_min, t_max, 2000)
            fine_q = np.clip(np.interp(fine_t, t, q), 0.0, None)
            if fine_q.sum() <= 0:
                return None
            ages = obs_jd - fine_t
            # Sort by age ascending: 0 d/recent → older dust.
            s = np.argsort(ages)
            ages_s = ages[s]
            q_s = fine_q[s]
            cum = np.cumsum(q_s)
            cum = cum / cum[-1]
            mean_age = float(np.sum(ages_s * q_s) / np.sum(q_s))
            def qtile(frac):
                return float(np.interp(frac, cum, ages_s))
            return dict(mean=mean_age,
                        w50=qtile(0.50), w80=qtile(0.80),
                        w90=qtile(0.90), w95=qtile(0.95),
                        total=float(np.trapezoid(fine_q, fine_t)),
                        t_min=t_min, t_max=t_max)

        def _plot_q_curve(t_jd, q_w, label, source_kind='Q(t)'):
            nonlocal has_data
            t_jd = np.asarray(t_jd, dtype=float)
            q_w  = np.asarray(q_w, dtype=float)
            m = np.isfinite(t_jd) & np.isfinite(q_w) & (q_w > 0)
            if m.sum() < 2:
                return None
            t = t_jd[m]
            q = q_w[m]
            order = np.argsort(t)
            t = t[order]
            q = q[order]
            dt = t - self.obs_jd
            q_norm = q / np.nanmax(q)
            q_norm = np.clip(q_norm, 1e-6, None)
            ax.plot(dt, q_norm, color=T_colors['cobs'], lw=1.7,
                    label=label, zorder=3)
            ax.scatter(dt, q_norm, s=8, alpha=0.35,
                       color=T_colors['cobs'], zorder=2)
            has_data = True
            return dict(t=t, q=q, dt=dt, q_norm=q_norm, source=source_kind)

        # ── Select Q(t) source: COBS, manual table, or steady ─────────────
        q_plot = None
        if self.rad_dmdt_cobs.isChecked():
            anchors = self._read_anchor_table()
            if not self._cobs_obs_list:
                info_parts.append(
                    "COBS Q(t) mode selected, but no COBS light curve is loaded. "
                    "Fetch COBS data first.")
            elif not anchors:
                info_parts.append(
                    "COBS Q(t) mode selected, but no Afρ anchor is entered. "
                    "At least one anchor is required to calibrate the proxy curve.")
            else:
                try:
                    import comet_tail_analyzer as cta
                    qt_result = cta.estimate_qt_from_lightcurve(
                        self.comet_el, self._cobs_obs_list, anchors,
                        p_v=float(self.sp_pv.value()),
                        smooth_window_days=float(self.sp_qt_smooth.value()))
                    self._last_qt_result = qt_result
                    q_plot = _plot_q_curve(qt_result['t_jd'], qt_result['Q_kg_s'],
                                           'All COBS-derived Q(t) points',
                                           source_kind='COBS')

                    # Coverage/inference safety gate.  The recommendation uses
                    # only the most recent continuous segment inside the F-P
                    # maximum-age interval; older disconnected apparitions are
                    # plotted for context but never integrated into the result.
                    coverage = cta.assess_qt_coverage(
                        qt_result['t_jd'], qt_result['Q_kg_s'], self.obs_jd,
                        max_age_days=fp_outer,
                        dominant_age_days=fp_dom,
                        smooth_window_days=float(self.sp_qt_smooth.value()))
                    self._qt_coverage = coverage
                    self._qt_recommendation_quality = coverage.get('quality', '')

                    if q_plot is not None:
                        # Mark the accepted recent segment with a heavier line.
                        seg_t = np.asarray(coverage.get('segment_t_jd', []), dtype=float)
                        seg_q = np.asarray(coverage.get('segment_q', []), dtype=float)
                        if seg_t.size >= 2 and seg_q.size == seg_t.size:
                            norm = max(float(np.nanmax(q_plot['q'])), 1e-30)
                            seg_y = np.clip(seg_q / norm, 1e-6, None)
                            quality_color = (T_colors['qw'] if coverage.get('quality') == 'RELIABLE'
                                             else T_colors['warn'])
                            ax.plot(seg_t - self.obs_jd, seg_y,
                                    color=quality_color, lw=3.0, alpha=0.95,
                                    label=(f"Recent continuous COBS segment "
                                           f"({coverage.get('quality','')})"), zorder=5)

                        # Mark anchor dates.  Anchor values are Afρ labels, not
                        # independent Q(t) y-values.
                        for jd_a, af_a in zip(qt_result['anchor_jd'],
                                              qt_result['anchor_afrho']):
                            y_a = np.interp(jd_a, q_plot['t'], q_plot['q_norm'])
                            ax.scatter([jd_a - self.obs_jd], [y_a],
                                       s=80, marker='*', color=T_colors['anchor'],
                                       zorder=6, label='Afρ anchor' if jd_a == qt_result['anchor_jd'][0] else None)
                            ax.annotate(f'{af_a:.0f} cm',
                                        (jd_a - self.obs_jd, y_a),
                                        textcoords='offset points', xytext=(5, 5),
                                        fontsize=7, color=T_colors['anchor'])

                    info_parts.append(
                        "Q(t) source: COBS light curve scaled by your Afρ anchor(s).")
                    info_parts.append(
                        f"COBS coverage: {coverage.get('quality','INSUFFICIENT')} | "
                        f"recent lookback={coverage.get('reliable_lookback_days',0.0):.0f} d | "
                        f"coverage={100.0*coverage.get('coverage_fraction',0.0):.0f}% | "
                        f"unique dates={coverage.get('unique_points',0)} | "
                        f"largest accepted gap={coverage.get('largest_gap_days',0.0):.1f} d.")
                    info_parts.append(coverage.get('reason', ''))
                except Exception as exc:
                    info_parts.append(f"Q(t) estimation error: {exc}")

        elif self.rad_dmdt_manual.isChecked():
            try:
                import comet_tail_analyzer as cta
                days_arr, dmdt_arr = self._read_dmdt_manual_table()
                per_jd = self.comet_el.get('T_jd') or cta.date_to_jd(self.comet_el['T'])
                t_jd = per_jd + np.asarray(days_arr, dtype=float)
                q_plot = _plot_q_curve(t_jd, np.asarray(dmdt_arr, dtype=float),
                                       'Manual dM/dt table', source_kind='manual')
                info_parts.append("Q(t) source: manual dM/dt table.")
            except Exception as exc:
                info_parts.append(f"Manual dM/dt table error: {exc}")

        else:
            # Steady mode: show the distribution actually implied by the MC
            # age sampler, but do not pretend there is a photometric Q(t)
            # recommendation.
            t_jd = np.array([self.obs_jd - window_d, self.obs_jd], dtype=float)
            q_w = np.array([1.0, 1.0], dtype=float)
            q_plot = _plot_q_curve(t_jd, q_w, 'Steady Q(t) over MC window',
                                   source_kind='steady')
            info_parts.append(
                "Q(t) source: steady/uniform. No light-curve-based release-window "
                "suggestion is available in this mode.")

        # ── MC release window shading ────────────────────────────────────
        ax.axvspan(-window_d, 0, alpha=0.12, color=T_colors['window'],
                   label=f'Current MC window ({window_d:.0f} d)')
        ax.axvline(0, color='#888', lw=0.8, ls='--', alpha=0.6)
        ax.text(0.985, 0.95, 'obs date', transform=ax.transAxes,
                fontsize=7, color='#888', ha='right', va='top')

        # ── Automatic F-P morphology-age overlay ────────────────────────
        if use_fp_age:
            if fp_dom > 0.0:
                ax.axvline(-fp_dom, color=T_colors['fp_band'], lw=1.3, ls=':',
                           alpha=0.95,
                           label=f'F-P dominant dust age ({fp_dom:.0f} d)')
            ax.axvline(-fp_outer, color=T_colors['fp'], lw=1.2, ls='--',
                       alpha=0.95,
                       label=f'F-P maximum dust age ({fp_outer:.0f} d)')
            dom_desc = (f"dominant dust age ≈ {fp_dom:.0f} d; "
                        if fp_dom > 0 else "")
            info_parts.append(
                f"F-P morphology guidance: {dom_desc}maximum dust age ≈ "
                f"{fp_outer:.0f} d. The maximum age is derived from the largest "
                f"synchrone and is used as the upper MC-window limit.")

        # ── Q-weighted dust-age statistics + suggested window ────────────
        if q_plot is not None and q_plot.get('source') != 'steady':
            source_kind = q_plot.get('source')
            recommendation_allowed = True

            if source_kind == 'COBS':
                coverage = getattr(self, '_qt_coverage', None)
                if isinstance(coverage, dict):
                    t = np.asarray(coverage.get('segment_t_jd', []), dtype=float)
                    q = np.asarray(coverage.get('segment_q', []), dtype=float)
                    recommendation_allowed = bool(
                        coverage.get('recommendation_allowed', False))
                else:
                    t = np.array([], dtype=float)
                    q = np.array([], dtype=float)
                    recommendation_allowed = False
            else:
                # A manual table is explicit user input, not a heterogeneous
                # archive.  It still obeys the F-P maximum-age interval.
                t = np.asarray(q_plot['t'], dtype=float)
                q = np.asarray(q_plot['q'], dtype=float)
                if use_fp_age:
                    m = ((t >= self.obs_jd - fp_outer) & (t <= self.obs_jd))
                    t = t[m]; q = q[m]

            stats_all = None
            stats_win = None
            if t.size >= 2 and q.size == t.size:
                # No extrapolation for inference: integrate only between the
                # first and last points of the accepted segment/table.
                t_min = float(np.nanmin(t))
                t_max = float(np.nanmax(t))
                stats_all = _weighted_age_stats(
                    t, q, self.obs_jd, t_min=t_min, t_max=t_max)
                win_min = max(t_min, self.obs_jd - window_d)
                win_max = min(t_max, self.obs_jd)
                if win_min < win_max:
                    stats_win = _weighted_age_stats(
                        t, q, self.obs_jd, t_min=win_min, t_max=win_max)

            if stats_all is not None:
                mean_all = stats_all['mean']
                stat_color = (T_colors['qw'] if recommendation_allowed
                              else T_colors['warn'])
                stat_prefix = '' if recommendation_allowed else 'provisional '
                ax.axvline(-mean_all, color=stat_color, lw=1.5,
                           ls=':', alpha=0.95,
                           label=f'{stat_prefix}Q-weighted mean ({mean_all:.0f} d)')
                ax.axvline(-stats_all['w80'], color=T_colors['warn'], lw=1.0,
                           ls='--', alpha=0.8,
                           label=f'{stat_prefix}80% cumulative Q ({stats_all["w80"]:.0f} d)')

                # Q(t)-only candidate from the accepted support.  F-P then
                # supplies a dominant-age floor and maximum-age cap.
                q_lo  = max(1.3 * mean_all, stats_all['w50'])
                q_mid = max(1.4 * mean_all, stats_all['w80'])
                q_hi  = max(1.6 * mean_all, stats_all['w90'])
                floor = fp_dom if fp_dom > 0.0 else 0.0
                cap = fp_outer if use_fp_age else float('inf')
                sug_lo  = min(cap, max(q_lo, floor))
                sug_mid = min(cap, max(q_mid, floor))
                sug_hi  = min(cap, max(q_hi, floor))
                self._qt_qonly_window = float(q_mid)

                if floor > 0 and q_mid < floor:
                    self._qt_fp_adjustment = f"raised to dominant age {floor:.0f} d"
                elif use_fp_age and q_mid > cap:
                    self._qt_fp_adjustment = f"capped at maximum age {cap:.0f} d"
                else:
                    self._qt_fp_adjustment = "none required"

                suggestion_rule = (
                    "Q(t): max(1.3–1.6 × weighted mean, W50/W80/W90), "
                    "computed only from supported dates; then apply the "
                    "dominant-age floor and maximum-age cap")

                if recommendation_allowed:
                    span_lo = min(sug_lo, sug_hi)
                    span_hi = max(sug_lo, sug_hi)
                    if span_hi > span_lo:
                        ax.axvspan(-span_hi, -span_lo,
                                   color=T_colors['qw'], alpha=0.08,
                                   label=(f'F-P+Q suggested {span_lo:.0f}–'
                                          f'{span_hi:.0f} d'))
                    ax.axvline(-sug_mid, color=T_colors['qw'], lw=1.2,
                               ls='-', alpha=0.85,
                               label=f'recommended MC window ({sug_mid:.0f} d)')
                    self._qt_suggested_window = float(sug_mid)
                    self._qt_suggested_window_range = (
                        float(span_lo), float(span_hi))
                    self._qt_suggested_window_note = (
                        f"Suggested effective MC window ≈ {sug_mid:.0f} d. "
                        f"Rule: {suggestion_rule}.")
                    self._qt_recommendation_quality = "RELIABLE"
                else:
                    # A partial/insufficient segment can still be informative,
                    # but must not be presented as an actionable result.
                    self._qt_provisional_window = float(sug_mid)
                    self._qt_suggested_window = None
                    self._qt_suggested_window_range = None
                    self._qt_suggested_window_note = (
                        "Q(t) statistics are provisional because continuous "
                        "COBS coverage is incomplete. Automatic Apply is disabled.")

                info_parts.append(
                    f"Q-weighted mean age: {mean_all:.1f} d | "
                    f"Q-only candidate: {q_mid:.0f} d | "
                    f"F-P-adjusted candidate: {sug_mid:.0f} d.")
                info_parts.append(f"Suggestion rule: {suggestion_rule}.")
                info_parts.append(
                    f"Cumulative-Q age scales: W50={stats_all['w50']:.0f} d, "
                    f"W80={stats_all['w80']:.0f} d, "
                    f"W90={stats_all['w90']:.0f} d.")
                if stats_win is not None:
                    info_parts.append(
                        f"Within supported dates inside the current {window_d:.0f} d "
                        f"window: mean={stats_win['mean']:.1f} d, "
                        f"W80={stats_win['w80']:.0f} d.")

                assessment = ""
                if recommendation_allowed:
                    if use_fp_age and window_d > fp_outer:
                        assessment = "Current value exceeds the F-P maximum dust age."
                    elif abs(window_d - sug_mid) <= max(3.0, 0.15 * sug_mid):
                        assessment = "Current value is close to the recommended window."
                    elif window_d > sug_mid:
                        assessment = "Current value is longer than the final recommendation."
                    else:
                        assessment = (
                            "Current value is shorter than the final recommendation; "
                            "older visible dust may be underrepresented.")
                else:
                    assessment = (
                        "No automatic recommendation: the displayed statistics "
                        "do not have sufficient continuous COBS support.")
                self._update_qt_suggestion_summary_labels(assessment)
                self._update_run_button_state()
            elif source_kind == 'COBS':
                self._qt_suggested_window_note = (
                    "No valid recent continuous Q(t) segment is available for inference.")
                self._update_qt_suggestion_summary_labels()
                self._update_run_button_state()

        # ── Axis styling ─────────────────────────────────────────────────
        ax.set_facecolor('#1a1a1a')
        for spine in ax.spines.values():
            spine.set_color('#444')
        ax.tick_params(colors='#aaa', labelsize=8)
        ax.set_xlabel('Days from observation (negative = before)',
                      color='#aaa', fontsize=8)
        ax.set_ylabel('Relative Q(t)',
                      color='#aaa', fontsize=8)
        ax.set_title('Dust-release weighting Q(t)', color='#ccc', fontsize=9)
        ax.grid(True, color=T_colors['grid'], lw=0.5, alpha=0.5)
        if has_data:
            ax.set_yscale('log')
            ax.set_ylim(bottom=1e-4, top=1.8)
            # Keep the diagnostic focused on the physically relevant F-P/MC
            # interval.  Historical COBS apparitions outside this span may be
            # present in memory, but should not compress the recent curve into
            # a few pixels or imply that they contribute to the recommendation.
            axis_lookback = max(window_d, fp_outer if use_fp_age else 0.0, 10.0)
            ax.set_xlim(-1.05 * axis_lookback, max(3.0, 0.03 * axis_lookback))
        handles, labels = ax.get_legend_handles_labels()
        # Deduplicate repeated labels, especially anchor markers.
        if handles:
            seen = set(); h2 = []; l2 = []
            for h, l in zip(handles, labels):
                if not l or l in seen:
                    continue
                seen.add(l); h2.append(h); l2.append(l)
            ax.legend(h2, l2, fontsize=7, facecolor='#222', edgecolor='#444',
                      labelcolor='#ccc', loc='upper left')

        if not has_data:
            ax.text(0.5, 0.5,
                    'No Q(t) curve available.\nFetch COBS + enter Afρ anchor,\nor use a manual dM/dt table.',
                    transform=ax.transAxes, ha='center', va='center',
                    color='#666', fontsize=9)

        fig.tight_layout(pad=0.8)
        self._qt_canvas.draw()
        self._qt_info_lbl.setText('\n'.join(info_parts) if info_parts else
                                   "No Q(t) statistics available.")

    def _apply_qt_suggested_window(self):
        """Apply the latest Q(t)/F-P suggested effective MC release window."""
        val = getattr(self, '_qt_suggested_window', None)
        if val is None:
            QMessageBox.information(
                self, "Apply suggested window",
                "No suggestion is available yet. Refresh the Q(t) plot first.")
            return
        try:
            val = float(val)
            self.sp_max_age.setValue(max(self.sp_max_age.minimum(),
                                         min(self.sp_max_age.maximum(), val)))
            note = getattr(self, '_qt_suggested_window_note', '')
            self.lbl_status.setText(note or f"Applied suggested MC window: {val:.0f} d")
            self._refresh_qt_plot()
            self._update_qt_suggestion_summary_labels()
            # Applying a suggestion changes only the Release-window value.
            # Keep the user on the current tab; manual workflows still require
            # the user to review every tab explicitly before Run is enabled.
            # Also make sure Run never becomes the dialog's default button.
            if hasattr(self, 'btn_run'):
                self.btn_run.setAutoDefault(False)
                self.btn_run.setDefault(False)
            self._update_run_button_state()
        except Exception as exc:
            QMessageBox.warning(self, "Apply suggested window", f"Could not apply:\n{exc}")

    def _reset_to_defaults(self):
        self._apply_state(self._default_state)
        self.lbl_status.setText("Fields reset to defaults.")

    def closeEvent(self, event):
        """Remember field values for next time this window is opened —
        a plain class attribute (not instance), since each open/close
        creates a brand-new MCWindow (WA_DeleteOnClose) that would
        otherwise lose everything typed in if the window is closed by
        accident before clicking Run."""
        MCWindow._saved_state = self._collect_state()
        super().closeEvent(event)

    def _collect_state(self) -> dict:
        """Generic snapshot of every input widget on this window, keyed
        by its self.xxx attribute name. Generic (rather than a hand-
        written field list) so newly added inputs are picked up for
        free and nothing gets silently missed."""
        state = {}
        for name, w in vars(self).items():
            if isinstance(w, QDoubleSpinBox):
                state[name] = ("dspin", w.value())
            elif isinstance(w, QSpinBox):
                state[name] = ("spin", w.value())
            elif isinstance(w, QCheckBox):
                state[name] = ("check", w.isChecked())
            elif isinstance(w, QRadioButton):
                state[name] = ("radio", w.isChecked())
            elif isinstance(w, QComboBox):
                state[name] = ("combo", w.currentIndex())
            elif isinstance(w, QLineEdit):
                state[name] = ("line", w.text())
            elif isinstance(w, QSlider):
                state[name] = ("slider", w.value())
            elif isinstance(w, QTableWidget):
                rows = []
                for r in range(w.rowCount()):
                    rows.append([
                        (w.item(r, c).text() if w.item(r, c) else "")
                        for c in range(w.columnCount())
                    ])
                state[name] = ("table", rows)
        return state

    def _apply_state(self, state: dict):
        """Inverse of _collect_state(). Silently skips any widget no
        longer present (e.g. after a future edit) instead of raising."""
        for name, (kind, val) in state.items():
            w = getattr(self, name, None)
            if w is None:
                continue
            if kind in ("dspin", "spin"):
                w.setValue(val)
            elif kind in ("check", "radio"):
                w.setChecked(val)
            elif kind == "combo":
                w.setCurrentIndex(val)
            elif kind == "line":
                w.setText(val)
            elif kind == "slider":
                w.setValue(val)
            elif kind == "table":
                # Keep Afρ anchor entry user-friendly after loading older or
                # sparse templates: always show at least three editable rows.
                min_rows = 3 if name == "tbl_anchors" else len(val)
                w.setRowCount(max(min_rows, len(val)))
                for r, row in enumerate(val):
                    for c, text in enumerate(row):
                        w.setItem(r, c, QTableWidgetItem(text))
                if name == "tbl_anchors":
                    for r in range(len(val), w.rowCount()):
                        for c in range(w.columnCount()):
                            if w.item(r, c) is None:
                                w.setItem(r, c, QTableWidgetItem(""))

        # Backward compatibility: older .mcin files may store the legacy
        # observation-time Sun reference. v3.1.1 always restores the physical
        # emission-time reference while preserving the rest of the file.
        if hasattr(self, "cmb_sunward_reference"):
            self.cmb_sunward_reference.setCurrentIndex(0)
        if hasattr(self, "grp_projection_diag") and hasattr(self, "chk_projected_sunward"):
            self.grp_projection_diag.setChecked(self.chk_projected_sunward.isChecked())

    def _save_input_file(self):
        """Save the current visible MC inputs as a NEW .mcin JSON file.

        This deliberately supports an edit-before-run workflow:

            Load .mcin  →  edit any GUI field  →  Save inputs as…  →  Run MC

        Loading a file never locks the controls and saving never depends on
        whether MC has already been run.  The saved file contains the current
        visible state at the moment the user clicks Save.
        """
        import json, os, re
        from datetime import datetime, timezone

        def _slug(text, fallback="MC_input"):
            text = (text or fallback).strip()
            text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
            text = text.strip("_")
            return text or fallback

        comet_name = self.comet_el.get("name", "") or "comet"
        obs_date = cta.jd_to_str(self.obs_jd)[:10]

        loaded_profile = getattr(self, "_loaded_input_profile", {}) or {}
        loaded_name = loaded_profile.get("name", "") if isinstance(loaded_profile, dict) else ""
        loaded_path = getattr(self, "_loaded_input_path", None)
        loaded_base = os.path.basename(loaded_path) if loaded_path else ""

        if loaded_name:
            default_name = _slug(loaded_name + "_edited") + ".mcin"
            profile_name = "Edited copy — " + loaded_name
        elif loaded_base:
            root, _ext = os.path.splitext(loaded_base)
            default_name = _slug(root + "_edited") + ".mcin"
            profile_name = f"Edited copy — {root}"
        else:
            default_name = f"MC_{comet_name.replace(' ','_').replace('/','')}_{obs_date}.mcin"
            profile_name = f"User saved MC input — {comet_name} {obs_date}"

        path, _ = QFileDialog.getSaveFileName(
            self, "Save edited MC input as", default_name,
            "CTA MC input (*.mcin);;JSON file (*.json);;All files (*)")
        if not path:
            return
        if not (path.lower().endswith(".mcin") or path.lower().endswith(".json")):
            path += ".mcin"

        notes = [
            "All model choices are stored as visible GUI inputs.",
            "This file does not activate any hidden object-specific model.",
            "This file can be loaded later, edited again, and saved as another input file."
        ]
        if loaded_base:
            notes.append(f"Created from editable input template: {loaded_base}")

        profile = {
            "name": profile_name,
            "target": comet_name,
            "observation_date_utc": obs_date,
            "mode": "editable CTA MC input",
            "intended_use": "Monte Carlo morphology / surface-brightness configuration",
            "created_by": "Comet Tail Analyzer",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "notes": notes,
        }
        if loaded_base:
            profile["parent_input_file"] = loaded_base
        if loaded_profile and isinstance(loaded_profile, dict):
            profile["parent_profile_name"] = loaded_profile.get("name", "")

        payload = {
            "schema": "CTA_MC_INPUT",
            "profile": profile,
            "comet": comet_name,
            "obs_jd": self.obs_jd,
            "fp_guidance": {
                "synchrone_ages_days": list(getattr(self, "fp_sync_ages", [])),
                "dominant_dust_age_days": float(getattr(self, "fp_dominant_age", 0.0) or 0.0),
                "maximum_dust_age_days": float(getattr(self, "fp_max_age", 0.0) or 0.0),
            },
            "state": self._collect_state(),
        }
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, ensure_ascii=False)
            # Treat the new file as the current editable template from now on.
            self._loaded_input_path = path
            self._loaded_input_profile = profile
            self._loaded_input_comet = comet_name
            self.lbl_status.setText(
                f"Edited inputs saved → {path}. You can continue editing or run MC.")
        except Exception as exc:
            QMessageBox.warning(self, "Save inputs", f"Could not save:\n{exc}")

    def _load_input_file(self):
        """Load a CTA .mcin input file.

        Supports CTA_MC_INPUT schema files only.

        Loading a file restores visible inputs only. It does not run MC, does
        not lock the controls, and does not switch to a hidden model. The user
        may edit any field after loading and save the edited state as a new
        .mcin file. Files without schema='CTA_MC_INPUT' are rejected to avoid
        ambiguous behaviour.
        """
        import json
        path, _ = QFileDialog.getOpenFileName(
            self, "Load MC input file", "",
            "CTA MC input (*.mcin *.json);;All files (*)")
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as f:
                payload = json.load(f)
        except Exception as exc:
            QMessageBox.warning(self, "Load inputs", f"Could not read file:\n{exc}")
            return

        schema = payload.get("schema", "")
        if schema != "CTA_MC_INPUT":
            QMessageBox.warning(
                self, "Load inputs",
                "Unsupported CTA input file.\n\n"
                "This build accepts only files with schema='CTA_MC_INPUT'.")
            return

        profile = payload.get("profile", {}) if isinstance(payload.get("profile", {}), dict) else {}
        saved_comet = payload.get("comet", "") or profile.get("target", "")
        this_comet  = self.comet_el.get("name", "")
        profile_name = profile.get("name", "")
        intended_use = profile.get("intended_use", "")
        notes = profile.get("notes", [])

        if saved_comet and this_comet and saved_comet != this_comet:
            from PyQt6.QtWidgets import QMessageBox as _MB
            msg = (
                f"This input file was saved for <b>{saved_comet}</b>, but the\n"
                f"current comet is <b>{this_comet}</b>.\n\n"
            )
            if profile_name:
                msg += f"Profile: {profile_name}\n"
            if intended_use:
                msg += f"Use: {intended_use}\n"
            msg += "\nLoad anyway?"
            ans = _MB.question(
                self, "Load inputs — comet mismatch", msg,
                _MB.StandardButton.Yes | _MB.StandardButton.Cancel)
            if ans != _MB.StandardButton.Yes:
                return

        fp_guidance = payload.get("fp_guidance", {})
        if isinstance(fp_guidance, dict):
            ages_saved = fp_guidance.get("synchrone_ages_days", [])
            dom_saved = fp_guidance.get("dominant_dust_age_days", 0.0)
            try:
                self.update_fp_guidance(ages_saved, dom_saved)
            except Exception:
                pass

        state = payload.get("state", {})
        if not isinstance(state, dict):
            QMessageBox.warning(self, "Load inputs", "This file has no valid 'state' object.")
            return
        self._apply_state(state)
        self._smooth_au_initialized = True
        self._loaded_input_path = path
        self._loaded_input_profile = profile
        self._loaded_input_comet = saved_comet or this_comet

        suffix = ""
        if profile_name:
            suffix += f"  [{profile_name}]"
        elif saved_comet and saved_comet != this_comet:
            suffix += f"  [saved for {saved_comet}]"
        if notes:
            suffix += "  [notes in file]"
        self.lbl_status.setText(
            f"Loaded editable inputs from {path}" + suffix +
            " — edit any field, then Run MC or Save inputs as a new .mcin file.")
        # A loaded .mcin file is already a structured setup, so move directly
        # to Simulation. It remains fully editable and is never run automatically.
        self._visited_tabs = {0, 1, 2}
        self._guided_first_use = False
        self.guide_banner.hide()
        self.tabs.setCurrentIndex(1)
        self._update_guided_navigation()
        self._update_run_button_state()

        # Loading a file should never trigger MC generation.  It only restores
        # visible values.  If the restored template uses COBS+Afρ Q(t), fetch
        # COBS automatically so the user can edit parameters and run MC without
        # first visiting the Q(t) tab manually.
        if self.rad_dmdt_cobs.isChecked():
            self._auto_prepare_qt_after_input_load()
        else:
            try:
                self._refresh_qt_plot()
            except Exception:
                pass

    def _auto_prepare_qt_after_input_load(self):
        """After loading an editable .mcin template, automatically prepare
        COBS-based Q(t) if that source is selected.

        This deliberately does NOT run MC and does NOT overwrite the loaded
        release window.  It only populates _cobs_obs_list, updates the COBS
        status label, and redraws the Q(t) preview/suggestion so the loaded
        template is ready for review and editing.
        """
        if self.main_window is None:
            try:
                self._refresh_qt_plot()
            except Exception:
                pass
            return

        anchors = []
        try:
            anchors = self._read_anchor_table()
        except Exception:
            anchors = []

        self._cobs_obs_list = None
        self.btn_fetch_cobs.setEnabled(False)
        if hasattr(self, 'lbl_cobs_status'):
            self.lbl_cobs_status.setText("Auto-fetching COBS for loaded input…")
        self.lbl_status.setText(
            "Loaded editable inputs. Auto-fetching COBS light curve for Q(t)…")

        def _on_cobs_ready(cobs_data):
            self.btn_fetch_cobs.setEnabled(True)
            if not cobs_data or not cobs_data.get("obs_list"):
                self._cobs_obs_list = None
                self._cobs_ready_cache = None
                self._pending_run_after_cobs = False
                if hasattr(self, 'lbl_cobs_status'):
                    self.lbl_cobs_status.setText(
                        "COBS auto-fetch failed or no usable observations. "
                        "You can edit inputs, retry Fetch COBS, or choose another Q(t) source.")
                self.lbl_status.setText(
                    "Loaded editable inputs, but COBS Q(t) is not ready. "
                    "Retry Fetch COBS before running a COBS-weighted MC model.")
                try:
                    self._refresh_qt_plot()
                except Exception:
                    pass
                return

            self._cobs_obs_list = cobs_data["obs_list"]
            self._cobs_ready_cache = list(self._cobs_obs_list)
            counts = cta.summarize_obs_methods(self._cobs_obs_list)
            method_str = ", ".join(
                f"{name or key}×{n}" for key, (n, name) in list(counts.items())[:4]
            ) if counts else "method info n/a"
            if hasattr(self, 'lbl_cobs_status'):
                self.lbl_cobs_status.setText(
                    f"Ready — {len(self._cobs_obs_list)} COBS points "
                    f"(source: {cobs_data.get('source','?')}). "
                    f"Methods: {method_str}")

            if anchors:
                self.lbl_status.setText(
                    f"Loaded editable inputs; Q(t) ready from COBS + "
                    f"{len(anchors)} Afρ anchor(s). Edit any field or Run MC.")
            else:
                self.lbl_status.setText(
                    "Loaded editable inputs; COBS fetched, but no valid Afρ "
                    "anchor was found. Add an anchor before COBS-weighted MC.")

            try:
                self._refresh_qt_plot()
            except Exception:
                pass
            self._update_run_button_state()

            # If the user clicked Run MC while COBS was not yet cached, run
            # automatically once the auto-fetch finishes.  This keeps the
            # Load → Edit → Run workflow one-click, without requiring the user
            # to refresh the embedded Q(t) preview or press Run a second time.
            if getattr(self, '_pending_run_after_cobs', False):
                self._pending_run_after_cobs = False
                QTimer.singleShot(0, self._run)

        self.main_window._ensure_cobs_fetched(_on_cobs_ready)

    def _load_tycho_photometry(self):
        """Open a Tycho Tracker ICQ-format .txt export, parse its Afρ
        summary and radius-vs-magnitude profile table, and auto-fill
        Observed Afρ + Faintest/Brightest level fields.

        Falls back gracefully: on any parse failure, or if the user
        cancels the file dialog, existing field values are left exactly
        as they were — this is a convenience auto-fill, not a required
        step, so failure here should never block the rest of the
        workflow."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Tycho Tracker photometry file", "",
            "Text file (*.txt);;All files (*)")
        if not path:
            return
        try:
            data = cta.parse_tycho_photometry_file(path)
        except Exception as exc:
            self.lbl_tycho_status.setText(f"⚠ Could not parse file: {exc}")
            self.lbl_tycho_status.setStyleSheet("color:#ff9955; font-size:9px;")
            return

        try:
            sb = cta.suggest_sb_contour_levels(
                data['radii_arcsec'], data['cumulative_mag'])
        except Exception as exc:
            self.lbl_tycho_status.setText(
                f"⚠ Parsed Afρ table OK but could not derive surface "
                f"brightness levels: {exc}")
            self.lbl_tycho_status.setStyleSheet("color:#ff9955; font-size:9px;")
            # Still fill in what we DID get (Afρ) even if the SB-level
            # calculation failed — partial success is still useful.
            if data.get('afrho_cm'):
                self.sp_afrho_cal.setValue(data['afrho_cm'])
                self.lbl_afrho_source.setText(
                    f"from Tycho file ({data.get('obs_date','?')})")
                self.lbl_afrho_source.setStyleSheet("color:#60e0a0; font-size:9px;")
            return

        # Clamp to the spinboxes' own [10, 35] range so an unusually
        # bright/faint comet can't silently push a value out of range
        # (setValue() would clamp anyway, but this keeps the displayed
        # numbers matching what was actually derived, for the status line).
        faint  = max(10.0, min(35.0, sb['suggested_faint']))
        bright = max(10.0, min(35.0, sb['suggested_bright']))
        self.sp_mag_faint.setValue(faint)
        self.sp_mag_bright.setValue(bright)

        if data.get('afrho_cm'):
            self.sp_afrho_cal.setValue(data['afrho_cm'])
            self.lbl_afrho_source.setText(
                f"from Tycho file ({data.get('obs_date','?')})")
            self.lbl_afrho_source.setStyleSheet("color:#60e0a0; font-size:9px;")

        n_pts = len(data['radii_arcsec'])
        self.lbl_tycho_status.setText(
            f"✓ Loaded {n_pts} radius points from "
            f"{data.get('mpc_code','?')} — "
            f"measured SB range {sb['sb_brightest_raw']:.1f}–"
            f"{sb['sb_faintest_raw']:.1f} mag/arcsec², "
            f"suggested levels applied with 0.3 mag margin.")
        self.lbl_tycho_status.setStyleSheet("color:#60e0a0; font-size:9px;")

    def _autofill_afrho(self):
        """Auto-populate sp_afrho_cal from the COBS anchor table when
        the source is COBS and at least one anchor date is close to the
        observation JD. Falls back to the main Analysis panel's last
        computed Afρ if available, otherwise leaves the field for the
        user to fill manually."""
        best_afrho = None
        best_label = None

        # ── Priority 1: COBS anchor table row closest to obs_jd ─────────
        if self.rad_dmdt_cobs.isChecked():
            best_dt = None
            for r in range(self.tbl_anchors.rowCount()):
                d_item = self.tbl_anchors.item(r, 0)
                v_item = self.tbl_anchors.item(r, 1)
                if not (d_item and v_item):
                    continue
                try:
                    t_jd  = cta.date_to_jd(d_item.text().strip())
                    afrho = float(v_item.text().strip())
                    dt    = abs(t_jd - self.obs_jd)
                    if best_dt is None or dt < best_dt:
                        best_dt    = dt
                        best_afrho = afrho
                        days_diff  = (t_jd - self.obs_jd)
                        best_label = (f"from COBS anchor {d_item.text().strip()} "
                                      f"({'obs date' if abs(days_diff)<0.6 else f'{days_diff:+.0f}d from obs'})")
                except Exception:
                    continue

        # ── Priority 2: main-window Analysis panel last Afρ value ───────
        if best_afrho is None and self.main_window is not None:
            try:
                afrho_analysis = self.main_window.last_afrho_cm
                if afrho_analysis and afrho_analysis > 0:
                    best_afrho = afrho_analysis
                    best_label = "from Analysis panel"
            except AttributeError:
                pass

        if best_afrho is not None:
            self.sp_afrho_cal.setValue(best_afrho)
            self.lbl_afrho_source.setText(best_label)
            self.lbl_afrho_source.setStyleSheet("color:#60e0a0; font-size:9px;")
        else:
            self.lbl_afrho_source.setText("enter manually (no anchor found)")
            self.lbl_afrho_source.setStyleSheet("color:#ff9955; font-size:9px;")

    @staticmethod
    def _rline(label: str, value: str, width: int = 28) -> str:
        """One aligned 'label value' report line — single column-width
        helper so every section lines up the same way without each call
        site hand-padding its own f-string (v3.1 report redesign)."""
        return f"{label:<{width}}{value}"

    def _build_report_lines(self) -> list[str]:
        """Build a compact, publication-style Monte Carlo model report.

        v3.1 report format: short section blocks intended for sharing in
        project notes, Facebook technical captions, and manuscript support
        material.  The report deliberately stays input-transparent: it lists
        the values actually visible in the GUI rather than applying any hidden
        object-specific preset.
        """
        import math

        R = self._rline
        lines = []
        dash = "-" * 60

        def _fmt_um(x: float) -> str:
            try:
                xf = float(x)
                if abs(xf - round(xf)) < 1e-6:
                    return f"{int(round(xf)):,}"
                return f"{xf:,.5g}"
            except Exception:
                return str(x)

        def _fmt_beta(x: float) -> str:
            try:
                xf = float(x)
                return f"{xf:.3g}" if abs(xf) < 1e-3 else f"{xf:.4g}"
            except Exception:
                return str(x)

        def _phase_law_label() -> str:
            idx = self.cmb_phase_law.currentIndex() if hasattr(self, 'cmb_phase_law') else 0
            if idx == 1:
                return "Linear–exponential custom [frame-wide scalar]"
            if idx == 2:
                return "None / relative morphology only"
            return "Schleicher composite dust phase function [frame-wide scalar]"

        # ── Header / observation context ─────────────────────────────────
        lines.append("Monte Carlo model report")
        lines.append(dash)
        lines.append(R("Comet:", self.comet_el.get('name', '?')))
        lines.append(R("Observation date (UT):", cta.jd_to_str(self.obs_jd)))

        r_h = None
        try:
            r_C, _ = cta.elem_to_state(self.comet_el, self.obs_jd)
            r_h = float(cta.vmag(r_C))
            lines.append(R("r_H at obs:", f"{r_h:.5f} AU"))
        except Exception:
            pass

        lines.append(f"q={self.comet_el.get('q',0):.6f} AU   "
                     f"e={self.comet_el.get('e',0):.6f}   "
                     f"i={self.comet_el.get('i',0):.4f}°")
        lines.append(R("Orbital element source:", self.comet_el.get('source', 'JPL Horizons')))

        # PsAng here follows the CTA label used on the canvas: apparent Sun
        # direction projected on the sky, i.e. comet → Sun, shown as
        # "sun→comet" in the report convention requested by the user.
        if self.main_window is not None:
            try:
                m = self.main_window._model
                if m:
                    sun_xi, sun_eta = m.get('sun_dir', (0.0, 0.0))
                    psang = math.degrees(math.atan2(sun_xi, sun_eta)) % 360.0
                    lines.append(R("PsAng (sun→comet):", f"{psang:.1f}°"))
            except Exception:
                pass
        lines.append("")

        # ── Grain size ───────────────────────────────────────────────────
        lines.append("GRAIN SIZE")
        lines.append(dash)
        rho = self.sp_rho.value()
        qpr = 1.0
        a_min, a_max = self.sp_r_min.value(), self.sp_r_max.value()
        beta_at_min = cta._radius_um_to_beta(a_min, rho)
        beta_at_max = cta._radius_um_to_beta(a_max, rho)
        lines.append(R("Grain density ρ:", f"{rho:.3g} g/cm³"))
        lines.append(R("Radiation pressure Qpr:", f"{qpr:.1f}"))
        lines.append(R("Grain radius (min):", f"{_fmt_um(a_min)} µm  (β {_fmt_beta(beta_at_min)})"))
        lines.append(R("Grain radius (max):", f"{_fmt_um(a_max)} µm  (β {_fmt_beta(beta_at_max)})"))
        lines.append(R("Size-distribution slope κ:", f"{self.sp_gamma_size.value():.3g}"))
        lines.append(R("β formula:", "β = 0.574·Qpr / (ρ·a)  [Burns, Lamy & Soter 1979 Eq.19]"))
        lines.append("")

        # ── Simulation ───────────────────────────────────────────────────
        lines.append("SIMULATION")
        lines.append(dash)
        lines.append(R("Particles:", f"{self.sp_n.value():,}"))
        lines.append(R("Grid resolution:", f"{self._grid_npix} × {self._grid_npix} px (fixed)"))
        lines.append(R("Dust age (release window):", f"{self.sp_max_age.value():.3g} days [approximate, ±2-3 days]"))

        if self.sp_v0.value() != 0:
            try:
                if r_h is None:
                    r_C, _ = cta.elem_to_state(self.comet_el, self.obs_jd)
                    r_h = float(cta.vmag(r_C))
                v_at_min = cta.real_ejection_speed_ms(
                    self.sp_v0.value(), beta_at_min, r_h,
                    self.sp_gamma.value(), self.sp_mexp.value())
                v_at_max = cta.real_ejection_speed_ms(
                    self.sp_v0.value(), beta_at_max, r_h,
                    self.sp_gamma.value(), self.sp_mexp.value())
                lines.append(R(f"Ejection speed ({_fmt_um(a_min)} µm):", f"{float(v_at_min):.2f} m/s (at r_H={r_h:.3f} AU)"))
                lines.append(R(f"Ejection speed ({_fmt_um(a_max)} µm):", f"{float(v_at_max):.2f} m/s"))
            except Exception:
                pass
        lines.append(R("V0 (β=1 reference):", f"{self.sp_v0.value():.3g} m/s [assumed — not independently fit]"))
        lines.append(R("Velocity law:", f"V = V0·β^{self.sp_gamma.value():.2g}·r_H^-{self.sp_mexp.value():.2g} [Whipple 1951 / Fulle 1987]"))
        if self.chk_seed.isChecked():
            lines.append(R("Random seed:", f"{self.sp_seed.value()} (fixed)"))
        lines.append("")

        # ── Ejection direction ────────────────────────────────────────────
        lines.append("EJECTION DIRECTION")
        lines.append(dash)
        if self.rad_mode_isotropic.isChecked():
            lines.append(R("Mode:", "Isotropic"))
        elif self.rad_mode_sunward.isChecked():
            cone = self.sp_sunward_cone.value()
            if cone >= 89.999:
                mode = "Sunward hemisphere [cos(z) modulation]"
            else:
                mode = f"Sunward cone {cone:.1f}° half-angle [cos(z) modulation]"
            lines.append(R("Mode:", mode))
            lines.append(R("Sunward reference:", "Emission time — physical"))
            lines.append(R("cos(z) exponent ε:", f"{self.sp_sunward_expocos.value():.3g}"))
            if (self.grp_projection_diag.isChecked() and
                    self.chk_projected_sunward.isChecked()):
                lines.append(R("Sky-plane apparent-Sun gate:", "ON [diagnostic; not used as physical default]"))
        else:
            lines.append(R("Mode:", "Active area (rotating nucleus)"))
            lines.append(R("Obliquity I:", f"{self.sp_nuc_inc.value():.3g}°"))
            lines.append(R("Subsolar phase Φ:", f"{self.sp_nuc_phi.value():.3g}°"))
            lines.append(R("Rotation period:", f"{self.sp_period.value():.4g} d"))
            lines.append(R("Active latitude:", f"{self.sp_lat_min.value():.3g}° to {self.sp_lat_max.value():.3g}°"))
            lines.append(R("Active longitude:", f"{self.sp_lon_min.value():.3g}° to {self.sp_lon_max.value():.3g}°"))
            lines.append(R("Sunlit-ground only:", str(self.chk_isun.isChecked())))
            lines.append(R("cos(z) exponent ε:", f"{self.sp_expocos.value():.3g}"))
        lines.append("")

        # ── Dust production / phase law ──────────────────────────────────
        lines.append("DUST PRODUCTION")
        lines.append(dash)
        if self.rad_dmdt_steady.isChecked():
            lines.append(R("Q(t) source:", "Steady constant production"))
        elif self.rad_dmdt_cobs.isChecked():
            lines.append(R("Q(t) source:", "COBS light curve (Afρ proxy)"))
            lines.append("Note: Afρ proxy assumes constant aperture & phase correction")
            for r in range(self.tbl_anchors.rowCount()):
                d = self.tbl_anchors.item(r, 0)
                v = self.tbl_anchors.item(r, 1)
                if d and v and d.text().strip() and v.text().strip():
                    lines.append(R(f"  Afρ anchor {d.text().strip()}:", f"{v.text().strip()} cm"))
            coverage = getattr(self, '_qt_coverage', None)
            if isinstance(coverage, dict):
                lines.append(R("COBS coverage quality:", coverage.get('quality', 'INSUFFICIENT')))
                lines.append(R(
                    "Continuous supported lookback:",
                    f"{float(coverage.get('reliable_lookback_days',0.0)):.1f} days "
                    f"({100.0*float(coverage.get('coverage_fraction',0.0)):.0f}% of F-P maximum age)"))
                lines.append(R("Coverage assessment:", coverage.get('reason', '')))
                extra_days = float(getattr(self, '_qt_run_extrapolation_days', 0.0) or 0.0)
                extrap_method = getattr(self, '_qt_run_extrapolation_method', 'none')
                if extra_days > 1.0:
                    lines.append(R(
                        "Q(t) extrapolated interval:",
                        f"{extra_days:.1f} days [method: {extrap_method}]"))
                    lines.append(
                        "Note: this extrapolation is a user-accepted MC assumption; "
                        "it is not used to validate the automatic window recommendation.")
        else:
            lines.append(R("Q(t) source:", "Manual dM/dt table"))
            for r in range(self.tbl_dmdt_manual.rowCount()):
                d = self.tbl_dmdt_manual.item(r, 0)
                v = self.tbl_dmdt_manual.item(r, 1)
                if d and v and d.text().strip() and v.text().strip():
                    lines.append(R(f"  days={d.text().strip()}:", f"rel.dM/dt={v.text().strip()}"))
        lines.append(R("Grain albedo p_v:", f"{self.sp_pv.value():.3g} [assumed — Hanner 1981]"))
        lines.append(R("Phase law:", _phase_law_label()))
        if hasattr(self, 'cmb_phase_law') and self.cmb_phase_law.currentIndex() == 1:
            lines.append(R("  β_α, m_oe, w_oe:", f"{self.sp_phase_beta.value():.3g} mag/deg, {self.sp_phase_moe.value():.3g} mag, {self.sp_phase_woe.value():.3g}°"))
        lines.append(R("Phase-law contour-shape effect:", "None for one epoch; brightness scale only"))
        lines.append("")

        # ── Contour display ──────────────────────────────────────────────
        lines.append("CONTOUR DISPLAY")
        lines.append(dash)
        if self.chk_sb_mode.isChecked():
            lines.append(R("Mode:", "Calibrated surface brightness [mag/arcsec²]"))
            lines.append(R("  Observed Afρ:", f"{self.sp_afrho_cal.value():.1f} cm"))
            lines.append(R("  Faintest level:", f"{self.sp_mag_faint.value():.1f} mag/arcsec²"))
            lines.append(R("  Brightest level:", f"{self.sp_mag_bright.value():.1f} mag/arcsec²"))
        else:
            lines.append(R("Mode:", "Relative percentile morphology contour"))
            lines.append(R("  Sensitivity floor:", f"{self.sp_floor.value():.3g} %ile"))
        lines.append(R("Contour rings:", f"{self.sp_nlevels.value()}"))
        try:
            lines.append(R("Presentation style:", self.cmb_mc_view_style.currentText()))
            lines.append(R("Background:", self.cmb_mc_background.currentText()))
            lines.append(R("Coordinates:", self.cmb_mc_coordinates.currentText()))
            lines.append(R("Observed isophotes:",
                           f"{self.cmb_observed_color.currentText()}, "
                           f"{self.cmb_observed_ls.currentText()}, "
                           f"{self.sp_observed_lw.value():.1f} pt"))
            lines.append(R("MC model contours:",
                           f"{self.cmb_model_color.currentText()}, "
                           f"{self.cmb_model_ls.currentText()}, "
                           f"{self.mc_lw_slider.value()/10:.1f} pt"))
        except Exception:
            pass
        lines.append("")

        # ── F-P/Q(t) bridge notes ─────────────────────────────────────────
        try:
            lines.append("F-P / Q(t) WINDOW GUIDANCE")
            lines.append(dash)
            dom = float(getattr(self, 'fp_dominant_age', 0.0) or 0.0)
            max_age = float(getattr(self, 'fp_max_age', 0.0) or 0.0)
            lines.append(R(
                "F-P dominant dust age:",
                f"{dom:.3g} days" if dom > 0 else "Not specified"))
            lines.append(R(
                "F-P maximum dust age:",
                f"{max_age:.3g} days [largest synchrone age]"
                if max_age > 0 else "Unavailable"))
            if getattr(self, '_qt_qonly_window', None) is not None:
                lines.append(R(
                    "Q(t)-only recommendation:",
                    f"{self._qt_qonly_window:.3g} days"))
            if getattr(self, '_qt_suggested_window', None) is not None:
                lines.append(R(
                    "Final recommended window:",
                    f"{self._qt_suggested_window:.3g} days"))
            elif getattr(self, '_qt_provisional_window', None) is not None:
                lines.append(R(
                    "Provisional Q(t) value:",
                    f"{self._qt_provisional_window:.3g} days [not automatically applicable]"))
            lines.append(R(
                "Interpretation:",
                "Effective morphology window, not unique activity onset"))
            lines.append("")
        except Exception:
            pass

        # ── Software citation ────────────────────────────────────────────
        lines.append("SOFTWARE")
        lines.append(dash)
        lines.append("Thaluang, T. (2026). Comet Tail Analyzer (CTA) v3.1.1.")
        lines.append("RNAAS, doi:10.3847/2515-5172/ae6f90")
        lines.append("Portions ported from py_COMTAILS (Moreno 2025, A&A 695, A263).")
        return lines

    def _preview_report(self):
        """Show the report exactly as it would be written, before
        committing to a file — a read-only monospace preview with its
        own Save button, so a typo or wrong tab doesn't get baked into
        a saved .txt with no chance to check first."""
        text = "\n".join(self._build_report_lines())
        dlg = QDialog(self)
        dlg.setWindowTitle("Monte Carlo report — preview")
        dlg.resize(560, 640)
        v = QVBoxLayout(dlg)
        edit = QTextEdit()
        edit.setReadOnly(True)
        edit.setPlainText(text)
        edit.setFont(QFont("Consolas, Courier New, monospace", 10))
        v.addWidget(edit)
        bot = QHBoxLayout()
        btn_save = QPushButton("💾 Save as…")
        btn_save.clicked.connect(lambda: self._save_report(text))
        btn_close = QPushButton("Close"); btn_close.setMaximumWidth(90)
        btn_close.clicked.connect(dlg.close)
        bot.addWidget(btn_save); bot.addStretch(); bot.addWidget(btn_close)
        v.addLayout(bot)
        dlg.exec()

    def _save_report(self, text: str | None = None):
        """Write the report to a plain-text file. If called from the
        preview dialog, `text` is the exact text already on screen there
        (no re-build, so what you previewed is what gets saved even if
        a field changed in the background somehow); called directly
        from the button it rebuilds fresh."""
        if text is None:
            text = "\n".join(self._build_report_lines())
        default_name = (
            f"MC_Report_{self.comet_el.get('name','comet').replace(' ', '_').replace('/', '')}"
            f"_{cta.jd_to_str(self.obs_jd)[:10]}.txt")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Monte Carlo report", default_name, "Text file (*.txt)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text + "\n")
            self.lbl_status.setText(f"Report saved to {path}")
        except Exception as exc:
            QMessageBox.warning(self, "Save Report", f"Could not save report:\n{exc}")

    def _assess_current_qt_coverage(self, qt_result: dict) -> dict:
        """Validate COBS support for inference and for the selected MC window."""
        coverage = cta.assess_qt_coverage(
            qt_result.get("t_jd", []), qt_result.get("Q_kg_s", []),
            self.obs_jd,
            max_age_days=float(getattr(self, 'fp_max_age', 0.0) or 0.0),
            dominant_age_days=float(getattr(self, 'fp_dominant_age', 0.0) or 0.0),
            smooth_window_days=float(self.sp_qt_smooth.value()))
        self._qt_coverage = coverage
        self._qt_recommendation_quality = coverage.get('quality', '')
        return coverage

    def _validated_cobs_weights(self, qt_result: dict):
        """Return the recent continuous COBS weights for MC sampling.

        The coverage guard remains strict for *automatic inference*: PARTIAL
        or INSUFFICIENT coverage cannot produce an automatically applicable
        MC-window recommendation.  A user-selected release window may still be
        simulated when a usable recent COBS segment exists.  If that window is
        longer than the segment, CTA asks for confirmation and then relies on
        the sampler's documented flat edge extrapolation: the earliest accepted
        Q(t) value is held constant over the unsupported older interval.
        """
        coverage = self._assess_current_qt_coverage(qt_result)
        t, q, source = self._usable_qt_sampling_arrays(qt_result)
        if t.size < 2 or q.size != t.size:
            QMessageBox.warning(
                self, "Invalid COBS Q(t)",
                "The calibrated COBS Q(t) curve does not contain at least two "
                "valid points at or before the observation date. Choose "
                "Steady/Manual production or obtain additional COBS data.")
            return None

        lookback = max(0.0, float(self.obs_jd) - float(np.min(t)))
        current_window = float(self.sp_max_age.value())
        extra = max(0.0, current_window - lookback)
        self._qt_run_extrapolation_days = extra
        self._qt_run_extrapolation_method = "none"

        if extra > 1.0:
            answer = QMessageBox.question(
                self, "Partial COBS coverage",
                f"Current MC release window: {current_window:.0f} d\n"
                f"Available calibrated Q(t) lookback: {lookback:.0f} d\n"
                f"Unsupported older interval: {extra:.0f} d\n\n"
                "The automatic MC-window recommendation remains disabled. "
                "To run the user-selected model, CTA will hold the earliest "
                "available Q(t) value constant across the unsupported older "
                "interval (flat edge extrapolation).\n\n"
                "Continue with this explicit assumption?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return None
            self._qt_run_extrapolation_method = "flat earliest-edge Q(t)"
        elif source == 'raw-qt-fallback':
            self._qt_run_extrapolation_method = (
                "none; raw calibrated Q(t) used because automatic F-P/coverage "
                "guidance was unavailable")
        elif coverage.get('quality') != 'RELIABLE':
            self._qt_run_extrapolation_method = "none; coverage provisional"

        return dict(t_jd=t, Q_kg_s=q)

    def _run(self):
        r_min_um = self.sp_r_min.value()
        r_max_um = self.sp_r_max.value()
        if r_min_um >= r_max_um:
            QMessageBox.warning(self, "Bad Range", "r_min must be less than r_max.")
            return
        # β ∝ 1/radius: smallest grain → largest β
        beta_max = cta._radius_um_to_beta(r_min_um, self.sp_rho.value())
        beta_min = cta._radius_um_to_beta(r_max_um, self.sp_rho.value())

        qt_weights = None
        self._last_qt_dips = []
        if self.rad_dmdt_cobs.isChecked():
            anchors = self._read_anchor_table()
            if not anchors:
                QMessageBox.warning(self, "No Afρ Anchor",
                    "Enter at least one (Date, Afρ) row — your own measured "
                    "value, e.g. from Tycho Tracker — or choose a different "
                    "dM/dt source.")
                return

            obs_for_qt = self._cobs_obs_list
            if not obs_for_qt:
                # Use the persistent cache filled by auto-fetch / explicit fetch.
                obs_for_qt = getattr(self, '_cobs_ready_cache', None)
                if obs_for_qt:
                    self._cobs_obs_list = obs_for_qt

            if not obs_for_qt:
                # As a last resort, if the Q(t) preview already computed a valid
                # curve, use that cached Q(t) directly rather than blocking Run.
                cached_qt = getattr(self, '_last_qt_result', None)
                if cached_qt is not None and len(cached_qt.get('t_jd', [])) >= 2:
                    qt_weights = self._validated_cobs_weights(cached_qt)
                    if qt_weights is None:
                        return
                    self._last_qt_dips = cta.find_qt_dips(
                        qt_weights["t_jd"], qt_weights["Q_kg_s"])
                else:
                    # Editable-input workflow: pressing Run should prepare COBS
                    # automatically if the input selected COBS mode but the cache
                    # is empty for any reason.
                    if self.main_window is not None:
                        self._pending_run_after_cobs = True
                        self.lbl_status.setText(
                            "COBS Q(t) was not cached. Auto-fetching now; MC will start when ready.")
                        self._auto_prepare_qt_after_input_load()
                        return
                    QMessageBox.warning(self, "No COBS Data",
                        "COBS Q(t) source is selected, but no COBS data is cached. "
                        "Fetch COBS light curve first, or choose a different dM/dt source.")
                    return
            else:
                try:
                    qt_result = cta.estimate_qt_from_lightcurve(
                        self.comet_el, obs_for_qt, anchors,
                        p_v=float(self.sp_pv.value()),
                        smooth_window_days=self.sp_qt_smooth.value())
                    self._last_qt_result = qt_result
                    qt_weights = self._validated_cobs_weights(qt_result)
                    if qt_weights is None:
                        return
                    self._last_qt_dips = cta.find_qt_dips(
                        qt_weights["t_jd"], qt_weights["Q_kg_s"])
                except Exception as e:
                    QMessageBox.warning(self, "Q(t) Estimation Error", str(e))
                    return
        elif self.rad_dmdt_manual.isChecked():
            try:
                days_arr, dmdt_arr = self._read_dmdt_manual_table()
            except ValueError as e:
                QMessageBox.warning(self, "Bad dM/dt Table", str(e))
                return
            per_jd_dmdt = self.comet_el.get("T_jd") or cta.date_to_jd(self.comet_el["T"])
            qt_weights = dict(t_jd=per_jd_dmdt + days_arr, Q_kg_s=dmdt_arr)
        # else: rad_dmdt_steady — qt_weights stays None, unchanged default behaviour

        size_dist_table = None
        if self.chk_size_table.isChecked():
            try:
                size_dist_table = self._read_size_over_time_table()
            except ValueError as e:
                QMessageBox.warning(self, "Bad Grain-Size Table", str(e))
                return

        self.btn_run.setEnabled(False)
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setFormat("0%")
        self.lbl_status.setText("Preparing Monte Carlo run…")

        active_area = None
        sunward = False
        if self.rad_mode_active.isChecked():
            if self.sp_lat_min.value() >= self.sp_lat_max.value():
                QMessageBox.warning(self, "Bad Range", "Latitude min must be less than max.")
                self.btn_run.setEnabled(True); self.progress.setRange(0, 1)
                return
            if self.sp_lon_min.value() >= self.sp_lon_max.value():
                QMessageBox.warning(self, "Bad Range", "Longitude min must be less than max.")
                self.btn_run.setEnabled(True); self.progress.setRange(0, 1)
                return
            per_jd = self.comet_el.get("T_jd") or cta.date_to_jd(self.comet_el["T"])
            active_area = dict(
                nuc_inc_deg=self.sp_nuc_inc.value(), nuc_phi_deg=self.sp_nuc_phi.value(),
                period_d=self.sp_period.value(), per_jd=per_jd,
                lat_min_deg=self.sp_lat_min.value(), lat_max_deg=self.sp_lat_max.value(),
                lon_min_deg=self.sp_lon_min.value(), lon_max_deg=self.sp_lon_max.value(),
                isun=self.chk_isun.isChecked(), expocos=self.sp_expocos.value())
        elif self.rad_mode_sunward.isChecked():
            sunward = True

        seed = self.sp_seed.value() if self.chk_seed.isChecked() else None
        self._worker = MCWorker(
            self.comet_el, self.obs_jd, (beta_min, beta_max),
            self.sp_gamma_size.value(), self.sp_max_age.value(), self.sp_n.value(),
            self.sp_v0.value(), self.sp_gamma.value(), self.sp_mexp.value(), seed,
            qt_weights=qt_weights, active_area=active_area, rho_g_cm3=self.sp_rho.value(),
            sunward=sunward, sunward_expocos=self.sp_sunward_expocos.value(),
            sunward_reference="emission",
            sunward_cone_half_angle_deg=self.sp_sunward_cone.value(),
            require_projected_sunward=(self.grp_projection_diag.isChecked() and
                                         self.chk_projected_sunward.isChecked()),
            phase_law=("linear_exponential" if self.cmb_phase_law.currentIndex() == 1 else
                       "none" if self.cmb_phase_law.currentIndex() == 2 else "schleicher"),
            phase_linear_beta=self.sp_phase_beta.value(),
            phase_linear_m_oe=self.sp_phase_moe.value(),
            phase_linear_w_oe=self.sp_phase_woe.value(),
            size_dist_table=size_dist_table, p_v=self.sp_pv.value(),
            grid_npix=self._grid_npix)
        self._run_generation += 1
        run_generation = self._run_generation

        def _mc_progress(pct, message, generation=run_generation):
            if generation == self._run_generation:
                self._on_progress(pct, message)

        def _mc_done(result, generation=run_generation):
            if generation == self._run_generation:
                self._on_done(result)

        def _mc_error(message, generation=run_generation):
            if generation == self._run_generation:
                self._on_error(message)

        self._worker.progress.connect(_mc_progress)
        self._worker.finished.connect(_mc_done)
        self._worker.error.connect(_mc_error)
        self._worker.start()

    def _read_anchor_table(self) -> list:
        import comet_tail_analyzer as cta
        """Parse (date_str, afrho_cm) rows out of tbl_anchors, silently
        skipping incomplete/unparseable rows (a blank trailing row from
        the default single empty row, or a typo) rather than erroring —
        only rows with BOTH a valid date and a positive number count."""
        anchors = []
        for row in range(self.tbl_anchors.rowCount()):
            date_item = self.tbl_anchors.item(row, 0)
            afrho_item = self.tbl_anchors.item(row, 1)
            if date_item is None or afrho_item is None:
                continue
            date_str = date_item.text().strip()
            afrho_str = afrho_item.text().strip()
            if not date_str or not afrho_str:
                continue
            try:
                cta.date_to_jd(date_str)   # validate parseable, value unused here
                afrho_val = float(afrho_str)
                if afrho_val > 0:
                    anchors.append((date_str, afrho_val))
            except Exception:
                continue
        return anchors

    def _read_dmdt_manual_table(self):
        """Parse (days_to_perihelion, relative_dM/dt) rows out of
        tbl_dmdt_manual, matching py_COMTAILS's own table-input
        convention. Raises ValueError (caught by the caller, shown as a
        QMessageBox) rather than silently skipping bad rows here — unlike
        the Afρ anchor table, this table IS the entire dM/dt(t) shape, so
        a silently-dropped row would change the result without any
        warning at all, instead of just being one option among several."""
        import numpy as np
        days, vals = [], []
        for row in range(self.tbl_dmdt_manual.rowCount()):
            d_item = self.tbl_dmdt_manual.item(row, 0)
            v_item = self.tbl_dmdt_manual.item(row, 1)
            if d_item is None or v_item is None or not d_item.text().strip():
                continue
            try:
                d = float(d_item.text().strip())
                v = float(v_item.text().strip())
            except ValueError:
                raise ValueError(f"Row {row+1}: \"{d_item.text()}\" / \"{v_item.text()}\" "
                                 f"isn't a valid number pair.")
            if v <= 0:
                raise ValueError(f"Row {row+1}: relative dM/dt must be positive, got {v}.")
            days.append(d); vals.append(v)
        if len(days) < 2:
            raise ValueError("Need at least 2 rows to define a dM/dt(t) shape.")
        order = np.argsort(days)
        return np.asarray(days)[order], np.asarray(vals)[order]

    def _read_size_over_time_table(self) -> dict:
        """Parse (days_to_perihelion, power, r_min_um, r_max_um) rows out
        of tbl_size_over_time into compute_morphology_mc()'s
        size_dist_table format. Raises ValueError (caught by the caller)
        on any bad/incomplete row, same reasoning as _read_dmdt_manual_
        table() — this table fully determines the simulated size
        distribution when active, so a silently-dropped row would be a
        silent change in physics, not a minor omission."""
        per_jd = self.comet_el.get("T_jd") or cta.date_to_jd(self.comet_el["T"])
        days, power, rmin, rmax = [], [], [], []
        for row in range(self.tbl_size_over_time.rowCount()):
            items = [self.tbl_size_over_time.item(row, c) for c in range(4)]
            if any(it is None or not it.text().strip() for it in items):
                continue
            try:
                d, p, rn, rx = (float(it.text().strip()) for it in items)
            except ValueError:
                raise ValueError(f"Row {row+1} contains a non-numeric value.")
            if not (0 < rn < rx):
                raise ValueError(f"Row {row+1}: need 0 < r_min < r_max (got {rn}, {rx}).")
            days.append(d); power.append(p); rmin.append(rn); rmax.append(rx)
        if len(days) < 2:
            raise ValueError("Need at least 2 rows to define how grain size changes over time.")
        return dict(per_jd=per_jd, days_to_per=days, power=power,
                   r_min_um=rmin, r_max_um=rmax)

    def _remove_anchor_row(self):
        row = self.tbl_anchors.currentRow()
        if row < 0:
            row = self.tbl_anchors.rowCount() - 1
        if row >= 0:
            self.tbl_anchors.removeRow(row)

    def _fetch_cobs_for_qt(self):
        if self.main_window is None:
            return
        self.btn_fetch_cobs.setEnabled(False)
        self.lbl_cobs_status.setText("Fetching…")

        def _on_cobs_ready(cobs_data):
            self.btn_fetch_cobs.setEnabled(True)
            if not cobs_data or not cobs_data.get("obs_list"):
                self._cobs_obs_list = None
                self._cobs_ready_cache = None
                self.lbl_cobs_status.setText(
                    "No usable COBS data (no light curve found, or "
                    "ephemeris matching failed).")
                return
            self._cobs_obs_list = cobs_data["obs_list"]
            self._cobs_ready_cache = list(self._cobs_obs_list)
            counts = cta.summarize_obs_methods(self._cobs_obs_list)
            method_str = ", ".join(f"{name or key}×{n}" for key, (n, name) in
                                   list(counts.items())[:4]) if counts else "method info n/a"
            self.lbl_cobs_status.setText(
                f"{len(self._cobs_obs_list)} COBS points with ephemeris "
                f"(source: {cobs_data.get('source','?')}). Methods: {method_str}")
            # Keep the Q(t) preview/suggestion synchronized whenever COBS
            # is fetched explicitly from the Simulate tab.
            try:
                self._refresh_qt_plot()
            except Exception:
                pass
            self._update_run_button_state()

        self.main_window._ensure_cobs_fetched(_on_cobs_ready)

    def _on_progress(self, pct: int, message: str):
        pct = max(0, min(100, int(pct)))
        self.progress.setRange(0, 100)
        self.progress.setValue(pct)
        self.progress.setFormat(f"{pct}%")
        self.lbl_status.setText(message)

    def _on_done(self, result: dict):
        self._mc_result = result
        self._mark_guided_setup_complete()
        self.progress.setRange(0, 100)
        self.progress.setValue(100)
        self.progress.setFormat("100%")
        self.btn_run.setEnabled(True)
        self.btn_reextract.setEnabled(True)
        info = result['info']

        # Auto-fill Smoothing (AU) to match the OLD pixel-based default
        # (1.5px) for THIS run's grid — but ONLY the first time. Doing
        # this on every run would silently reproduce the exact bug this
        # control exists to fix: re-syncing to "1.5 pixels of THIS run's
        # grid" every time means a smaller r_min (bigger grid) still gets
        # a bigger absolute smoothing radius automatically, same as
        # before. Filling it once gives a sensible starting point; from
        # then on the value the user sees is the value actually used,
        # stable across later r_min/r_max changes unless they edit it
        # themselves.
        # v3.1 — update smooth floor on EVERY run (grid size may change
        # if grain range changed). Only move the slider if the user hasn't
        # manually edited it — if they have, preserve their AU value but
        # clamp upward to the new floor if needed.
        self.sp_smooth_au.setValue(0.0)
        self._smooth_au_initialized = True

        # Sync ejection velocity into the main panel + recompute the F-P
        # model with it, so both windows agree (v3.1 — see class docstring
        # for why v0 maps to v_N0 specifically; edit here to change that).
        # BUG FIX (v3.1): must WAIT for that recompute to actually finish
        # before drawing the MC contour — recompute_and_then() handles
        # this; calling ctrl._emit_compute() and immediately continuing
        # (the old code) drew the MC contour against whatever F-P model
        # was left over from a PREVIOUS run, since _emit_compute() only
        # starts a background computation and returns right away. See
        # recompute_and_then()'s own docstring for the full symptom this
        # caused (F-P appearing to ignore a just-changed v0=0).
        if self.main_window is not None:
            self.main_window.ctrl.set_ejection_params(dict(
                v_R0=0.0, v_T0=0.0, v_N0=info['v0_coeff'],
                gamma=info['gamma'], m_exp=info['m_exp']))
            self.main_window.recompute_and_then(lambda: self._extract_and_send(info))
        else:
            self._extract_and_send(info)

    def _extract_and_send(self, info):
        """Shared by _on_done() and _reextract_contours()."""
        try:
            if self.chk_sb_mode.isChecked():
                afrho_cm  = self.sp_afrho_cal.value()
                delta_au  = info.get('r_geo_au', 1.5)
                r_h_au    = info.get('r_helio_au', 2.0)
                au_per_px = (self.main_window.ctrl.au_px_spin.value()
                             if self.main_window else 0.003)
                sb_result = cta.calibrate_mc_to_surface_brightness(
                    self._mc_result, afrho_cm=afrho_cm,
                    delta_au=delta_au, r_h_au=r_h_au, au_per_px=au_per_px)
                n_rings    = self.sp_nlevels.value()
                mag_faint  = self.sp_mag_faint.value()
                mag_bright = self.sp_mag_bright.value()
                sb_peak    = sb_result.get('sb_peak',  float('nan'))
                sb_median  = sb_result.get('sb_median', float('nan'))
                # Clamp levels to physically meaningful range
                import math
                if math.isfinite(sb_peak) and math.isfinite(sb_median):
                    mag_bright = max(mag_bright, sb_peak + 0.5)
                    mag_faint  = min(mag_faint,  sb_median - 1.0)
                mag_levels = np.linspace(mag_faint, mag_bright, n_rings).tolist()
                paths = cta.extract_contours_at_magnitude_levels(
                    sb_result, mag_levels=mag_levels)
                sb_info = (f"  model peak={sb_peak:.1f}, median={sb_median:.1f} mag/arcsec²"
                           if math.isfinite(sb_peak) else "")
                mode_note = (f"SB-calibrated  Afρ={afrho_cm:.0f} cm  "
                             f"levels={mag_bright:.1f}–{mag_faint:.1f} mag/arcsec²{sb_info}")
            else:
                paths = cta.extract_morphology_contours(
                    self._mc_result,
                    n_levels=self.sp_nlevels.value(),
                    percentile_floor=self.sp_floor.value(),
                    smooth_sigma_au=self.sp_smooth_au.value())
                mode_note = (f"floor={self.sp_floor.value():.0f}%ile, "
                             f"{self.sp_nlevels.value()} levels")
        except Exception as e:
            self.lbl_status.setText(f"Contour extraction error: {e}")
            return

        if not paths:
            self.lbl_status.setText(
                f"n_used={info['n_used']}/{info['n_particles']} — not enough "
                f"signal. Try adjusting thresholds or increasing particles.")
            return

        if self.main_window is not None:
            self.main_window.set_mc_contours(paths, info=info)
            qt_note = "  ·  Q(t)-weighted release times" if info.get('qt_weighted') else ""
            if info.get('active_area_used'):
                qt_note += "  ·  active-area ejection"
            elif info.get('sunward_used'):
                cone = float(info.get('sunward_cone_half_angle_deg', 90.0))
                ref = info.get('sunward_reference', 'emission')
                qt_note += ("  ·  sunward hemisphere" if cone >= 89.999
                            else f"  ·  sunward cone {cone:.0f}°")
                qt_note += f" ({ref} reference)"
                if info.get('require_projected_sunward'):
                    qt_note += "  ·  apparent-Sun gate"
            self.lbl_status.setText(
                f"n_used={info['n_used']}/{info['n_particles']}  ·  "
                f"{len(paths)} contour path(s) ({mode_note}) sent to main canvas  ·  "
                f"main F-P model recomputed with v_N0={info['v0_coeff']:.2f} m/s, "
                f"γ={info['gamma']:.2g}, m={info['m_exp']:.2g}.{qt_note}")
        else:
            self.lbl_status.setText(
                f"n_used={info['n_used']}/{info['n_particles']}  ·  "
                f"{len(paths)} contour path(s) ({mode_note})")

    def _reextract_contours(self):
        """Re-threshold the CACHED Monte Carlo result with the current
        Floor/N levels values — no new sampling, so this is fast even
        though Run Monte Carlo itself (N=tens of thousands of particles)
        is not. Use this to dial in Floor/N levels without waiting for a
        fresh run each time."""
        if self._mc_result is None:
            return
        self._extract_and_send(self._mc_result['info'])

    def _export_current_mc_figure(self):
        """Export the main-canvas Matplotlib figure without any Qt chrome.

        The current Display-tab presentation state is used exactly as shown.
        Raster formats are written at 300 dpi; PDF/SVG retain vector contour
        paths, which is preferable for journal figures and later layout work.
        """
        if self.main_window is None or self.main_window._model is None:
            QMessageBox.information(
                self, "No Figure",
                "Run the F-P model and Monte Carlo model before exporting.")
            return

        import os
        import re

        name = self.comet_el.get("name", "Comet")
        name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "Comet"
        try:
            obs = cta.jd_to_str(self.obs_jd).replace(" UT", "")[:10]
        except Exception:
            obs = "epoch"
        mode = self.cmb_mc_view_style.currentData() or "figure"
        default_name = f"{name}_{obs}_MC_{mode}.png"

        path, selected_filter = QFileDialog.getSaveFileName(
            self, "Export Monte Carlo figure", default_name,
            "PNG image (*.png);;TIFF image (*.tif *.tiff);;PDF vector (*.pdf);;SVG vector (*.svg)")
        if not path:
            return

        # Add a sensible suffix if the user omitted one.
        root, ext = os.path.splitext(path)
        if not ext:
            if "TIFF" in selected_filter:
                ext = ".tif"
            elif "PDF" in selected_filter:
                ext = ".pdf"
            elif "SVG" in selected_filter:
                ext = ".svg"
            else:
                ext = ".png"
            path = root + ext

        try:
            fig = self.main_window.canvas.fig
            fig.savefig(
                path,
                dpi=300,
                bbox_inches="tight",
                facecolor=fig.get_facecolor(),
                edgecolor="none")
            self.lbl_status.setText(f"Figure exported: {path}")
        except Exception as exc:
            QMessageBox.warning(
                self, "Export Failed",
                f"Could not export the figure:\n{exc}")

    def _on_error(self, msg: str):
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setFormat("Error")
        self.btn_run.setEnabled(True)
        self.lbl_status.setText(f"Error: {msg}")
        QMessageBox.warning(self, "Monte Carlo Error", msg)


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN WINDOW
# ─────────────────────────────────────────────────────────────────────────────
def _model_extent_au(model) -> tuple[float, float]:
    """Max |xi|, |eta| (AU) across every syndyne+synchrone in one
    compute_model() frame, ignoring NaNs. Used by the Animator's
    Auto-suggest button to size the frame to fit the widest tail in the
    selected date range."""
    max_xi = max_eta = 0.0
    for s in model.get("syndynes", []) + model.get("synchrones", []):
        xi = s.get("xi"); eta = s.get("eta")
        if xi is not None and np.any(np.isfinite(xi)):
            max_xi = max(max_xi, float(np.nanmax(np.abs(xi))))
        if eta is not None and np.any(np.isfinite(eta)):
            max_eta = max(max_eta, float(np.nanmax(np.abs(eta))))
    return max_xi, max_eta


# ─────────────────────────────────────────────────────────────────────────────
#  DUST PARTICLE RADIUS CALCULATOR  (v3.0) — standalone β→radius tool,
#  Calculation > Dust particle radius… Needs only β/ρ_d/Qpr, no comet or
#  model required, so it can be opened any time. Always delegates to
#  beta_to_size() — never recomputes the formula itself — so this can
#  never drift out of sync with any other place that reports a grain
#  radius, the way the v2.4 duplicated-constant bug did.
class GrainRadiusDialog(QDialog):
    def __init__(self, parent, default_rho: float = 0.5):
        super().__init__(parent)
        self.setWindowTitle("Dust Particle Radius Calculator")
        self.setMinimumWidth(440)
        v = QVBoxLayout(self)

        formula = QLabel(
            "<b>r&nbsp;[µm]&nbsp;=&nbsp;0.574&nbsp;·&nbsp;Qpr&nbsp;/&nbsp;"
            "(ρ&nbsp;[g/cm³]&nbsp;·&nbsp;β)</b><br>"
            "<span style='color:#7090b0; font-size:10px;'>"
            "Burns, Lamy &amp; Soter (1979), Icarus 40, 1, Eq. 19 — r is "
            "the particle's RADIUS (not diameter).</span>")
        formula.setWordWrap(True)
        formula.setTextFormat(Qt.TextFormat.RichText)
        formula.setStyleSheet(
            "padding:10px; background:rgba(255,255,255,12); border-radius:4px;")
        v.addWidget(formula)

        # ── Mode toggle ──────────────────────────────────────────────
        mode_row = QHBoxLayout()
        self.rb_single = QRadioButton("Single β")
        self.rb_range  = QRadioButton("β range")
        self.rb_single.setChecked(True)
        mode_row.addWidget(self.rb_single)
        mode_row.addWidget(self.rb_range)
        v.addLayout(mode_row)

        self.single_box = QWidget()
        sf = QFormLayout(self.single_box); sf.setContentsMargins(0,4,0,4)
        self.sp_beta = QDoubleSpinBox()
        self.sp_beta.setRange(1e-6, 100.0)
        self.sp_beta.setDecimals(6)
        self.sp_beta.setValue(0.01)
        sf.addRow("β:", self.sp_beta)
        v.addWidget(self.single_box)

        self.range_box = QWidget()
        rf = QFormLayout(self.range_box); rf.setContentsMargins(0,4,0,4)
        self.sp_beta_min = QDoubleSpinBox()
        self.sp_beta_min.setRange(1e-6, 100.0)
        self.sp_beta_min.setDecimals(6)
        self.sp_beta_min.setValue(0.001)
        self.sp_beta_max = QDoubleSpinBox()
        self.sp_beta_max.setRange(1e-6, 100.0)
        self.sp_beta_max.setDecimals(6)
        self.sp_beta_max.setValue(0.1)
        rf.addRow("β min:", self.sp_beta_min)
        rf.addRow("β max:", self.sp_beta_max)
        v.addWidget(self.range_box)
        self.range_box.setVisible(False)

        self.rb_single.toggled.connect(self._on_mode_toggle)
        for w in (self.sp_beta, self.sp_beta_min, self.sp_beta_max):
            w.valueChanged.connect(self._recompute)

        # ── Other parameters (defaulted, all editable) ──────────────
        form2 = QFormLayout()
        self.sp_rho = QDoubleSpinBox()
        self.sp_rho.setRange(0.05, 5.0)
        self.sp_rho.setDecimals(3)
        self.sp_rho.setSingleStep(0.05)
        self.sp_rho.setSuffix("  g/cm³")
        self.sp_rho.setValue(default_rho)
        self.sp_rho.valueChanged.connect(self._recompute)
        form2.addRow("Grain density ρ:", self.sp_rho)
        form2.addRow(QLabel(
            "<i>default 0.5 g/cm³ — Fulle et al. (2016), in-situ Rosetta/67P</i>"))

        self.sp_qpr = QDoubleSpinBox()
        self.sp_qpr.setRange(0.01, 5.0)
        self.sp_qpr.setDecimals(3)
        self.sp_qpr.setSingleStep(0.1)
        self.sp_qpr.setValue(1.0)
        self.sp_qpr.valueChanged.connect(self._recompute)
        form2.addRow("Scattering efficiency Qpr:", self.sp_qpr)
        form2.addRow(QLabel(
            "<i>default 1 — appropriate for grains large relative to the "
            "wavelength (Burns et al. 1979)</i>"))
        v.addLayout(form2)

        # ── Output ───────────────────────────────────────────────────
        self.lbl_result = QLabel()
        self.lbl_result.setStyleSheet(
            "font-size:15px; font-weight:bold; padding:10px; "
            "background:rgba(60,120,200,30); border-radius:4px;")
        self.lbl_result.setWordWrap(True)
        v.addWidget(self.lbl_result)

        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        btns.rejected.connect(self.accept)
        v.addWidget(btns)

        self._recompute()

    def _on_mode_toggle(self, _checked):
        self.single_box.setVisible(self.rb_single.isChecked())
        self.range_box.setVisible(not self.rb_single.isChecked())
        self._recompute()

    def _recompute(self):
        rho = self.sp_rho.value()
        qpr = self.sp_qpr.value()
        if self.rb_single.isChecked():
            b = self.sp_beta.value()
            r = cta.beta_to_size(b, rho, qpr)
            self.lbl_result.setText(f"β = {b:g}  →  grain radius = {r}")
        else:
            b_lo = min(self.sp_beta_min.value(), self.sp_beta_max.value())
            b_hi = max(self.sp_beta_min.value(), self.sp_beta_max.value())
            r_at_lo = cta.beta_to_size(b_lo, rho, qpr)   # smaller β → larger grain
            r_at_hi = cta.beta_to_size(b_hi, rho, qpr)
            self.lbl_result.setText(
                f"β = {b_lo:g}  →  grain radius = {r_at_lo}\n"
                f"β = {b_hi:g}  →  grain radius = {r_at_hi}")


# ─────────────────────────────────────────────────────────────────────────────
#  DUST PRODUCTION RATE CALCULATOR  (v3.0) — standalone Afρ→Q_d tool,
#  Calculation > Dust production rate… Requires an active comet (to get
#  r_h from Horizons for the date entered); delegates to
#  compute_dust_production_rate() — the same function generate_dust_
#  analysis()'s inline Afρ section uses — so the two can never drift apart.
class DustProductionDialog(QDialog):
    def __init__(self, parent, comet_el: dict, obs_jd: float,
                v_dust_default: float | None, p_v_default: float):
        super().__init__(parent)
        self.comet_el = comet_el
        self.setWindowTitle(f"Dust Production Rate Calculator — {comet_el.get('name','Comet')}")
        self.setMinimumWidth(460)
        v = QVBoxLayout(self)

        # Compute r_h up front (needed both to display it and, if no
        # v_dust override was supplied, to seed v_dust from the r^-0.5 law).
        try:
            r_C, _ = cta.elem_to_state(comet_el, obs_jd)
            self._r_h = float(cta.vmag(r_C))
        except Exception:
            self._r_h = 1.0

        formula = QLabel(
            "<b>Afρ_norm(1AU) = Afρ_obs · r_h²</b><br>"
            "<b>Q_d ≈ Afρ_norm · v_dust / (2·p_v)</b>&nbsp;&nbsp;[kg/s]<br>"
            "<span style='color:#7090b0; font-size:10px;'>"
            "Order-of-magnitude estimate only (A'Hearn et al. 1984 Afρ "
            "convention) — not a substitute for full coma photometry/Mie "
            "modeling.</span>")
        formula.setWordWrap(True)
        formula.setTextFormat(Qt.TextFormat.RichText)
        formula.setStyleSheet(
            "padding:10px; background:rgba(255,255,255,12); border-radius:4px;")
        v.addWidget(formula)

        form = QFormLayout()

        # Date — editable; r_h is re-derived from it via Horizons/orbital
        # elements below and is NOT itself editable.
        self.ed_date = QLineEdit(cta.jd_to_str(obs_jd)[:10])
        self.ed_date.editingFinished.connect(self._refresh_rh)
        form.addRow("Date (YYYY-MM-DD):", self.ed_date)

        self.lbl_rh = QLabel(f"{self._r_h:.4f} AU")
        self.lbl_rh.setStyleSheet("color:#a0c0e0;")
        self.lbl_rh.setToolTip(
            "Heliocentric distance for the date above, propagated from "
            "the comet's orbital elements (Horizons) — read-only.")
        form.addRow("r_h  (from Horizons):", self.lbl_rh)

        self.sp_afrho = QDoubleSpinBox()
        self.sp_afrho.setRange(0.0, 1_000_000.0)
        self.sp_afrho.setDecimals(1)
        self.sp_afrho.setValue(100.0)
        self.sp_afrho.setSuffix("  cm")
        self.sp_afrho.valueChanged.connect(self._recompute)
        form.addRow("Afρ (observed):", self.sp_afrho)

        v_dust_seed = (v_dust_default if v_dust_default is not None
                      else round(0.1 * self._r_h ** -0.5, 4))
        self.sp_vdust = QDoubleSpinBox()
        self.sp_vdust.setRange(0.0001, 10.0)
        self.sp_vdust.setDecimals(4)
        self.sp_vdust.setValue(v_dust_seed)
        self.sp_vdust.setSuffix("  km/s")
        self.sp_vdust.valueChanged.connect(self._recompute)
        vdust_row = QHBoxLayout()
        vdust_row.addWidget(self.sp_vdust)
        self.btn_auto_v = QPushButton("Reset to r_h-based default")
        self.btn_auto_v.clicked.connect(
            lambda: self.sp_vdust.setValue(round(0.1 * self._r_h ** -0.5, 4)))
        vdust_row.addWidget(self.btn_auto_v)
        form.addRow("Dust velocity v_dust:", vdust_row)
        form.addRow(QLabel(
            "<i>default from 0.1·r_h⁻⁰·⁵ km/s (empirical, no fixed "
            "reference) — override with a measured value if available, "
            "e.g. 29P quiescent ~0.01-0.05, outburst ~0.1-0.3 km/s "
            "(Schambeau et al. 2017, 2019)</i>"))

        self.sp_pv = QDoubleSpinBox()
        self.sp_pv.setRange(0.01, 1.0)
        self.sp_pv.setDecimals(3)
        self.sp_pv.setValue(p_v_default)
        self.sp_pv.valueChanged.connect(self._recompute)
        form.addRow("Geometric albedo p_v:", self.sp_pv)
        v.addLayout(form)

        self.lbl_result = QLabel()
        self.lbl_result.setStyleSheet(
            "font-size:14px; font-weight:bold; padding:10px; "
            "background:rgba(60,120,200,30); border-radius:4px;")
        self.lbl_result.setWordWrap(True)
        v.addWidget(self.lbl_result)

        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        btns.rejected.connect(self.accept)
        v.addWidget(btns)

        self._recompute()

    def _refresh_rh(self):
        try:
            jd = cta.date_to_jd(self.ed_date.text().strip())
            r_C, _ = cta.elem_to_state(self.comet_el, jd)
            self._r_h = float(cta.vmag(r_C))
            self.lbl_rh.setText(f"{self._r_h:.4f} AU")
        except Exception as e:
            self._r_h = None
            self.lbl_rh.setText("— (bad date / elements)")
        self._recompute()

    def _recompute(self):
        if self._r_h is None or self._r_h <= 0:
            self.lbl_result.setText("Enter a valid date to get r_h first.")
            return
        qd = cta.compute_dust_production_rate(
            self.sp_afrho.value(), self._r_h,
            self.sp_vdust.value(), self.sp_pv.value())
        self.lbl_result.setText(
            f"Afρ_norm(1AU) = {qd['afrho_norm']:.0f} cm\n"
            f"Q_d ≈ {qd['Qd_rough']:.1e} kg/s")


# ─────────────────────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("☄ Comet Tail Analyzer  —  Finson–Probstein Model  ·  v3.1.1")
        # Keep a practical lower bound for compact laptop displays, while
        # opening maximized at startup (see main()). The left and right
        # panels are scrollable, so a 960 px-wide workspace remains usable.
        self.setMinimumSize(960, 640)
        self.resize(1366, 860)
        self._startup_layout_applied = False

        self._model    = None
        self._comet_el = None
        self._worker   = None
        # Visibility before MC auto-hides F-P curves; restored by Clear All.
        self._pre_mc_fp_visibility = None
        # Incremented whenever results are cleared.  Async callbacks captured
        # before the increment are ignored so a late worker cannot repopulate
        # a canvas that the user has explicitly cleared.
        self._model_generation = 0
        self._cobs_data     = None
        self._cobs_data_for = None # comet name self._cobs_data was fetched for
        self._pending_cobs_callback  = None

        # Update check (v3.1) — see _check_for_update_silent/_manual below.
        self._update_worker = None

        # Embedded Animator (v3.0) state — see _anim_* methods below.
        self._anim_frames    = []
        self._anim_frame_idx = 0
        self._anim_worker    = None
        self._anim_timer     = QTimer(self)
        self._anim_timer.timeout.connect(self._anim_advance_frame)

        # Static literature defaults — seed the standalone calculators
        # (Dust particle radius…, Dust production rate…) and the β-hint
        # next to the β input field. Each calculator's own inputs are
        # self-contained (v3.0) — nothing writes back here, so these stay
        # fixed at their literature values for the whole session.
        self.phys_params = {
            'rho_d':  0.5,   # g/cm³, Fulle et al. (2016)
            'v_dust': None,  # km/s, None = auto (0.1·r_h^-0.5 empirical law)
            'p_v':    0.04,  # geometric albedo, SEPPCoN/Schambeau et al. 2021
        }

        self._build_menu()
        self._build_ui()
        self._build_status()
        self._wire_animator_controls()

        # The maximized client size is only final after the window has been
        # shown. Apply the startup splitter proportions on the first event-
        # loop turn so they are derived from the user's actual display.
        QTimer.singleShot(0, self._apply_responsive_layout)

        # Silent update check (v3.1), a few seconds after the window is
        # actually usable — never blocks startup, never shows anything
        # if there's no update, no network, or GitHub rate-limits us.
        QTimer.singleShot(3000, self._check_for_update_silent)

    # ── Menu ──────────────────────────────────────────────────────────────
    def _build_menu(self):
        mb = self.menuBar()

        # File — opening, saving, exporting, and application lifecycle only.
        file_m = mb.addMenu("File")

        act_image_setup = QAction("Open / Setup Image…", self)
        act_image_setup.setShortcut("Ctrl+O")
        act_image_setup.setToolTip(
            "Open the Image Setup & Calibration dialog to load or configure "
            "the observed image overlay.")
        act_image_setup.triggered.connect(lambda: self.ctrl._open_image_dialog())
        file_m.addAction(act_image_setup)
        file_m.addSeparator()

        act_save_png = QAction("Save plot as PNG…", self)
        act_save_png.setShortcut("Ctrl+S")
        act_save_png.triggered.connect(self._save_png)
        act_export_csv = QAction("Export data as CSV…", self)
        act_export_csv.setShortcut("Ctrl+E")
        act_export_csv.triggered.connect(self._export_csv)
        file_m.addAction(act_save_png)
        file_m.addAction(act_export_csv)
        file_m.addSeparator()

        act_quit = QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        file_m.addAction(act_quit)

        # Model — every action that creates, refines, or clears a model result.
        model_m = mb.addMenu("Model")

        self._act_run_fp = QAction("▶  Run Finson–Probstein Model", self)
        self._act_run_fp.setShortcut("F5")
        self._act_run_fp.setToolTip(
            "Run the Finson–Probstein model using the current inputs "
            "(same action as COMPUTE MODEL).")
        self._act_run_fp.setStatusTip(
            "Run the Finson–Probstein model with the current inputs (F5).")
        self._act_run_fp.triggered.connect(self._run_fp_model_from_menu)
        model_m.addAction(self._act_run_fp)

        self._act_open_mc = QAction("🎲  Open Monte Carlo Model…", self)
        self._act_open_mc.setShortcut("F6")
        self._act_open_mc.setToolTip(
            "Open the Monte Carlo setup window. F6 opens the inputs; it does "
            "not start a long particle simulation automatically.")
        self._act_open_mc.setStatusTip("Open Monte Carlo morphology setup (F6).")
        self._act_open_mc.triggered.connect(self._open_mc_window)
        model_m.addAction(self._act_open_mc)

        self._act_reextract_mc = QAction("↻  Re-extract MC Contours", self)
        self._act_reextract_mc.setToolTip(
            "Rebuild contour paths from the cached MC density grid without "
            "resampling or propagating particles.")
        self._act_reextract_mc.setEnabled(False)
        self._act_reextract_mc.triggered.connect(self._reextract_mc_from_menu)
        model_m.addAction(self._act_reextract_mc)

        model_m.addSeparator()
        self._act_clear_mc = QAction("✕  Clear All Models", self)
        self._act_clear_mc.setToolTip(
            "Remove F-P and MC results, restore the zero-ejection F-P baseline, "
            "and preserve the loaded image and editable MC inputs.")
        self._act_clear_mc.setEnabled(False)
        self._act_clear_mc.triggered.connect(self._confirm_clear_all_models)
        model_m.addAction(self._act_clear_mc)

        # Compact main-window model toolbar for discoverability.
        self._model_toolbar = QToolBar("Model", self)
        self._model_toolbar.setObjectName("ModelToolbar")
        self._model_toolbar.setMovable(False)
        self._model_toolbar.setFloatable(False)
        self._model_toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self._model_toolbar.addAction(self._act_run_fp)
        self._model_toolbar.addAction(self._act_open_mc)
        self._model_toolbar.addSeparator()
        self._model_toolbar.addAction(self._act_clear_mc)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self._model_toolbar)

        # View — visualization windows and appearance only.
        view_m = mb.addMenu("View")
        act_reset = QAction("Reset zoom", self)
        act_reset.triggered.connect(lambda: self.canvas.canvas.toolbar.home())
        view_m.addAction(act_reset)
        view_m.addSeparator()

        act_orbit_view = QAction("🪐  Orbit position diagram…", self)
        act_orbit_view.setShortcut("Ctrl+Shift+O")
        act_orbit_view.triggered.connect(self._open_orbit_window)
        view_m.addAction(act_orbit_view)

        act_lc_view = QAction("📈  Light curve…", self)
        act_lc_view.setShortcut("Ctrl+L")
        act_lc_view.setToolTip(
            "Fetches COBS photometry automatically if not already fetched "
            "for the selected comet.")
        act_lc_view.triggered.connect(self._menu_open_lc)
        view_m.addAction(act_lc_view)
        view_m.addSeparator()

        theme_m = view_m.addMenu("🎨  Theme")
        self._act_dark  = QAction("☾  Dark  (Space)",  self, checkable=True)
        self._act_light = QAction("☀  Light (Observatory)", self, checkable=True)
        self._act_dark.setChecked(True)
        self._act_dark.triggered.connect(lambda: self._set_theme("dark"))
        self._act_light.triggered.connect(lambda: self._set_theme("light"))
        theme_m.addAction(self._act_dark)
        theme_m.addAction(self._act_light)

        # Calculation — standalone physical calculators.
        calc_m = mb.addMenu("Calculation")
        act_radius = QAction("🪨  Dust particle radius…", self)
        act_radius.setToolTip(
            "β → grain radius calculator. Needs only β, ρ, Qpr — no "
            "comet/model required.")
        act_radius.triggered.connect(self._open_grain_radius_dialog)
        calc_m.addAction(act_radius)

        act_qd = QAction("💨  Dust production rate…", self)
        act_qd.setToolTip(
            "Afρ → Q_d calculator. r_h comes from the selected comet's "
            "Horizons elements for a chosen date.")
        act_qd.triggered.connect(self._open_dust_production_dialog)
        calc_m.addAction(act_qd)

        # Help
        help_m = mb.addMenu("Help")
        act_about = QAction("About…", self)
        act_about.triggered.connect(self._about)
        help_m.addAction(act_about)
        help_m.addSeparator()
        act_update = QAction("🔔  Check for Updates…", self)
        act_update.triggered.connect(self._check_for_update_manual)
        help_m.addAction(act_update)

    # ── UI ────────────────────────────────────────────────────────────────
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        hbox = QHBoxLayout(central)
        hbox.setContentsMargins(0, 0, 0, 0)
        hbox.setSpacing(0)

        # Splitter: left | center | right
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(2)
        self.splitter = splitter

        # Left
        self.ctrl = ControlPanel()
        self.ctrl.compute_requested.connect(self._on_compute)
        self.ctrl.fetch_requested.connect(self._on_fetch)
        self.ctrl.image_loaded.connect(self._on_image_loaded)
        self.ctrl.comet_ready.connect(self._on_comet_ready)
        splitter.addWidget(self.ctrl)

        # Center (canvas)
        self.canvas = PlotCanvas()
        self.canvas.nucleus_clicked.connect(self._on_nucleus_clicked)
        self.ctrl.btn_pick_nuc.clicked.connect(
            lambda: self.canvas.set_picking_nucleus(True))
        splitter.addWidget(self.canvas)

        # Right
        self.info = InfoPanel()
        splitter.addWidget(self.info)

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 0)
        self.info.setMinimumWidth(240)
        # Temporary pre-show sizes. _apply_responsive_layout() replaces
        # these once the maximized client width is known.
        splitter.setSizes([320, 780, 260])

        hbox.addWidget(splitter)

    def _apply_responsive_layout(self):
        """Set sensible startup splitter widths for the current display.

        This runs once after the maximized window is shown. It intentionally
        does not keep forcing proportions during later resizes, so users can
        drag splitter handles and retain their preferred workspace.
        """
        if self._startup_layout_applied or not hasattr(self, "splitter"):
            return

        total = self.centralWidget().width() if self.centralWidget() else self.width()
        if total <= 0:
            # Window geometry may not yet be committed on some platforms.
            QTimer.singleShot(50, self._apply_responsive_layout)
            return

        # Side panels stay within their widget constraints while scaling
        # gently with screen width. The plot receives all remaining space.
        left = max(320, min(400, int(total * 0.24)))
        right = max(240, min(300, int(total * 0.19)))
        handles = self.splitter.handleWidth() * 2
        center = max(360, total - left - right - handles)
        self.splitter.setSizes([left, center, right])
        self._startup_layout_applied = True

    def _run_fp_model_from_menu(self):
        """Run F-P through the exact same path as the sidebar button."""
        if not hasattr(self, "ctrl"):
            return
        if not self.ctrl.btn_compute.isEnabled():
            self.status.showMessage(
                "Finson–Probstein computation is already running…", 3000)
            return
        self.ctrl.btn_compute.click()

    def _build_status(self):
        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self.status_comet  = QLabel("No comet selected")
        self.status_helio  = QLabel("")
        self.status_geo    = QLabel("")
        self.status_phase  = QLabel("")
        self.status_mode   = QLabel("Standalone mode")
        for w in [self.status_comet, self.status_helio, self.status_geo,
                  self.status_phase, self.status_mode]:
            w.setStyleSheet("padding: 0 12px; color:#2a4060;")
            self.status.addPermanentWidget(w)

    def _wire_animator_controls(self):
        """Connect the embedded ANIMATOR group box (InfoPanel, INFO tab,
        below ORBITAL ELEMENTS) to the _anim_* orchestration methods.
        Called once after _build_ui() so self.info already exists."""
        self.info.anim_btn_auto.clicked.connect(self._anim_auto_suggest)
        self.info.anim_btn_compute.clicked.connect(self._anim_compute_frames)
        self.info.anim_btn_play.clicked.connect(self._anim_toggle_play)
        self.info.anim_slider.valueChanged.connect(self._anim_show_frame)


    # ── v3.0: comet lookup + auto-chain helpers ─────────────────────────
    def _active_comet_el(self) -> dict | None:
        """The comet currently in play. ControlPanel's _comet_el is the
        freshest signal — it updates the instant "Use this comet" is
        clicked. self._comet_el (this window's own) only updates lazily
        after a Compute Model run finishes, so if the user switches to a
        different comet WITHOUT recomputing, self._comet_el would still
        be holding the *previous* comet — which is exactly the bug where
        Light curve…/Orbit… kept showing the old comet after switching.
        self._comet_el is kept only as a last-resort fallback (e.g. if
        ControlPanel somehow has none), never as the preferred source."""
        return getattr(self.ctrl, "_comet_el", None) or self._comet_el

    def _active_comet_name(self) -> str:
        el = self._active_comet_el()
        return el.get("name", "") if el else ""

    def _on_comet_ready(self, el: dict):
        """Fired by ControlPanel.comet_ready whenever a comet is freshly
        selected/fetched (preset, FETCH JPL, or manual entry) — re-links
        the embedded Animator to it: Start = comet's obs date, End =
        Start+360d. Any frames already computed belonged to the PREVIOUS
        comet, so they're invalidated (play/slider disabled) rather than
        left showing stale data — Compute frames must be run again."""
        obs_jd = el.get("obs_jd")
        if obs_jd is not None:
            self.ctrl.obs_date.setText(cta.jd_to_str(obs_jd)[:10])
            self.info.anim_start.setText(cta.jd_to_str(obs_jd)[:10])
            self.info.anim_end.setText(cta.jd_to_str(obs_jd + 360)[:10])

        self._anim_timer.stop()
        self._anim_frames    = []
        self._anim_frame_idx = 0
        self.info.anim_btn_play.setText("▶")
        self.info.anim_btn_play.setEnabled(False)
        self.info.anim_slider.setEnabled(False)
        self.info.anim_slider.setRange(0, 0)
        self.info.anim_lbl_frame.setText(
            f"Linked to {el.get('name','this comet')} — "
            f"click COMPUTE FRAMES to (re)build the animation.")

        # BUG FIX (v3.1): an MC contour computed for the PREVIOUS comet was
        # being left on the canvas (self.canvas.mc_contours/mc_info) after
        # switching to a new one — same "stale data from before" issue as
        # the Animator frames above, just not yet handled for this piece of
        # state. Clear it here too, for the same reason.
        self.canvas.mc_contours = []
        self.canvas.mc_info = None
        self.ctrl._mc_contour_visible = False
        try:
            self.ctrl.btn_clear_mc.setEnabled(False)
            self._act_clear_mc.setEnabled(False)
        except Exception:
            pass
        # Clear MCWindow's remembered state so next open starts fresh
        # for the new comet — old grain size range, release window,
        # COBS anchors etc. from the previous comet should not persist.
        MCWindow._saved_state = None

    def _ensure_cobs_fetched(self, then_callback):
        """Call then_callback(cobs_data) once COBS data exists for the
        *current* comet — reusing the cache if it's still for the same
        comet, otherwise fetching first. then_callback always eventually
        fires (with cobs_data=None on failure/no comet) so callers never
        hang waiting for it."""
        name = self._active_comet_name()
        if not name:
            then_callback(None)
            return
        if self._cobs_data is not None and self._cobs_data_for == name:
            then_callback(self._cobs_data)
            return
        self._pending_cobs_callback = then_callback
        self._fetch_cobs(name)

    def _open_grain_radius_dialog(self):
        # No comet/model needed — the β→radius formula is comet-independent.
        dlg = GrainRadiusDialog(self, self.phys_params['rho_d'])
        dlg.exec()

    def _open_dust_production_dialog(self):
        el = self._active_comet_el()
        if not el:
            QMessageBox.warning(self, "No Comet",
                                "Select or fetch a comet first — r_h is "
                                "propagated from its orbital elements.")
            return
        obs_jd = self._current_obs_jd()
        dlg = DustProductionDialog(
            self, el, obs_jd,
            self.phys_params['v_dust'], self.phys_params['p_v'])
        dlg.exec()

    def _menu_open_lc(self):
        if not self._active_comet_name():
            QMessageBox.warning(self, "No Comet",
                                "Select or fetch a comet first.")
            return

        # Light Curve is an auto-chain action: fetch COBS first when the
        # current comet has no cached data, then open the window only after
        # the worker succeeds.  On a real fetch error, _fetch_cobs() already
        # reports the error; do not continue into _open_lc_window(), which
        # would otherwise show the misleading "Fetch COBS light curve first"
        # message even though an automatic fetch was just attempted.
        self._ensure_cobs_fetched(
            lambda d: self._open_lc_window() if d is not None else None)

    def _anim_parse_range(self):
        try:
            start_jd = cta.date_to_jd(self.info.anim_start.text().strip())
            end_jd   = cta.date_to_jd(self.info.anim_end.text().strip())
        except Exception:
            QMessageBox.warning(self, "Bad Date", "Enter dates as YYYY-MM-DD.")
            return None, None
        if end_jd <= start_jd:
            QMessageBox.warning(self, "Bad Range",
                                "End date must be after start date.")
            return None, None
        return start_jd, end_jd

    def _anim_gather_inputs(self):
        """(comet_el, betas, ages, max_age, n_pts, ejection) for the CURRENT
        comet/model-panel settings, or None (with a warning already shown)
        if anything's missing — gathered fresh on every Compute/Auto-suggest
        click rather than once up front, since unlike the old popup this
        panel stays open across comet/setting changes. ejection is the
        v3.1 ejection-velocity dict (all-zero by default == v3.0 model)."""
        el = self._active_comet_el()
        if not el:
            QMessageBox.warning(self, "No Comet",
                                "Select or fetch a comet first.")
            return None
        betas = self.ctrl._parse_floats(self.ctrl.beta_str.text())
        ages  = self.ctrl._parse_ints(self.ctrl.age_str.text())
        if not betas or not ages:
            QMessageBox.warning(self, "Bad Input",
                "Enter at least one valid β value and synchrone age in "
                "the model panel first.")
            return None
        return (el, betas, ages, self.ctrl.max_age.value(), self.ctrl.n_pts.value(),
                self.ctrl.get_ejection_params())

    def _anim_auto_suggest(self):
        start_jd, end_jd = self._anim_parse_range()
        if start_jd is None:
            return
        inputs = self._anim_gather_inputs()
        if inputs is None:
            return
        el, betas, ages, max_age, n_pts, ejection = inputs
        max_au = max_deg = 0.0
        for jd in np.linspace(start_jd, end_jd, 5):
            try:
                m = cta.compute_model(el, float(jd), betas, ages, max_age,
                                      min(n_pts, 60), ejection=ejection)
            except Exception:
                continue
            ex, ey = _model_extent_au(m)
            ext   = max(ex, ey)
            r_geo = m["info"].get("r_geo", 1.0)
            K     = (180.0 / np.pi) / r_geo
            max_au  = max(max_au, ext)
            max_deg = max(max_deg, ext * K)
        if max_au <= 0:
            QMessageBox.information(self, "Auto-suggest",
                "Couldn't sample the model — check the comet/date range.")
            return
        margin = 1.2
        if self.info.anim_rb_fov.isChecked():
            self.info.anim_size.setValue(round(max_deg * 60 * 2 * margin, 2))
        else:
            self.info.anim_size.setValue(round(max_au * 2 * margin, 4))

    def _anim_sync_fov_to_overlay(self, el: dict):
        """If an image is currently loaded for overlay, set the Animator's
        Fixed FOV to that image's own angular width (using its au_per_px
        scale and Δ at the current working date) — so animation frames
        are shown at the same scale as the overlay was built at, instead
        of whatever Fixed FOV happened to be left over from a previous
        comet/session. No-op if no image is loaded."""
        img = self.ctrl._img_arr
        if img is None:
            return
        au_per_px = self.ctrl.au_px_spin.value()
        if au_per_px <= 0:
            return
        w_img = img.shape[1]
        width_au = w_img * au_per_px
        try:
            jd = self._current_obs_jd()
            r_C, _ = cta.elem_to_state(el, jd)
            r_E    = cta.earth_pos(jd)
            r_geo  = float(cta.vmag(r_C - r_E))
        except Exception:
            return
        if r_geo <= 0:
            return
        K = (180.0 / np.pi) / r_geo   # deg per AU
        fov_arcmin = width_au * K * 60.0
        self.info.anim_rb_fov.setChecked(True)
        self.info.anim_size.setValue(round(fov_arcmin, 3))
        self.status.showMessage(
            f"Animator FOV synced to overlay image: {fov_arcmin:.1f} arcmin "
            f"({w_img}px × {au_per_px:.6f} AU/px at Δ={r_geo:.3f} AU)", 6000)

    def _anim_compute_frames(self):
        start_jd, end_jd = self._anim_parse_range()
        if start_jd is None:
            return
        inputs = self._anim_gather_inputs()
        if inputs is None:
            return
        el, betas, ages, max_age, n_pts, ejection = inputs
        self._anim_sync_fov_to_overlay(el)
        step = self.info.anim_step.value()
        n_frames = int(round((end_jd - start_jd) / step)) + 1
        if n_frames < 2:
            QMessageBox.warning(self, "Bad Range",
                "Range/step combination produces fewer than 2 frames.")
            return
        if n_frames > 300:
            QMessageBox.warning(self, "Too Many Frames",
                f"{n_frames} frames requested (max 300) — increase the "
                f"step or shorten the date range.")
            return

        self._anim_timer.stop()
        self.info.anim_btn_play.setText("▶")
        self._anim_frames = []
        self.info.anim_btn_compute.setEnabled(False)
        self.info.anim_progress.setMaximum(n_frames)
        self.info.anim_progress.setValue(0)
        self.info.anim_lbl_frame.setText(f"Computing {n_frames} frames…")

        class AnimWorker(QThread):
            frame_done   = pyqtSignal(dict)
            progress     = pyqtSignal(int)
            finished_all = pyqtSignal()
            error        = pyqtSignal(str)

            def __init__(self, el, betas, ages, max_age, n_pts,
                        start_jd, step, n_frames, ejection=None):
                super().__init__()
                self.el, self.betas, self.ages = el, betas, ages
                self.max_age, self.n_pts = max_age, n_pts
                self.start_jd, self.step, self.n_frames = start_jd, step, n_frames
                self.ejection = ejection

            def run(self):
                try:
                    for i in range(self.n_frames):
                        jd = self.start_jd + i * self.step
                        m  = cta.compute_model(self.el, jd, self.betas,
                                               self.ages, self.max_age,
                                               self.n_pts, ejection=self.ejection)
                        # compute_model() doesn't set info['name'] itself
                        # (only _on_model_ready does, for the regular
                        # Compute Model flow) — without this, the title
                        # falls back to the literal default "Comet".
                        m["info"]["name"] = self.el.get("name", "Comet")
                        self.frame_done.emit(m)
                        self.progress.emit(i + 1)
                    self.finished_all.emit()
                except Exception as e:
                    self.error.emit(str(e))

        self._anim_worker = AnimWorker(el, betas, ages, max_age, n_pts,
                                       start_jd, step, n_frames, ejection=ejection)
        self._anim_worker.frame_done.connect(self._anim_frames.append)
        self._anim_worker.progress.connect(self.info.anim_progress.setValue)
        self._anim_worker.finished_all.connect(self._anim_on_compute_done)
        self._anim_worker.error.connect(self._anim_on_compute_error)
        self._anim_worker.start()

    def _anim_on_compute_done(self):
        self.info.anim_btn_compute.setEnabled(True)
        if not self._anim_frames:
            QMessageBox.warning(self, "No Frames", "No frames were computed.")
            return
        self.info.anim_slider.setEnabled(True)
        self.info.anim_slider.setRange(0, len(self._anim_frames) - 1)
        self.info.anim_btn_play.setEnabled(True)
        self.info.anim_slider.setValue(0)
        self._anim_show_frame(0)

    def _anim_on_compute_error(self, msg):
        self.info.anim_btn_compute.setEnabled(True)
        self.info.anim_lbl_frame.setText("Compute failed.")
        QMessageBox.warning(self, "Compute Error", msg)

    def _anim_frame_limits(self, model):
        r_geo = model["info"].get("r_geo", 1.0)
        K = (180.0 / np.pi) / r_geo
        if self.info.anim_rb_fov.isChecked():
            half_deg = (self.info.anim_size.value() / 60.0) / 2.0
        else:
            half_deg = (self.info.anim_size.value() / 2.0) * K
        # East LEFT = inverted x-axis, matching PlotCanvas's own convention
        return (half_deg, -half_deg), (-half_deg, half_deg)

    def _anim_show_frame(self, idx: int):
        if not (0 <= idx < len(self._anim_frames)):
            return
        self._anim_frame_idx = idx
        model = self._anim_frames[idx]
        xlim, ylim = self._anim_frame_limits(model)
        self.canvas.draw_model(model, img_arr=None,
                               fixed_xlim=xlim, fixed_ylim=ylim)
        info = model["info"]
        # v3.0: keep obs_date AND the EPHEMERIS display in lockstep with
        # whatever frame is showing — consistent with the Orbit view link
        # above (orbital elements themselves don't change with date, so
        # only the ephemeris half of update_info() actually changes here).
        self.ctrl.obs_date.setText(info["obs_str"][:10])
        el = self._active_comet_el()
        if el:
            self.info.update_info(model, el)
        mode = "FOV" if self.info.anim_rb_fov.isChecked() else "dist"
        self.info.anim_lbl_frame.setText(
            f"Frame {idx+1}/{len(self._anim_frames)} · {info['obs_str'][:10]} · "
            f"r☉={info['r_helio']:.3f} Δ={info['r_geo']:.3f} · Fixed {mode}")

    def _anim_toggle_play(self):
        if self._anim_timer.isActive():
            self._anim_timer.stop()
            self.info.anim_btn_play.setText("▶")
        else:
            self._anim_timer.start(200)   # ms/frame
            self.info.anim_btn_play.setText("⏸")

    def _anim_advance_frame(self):
        nxt = self._anim_frame_idx + 1
        if nxt >= len(self._anim_frames):
            nxt = 0
        self.info.anim_slider.setValue(nxt)   # triggers _anim_show_frame via signal

    # ── LC popup window ────────────────────────────────────────────────
    def _open_lc_window(self):
        """Open a larger standalone light curve window.

        A fitted H₀/n pair is optional.  The COBS viewer can still plot the
        raw date/magnitude observations when ephemeris enrichment or the
        H₀/n fit is unavailable.  Requiring H₀ here caused a successful COBS
        download with raw observations to be rejected as "No Data".
        """
        d = getattr(self, "_cobs_data", None)
        has_raw = bool(d and d.get("raw_obs"))
        has_eph = bool(d and d.get("obs_list"))
        has_fit = bool(d and d.get("H0") is not None and d.get("n") is not None)
        if not (has_raw or has_eph or has_fit):
            QMessageBox.information(
                self, "No COBS Data",
                "No COBS light-curve observations are available for the "
                "selected comet.")
            return
        try:
            # Compute today's actual r_helio and delta from orbital elements.
            # info.r_helio / r_geo are for the *model obs date*, not today,
            # so using them gives a wrong "Now" magnitude on the light curve.
            import numpy as np
            r_now = None; delta_now = None
            if getattr(self, '_comet_el', None):
                today      = cta.today_jd()
                r_C, _     = cta.elem_to_state(self._comet_el, today)
                r_E        = cta.earth_pos(today)
                r_now      = float(np.linalg.norm(r_C))
                delta_now  = float(np.linalg.norm(r_C - r_E))
            else:
                info       = self._model["info"] if self._model else {}
                r_now      = info.get("r_helio")
                delta_now  = info.get("r_geo")
            # obs_jd: user's observation date; T_jd: perihelion from orbital elements
            obs_jd = getattr(self._comet_el, 'get', lambda k, v=None: None)('obs_jd') \
                     if self._comet_el else None
            if obs_jd is None and self._comet_el:
                obs_jd = self._comet_el.get('obs_jd')
            T_jd = self._comet_el.get('T_jd') if self._comet_el else None
            # Store reference — prevents garbage collection crash
            self._lc_win = LCWindow(d, r_now, delta_now,
                                    obs_jd=obs_jd, T_jd=T_jd, parent=None)
            self._lc_win.show()
            self._lc_win.raise_()
            self._lc_win.activateWindow()
        except Exception as e:
            QMessageBox.warning(self, "Plot Error", str(e))

    # ── COBS light curve ───────────────────────────────────────────────
    def _fetch_cobs(self, name: str):
        if not name:
            self.status.showMessage("No comet loaded.", 2000); return
        self.status.showMessage(f"Fetching COBS data for {name}…")

        class COBSWorker(QThread):
            done  = pyqtSignal(dict)
            error = pyqtSignal(str)
            def __init__(self, n, el): super().__init__(); self.n=n; self.el=el
            def run(self):
                try: self.done.emit(cta.fetch_from_cobs(self.n, self.el))
                except Exception as e: self.error.emit(str(e))

        self._cobs_worker = COBSWorker(name, self._active_comet_el() or {})
        def _on_done(d):
            self._cobs_data     = d
            self._cobs_data_for = name
            H0 = d.get("H0"); n_val = d.get("n")
            src = d.get("source","")
            msg = f"H₀={H0:.2f}  n={n_val:.2f}" if H0 else "COBS: no fit"
            self.status.showMessage(msg + f"  ({src[:35]})", 8000)
            cb, self._pending_cobs_callback = self._pending_cobs_callback, None
            if cb: cb(d)
        def _on_err(msg):
            logging.warning("COBS fetch failed for %s: %s", name, msg)
            short_msg = msg.replace("\n", " ")
            if len(short_msg) > 500:
                short_msg = short_msg[:497] + "…"
            self.status.showMessage(f"COBS fetch failed: {short_msg}", 12000)

            # A windowed PyInstaller build has no visible console.  Show the
            # actual network/API error instead of allowing the Light Curve
            # action to fall through to a generic "fetch first" message.
            QMessageBox.warning(
                self, "COBS Fetch Failed",
                f"Could not retrieve the COBS light curve for {name}.\n\n"
                f"{short_msg}")

            cb, self._pending_cobs_callback = self._pending_cobs_callback, None
            if cb: cb(None)
        self._cobs_worker.done.connect(_on_done)
        self._cobs_worker.error.connect(_on_err)
        self._cobs_worker.start()

    def _set_theme(self, theme_name: str):
        """Switch between dark and light themes."""
        app = QApplication.instance()
        apply_theme(theme_name, app)

        # Update matplotlib canvas colors
        is_dark = (theme_name == "dark")
        bg = T["mpl_bg"]
        self.canvas.fig.set_facecolor(bg)
        self.canvas.ax.set_facecolor(bg)
        for sp in self.canvas.ax.spines.values():
            sp.set_edgecolor(T["border"])
        self.canvas.ax.tick_params(labelcolor=T["mpl_tick"], color=T["border"])
        self.canvas.ax.xaxis.label.set_color(T["mpl_label"])
        self.canvas.ax.yaxis.label.set_color(T["mpl_label"])
        self.canvas.ax.title.set_color(T["mpl_title"])
        self.canvas.canvas.figure.set_facecolor(bg)
        self.canvas.canvas.draw_idle()

        # Update toolbar color
        tb = self.canvas.toolbar
        tb_style = f"background:{T['panel_bg']}; border-bottom:1px solid {T['border']};"
        tb.setStyleSheet(tb_style)

        # Update button icon and menu checkmarks
        if hasattr(self.ctrl, 'btn_theme'):
            self.ctrl.btn_theme.setText("☾" if is_dark else "☀")
            self.ctrl.btn_theme.setToolTip(
                "Switch to Light mode" if is_dark else "Switch to Dark mode")
        if hasattr(self, '_act_dark'):
            self._act_dark.setChecked(is_dark)
            self._act_light.setChecked(not is_dark)

        # Re-apply all hardcoded styles in the right info panel
        self.info.update_theme()

        # Redraw model if present
        if self._model:
            ov = self.ctrl.get_overlay()
            self.canvas.draw_model(
                self._model, ov["img_arr"] if self.ctrl._img_arr is not None else None)

        self.status.showMessage(
            f"Theme: {'Dark ☾' if is_dark else 'Light ☀'}", 2500)

    def recompute_and_then(self, callback):
        """
        v3.1 bug fix — MCWindow's "Run Monte Carlo" syncs ejection params
        then needs the F-P model recomputed with them BEFORE drawing the
        MC contour alongside it. The obvious way — call ctrl._emit_compute()
        then immediately draw — is WRONG: _emit_compute() only *starts* a
        background ComputeWorker (see _on_compute()) and returns right
        away, so the immediate draw used whatever self._model was left
        over from a PREVIOUS run (e.g. an earlier nonzero-v0 test) — not
        the new one. This produced a real, reproducible symptom: setting
        v0 back to 0 and re-running still showed the F-P curves diverging
        from the MC contour, because the visible F-P curves weren't
        actually at v0=0 yet, just stale.
        This method does the SAME recompute but calls `callback()` only
        once the new model has actually arrived (chained onto the same
        ComputeWorker.finished signal _on_model_ready() itself uses), so
        a caller can safely redraw using the new self._model right after.
        """
        betas = self.ctrl._parse_floats(self.ctrl.beta_str.text())
        ages  = self.ctrl._parse_ints(self.ctrl.age_str.text())
        if not betas or not ages or self.ctrl._comet_el is None:
            callback()   # nothing sensible to recompute; don't block the caller
            return

        comet_el = self.ctrl._comet_el
        obs_jd = comet_el.get("obs_jd", cta.today_jd())
        ov = self.ctrl.get_overlay()
        ejection = self.ctrl.get_ejection_params()

        self._comet_el = comet_el
        self.canvas._vis = self.ctrl.get_vis()
        self.canvas.nuc_x   = ov["nuc_x"]
        self.canvas.nuc_y   = ov["nuc_y"]
        self.canvas.au_per_px = ov["au_per_px"]
        self.canvas.north_pa  = ov["north_pa"]
        comet_el["obs_jd"] = obs_jd

        self.ctrl.set_computing(True)
        self._act_run_fp.setEnabled(False)
        self.status.showMessage("Computing Finson–Probstein model…")

        self._model_generation += 1
        generation = self._model_generation
        worker = ComputeWorker(comet_el, obs_jd, betas, ages,
                               self.ctrl.max_age.value(), self.ctrl.n_pts.value(),
                               ejection=ejection)
        self._worker = worker   # keep a reference, same as _on_compute() does

        def _done(model, expected=generation):
            if expected != self._model_generation:
                return
            self._on_model_ready(model, betas, ov)
            callback()

        def _progress(v, msg, expected=generation):
            if expected == self._model_generation:
                self.ctrl.progress_bar.setValue(v)
                self.status.showMessage(msg)

        def _error(msg, expected=generation):
            if expected == self._model_generation:
                self._on_compute_error(msg)

        worker.progress.connect(_progress)
        worker.finished.connect(_done)
        worker.error.connect(_error)
        worker.start()

    def set_mc_contours(self, paths, info=None):
        """
        v3.1, Phase 2 — receives contour paths (from extract_morphology_
        contours()) pushed by the Monte Carlo window's Run button, via an
        explicit self.main_window reference held by MCWindow (that window
        uses parent=None for the same reason OrbitWindow/LCWindow do — see
        their comments — so this can't go through Qt parent/child
        signaling). Redraws immediately using the already-cached F-P
        model/image, same "redraw without recompute" pattern as
        _set_theme() above. info (the MC result's info dict) is stored
        for the optional on-canvas parameter box — see draw_model()'s
        mc_info_box handling.

        BUG FIX (v3.1): _on_model_ready() rotates the F-P model's
        syndynes/synchrones/sun_dir/antivel_dir by the "Grid rotation
        (match observed tail)" slider, but extract_morphology_contours()
        deliberately returns PURE, unrotated sky-plane coordinates (see
        its docstring) so it can go through the main canvas's own to_px()
        — which only applies north_pa, not this separate grid-rotation
        offset. Net effect: with a nonzero Grid Rotation, the F-P curves
        rotated correctly but the MC contour silently stayed at its
        original, un-rotated orientation — exactly the kind of "model
        points one way, should point another" mismatch this was caught
        from. Apply the SAME rotation here, the same formula as _on_
        model_ready(), so both stay locked together. Like self._model's
        own rotation, this is baked in at the time this method runs, not
        live — if Grid Rotation is changed afterward without re-running
        Monte Carlo, re-run it to pick up the new angle (same convention
        as the F-P model itself).
        """
        rot_deg = self.ctrl.get_rotation_offset()
        if abs(rot_deg) > 1e-3 and paths:
            theta = np.radians(-rot_deg)
            cos_t, sin_t = np.cos(theta), np.sin(theta)
            paths = [
                np.column_stack([p[:, 0]*cos_t - p[:, 1]*sin_t,
                                 p[:, 0]*sin_t + p[:, 1]*cos_t])
                for p in paths
            ]
        # A Monte Carlo population contour is intended as the final morphology
        # comparison layer. Hide the discrete F-P diagnostic curves automatically
        # so the canvas is not cluttered by syndynes/synchrones. The checkboxes are
        # updated visibly; users may turn either layer back on manually afterward.
        if paths:
            try:
                if self._pre_mc_fp_visibility is None:
                    self._pre_mc_fp_visibility = (
                        self.ctrl.chk_synd.isChecked(),
                        self.ctrl.chk_sync.isChecked())
                self.ctrl.chk_synd.setChecked(False)
                self.ctrl.chk_sync.setChecked(False)
            except Exception:
                pass

        self.canvas.mc_contours = paths
        self.canvas.mc_info = info
        self.ctrl._mc_contour_visible = True
        try:
            self.ctrl.btn_clear_mc.setEnabled(True)
            self._act_clear_mc.setEnabled(True)
        except Exception:
            pass
        self.canvas._vis = self.ctrl.get_vis()
        if self._model:
            ov = self.ctrl.get_overlay()
            self.canvas.draw_model(
                self._model, ov["img_arr"] if self.ctrl._img_arr is not None else None)
        try:
            self._act_reextract_mc.setEnabled(bool(paths) and
                getattr(getattr(self, "_mc_win", None), "_mc_result", None) is not None)
        except Exception:
            pass
        self.status.showMessage(
            f"Monte Carlo contour ({len(paths)} path segment(s)) shown on main canvas.", 5000)

    def _reextract_mc_from_menu(self):
        """Re-extract contours from the open MC window's cached density grid."""
        win = getattr(self, "_mc_win", None)
        if win is None or getattr(win, "_mc_result", None) is None:
            QMessageBox.information(
                self, "Re-extract MC Contours",
                "No cached Monte Carlo density grid is available.\n\n"
                "Open the Monte Carlo window and run the model first.")
            self._act_reextract_mc.setEnabled(False)
            return
        try:
            win.show()
            win.raise_()
            win.activateWindow()
            win._reextract_contours()
        except RuntimeError:
            self._mc_win = None
            self._act_reextract_mc.setEnabled(False)

    def _confirm_clear_all_models(self):
        """Ask before clearing displayed results; saved files remain untouched."""
        has_result = bool(self._model is not None or self.canvas._model is not None or
                          self.canvas.mc_contours or self.canvas.mc_info)
        if not has_result:
            self.status.showMessage("No computed model is currently displayed.", 4000)
            return
        answer = QMessageBox.question(
            self, "Clear All Models",
            "Clear all Finson–Probstein and Monte Carlo results?\n\n"
            "This also resets the active F-P ejection velocity to zero. "
            "The loaded image, editable MC inputs, and saved .mcin files are preserved.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel)
        if answer == QMessageBox.StandardButton.Yes:
            self._clear_all_models()

    def _clear_all_models(self):
        """Clear all computed results and return F-P to a clean baseline.

        Editable MC fields remain available in the MC window/.mcin workflow,
        but hidden state that can silently alter a later F-P result is reset:
        the MC-coupled ejection velocity, contour-only presentation mode, and
        auto-hidden F-P visibility.
        """
        had_fp = bool(self._model is not None or self.canvas._model is not None)
        had_mc = bool(self.canvas.mc_contours or self.canvas.mc_info)

        # Invalidate any result callback that was already in flight when the
        # user pressed Clear.  QThread computation may finish in the background,
        # but its captured generation no longer matches and it is ignored.
        self._model_generation += 1
        try:
            if self._worker is not None and self._worker.isRunning():
                self._worker.requestInterruption()
        except Exception:
            pass

        self._model = None
        self.canvas._model = None
        self.canvas.mc_contours = []
        self.canvas.mc_info = None
        self.ctrl._mc_contour_visible = False

        # MC Run writes its velocity law into a hidden main-panel state so F-P
        # and MC are physically synchronized.  That state must not survive an
        # explicit Clear All command; otherwise the next F-P model silently
        # remains a non-zero-velocity model.
        self.ctrl.reset_ejection_params()

        # Restore the user's pre-MC F-P visibility.  If no snapshot exists,
        # use the normal diagnostic defaults.
        try:
            if self._pre_mc_fp_visibility is not None:
                show_synd, show_sync = self._pre_mc_fp_visibility
            else:
                show_synd, show_sync = True, True
            self.ctrl.chk_synd.setChecked(bool(show_synd))
            self.ctrl.chk_sync.setChecked(bool(show_sync))
        except Exception:
            pass
        self._pre_mc_fp_visibility = None

        # Publication/contour-only modes belong to an MC comparison.  Reset to
        # Analysis Overlay so the next standalone F-P computation is visible.
        self.ctrl.reset_mc_presentation()

        mc_win = getattr(self, "_mc_win", None)
        if mc_win is not None:
            try:
                mc_win._mc_result = None
                mc_win._run_generation = getattr(mc_win, "_run_generation", 0) + 1
                if mc_win._worker is not None and mc_win._worker.isRunning():
                    mc_win._worker.requestInterruption()
                mc_win.btn_reextract.setEnabled(False)
                mc_win.progress.setRange(0, 100)
                mc_win.progress.setValue(0)
                mc_win.progress.setFormat("0%")
                # Keep the open dialog consistent with the main-canvas reset.
                idx = mc_win.cmb_mc_view_style.findData("analysis")
                if idx >= 0:
                    mc_win.cmb_mc_view_style.setCurrentIndex(idx)
                mc_win.lbl_status.setText(
                    "All model results cleared; F-P reset to zero ejection velocity.")
            except Exception:
                pass

        try:
            self.ctrl.btn_clear_mc.setEnabled(False)
            self._act_clear_mc.setEnabled(False)
            self._act_reextract_mc.setEnabled(False)
            self.ctrl.set_computing(False)
            self._act_run_fp.setEnabled(True)
        except Exception:
            pass

        # Redraw observation only.  Model-derived contours and annotations are
        # absent; image-derived isophotes remain governed by the main checkbox.
        self.canvas._vis = self.ctrl.get_vis()
        if self.ctrl._img_arr is not None:
            ov = self.ctrl.get_overlay()
            self.canvas.nuc_x = ov.get("nuc_x", self.canvas.nuc_x)
            self.canvas.nuc_y = ov.get("nuc_y", self.canvas.nuc_y)
            self.canvas.north_pa = ov.get("north_pa", self.canvas.north_pa)
            self.canvas.draw_image_preview(self.ctrl._img_arr)
        else:
            self.canvas._imgArr = None
            self.canvas._draw_empty()

        if had_fp and had_mc:
            msg = "F-P and MC models cleared; F-P reset to zero ejection velocity."
        elif had_fp:
            msg = "F-P model cleared; zero-ejection baseline restored."
        elif had_mc:
            msg = "MC model cleared; zero-ejection F-P baseline restored."
        else:
            msg = "No computed model is currently displayed."
        self.status.showMessage(msg, 6000)

    # Backward-compatible internal alias for any older signal wiring.
    def _clear_mc_model(self):
        self._clear_all_models()

    def _on_compute(self, comet_el, obs_jd, betas, ages, max_age, n_pts, ejection=None):
        self._comet_el = comet_el
        vis = self.ctrl.get_vis()
        ov  = self.ctrl.get_overlay()
        # Apply vis to canvas
        self.canvas._vis = vis
        self.canvas.nuc_x   = ov["nuc_x"]
        self.canvas.nuc_y   = ov["nuc_y"]
        self.canvas.au_per_px = ov["au_per_px"]
        self.canvas.north_pa  = ov["north_pa"]

        comet_el["obs_jd"] = obs_jd

        self.ctrl.set_computing(True)
        self._act_run_fp.setEnabled(False)
        self.status.showMessage("Computing Finson–Probstein model…")

        self._model_generation += 1
        generation = self._model_generation
        self._worker = ComputeWorker(comet_el, obs_jd, betas, ages, max_age, n_pts,
                                     ejection=ejection)

        def _progress(v, msg, expected=generation):
            if expected == self._model_generation:
                self.ctrl.progress_bar.setValue(v)
                self.status.showMessage(msg)

        def _done(model, expected=generation):
            if expected == self._model_generation:
                self._on_model_ready(model, betas, ov)

        def _error(msg, expected=generation):
            if expected == self._model_generation:
                self._on_compute_error(msg)

        self._worker.progress.connect(_progress)
        self._worker.finished.connect(_done)
        self._worker.error.connect(_error)
        self._worker.start()

    def _on_model_ready(self, model, betas, ov):
        self._model = model
        # A Finson–Probstein result is now displayed, so CLEAR ALL MODELS
        # must be available even when no Monte Carlo result exists yet.
        try:
            self.ctrl.btn_clear_mc.setEnabled(True)
            self._act_clear_mc.setEnabled(True)
        except Exception:
            pass
        model["info"]["name"] = self._comet_el.get("name","Comet")
        info = model["info"]

        # ── Apply PA rotation offset to all xi/eta coordinates ───────────
        rot_deg = self.ctrl.get_rotation_offset()
        if abs(rot_deg) > 1e-3:
            # Rotation in sky plane: xi' = xi*cos(θ) - eta*sin(θ)
            #                        eta'= xi*sin(θ) + eta*cos(θ)
            # PA rotation convention: positive rot_deg rotates grid CW on sky
            # = adds to PA of each point → negative θ in standard math rotation
            θ = np.radians(-rot_deg)   # negative: CW on sky = CCW in standard math
            cos_t, sin_t = np.cos(θ), np.sin(θ)
            for s in model["syndynes"]:
                xi, eta = s["xi"].copy(), s["eta"].copy()
                s["xi"]  =  xi*cos_t - eta*sin_t
                s["eta"] =  xi*sin_t + eta*cos_t
            for s in model["synchrones"]:
                xi, eta = s["xi"].copy(), s["eta"].copy()
                s["xi"]  =  xi*cos_t - eta*sin_t
                s["eta"] =  xi*sin_t + eta*cos_t
            # Rotate sun_dir
            sxi, seta = model["sun_dir"]
            model["sun_dir"] = (sxi*cos_t - seta*sin_t,
                                sxi*sin_t + seta*cos_t)
            # Rotate antivel_dir (v3.0) — must stay locked to sun_dir's frame
            avxi, aveta = model.get("antivel_dir", (0.0, 0.0))
            model["antivel_dir"] = (avxi*cos_t - aveta*sin_t,
                                    avxi*sin_t + aveta*cos_t)
            # Rotate orbit points
            model["orbit"] = [(x*cos_t - y*sin_t, x*sin_t + y*cos_t, dt)
                               for x, y, dt in model["orbit"]]

        # ── Update rot_info label ─────────────────────────────────────────
        try:
            _sxi, _set = model["sun_dir"]
            _sl = np.sqrt(_sxi**2+_set**2)
            _PsAng_post = ((np.degrees(np.arctan2(_sxi/_sl,_set/_sl))+180.0)%360.0) if _sl>1e-10 else 0.0
            _PsAng_pre  = _PsAng_post + rot_deg   # before rotation
            self.ctrl.rot_info_lbl.setText(
                f"PsAng: {_PsAng_pre:.1f}°  →  {_PsAng_post:.1f}°  (after {rot_deg:+.1f}° grid rotation)")
        except Exception:
            pass
        obs_ovr = self.ctrl.get_obs_override()
        if obs_ovr:
            info["r_helio"]     = obs_ovr["r_helio"]
            info["r_geo"]       = obs_ovr["r_geo"]
            info["phase_angle"] = obs_ovr["phase_angle"]
            # Override Sun direction from PsAng (position angle of Sun from Horizons)
            # PsAng = PA of Sun as seen from comet (from Earth's perspective)
            # Anti-solar direction = PsAng + 180° (dust tail points away from Sun)
            # We need to convert PA to sky-plane offset (xi, eta)
            # xi = East offset, eta = North offset on sky
            # PA measured from North toward East: xi = sin(PA), eta = cos(PA)
            psang_rad = np.radians(obs_ovr["psang"])
            # Sun is in direction of PsAng from comet, so (xi_sun,eta_sun) points toward Sun
            model["sun_dir"] = (np.sin(psang_rad), np.cos(psang_rad))
            info["obs_override"] = True
            # v3.0: antivel_dir was computed from the fetched orbital elements,
            # which may not match a manually-overridden geometry. There is no
            # manual PsAMV (position-angle-of-motion-vector) input yet to keep
            # it consistent, so suppress the arrow rather than show a
            # potentially mismatched angle relative to the overridden Sun
            # direction. (A future version could add a psamv override field
            # alongside psang, mirroring JPL Horizons' own PsAMV column.)
            model["antivel_dir"] = (0.0, 0.0)

        # ── Update AU/pixel from effective plate scale ──
        # Use helper to get plate scale with correct priority (user > WCS > none)
        ps_arcsec = self._effective_plate_scale_arcsec()
        if ps_arcsec is not None and ps_arcsec > 0:
            r_geo_use = obs_ovr["r_geo"] if obs_ovr else info["r_geo"]
            # Convert arcsec/px → rad/px → AU/px
            au_per_px = ps_arcsec * (np.pi / 180.0 / 3600.0) * r_geo_use
            # Update au_px_spin ONLY if it differs significantly (avoid infinite loops)
            current_au_px = self.ctrl.au_px_spin.value()
            if abs(current_au_px - au_per_px) > 1e-9:
                self.ctrl.au_px_spin.setValue(round(au_per_px, 9))
            ov = self.ctrl.get_overlay()
            self.canvas.au_per_px = ov["au_per_px"]

        self.canvas._vis = self.ctrl.get_vis()
        img_arr = ov["img_arr"]
        self.canvas.nuc_x    = ov["nuc_x"]
        self.canvas.nuc_y    = ov["nuc_y"]
        self.canvas.au_per_px= ov["au_per_px"]
        self.canvas.north_pa = ov["north_pa"]
        if img_arr is None:
            # Standalone (no image overlay): use a 600 arcmin Fixed FOV as
            # the default view instead of the old auto-fit percentile box,
            # so Compute Model's result is the same starting frame size
            # the embedded Animator below would continue from — no jump
            # in apparent zoom between "look at the static model" and
            # "now play the animation". The Animator's own size field
            # defaults to the same 600 arcmin for the same reason; both
            # are still freely editable afterward.
            half_deg = (600.0 / 60.0) / 2.0
            self.canvas.draw_model(model, img_arr,
                                   fixed_xlim=(half_deg, -half_deg),
                                   fixed_ylim=(-half_deg, half_deg))
        else:
            # Overlay (image loaded): always fit the plot area to the
            # image's own native pixel extent — unchanged from the
            # earliest versions, never constrained by the Fixed-FOV
            # default above.
            self.canvas.draw_model(model, img_arr)
        self.info.update_info(model, self._comet_el)

        # v3.0: invalidate any COBS data cached for a *different* comet —
        # avoids silently showing a previous comet's light curve after
        # switching comets without an intervening explicit COBS re-fetch.
        _name = self._comet_el.get("name", "")
        if self._cobs_data_for is not None and self._cobs_data_for != _name:
            self._cobs_data     = None
            self._cobs_data_for = None

        r_helio_disp = obs_ovr["r_helio"]   if obs_ovr else info["r_helio"]
        r_geo_disp   = obs_ovr["r_geo"]     if obs_ovr else info["r_geo"]
        phase_disp   = obs_ovr["phase_angle"]if obs_ovr else info["phase_angle"]
        ovr_tag      = " ★obs" if obs_ovr else ""

        self.status_comet.setText(f"☄  {self._comet_el.get('name','')[:35]}")
        self.status_helio.setText(f"r☉ = {r_helio_disp:.4f} AU{ovr_tag}")
        self.status_geo.setText(f"Δ = {r_geo_disp:.4f} AU{ovr_tag}")
        self.status_phase.setText(f"Phase = {phase_disp:.1f}°{ovr_tag}")
        mode = "Overlay" if img_arr is not None else "Standalone"
        if img_arr is not None and hasattr(self.ctrl,"_wcs_ps_deg") and self.ctrl._wcs_ps_deg:
            ps_as = self.ctrl._wcs_ps_deg * 3600
            mode += f"  ·  {ps_as:.3f}\"/px  ·  {ov['au_per_px']:.2e} AU/px"
        if obs_ovr:
            mode += "  ·  ★ ephemeris override active"
        self.status_mode.setText(mode)
        self.status.showMessage("Model computed successfully", 3000)
        self.ctrl.set_computing(False)
        self._act_run_fp.setEnabled(True)

    def _effective_plate_scale_arcsec(self):
        """
        Return the effective plate scale in arcsec/px with priority:
          1. User's manual entry in the Image Setup dialog ps_spin
          2. FITS WCS-derived _wcs_ps_deg (if a FITS was loaded)
          3. None (no plate scale available)

        Both attributes live on ControlPanel (self.ctrl), not MainWindow.
        """
        ctrl = self.ctrl

        # Priority 1: Image Setup dialog — user's manually entered value
        if hasattr(ctrl, '_img_dialog') and ctrl._img_dialog is not None:
            ps_arcsec = ctrl._img_dialog.ps_spin.value()
            if ps_arcsec > 0:
                return ps_arcsec

        # Priority 2: FITS WCS header (auto-populated on FITS load)
        if hasattr(ctrl, '_wcs_ps_deg') and ctrl._wcs_ps_deg is not None:
            return ctrl._wcs_ps_deg * 3600.0  # deg/px → arcsec/px

        return None

    def _on_compute_error(self, msg):
        QMessageBox.critical(self, "Compute Error", msg)
        self.status.showMessage(f"Error: {msg}", 5000)
        self.ctrl.set_computing(False)
        self._act_run_fp.setEnabled(True)

    # ── Fetch ──────────────────────────────────────────────────────────────
    def _on_fetch(self, desig, date):
        self._fetch_worker = FetchWorker(desig, date)
        self._fetch_worker.finished.connect(self.ctrl.on_fetch_done)
        self._fetch_worker.error.connect(self.ctrl.on_fetch_error)
        self._fetch_worker.start()

    # ── Image ──────────────────────────────────────────────────────────────
    def _on_image_loaded(self, arr):
        if arr is None:
            self.canvas._imgArr = None
            if self._model:
                self.canvas.draw_model(self._model, None)
            else:
                self.canvas._draw_empty()
            self.status_mode.setText("Standalone mode")
            return

        # ── Show image immediately ────────────────────────────────────────
        self.canvas._imgArr = arr
        ov = self.ctrl.get_overlay()
        self.canvas.nuc_x    = ov["nuc_x"]
        self.canvas.nuc_y    = ov["nuc_y"]
        self.canvas.au_per_px= ov["au_per_px"]
        self.canvas.north_pa = ov["north_pa"]

        # Gather WCS info for overlay
        hdr = getattr(self.ctrl, "_fits_header", None)
        wcs_info = None
        if hdr is not None:
            cd12 = float(hdr.get("CD1_2", 0))
            cd22 = float(hdr.get("CD2_2", 1e-4))
            ps_as = np.sqrt(cd12**2 + cd22**2) * 3600
            wcs_info = dict(
                ps_arcsec = ps_as,
                npa       = ov["north_pa"],
                object    = hdr.get("OBJECT", ""),
                date_obs  = hdr.get("DATE-OBS", ""),
            )
        fname = ""
        if hasattr(self.ctrl, "_last_image_path"):
            fname = os.path.basename(self.ctrl._last_image_path)

        self.canvas.draw_image_preview(arr, filename=fname, wcs_info=wcs_info)
        h, w = arr.shape[:2]
        self.status_mode.setText(f"Image loaded  {w}×{h} px  ·  Adjust stretch → COMPUTE")

    def _on_nucleus_clicked(self, x, y):
        self.ctrl.nuc_x_spin.setValue(x)
        self.ctrl.nuc_y_spin.setValue(y)
        self.status.showMessage(f"Nucleus set to ({x:.0f}, {y:.0f}) px", 3000)

    # ── Export ─────────────────────────────────────────────────────────────
    def _current_obs_jd(self) -> float:
        """The 'current working date' for Orbit/calculators. Prefers the
        obs_date text field — which the embedded Animator keeps in sync
        with whichever frame is currently shown (v3.0: scrub/play the
        Animator, then open Orbit, and it shows THAT date) — falling back
        to the comet's static fetch-time obs_jd if obs_date is empty or
        unparseable, then today as a last resort."""
        try:
            txt = self.ctrl.obs_date.text().strip()
            if txt:
                return cta.date_to_jd(txt)
        except Exception:
            pass
        el = self._active_comet_el()
        if el and el.get("obs_jd") is not None:
            return el["obs_jd"]
        return cta.today_jd()

    def _open_orbit_window(self):
        el = self._active_comet_el()
        if not el:
            QMessageBox.warning(self, "No Comet",
                                "Select or fetch a comet first.")
            return
        obs_jd = self._current_obs_jd()
        # NOTE: parent=None (not parent=self) — a non-modal QDialog parented
        # to MainWindow is kept "transient for" it by the window manager,
        # which on most platforms means it's pinned ABOVE MainWindow at all
        # times, even after MainWindow regains focus (this was the "Orbit
        # window stays in front always" bug). LCWindow already avoids this
        # the same way. Stored as self._orbit_win (not a throwaway local)
        # to prevent a garbage-collection crash, same reasoning as LCWindow.
        self._orbit_win = OrbitWindow(el, obs_jd, parent=None)
        self._orbit_win.show()
        self._orbit_win.raise_()
        self._orbit_win.activateWindow()

    def _open_mc_window(self):
        # Keep one editable MC window per main-window session.  Repeated clicks
        # on Model → Open Monte Carlo Model should reveal/raise the existing
        # window instead of creating another copy with diverging inputs.
        existing = getattr(self, '_mc_win', None)
        if existing is not None:
            try:
                ages_now = self.ctrl._parse_floats(self.ctrl.age_str.text())
                existing.update_fp_guidance(
                    ages_now, float(self.ctrl.dominant_age.value()))
                # Restore a minimized window before bringing it to front.
                if existing.windowState() & Qt.WindowState.WindowMinimized:
                    existing.setWindowState(
                        existing.windowState() & ~Qt.WindowState.WindowMinimized)
                existing.show()
                existing.raise_()
                existing.activateWindow()
                return
            except RuntimeError:
                # WA_DeleteOnClose can leave a stale Python wrapper briefly.
                self._mc_win = None

        el = self._active_comet_el()
        if not el:
            QMessageBox.warning(self, "No Comet",
                                "Select or fetch a comet first.")
            return
        obs_jd  = self._current_obs_jd()
        betas   = self.ctrl._parse_floats(self.ctrl.beta_str.text())
        ages    = self.ctrl._parse_floats(self.ctrl.age_str.text())
        # Keep the MC/F-P upper age synchronized with the actual synchrone
        # list.  The hidden main-panel widget is only a compatibility carrier.
        max_age = max(ages) if ages else self.ctrl.max_age.value()
        # Same parent=None reasoning as OrbitWindow/LCWindow (see their
        # comments) — keeps this non-modal window from getting perma-
        # pinned above MainWindow by the window manager. Stored on self
        # to avoid a garbage-collection crash, same as those two.
        win = MCWindow(el, obs_jd, beta_values=betas, max_age=max_age,
                       fp_sync_ages=ages,
                       fp_dominant_age=float(self.ctrl.dominant_age.value()),
                       fp_model=self._model, parent=None,
                       main_window=self)
        self._mc_win = win

        # Clear the reference only when this exact window is destroyed.  This
        # makes the next menu click create a fresh instance after the user has
        # genuinely closed the previous one.
        def _forget_mc_window(*_args, _win=win):
            if getattr(self, '_mc_win', None) is _win:
                self._mc_win = None
                try:
                    self._act_reextract_mc.setEnabled(False)
                except Exception:
                    pass
        win.destroyed.connect(_forget_mc_window)

        win.show()
        win.raise_()
        win.activateWindow()

    def _save_png(self):
        if not self._model:
            QMessageBox.information(self, "No Model", "Compute a model first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "", "PNG Image (*.png);;PDF (*.pdf)")
        if path:
            # 300 dpi = standard print/publication minimum for the PNG path;
            # the PDF option is vector and prints crisp at any size regardless.
            self.canvas.fig.savefig(path, dpi=300, bbox_inches="tight",
                                    facecolor=self.canvas.fig.get_facecolor())
            self.status.showMessage(f"Saved → {path}", 4000)

    def _export_csv(self):
        if not self._model:
            QMessageBox.information(self, "No Model", "Compute a model first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export CSV", "", "CSV File (*.csv)")
        if not path: return
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["type","parameter","xi_AU","eta_AU"])
            for s in self._model["syndynes"]:
                for xi,eta,age in zip(s["xi"],s["eta"],s["age"]):
                    if np.isfinite(xi) and np.isfinite(eta):
                        w.writerow(["syndyne",s["beta"],f"{xi:.8f}",f"{eta:.8f}"])
            for s in self._model["synchrones"]:
                for xi,eta,b in zip(s["xi"],s["eta"],s["beta"]):
                    if np.isfinite(xi) and np.isfinite(eta):
                        w.writerow(["synchrone",s["age"],f"{xi:.8f}",f"{eta:.8f}"])
        self.status.showMessage(f"CSV exported → {path}", 4000)

    # ── Update check (v3.1.1) ────────────────────────────────────────────
    # Startup checks stay silent unless a genuine newer stable release exists.
    # Manual checks distinguish "current" from a network/API/SSL failure.
    def _update_settings(self) -> QSettings:
        return QSettings("MPC-O58", "CometTailAnalyzer")

    def _check_for_update_silent(self):
        if self._update_worker is not None and self._update_worker.isRunning():
            return
        self._update_worker = UpdateCheckWorker()
        self._update_worker.finished.connect(self._on_update_check_silent)
        self._update_worker.start()

    def _on_update_check_silent(self, info: dict | None):
        if not isinstance(info, dict) or info.get("status") != "update":
            return
        skipped = self._update_settings().value("skipped_update_version", "")
        if skipped == info.get("latest", ""):
            return
        self._show_update_dialog(info, manual=False)

    def _check_for_update_manual(self):
        if self._update_worker is not None and self._update_worker.isRunning():
            self.status.showMessage("An update check is already in progress.", 3000)
            return
        self.status.showMessage("Checking for updates…")
        self._update_worker = UpdateCheckWorker()
        self._update_worker.finished.connect(self._on_update_check_manual)
        self._update_worker.start()

    def _on_update_check_manual(self, info: dict | None):
        self.status.showMessage("", 0)
        if not isinstance(info, dict):
            info = {"status": "error", "reason": "No response from the update service."}

        status = info.get("status")
        if status == "update":
            self._show_update_dialog(info, manual=True)
            return
        if status == "current":
            latest = info.get("tag") or info.get("latest") or cta.__version__
            QMessageBox.information(
                self, "Check for Updates",
                f"Comet Tail Analyzer is up to date.\n\n"
                f"Installed version: {cta.__version__}\n"
                f"Latest stable release: {latest}")
            return

        reason = info.get("reason") or "The GitHub Releases service could not be reached."
        box = QMessageBox(self)
        box.setWindowTitle("Unable to Check for Updates")
        box.setIcon(QMessageBox.Icon.Warning)
        box.setText("CTA could not verify the latest release.")
        box.setInformativeText(
            f"{reason}\n\nInstalled version: {cta.__version__}\n"
            "The application can continue to be used normally.")
        btn_open = box.addButton("Open Releases Page", QMessageBox.ButtonRole.AcceptRole)
        box.addButton("Close", QMessageBox.ButtonRole.RejectRole)
        box.exec()
        if box.clickedButton() is btn_open:
            webbrowser.open(f"https://github.com/{cta.GITHUB_REPO}/releases/latest")

    def _show_update_dialog(self, info: dict, manual: bool):
        box = QMessageBox(self)
        box.setWindowTitle("Update Available")
        box.setIcon(QMessageBox.Icon.Information)
        notes = info.get("notes") or ""
        if len(notes) >= 500:
            notes = notes.rstrip() + "…"
        box.setText(
            f"<b>Comet Tail Analyzer {info['tag']}</b> is available "
            f"(you have {cta.__version__}).")
        if notes:
            box.setInformativeText(notes)
        btn_open = box.addButton("Open Release Page", QMessageBox.ButtonRole.AcceptRole)
        box.addButton("Later", QMessageBox.ButtonRole.RejectRole)
        btn_skip = None
        if not manual:
            btn_skip = box.addButton("Skip This Version", QMessageBox.ButtonRole.DestructiveRole)
        box.exec()
        clicked = box.clickedButton()
        if clicked is btn_open:
            webbrowser.open(info["url"])
        elif btn_skip is not None and clicked is btn_skip:
            self._update_settings().setValue("skipped_update_version", info["latest"])

    def _about(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("About Comet Tail Analyzer")
        dlg.setMinimumWidth(520)
        dlg.setMaximumWidth(560)
        dlg.setStyleSheet("background:#07090f;")
        vb = QVBoxLayout(dlg); vb.setContentsMargins(0,0,0,0); vb.setSpacing(0)

        # Header banner
        banner = QLabel()
        banner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        banner.setStyleSheet(
            "background: qlineargradient(x1:0,y1:0,x2:1,y2:0,"
            "stop:0 #050a18, stop:0.5 #0a1835, stop:1 #050a18);"
            "padding: 22px 16px;")
        banner.setText(
            "<div style='font-family:'Segoe UI',Arial,sans-serif;text-align:center'>"
            "<div style='font-size:36px;margin-bottom:4px'>☄</div>"
            "<div style='font-size:18px;font-weight:bold;color:#7ac8ff;letter-spacing:2px'>"
            "COMET TAIL ANALYZER</div>"
            "<div style='font-size:10px;color:#2a5070;letter-spacing:3px;margin-top:4px'>"
            "FINSON–PROBSTEIN DUST TAIL MODEL  ·  1968</div>"
            "<div style='font-size:11px;color:#3a6090;margin-top:8px'>Version 3.1.1  ·  2026</div>"
            "</div>")
        vb.addWidget(banner)

        # Divider
        div = QFrame(); div.setFrameShape(QFrame.Shape.HLine)
        div.setStyleSheet("color:#1a2540; background:#1a2540; height:1px;")
        vb.addWidget(div)

        # Body
        body = QScrollArea(); body.setWidgetResizable(True)
        body.setStyleSheet("border:none; background:#07090f;")
        body.setMaximumHeight(340)
        inner = QWidget(); inner.setStyleSheet("background:#07090f;")
        bl = QVBoxLayout(inner); bl.setContentsMargins(24,18,24,16); bl.setSpacing(12)

        def section(title, content_html):
            t = QLabel(f"<span style='color:#6ab0d8;font-size:10px;"
                       f"letter-spacing:2px;font-family:'Segoe UI',Arial,sans-serif'>{title}</span>")
            c = QLabel(content_html)
            c.setWordWrap(True)
            c.setStyleSheet("color:#a0c0d8; font-size:11px; font-family:'Segoe UI',Arial,sans-serif; line-height:1.7;")
            c.setTextFormat(Qt.TextFormat.RichText)
            bl.addWidget(t); bl.addWidget(c)

        section("FINSON–PROBSTEIN MODEL (1968)",
            "<b style='color:#80c8e8'>β = F<sub>rad</sub> / F<sub>grav</sub></b>"
            "<span style='color:#5a90a8'> — radiation-to-gravity force ratio.<br>"
            "Larger β = smaller grains, pushed further into the tail.</span><br><br>"
            "<span style='color:#ff9090'>● Syndynes</span>"
            "<span style='color:#a0c0d8'> — same β (grain radius), different emission times.</span><br>"
            "<span style='color:#ffe080'>● Synchrones</span>"
            "<span style='color:#a0c0d8'> — same emission time, different β (grain radii).</span><br><br>"
            "<span style='color:#3a6080'>⚠ Zero ejection velocity assumed.</span>")

        section("DEVELOPER",
            "<b style='color:#80d8ff;font-size:13px'>Teerasak Thaluang</b><br>"
            "<span style='color:#5090a8'>MPC Station Code: O51 &nbsp;·&nbsp; O58</span><br>"
            "<span style='color:#3a6070'>Comet and Asteroid Observer </span>")

        section("REFERENCES",
            "Finson M.L. &amp; Probstein R.F. (1968a)<br>"
            "&nbsp;&nbsp;<i>ApJ</i> 154, 327–352 — Theory of dust comets I<br>"
            "Finson M.L. &amp; Probstein R.F. (1968b)<br>"
            "&nbsp;&nbsp;<i>ApJ</i> 154, 353–380 — Theory of dust comets II<br>"
            "Bredichin T.A. (1884) — Syndyne / Synchrone concepts")

        section("PHYSICS ENGINE",
            "Analytical Keplerian propagation under<br>"
            "&nbsp;&nbsp;μ<sub>eff</sub> = μ<sub>☉</sub>·(1−β) &nbsp;[AU³/day²]<br>"
            "RK4 numerical integrator (fallback for β &gt; 1)<br>"
            "Gnomonic sky-plane projection (RA/Dec tangent plane)<br>"
            "Earth ephemeris: astropy VSOP87 / analytic backup")

        section("SOFTWARE STACK",
            "Python 3  ·  NumPy  ·  SciPy  ·  Astropy<br>"
            "Matplotlib  ·  PyQt6  ·  Pillow  ·  astroquery<br>"
            "Data sources: JPL Horizons (sole orbital element source)  ·  COBS")

        section("LICENSE",
            "Open-source for scientific and educational use.<br>"
            "Please cite: <i>Finson &amp; Probstein (1968)</i> when publishing results.<br>"
            "Credit: <b>Teerasak Thaluang</b> for software.")

        bl.addStretch()
        body.setWidget(inner)
        vb.addWidget(body)

        # Footer buttons
        div2 = QFrame(); div2.setFrameShape(QFrame.Shape.HLine)
        div2.setStyleSheet("color:#1a2540; background:#1a2540; height:1px;")
        vb.addWidget(div2)
        foot = QWidget(); foot.setStyleSheet("background:#06080e; padding:8px;")
        fl = QHBoxLayout(foot); fl.setContentsMargins(16,8,16,8)
        lbl_ver = QLabel(
            "<span style='color:#1a3050;font-size:10px;font-family:'Segoe UI',Arial,sans-serif'>"
            "☄ Comet Tail Analyzer v3.0 · Teerasak Thaluang · MPC-O51/O58</span>")
        fl.addWidget(lbl_ver,1)
        ok_btn = QPushButton("  OK  ")
        ok_btn.clicked.connect(dlg.accept)
        ok_btn.setFixedWidth(80)
        fl.addWidget(ok_btn)
        vb.addWidget(foot)

        dlg.exec()

# ─────────────────────────────────────────────────────────────────────────────
#  ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────
def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    # Apply initial dark theme
    apply_theme("dark", app)

    win = MainWindow()
    # Maximized keeps the normal title bar/taskbar behavior while using the
    # full available desktop area. This is intentionally not borderless
    # fullscreen, so minimize/restore and multi-window work remain normal.
    win.showMaximized()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
