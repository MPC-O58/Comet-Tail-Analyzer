#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
  CometTailGUI.py  —  Finson–Probstein Comet Dust Tail Analyzer
  Version 2.4   ·   Teerasak Thaluang (MPC O51/O58)
  Native Desktop Application (PyQt6 + Matplotlib)
=============================================================================
  Changelog:
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

__version__ = "2.4"

import sys, os, warnings, csv
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
    QMenu,
)
from PyQt6.QtCore import (
    Qt, QThread, pyqtSignal, QTimer, QSize,
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
BG      = T["mpl_bg"]
PANEL_BG= T["panel_bg"]
MF      = "DejaVu Sans"   # clean sans-serif for plot labels

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

    def __init__(self, comet_el, obs_jd, betas, ages, max_age, n_pts):
        super().__init__()
        self.comet_el = comet_el
        self.obs_jd   = obs_jd
        self.betas    = betas
        self.ages     = ages
        self.max_age  = max_age
        self.n_pts    = n_pts

    def run(self):
        try:
            self.progress.emit(10, "Computing dust trajectories…")
            model = cta.compute_model(
                self.comet_el, self.obs_jd,
                self.betas, self.ages,
                self.max_age, self.n_pts)
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
        self._vis = dict(synd=True, sync=True, orbit=True, lw=1.5, alpha=0.85)
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
        self._imgArr = img_arr
        self.ax.clear()
        self.fig.patch.set_facecolor("black")
        self.ax.set_facecolor("black")
        for sp in self.ax.spines.values():
            sp.set_edgecolor("#1a2540")

        h, w = img_arr.shape[:2]
        self.ax.imshow(img_arr, origin="upper", aspect="equal",
                       extent=[0, w, h, 0], zorder=0)
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

    # ── Main draw ─────────────────────────────────────────────────────────
    def draw_model(self, model, img_arr=None):
        self._model  = model
        self._imgArr = img_arr
        self.ax.clear()
        self.ax.set_facecolor(BG)
        self.fig.patch.set_facecolor(BG)
        for sp in self.ax.spines.values():
            sp.set_edgecolor("#1a2540")
        self.ax.tick_params(labelcolor="#2a4060", labelsize=8, color="#1a2540")
        self.ax.grid(color="#0d1a2e", lw=0.4, zorder=0)

        overlay = img_arr is not None

        # ── AU → degree conversion factor (standalone only) ──────────────
        r_geo = model["info"].get("r_geo", 1.0)
        K = (180.0 / np.pi) / r_geo   # deg per AU  (small-angle, accurate for <5°)

        if overlay:
            # Draw image
            h_img, w_img = img_arr.shape[:2]
            self.ax.imshow(img_arr, origin="upper", aspect="equal",
                           extent=[0, w_img, h_img, 0], zorder=0)
            self.ax.set_xlim(0, w_img)
            self.ax.set_ylim(h_img, 0)
            self.ax.set_xlabel("X (pixels)", color="#3a6080", fontfamily=MF)
            self.ax.set_ylabel("Y (pixels)", color="#3a6080", fontfamily=MF)

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
            if len(all_xi) > 1:
                xr  = np.percentile(all_xi,  [2, 98])
                yr  = np.percentile(all_eta, [2, 98])
                pad = 0.12
                dx  = max((xr[1]-xr[0])*pad, 0.01)
                dy  = max((yr[1]-yr[0])*pad, 0.01)
                # East LEFT = inverted x-axis (E is positive xi, appears left)
                self.ax.set_xlim(xr[1]+dx, xr[0]-dx)
                self.ax.set_ylim(yr[0]-dy, yr[1]+dy)
            self.ax.axhline(0, color="#1a2a40", lw=0.7, zorder=1)
            self.ax.axvline(0, color="#1a2a40", lw=0.7, zorder=1)
            self.ax.set_xlabel("← East  ·  Δ (°)  ·  West →",
                               color="#3a6080", fontsize=9, fontfamily=MF)
            self.ax.set_ylabel("↑  Δ North (°)",
                               color="#3a6080", fontsize=9, fontfamily=MF)
            # Starfield (uses sorted xlim to avoid issues with inverted axis)
            rng = np.random.default_rng(42)
            xl  = sorted(self.ax.get_xlim()); yl = self.ax.get_ylim()
            sx  = rng.uniform(*xl, 200); sy = rng.uniform(*yl, 200)
            self.ax.scatter(sx, sy, s=0.3, c="white", alpha=0.2, zorder=0)

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
        if self._vis["orbit"] and model["orbit"]:
            op = np.array(model["orbit"])
            if overlay:
                xs, ys = zip(*[to_px(p[0], p[1]) for p in op])
            else:
                xs = op[:,0] * K
                ys = op[:,1] * K
            self.ax.plot(xs, ys, "--", color="#2050a0", lw=0.9,
                         alpha=0.5, zorder=2, label="Orbital path")

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
        if self._vis["synd"]:
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
        if self._vis["sync"]:
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

        # ── Sun direction arrow ───────────────────────────────────────────
        # sun_dir = (xi, eta) pointing FROM COMET TOWARD SUN
        # JPL Horizons PsAng = PA of anti-solar direction = PA toward Sun + 180°
        sun_xi, sun_eta = model["sun_dir"]
        slen = np.sqrt(sun_xi**2 + sun_eta**2)
        if slen > 1e-10:
            # Compute PsAng (anti-solar = dust tail direction)
            # PA measured from North (eta) toward East (xi): atan2(xi, eta)
            PA_to_sun  = np.degrees(np.arctan2(sun_xi/slen, sun_eta/slen)) % 360
            PsAng      = (PA_to_sun + 180.0) % 360   # anti-solar (Horizons PsAng)

            if overlay:
                sc = 60 * self.au_per_px
                ax_end, ay_end = to_px(sun_xi/slen*sc, sun_eta/slen*sc)
                ax_start, ay_start = self.nuc_x, self.nuc_y
            else:
                # ── FIX: abs() to prevent sign flip from inverted x-axis ──
                xl2 = self.ax.get_xlim()
                sc  = abs(xl2[1] - xl2[0]) * 0.13  # always positive ✓
                # In standalone, xi/eta already in degrees via K
                ax_end   = (sun_xi/slen) * (sc * K / K)   # same as sun_xi/slen*sc_deg
                ay_end   = (sun_eta/slen) * (sc * K / K)
                # Actually sc is in degrees since axis is in degrees:
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

            # PsAng label (standalone only)
            if not overlay:
                xl3 = self.ax.get_xlim(); yl3 = self.ax.get_ylim()
                # bottom-left in screen = West+South = xr[0]-dx (right data), yr[0]-dy
                lx = xl3[0] + (xl3[1]-xl3[0])*0.04   # right side (West) slightly inside
                ly = yl3[0] + (yl3[1]-yl3[0])*0.06
                self.ax.text(lx, ly,
                             f"PsAng = {PsAng:.1f}°\n(anti-solar / tail dir)",
                             color="#888830", fontsize=7.5, ha="right",
                             fontfamily=MF, va="bottom", zorder=8)

        # Nucleus
        nx = self.nuc_x if overlay else 0.0
        ny = self.nuc_y if overlay else 0.0
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
            # Overlay: image pixel coords, y increases downward
            h_img, w_img = self._imgArr.shape[:2] if self._imgArr is not None else (1000, 1000)
            aL_c = min(w_img, h_img) * 0.038
            cx   = w_img * 0.07
            cy   = h_img * 0.07
            # North: screen-up = negative py
            n_dx =  np.sin(npa_r) * aL_c
            n_dy = -np.cos(npa_r) * aL_c   # negative py = up in image ✓
            # East: screen-left = negative px
            e_dx = -np.cos(npa_r) * aL_c   # negative px = left ✓
            e_dy = -np.sin(npa_r) * aL_c

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
        if self._vis["synd"]:
            handles.append(Line2D([0],[0], color=SYNDYNE_COLORS[0], lw=1.5,
                                  label="Syndynes  (β = const)"))
        if self._vis["sync"]:
            handles.append(Line2D([0],[0], color=SYNC_COLORS[0], lw=1.2,
                                  ls="-.", label="Synchrones (age = const)"))
        if self._vis["orbit"]:
            handles.append(Line2D([0],[0], color="#2050a0", lw=0.9,
                                  ls="--", label="Orbital path"))
        handles.append(Line2D([0],[0], color="#ffe030", lw=1.5,
                              marker=">", ms=7, label="Sun direction"))
        if handles:
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
        self.ax.set_title(
            f'{info.get("name","Comet")}   ·   {info["obs_str"]}   ·   Finson–Probstein\n'
            f'r☉={info["r_helio"]:.4f} AU  ·  Δ={r_g:.4f} AU  ·  '
            f'Phase={info["phase_angle"]:.1f}°  ·  PsAng={_PsAng:.1f}°{ovr_tag}',
            color="#7ab8ff", fontsize=9, fontfamily=MF, pad=8)

        # Force axes back to exact image bounds after drawing
        # (matplotlib may auto-expand axes when curves extend beyond)
        if overlay:
            self.ax.set_xlim(0, w_img)
            self.ax.set_ylim(h_img, 0)

        self.canvas.draw_idle()

# ─────────────────────────────────────────────────────────────────────────────
#  LEFT CONTROL PANEL
# ─────────────────────────────────────────────────────────────────────────────
class ControlPanel(QScrollArea):
    compute_requested = pyqtSignal(dict, float, list, list, int, int)
    fetch_requested   = pyqtSignal(str, str)
    image_loaded      = pyqtSignal(object)   # object allows None or ndarray

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWidgetResizable(True)
        self.setMinimumWidth(320)
        self.setMaximumWidth(400)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

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
        gm.addWidget(self.age_str)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Max age (d)"))
        self.max_age = QSpinBox()
        self.max_age.setRange(10, 3650); self.max_age.setValue(200)
        row1.addWidget(self.max_age)
        gm.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Resolution"))
        self.n_pts = QSpinBox()
        self.n_pts.setRange(20, 200); self.n_pts.setValue(70)
        row2.addWidget(self.n_pts)
        gm.addLayout(row2)

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

        # ── Display options ───────────────────────────────────────────────
        grp_disp = QGroupBox("DISPLAY")
        gd = QVBoxLayout(grp_disp); gd.setSpacing(4)
        self.chk_synd  = QCheckBox("Syndynes (β = const)")
        self.chk_sync  = QCheckBox("Synchrones (age = const)")
        self.chk_orbit = QCheckBox("Orbital path")
        self.chk_synd.setChecked(True)
        self.chk_sync.setChecked(True)
        self.chk_orbit.setChecked(True)
        self.chk_synd.setStyleSheet("color:#ff7070;")
        self.chk_sync.setStyleSheet("color:#ffd060;")
        self.chk_orbit.setStyleSheet("color:#4a8adf;")
        for w in [self.chk_synd, self.chk_sync, self.chk_orbit]:
            gd.addWidget(w)

        row_lw = QHBoxLayout()
        row_lw.addWidget(QLabel("Line width"))
        self.lw_slider = QSlider(Qt.Orientation.Horizontal)
        self.lw_slider.setRange(5, 50); self.lw_slider.setValue(15)
        self.lw_label  = QLabel("1.5")
        self.lw_slider.valueChanged.connect(lambda v: self.lw_label.setText(f"{v/10:.1f}"))
        row_lw.addWidget(self.lw_slider); row_lw.addWidget(self.lw_label)
        gd.addLayout(row_lw)

        row_op = QHBoxLayout()
        row_op.addWidget(QLabel("Opacity  "))
        self.op_slider = QSlider(Qt.Orientation.Horizontal)
        self.op_slider.setRange(10, 100); self.op_slider.setValue(85)
        self.op_label  = QLabel("0.85")
        self.op_slider.valueChanged.connect(lambda v: self.op_label.setText(f"{v/100:.2f}"))
        row_op.addWidget(self.op_slider); row_op.addWidget(self.op_label)
        gd.addLayout(row_op)
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
        self.sl_lo      = QSlider(Qt.Orientation.Horizontal)
        self.sl_hi      = QSlider(Qt.Orientation.Horizontal)
        self.sl_k       = QSlider(Qt.Orientation.Horizontal)
        self.sl_gamma   = QSlider(Qt.Orientation.Horizontal)
        self.sl_lo.setRange(0,500);     self.sl_lo.setValue(450)   # 4.50%
        self.sl_hi.setRange(9000,10000); self.sl_hi.setValue(9850)  # 98.50%
        self.sl_k.setRange(1,200);      self.sl_k.setValue(200)    # 2.00
        self.sl_gamma.setRange(50,300); self.sl_gamma.setValue(300) # 3.00
        self.btn_pick_nuc  = QPushButton("⊕  CLICK NUCLEUS ON PLOT")
        self.btn_rot_nup   = QPushButton("↻  ROTATE N-UP / E-LEFT")
        self.rot_status    = QLabel("")
        self.btn_auto_stretch = QPushButton("✨  AUTO STRETCH")
        self.btn_auto_stretch.setToolTip("Automatically compute optimal Shadow/Highlight from image histogram")
        self.btn_restretch = QPushButton("↺  APPLY STRETCH")
        self.btn_auto_stretch.clicked.connect(self._auto_stretch)
        self.btn_restretch.clicked.connect(self._restretch)
        self.btn_rot_nup.clicked.connect(self._rotate_north_up)
        self.npa_spin.valueChanged.connect(self._update_npa_hint)

        # Give all stub widgets the hidden container as parent
        for w in [self.nuc_x_spin, self.nuc_y_spin, self.au_px_spin,
                  self.npa_spin, self.npa_hint, self.sl_lo, self.sl_hi,
                  self.sl_k, self.sl_gamma, self.btn_pick_nuc,
                  self.btn_rot_nup, self.rot_status, self.btn_restretch]:
            _stub_lay.addWidget(w)

        # btn_clear_img — must also have a parent, not be orphaned
        self.btn_clear_img = QPushButton("CLEAR", _stub)
        self.btn_clear_img.setProperty("class","danger")
        self.btn_clear_img.clicked.connect(self._clear_image)

        # overlay_params stub
        self.overlay_params = QWidget(_stub)
        self.btn_img = self.btn_img_open   # alias

        vbox.addWidget(grp_img)

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
        self._img_dialog.raise_()
        self._img_dialog.activateWindow()

    def _toggle_theme(self):
        """Toggle dark/light theme — called from the ☀/☾ button in the header."""
        mw = self.window()
        if hasattr(mw, '_set_theme'):
            new = "light" if CURRENT_THEME == "dark" else "dark"
            mw._set_theme(new)

    def _update_beta_hint(self):
        betas = self._parse_floats(self.beta_str.text())
        hints = [f"β{b}→{cta.beta_to_size(b)}" for b in betas[:4]]
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

            # Pre-fill observation date from FITS header
            date_obs = hdr.get("DATE-OBS","")
            if date_obs and hasattr(self, "obs_date"):
                d = date_obs[:10]
                self.obs_date.setText(d)

            # Pre-fill object name if recognised
            obj = hdr.get("OBJECT","")
            if obj:
                # Try to match against DB
                for i, name in enumerate(list(cta.COMET_DB.keys())):
                    if any(tok in name for tok in obj.split()):
                        self.combo_comet.setCurrentIndex(i)
                        break

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

    def _stretch_raw(self, data: np.ndarray) -> np.ndarray:
        """Apply current slider-based stretch to raw float32 FITS data → RGB uint8."""
        lo_pct = self.sl_lo.value()    / 100.0   # e.g. 0.30
        hi_pct = self.sl_hi.value()    / 100.0   # e.g. 99.98
        k      = self.sl_k.value()     / 100.0   # asinh softness  e.g. 0.05
        gamma  = self.sl_gamma.value() / 100.0   # gamma  e.g. 1.00
        lo  = np.percentile(data, lo_pct)
        hi  = np.percentile(data, hi_pct)
        d   = np.clip(data, lo, hi)
        d   = (d - lo) / (hi - lo + 1e-8)
        d   = np.arcsinh(d / k) / np.arcsinh(1.0 / k)
        if abs(gamma - 1.0) > 0.01:
            d = np.power(np.clip(d, 1e-8, 1.0), gamma)
        d   = np.clip(d, 0.0, 1.0)
        img8 = (d * 255).astype(np.uint8)
        from PIL import Image as PILImage
        return np.array(PILImage.fromarray(img8, "L").convert("RGB"))

    def _auto_stretch(self):
        """
        Apply optimal default stretch settings for comet images.
        Sets Shadow=4.50%, Highlight=98.50%, Softness=2.00, Gamma=3.00
        """
        if not hasattr(self, "_fits_raw") or self._fits_raw is None:
            QMessageBox.information(self, "Auto Stretch",
                "No FITS image loaded. Auto stretch only works with FITS files.")
            return
        
        # Set optimal default values for comet imaging
        self.sl_lo.setValue(450)      # Shadow: 4.50%
        self.sl_hi.setValue(9850)     # Highlight: 98.50%
        self.sl_k.setValue(200)       # Softness: 2.00
        self.sl_gamma.setValue(300)   # Gamma: 3.00
        
        # Auto-apply stretch
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
        if not betas or not ages:
            QMessageBox.warning(self, "Bad Input",
                                "Enter at least one valid β value and synchrone age.")
            return
        obs_jd = self._comet_el.get("obs_jd",
                    cta.date_to_jd(self.obs_date.text()) if self.comet_tabs.currentIndex()==0
                    else cta.today_jd())
        self.compute_requested.emit(
            self._comet_el, obs_jd, betas, ages,
            self.max_age.value(), self.n_pts.value())

    def get_rotation_offset(self) -> float:
        """Return grid rotation offset in degrees (positive = CW = more toward West)."""
        return self.rot_slider.value() / 10.0

    def get_vis(self):
        return dict(
            synd=self.chk_synd.isChecked(),
            sync=self.chk_sync.isChecked(),
            orbit=self.chk_orbit.isChecked(),
            lw=self.lw_slider.value()/10,
            alpha=self.op_slider.value()/100,
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
        self.setMinimumWidth(280)
        self.setMaximumWidth(320)

        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        # Thin header strip (stored so update_theme() can re-style it)
        self.hdr = QLabel("  ANALYSIS")
        self.hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_muted']}; font-size:10px;"
            f"letter-spacing:2px; font-weight:bold; padding:6px 8px;"
            f"border-bottom:1px solid {T['border']};")
        vbox.addWidget(self.hdr)

        self.tabs = QTabWidget()
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
        v1.addStretch()
        self._apply_info_label_styles()   # initial colour pass
        self.tabs.addTab(sa1, "INFO")

        # ── Tab 2 : β TABLE ─────────────────────────────────────────────
        sa2 = QScrollArea(); sa2.setWidgetResizable(True)
        sa2.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        w2 = QWidget(); sa2.setWidget(w2)
        v2 = QVBoxLayout(w2); v2.setContentsMargins(8,8,8,8); v2.setSpacing(6)

        hdr_lbl = QLabel("Select β values that match the observed tail:")
        hdr_lbl.setStyleSheet(f"color:{T['text_muted']}; font-size:11px;")
        hdr_lbl.setWordWrap(True)
        self._beta_hdr_lbl = hdr_lbl   # theme ref
        v2.addWidget(hdr_lbl)

        # Table: ✓ | β | Size | Effect
        self.beta_table = QTableWidget(0, 4)
        self.beta_table.setHorizontalHeaderLabels(["✓","β","Size","Effect"])
        hh = self.beta_table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self.beta_table.setColumnWidth(0, 28)
        self.beta_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.beta_table.setAlternatingRowColors(True)
        self.beta_table.verticalHeader().setVisible(False)
        self.beta_table.setStyleSheet(
            "QTableWidget{font-size:11px;}"
            "QTableWidget::item{padding:3px 4px;}")
        self.beta_table.setSelectionMode(
            QTableWidget.SelectionMode.NoSelection)
        v2.addWidget(self.beta_table)

        note = QLabel("☑ Tick β values whose syndyne aligns with the observed tail.\n"
                      "Selected β → grain sizes used in Analysis report.")
        note.setStyleSheet(f"color:{T['text_muted']}; font-size:10px; padding:4px 0;")
        note.setWordWrap(True)
        self._beta_note_lbl = note   # theme ref
        v2.addWidget(note)
        v2.addStretch()
        self.tabs.addTab(sa2, "β TABLE")

        # ── Tab 3 : LIGHT CURVE ──────────────────────────────────────────
        w4 = QWidget()
        v4 = QVBoxLayout(w4); v4.setContentsMargins(8,8,8,8); v4.setSpacing(8)

        # Buttons row
        fetch_row = QHBoxLayout(); fetch_row.setSpacing(6)
        self.btn_fetch_cobs = QPushButton("📡  FETCH COBS")
        self.btn_fetch_cobs.setProperty("class","success")
        self.btn_fetch_cobs.setEnabled(False)
        self.btn_fetch_cobs.setMinimumHeight(34)
        fetch_row.addWidget(self.btn_fetch_cobs)

        self.btn_open_lc = QPushButton("📈  PLOT LC")
        self.btn_open_lc.setProperty("class","primary")
        self.btn_open_lc.setEnabled(False)
        self.btn_open_lc.setMinimumHeight(34)
        self.btn_open_lc.setToolTip("Open full light curve window")
        fetch_row.addWidget(self.btn_open_lc)
        v4.addLayout(fetch_row)

        # H₀/n results display
        res_grp = QGroupBox("PHOTOMETRIC PARAMETERS")
        res_grid = QGridLayout(res_grp); res_grid.setSpacing(6)

        lbl_h0 = QLabel("H₀"); lbl_h0.setStyleSheet(f"color:{T['text_muted']};font-size:11px;")
        self._lbl_h0 = lbl_h0   # theme ref
        self.val_h0 = QLabel("—")
        self.val_h0.setStyleSheet(f"color:{T['text_value']};font-size:15px;font-weight:700;")
        lbl_n = QLabel("n"); lbl_n.setStyleSheet(f"color:{T['text_muted']};font-size:11px;")
        self._lbl_n = lbl_n     # theme ref
        self.val_n = QLabel("—")
        self.val_n.setStyleSheet(f"color:{T['text_value']};font-size:15px;font-weight:700;")
        res_grid.addWidget(lbl_h0, 0, 0)
        res_grid.addWidget(self.val_h0, 0, 1)
        res_grid.addWidget(lbl_n, 0, 2)
        res_grid.addWidget(self.val_n, 0, 3)
        v4.addWidget(res_grp)

        # Status/source label
        self.cobs_lbl = QLabel("Press FETCH COBS to load light curve data")
        self.cobs_lbl.setStyleSheet(f"color:{T['text_muted']}; font-size:10px;")
        self.cobs_lbl.setWordWrap(True)
        v4.addWidget(self.cobs_lbl)
        v4.addStretch()
        self.tabs.addTab(w4, "LIGHT CURVE")

        # ── Tab 5 : ANALYSIS ─────────────────────────────────────────────
        w5 = QWidget()
        v5 = QVBoxLayout(w5); v5.setContentsMargins(6,6,6,6); v5.setSpacing(6)

        # Afρ input row
        afrho_row = QHBoxLayout()
        afrho_row.addWidget(QLabel("Afρ (cm)"))
        self.afrho_spin = QDoubleSpinBox()
        self.afrho_spin.setRange(0, 1e6); self.afrho_spin.setDecimals(0)
        self.afrho_spin.setValue(0); self.afrho_spin.setSpecialValueText("—")
        self.afrho_spin.setSuffix(" cm")
        self.afrho_spin.setToolTip(
            "Afρ parameter measured from the image.\n"
            "A'Hearn et al. (1984)  [cm]")
        afrho_row.addWidget(self.afrho_spin)
        v5.addLayout(afrho_row)

        # Generate analysis button
        self.btn_analyse = QPushButton("🔬  GENERATE ANALYSIS")
        self.btn_analyse.setProperty("class","success")
        self.btn_analyse.setEnabled(False)
        v5.addWidget(self.btn_analyse)

        # Analysis text output
        self.analysis_text = QTextEdit()
        self.analysis_text.setReadOnly(True)
        self._apply_analysis_text_style()   # sets bg/border/color from T
        self.analysis_text.setPlaceholderText(
            "Compute model first, then click\n'GENERATE ANALYSIS'")
        v5.addWidget(self.analysis_text)

        # Copy button
        self.btn_copy_analysis = QPushButton("📋  COPY TEXT")
        self.btn_copy_analysis.clicked.connect(self._copy_analysis)
        v5.addWidget(self.btn_copy_analysis)
        self.tabs.addTab(w5, "ANALYSIS")


    def _copy_analysis(self):
        from PyQt6.QtWidgets import QApplication
        QApplication.clipboard().setText(self.analysis_text.toPlainText())

    def update_info(self, model, comet_el, betas):
        info = model["info"]
        self.eph_labels["r☉"].setText(f"{info['r_helio']:.5f} AU")
        self.eph_labels["Δ"].setText(f"{info['r_geo']:.5f} AU")
        self.eph_labels["Phase"].setText(f"{info['phase_angle']:.2f}°")
        self.eph_labels["RA"].setText(f"{info['RA']:.4f}°")
        self.eph_labels["Dec"].setText(f"{info['Dec']:.4f}°")
        self.eph_labels["Date"].setText(info["obs_str"][:16])

        # β TABLE with checkboxes
        self._betas = betas   # store for get_selected_betas()
        self.beta_table.setRowCount(len(betas))
        effects = {0.001:"~orbital trail",0.005:"barely pushed",0.01:"slightly pushed",
                   0.05:"moderately pushed",0.1:"pushed outward",
                   0.3:"strongly pushed",0.6:"very strongly",1.0:"no gravity"}
        for i, b in enumerate(betas):
            eff = next((v for k,v in effects.items()
                        if abs(b-k) < max(b*0.4, 0.002)), "pushed")
            # Col 0: checkbox
            chk = QTableWidgetItem()
            chk.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            chk.setCheckState(Qt.CheckState.Unchecked)
            chk.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.beta_table.setItem(i, 0, chk)
            # Col 1: β value (coloured)
            b_item = QTableWidgetItem(str(b))
            b_item.setForeground(QColor(SYNDYNE_COLORS[i % len(SYNDYNE_COLORS)]))
            b_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.beta_table.setItem(i, 1, b_item)
            # Col 2: grain size
            sz_item = QTableWidgetItem(cta.beta_to_size(b))
            sz_item.setForeground(QColor(T["text_sec"]))
            sz_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.beta_table.setItem(i, 2, sz_item)
            # Col 3: effect
            eff_item = QTableWidgetItem(eff)
            eff_item.setForeground(QColor(T["text_muted"]))
            eff_item.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            self.beta_table.setItem(i, 3, eff_item)
        self.beta_table.resizeRowsToContents()

        self.orb_labels["q"].setText(f"{comet_el['q']:.6f} AU")
        self.orb_labels["e"].setText(f"{comet_el['e']:.7f}")
        self.orb_labels["i"].setText(f"{comet_el['i']:.4f}°")
        self.orb_labels["Ω"].setText(f"{comet_el['Omega']:.4f}°")
        self.orb_labels["ω"].setText(f"{comet_el['omega']:.4f}°")
        self.orb_labels["T"].setText(comet_el.get("T","")[:10])
        self.orb_labels["Source"].setText(comet_el.get("source", "JPL Horizons"))
        self.tabs.setCurrentIndex(0)

    def get_selected_betas(self) -> list:
        """Return list of β values ticked in the β table."""
        sel = []
        for i in range(self.beta_table.rowCount()):
            chk = self.beta_table.item(i, 0)
            if chk and chk.checkState() == Qt.CheckState.Checked:
                try:
                    b = float(self.beta_table.item(i, 1).text())
                    sel.append(b)
                except Exception: pass
        return sel or getattr(self, '_betas', [])   # fallback: all betas

    # ── Theme helpers ──────────────────────────────────────────────────────
    def _apply_analysis_text_style(self):
        """Re-apply the analysis QTextEdit colours from the active theme."""
        self.analysis_text.setStyleSheet(
            f"background:{T['input_bg']}; border:1px solid {T['border']};"
            f"border-radius:6px; color:{T['text']}; font-size:10px;")

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
        """Re-apply β TABLE, LC, and header strip colours from the active theme."""
        # Header strip
        self.hdr.setStyleSheet(
            f"background:{T['panel_bg']}; color:{T['text_muted']}; font-size:10px;"
            f"letter-spacing:2px; font-weight:bold; padding:6px 8px;"
            f"border-bottom:1px solid {T['border']};")
        # β TABLE labels
        if hasattr(self, '_beta_hdr_lbl'):
            self._beta_hdr_lbl.setStyleSheet(
                f"color:{T['text_muted']}; font-size:11px;")
        if hasattr(self, '_beta_note_lbl'):
            self._beta_note_lbl.setStyleSheet(
                f"color:{T['text_muted']}; font-size:10px; padding:4px 0;")
        # LIGHT CURVE tab labels
        if hasattr(self, '_lbl_h0'):
            self._lbl_h0.setStyleSheet(f"color:{T['text_muted']}; font-size:11px;")
        if hasattr(self, '_lbl_n'):
            self._lbl_n.setStyleSheet(f"color:{T['text_muted']}; font-size:11px;")
        if hasattr(self, 'val_h0'):
            self.val_h0.setStyleSheet(
                f"color:{T['text_value']}; font-size:15px; font-weight:700;")
        if hasattr(self, 'val_n'):
            self.val_n.setStyleSheet(
                f"color:{T['text_value']}; font-size:15px; font-weight:700;")
        if hasattr(self, 'cobs_lbl'):
            self.cobs_lbl.setStyleSheet(
                f"color:{T['text_muted']}; font-size:10px;")

    def update_theme(self):
        """Re-apply ALL theme-sensitive styles in the right panel."""
        self._apply_analysis_text_style()
        self._apply_info_label_styles()
        self._apply_misc_label_styles()



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
      - FITS stretch (Shadow, Highlight, Softness, Gamma)
      - Observed Ephemeris Override (PsAng, r☉, Δ, Phase)
    """
    def __init__(self, ctrl: "ControlPanel"):
        super().__init__(ctrl.window() if ctrl.window() else ctrl,
                         Qt.WindowType.Tool |
                         Qt.WindowType.WindowStaysOnTopHint)
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

        vb.addLayout(sl_row("Shadow %",   ctrl.sl_lo,    0,500,   450,  100, "{:.2f}%"))
        vb.addLayout(sl_row("Highlight %",ctrl.sl_hi,    9000,10000,9850,100,"{:.2f}%"))
        vb.addLayout(sl_row("Softness",   ctrl.sl_k,     1,200,    200,  100, "{:.2f}"))
        vb.addLayout(sl_row("Gamma",      ctrl.sl_gamma, 50,300,  300,  100, "{:.2f}"))
        
        # Auto Stretch + Apply Stretch buttons
        stretch_btns = QHBoxLayout()
        stretch_btns.addWidget(ctrl.btn_auto_stretch)
        stretch_btns.addWidget(ctrl.btn_restretch)
        vb.addLayout(stretch_btns)

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
#  MAIN WINDOW
# ─────────────────────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("☄ Comet Tail Analyzer  —  Finson–Probstein Model  ·  v2.4")
        self.setMinimumSize(1300, 800)
        self.resize(1440, 900)

        self._model    = None
        self._comet_el = None
        self._worker   = None

        self._build_menu()
        self._build_ui()
        self._build_status()

    # ── Menu ──────────────────────────────────────────────────────────────
    def _build_menu(self):
        mb = self.menuBar()

        # File
        file_m = mb.addMenu("File")
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

        # View
        view_m = mb.addMenu("View")
        act_reset = QAction("Reset zoom", self)
        act_reset.triggered.connect(lambda: self.canvas.canvas.toolbar.home())
        view_m.addAction(act_reset)
        view_m.addSeparator()

        # Theme submenu
        theme_m = view_m.addMenu("🎨  Theme")
        self._act_dark  = QAction("☾  Dark  (Space)",  self, checkable=True)
        self._act_light = QAction("☀  Light (Observatory)", self, checkable=True)
        self._act_dark.setChecked(True)
        self._act_dark.triggered.connect(lambda: self._set_theme("dark"))
        self._act_light.triggered.connect(lambda: self._set_theme("light"))
        theme_m.addAction(self._act_dark)
        theme_m.addAction(self._act_light)

        # Help
        help_m = mb.addMenu("Help")
        act_about = QAction("About…", self)
        act_about.triggered.connect(self._about)
        help_m.addAction(act_about)

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

        # Left
        self.ctrl = ControlPanel()
        self.ctrl.compute_requested.connect(self._on_compute)
        self.ctrl.fetch_requested.connect(self._on_fetch)
        self.ctrl.image_loaded.connect(self._on_image_loaded)
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
        splitter.setSizes([370, 800, 300])

        hbox.addWidget(splitter)

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


    # ── LC popup window ────────────────────────────────────────────────
    def _open_lc_window(self):
        """Open a larger standalone light curve window."""
        d = getattr(self, "_cobs_data", None)
        if not d or d.get("H0") is None:
            QMessageBox.information(self, "No Data",
                "Fetch COBS light curve first."); return
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
        self.info.btn_fetch_cobs.setEnabled(False)
        self.info.btn_fetch_cobs.setText("⏳  Fetching…")

        class COBSWorker(QThread):
            done  = pyqtSignal(dict)
            error = pyqtSignal(str)
            def __init__(self, n, el): super().__init__(); self.n=n; self.el=el
            def run(self):
                try: self.done.emit(cta.fetch_from_cobs(self.n, self.el))
                except Exception as e: self.error.emit(str(e))

        self._cobs_worker = COBSWorker(name, self._comet_el or {})
        def _on_done(d):
            self._cobs_data = d
            H0 = d.get("H0"); n_val = d.get("n")
            # Show H₀/n prominently
            try:
                self.info.val_h0.setText(f"{H0:.2f}" if H0 else "—")
                self.info.val_n.setText(f"{n_val:.2f}" if n_val else "—")
            except Exception: pass
            raw  = len(d.get("raw_obs", []))
            nos  = len(d.get("obs_list", []))
            src  = d.get("source","")
            last = d.get("last_date","")
            parts = []
            if raw:  parts.append(f"{raw} COBS observations")
            if nos:  parts.append(f"{nos} with ephemeris for fit")
            if last: parts.append(f"Latest: m={d.get('last_mag','?')}  ({last})")
            parts.append(f"Source: {src}")
            self.info.cobs_lbl.setText("\n".join(parts))
            self.info.btn_fetch_cobs.setEnabled(True)
            self.info.btn_fetch_cobs.setText("📡  FETCH COBS")
            try: self.info.btn_open_lc.setEnabled(True)
            except Exception: pass
            msg = f"H₀={H0:.2f}  n={n_val:.2f}" if H0 else "COBS: no fit"
            self.status.showMessage(msg + f"  ({src[:35]})", 8000)
        def _on_err(msg):
            self.info.btn_fetch_cobs.setEnabled(True)
            self.info.btn_fetch_cobs.setText("📡  FETCH COBS")
            self.info.cobs_lbl.setText(f"⚠ COBS fetch failed:\n{msg[:120]}")
            self.status.showMessage("COBS fetch failed.", 3000)
        self._cobs_worker.done.connect(_on_done)
        self._cobs_worker.error.connect(_on_err)
        self._cobs_worker.start()

    # ── Analysis report ───────────────────────────────────────────────
    def _run_analysis(self, comet_el: dict, model: dict):
        if not model: return
        info = model["info"]
        afrho = self.info.afrho_spin.value() or None
        selected = self.info.get_selected_betas()
        info["best_betas"] = selected
        info["obs_jd"] = cta.date_to_jd(
            self.ctrl.obs_date.text() or cta.jd_to_str(cta.today_jd()))
        cobs = getattr(self, "_cobs_data", None)
        txt = cta.generate_dust_analysis(comet_el, info, cobs, afrho)
        self.info.analysis_text.setPlainText(txt)
        self.info.tabs.setCurrentIndex(3)
        self.status.showMessage(f"Analysis: {len(selected)} β selected", 2500)

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

    def _on_compute(self, comet_el, obs_jd, betas, ages, max_age, n_pts):
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
        self.status.showMessage("Computing Finson–Probstein model…")

        self._worker = ComputeWorker(comet_el, obs_jd, betas, ages, max_age, n_pts)
        self._worker.progress.connect(
            lambda v, msg: (self.ctrl.progress_bar.setValue(v),
                            self.status.showMessage(msg)))
        self._worker.finished.connect(lambda m: self._on_model_ready(m, betas, ov))
        self._worker.error.connect(self._on_compute_error)
        self._worker.start()

    def _on_model_ready(self, model, betas, ov):
        self._model = model
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
        self.canvas.draw_model(model, img_arr)
        self.info.update_info(model, self._comet_el, betas)

        # Enable COBS + analysis buttons
        try: self.info.btn_fetch_cobs.setEnabled(True)
        except Exception: pass
        try: self.info.btn_fetch_cobs.clicked.disconnect()
        except Exception: pass
        try:
            _name = self._comet_el.get("name","")
            self.info.btn_fetch_cobs.clicked.connect(
                lambda: self._fetch_cobs(_name))
        except Exception: pass


        try: self.info.btn_open_lc.setEnabled(False)  # enabled after COBS fetch
        except Exception: pass
        try: self.info.btn_open_lc.clicked.disconnect()
        except Exception: pass
        try:
            self.info.btn_open_lc.clicked.connect(self._open_lc_window)
        except Exception: pass

        try: self.info.btn_analyse.setEnabled(True)
        except Exception: pass
        try: self.info.btn_analyse.clicked.disconnect()
        except Exception: pass
        try:
            _el = self._comet_el; _mdl = model
            self.info.btn_analyse.clicked.connect(
                lambda: self._run_analysis(_el, _mdl))
        except Exception: pass

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
    def _save_png(self):
        if not self._model:
            QMessageBox.information(self, "No Model", "Compute a model first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "", "PNG Image (*.png);;PDF (*.pdf)")
        if path:
            self.canvas.fig.savefig(path, dpi=200, bbox_inches="tight",
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
            "<div style='font-size:11px;color:#3a6090;margin-top:8px'>Version 2.4  ·  2026</div>"
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
            "<span style='color:#a0c0d8'> — same β (grain size), different emission times.</span><br>"
            "<span style='color:#ffe080'>● Synchrones</span>"
            "<span style='color:#a0c0d8'> — same emission time, different β (grain sizes).</span><br><br>"
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
            "☄ Comet Tail Analyzer v2.4 · Teerasak Thaluang · MPC-O51/O58</span>")
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
    win.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
