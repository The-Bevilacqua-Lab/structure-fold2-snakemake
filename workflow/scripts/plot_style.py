"""Shared plot style: Arial in every figure, plotnine and matplotlib alike.

Import it before drawing:  import plot_style  (the scripts add workflow/scripts to sys.path first).
Arial comes from conda-forge's mscorefonts package (in the plotting envs), which installs into
<env>/fonts, a folder matplotlib doesn't search, so the fonts are registered here.
"""

import sys
from pathlib import Path

import matplotlib
from matplotlib import font_manager

FONT = "Arial"

for _ttf in Path(sys.prefix, "fonts").glob("*.ttf"):
    font_manager.fontManager.addfont(str(_ttf))

if not any(f.name == FONT for f in font_manager.fontManager.ttflist):
    print(f"plot_style: {FONT} not found (is mscorefonts in the env?); using DejaVu Sans", file=sys.stderr)

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [FONT, "DejaVu Sans"]

try:   # plotnine themes take their font from this option
    from plotnine import options
    options.base_family = FONT
except ImportError:   # matplotlib-only env (upsetplot)
    pass
