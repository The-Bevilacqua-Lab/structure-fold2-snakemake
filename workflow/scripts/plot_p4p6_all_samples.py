##########################################################################
# Colors p4p6 reactivity for multiple samples onto the known P4-P6 domain
# structure diagram, on a shared color scale. Only meaningful for the p4p6
# positive control -- see the note in plot_p4p6.py.
##########################################################################

import argparse
import os
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

RESOURCES_DIR = Path(__file__).resolve().parent.parent.parent / "resources" / "p4p6"
DEFAULT_PAIRING_FILE = str(RESOURCES_DIR / "p4p6_pairing.txt")
DEFAULT_IMAGE_FILE = str(RESOURCES_DIR / "p4p6_structure.png")

# Pixel coordinates for each A/C position on the P4-P6 structure image
SEQ_DICT = {
    "1": (1370, 1200), "2": (1330, 1200), "6": (1250, 965),
    "10": (1280, 790), "11": (1295, 725), "12": (1295, 660),
    "18": (1250, 370), "19": (1250, 327), "20": (510, 267),
    "21": (510, 325), "22": (510, 375), "24": (522, 465),
    "25": (522, 505), "29": (522, 648), "30": (522, 688),
    "33": (522, 940), "34": (522, 984), "35": (522, 1025),
    "36": (522, 1255), "37": (522, 1297), "40": (522, 1425),
    "42": (522, 1507), "43": (522, 1550), "48": (547, 1774),
    "49": (600, 1774), "50": (635, 1734), "51": (625, 1680),
    "56": (625, 1465), "58": (625, 1380), "62": (790, 1240),
    "63": (835, 1240), "67": (975, 1185), "68": (960, 1145),
    "69": (917, 1120), "70": (875, 1130), "75": (700, 1053),
    "80": (660, 920), "81": (680, 890), "83": (660, 810),
    "84": (630, 795), "86": (628, 730), "89": (628, 605),
    "90": (628, 550), "93": (626, 426), "94": (637, 375),
    "95": (637, 325), "100": (1145, 499), "101": (1145, 540),
    "103": (1102, 695), "104": (1112, 772), "105": (1147, 842),
    "106": (1147, 882), "107": (1107, 905), "108": (1147, 928),
    "110": (1146, 1013), "111": (1144, 1052), "113": (1147, 1201),
    "114": (1147, 1245), "115": (1130, 1300), "116": (1130, 1360),
    "119": (1147, 1505), "120": (1147, 1550), "122": (1097, 1637),
    "123": (1113, 1686), "126": (1145, 1807), "127": (1145, 1848),
    "128": (1145, 1892), "129": (1145, 1933), "130": (1145, 1975),
    "132": (1133, 2059), "134": (1200, 2120), "137": (1251, 2018),
    "143": (1251, 1765), "145": (1281, 1670), "149": (1250, 1461),
    "152": (1270, 1300), "153": (1250, 1240), "157": (1330, 970),
    "158": (1370, 970),
}


def load_sample(path, accepted):
    df = pd.read_csv(path)
    df.columns = ["transcript", "index", "seq", "reactivity"]
    df = df[df["transcript"] == "p4p6"].copy()
    df["index"] = df["index"] - 2
    merged = pd.merge(accepted, df, on=["index", "seq"])
    merged = merged[(merged["seq"] == "A") | (merged["seq"] == "C")]
    merged = merged[(merged["pairing"] == "0") | (merged["pairing"] == "1")]
    merged = merged[~merged["reactivity"].isna()]
    return merged


def main():
    parser = argparse.ArgumentParser(description="Plot p4p6 reactivity for all samples on a common scale")
    parser.add_argument("--input", nargs="+", required=True, help="Reactivity CSV files, one per sample")
    parser.add_argument("--samples", nargs="+", default=None, help="Sample labels (defaults to parent directory names)")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--pairing", default=DEFAULT_PAIRING_FILE, help="p4p6 pairing table")
    parser.add_argument("--image", default=DEFAULT_IMAGE_FILE, help="p4p6 structure image")
    parser.add_argument("--ncols", type=int, default=4, help="Max columns in subplot grid")
    args = parser.parse_args()

    samples = args.samples if args.samples else [os.path.basename(os.path.dirname(p)) for p in args.input]

    accepted = pd.read_csv(args.pairing)
    img = mpimg.imread(args.image)

    data_per_sample = [load_sample(path, accepted) for path in args.input]

    all_reactivities = pd.concat([d["reactivity"] for d in data_per_sample])
    vmin, vmax = all_reactivities.min(), all_reactivities.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("white_red", ["white", "red"])

    n = len(data_per_sample)
    ncols = min(n, args.ncols)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 8 * nrows), squeeze=False,
                              gridspec_kw={"wspace": 0.0, "hspace": 0})

    for idx, (data, sample) in enumerate(zip(data_per_sample, samples)):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        ax.imshow(img)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(sample, fontsize=10)
        for spine in ax.spines.values():
            spine.set_visible(False)

        for _, row_data in data.iterrows():
            key = str(int(row_data["index"]))
            if key in SEQ_DICT:
                color = cmap(norm(row_data["reactivity"]))
                ax.scatter(SEQ_DICT[key][0], SEQ_DICT[key][1], color=color, marker="o", s=50, zorder=-10)

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.02, pad=0.04, shrink=0.6)
    cbar.set_label("Reactivity", rotation=270, labelpad=15)

    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
