##########################################################################
# Colors p4p6 reactivity onto the known P4-P6 domain structure diagram.
#
# Only meaningful for the p4p6 positive control (Tetrahymena group I intron
# P4-P6 domain) -- the pixel coordinates in SEQ_DICT and the structure image
# are specific to that one construct. The pipeline only invokes this script
# when config['positive_control_name'] == "p4p6" (see
# workflow/rules/positive_control.smk); it is not applicable to any other
# organism or positive control.
##########################################################################

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import pandas as pd

RESOURCES_DIR = Path(__file__).resolve().parent.parent.parent / "resources" / "p4p6"
DEFAULT_PAIRING_FILE = RESOURCES_DIR / "p4p6_pairing.txt"
DEFAULT_IMAGE_FILE = RESOURCES_DIR / "p4p6_structure.png"

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


def main():
    parser = argparse.ArgumentParser(description="Plot p4p6 reactivity onto the known structure diagram")
    parser.add_argument("--input", type=str, required=True, help="reactivity.csv (transcript,position,base,reactivity)")
    parser.add_argument("--output", type=str, required=True, help="Path to save the output plot")
    parser.add_argument("--pairing", type=str, default=str(DEFAULT_PAIRING_FILE),
                         help="p4p6 paired(1)/unpaired(0) reference table")
    parser.add_argument("--image", type=str, default=str(DEFAULT_IMAGE_FILE),
                         help="p4p6 structure diagram PNG")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df.columns = ["transcript", "index", "seq", "reactivity"]
    df = df[df["transcript"] == "p4p6"].copy()
    df["index"] = df["index"] - 2

    accepted = pd.read_csv(args.pairing)
    merged = pd.merge(accepted, df, on=["index", "seq"])
    merged = merged[(merged["seq"] == "A") | (merged["seq"] == "C")]
    merged = merged[(merged["pairing"] == "0") | (merged["pairing"] == "1")]
    merged = merged[~merged["reactivity"].isna()]

    values = merged["reactivity"].tolist()
    norm = mcolors.Normalize(vmin=min(values), vmax=max(values))
    cmap = cm.Reds
    merged["color"] = [mcolors.to_hex(cmap(norm(v))) for v in values]

    fig, ax = plt.subplots(figsize=(5, 8))
    for spine in ax.spines.values():
        spine.set_visible(False)

    img = mpimg.imread(args.image)
    ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Reactivity", rotation=270, labelpad=15)

    for _, row in merged.iterrows():
        key = str(int(row["index"]))
        if key in SEQ_DICT:
            ax.scatter(SEQ_DICT[key][0], SEQ_DICT[key][1], color=row["color"], marker="o", s=50, zorder=-1)

    plt.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
