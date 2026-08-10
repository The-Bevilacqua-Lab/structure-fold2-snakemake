import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot mean 3' TP coverage vs. trim amount from end-coverage sweep TSV."
    )
    parser.add_argument("--input", required=True, help="TSV with columns trim, mean_tp_coverage")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--id", default=None, help="Sample ID for plot title")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    df = df[df["mean_tp_coverage"] != "NA"].copy()
    df["mean_tp_coverage"] = pd.to_numeric(df["mean_tp_coverage"], errors="coerce")
    df = df.dropna(subset=["mean_tp_coverage"])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df["trim"], df["mean_tp_coverage"], linewidth=1.5, color="#2c7bb6")
    ax.set_xlabel("3' Trim (nt)", fontsize=12)
    ax.set_ylabel("Mean 3' TP Coverage", fontsize=12)
    title = "3' End Coverage vs. Trim Amount"
    if args.id:
        title += f"\n{args.id}"
    ax.set_title(title, fontsize=13)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()


if __name__ == "__main__":
    main()
