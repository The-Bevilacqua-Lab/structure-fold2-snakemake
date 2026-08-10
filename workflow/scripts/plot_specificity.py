import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser(description="Plot per-sample base specificity as horizontal stacked bar graph")
parser.add_argument("--inputs", nargs="+", required=True)
parser.add_argument("--samples", nargs="+", required=True)
parser.add_argument("--conditions", nargs="+", required=True)
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

frames = []
for path, sample, condition in zip(args.inputs, args.samples, args.conditions):
    df = pd.read_csv(path)
    tmp = df[[df.columns[0], df.columns[2]]].copy()
    tmp.columns = ["base", "specificity"]
    tmp["sample"] = sample
    tmp["condition"] = condition
    frames.append(tmp)

long_df = pd.concat(frames, ignore_index=True)
data = long_df.pivot_table(index=["sample", "condition"], columns="base", values="specificity")

# minus samples first, then plus, each block in the order they were passed in
sample_order = [s for s, c in zip(args.samples, args.conditions) if c == "minus"] + \
               [s for s, c in zip(args.samples, args.conditions) if c == "plus"]
condition_map = dict(zip(args.samples, args.conditions))

data = data.reset_index().set_index("sample").loc[sample_order]

base_colors = {"A": "#4C72B0", "T": "#DD8452", "G": "#55A868", "C": "#C44E52"}
bases = [b for b in ["T", "G", "C", "A"] if b in data.columns]

fig, ax = plt.subplots(figsize=(4.5, 0.4 * len(sample_order) + 1))

left = [0.0] * len(data)
for base in bases:
    values = data[base].fillna(0).tolist()
    ax.barh(range(len(data)), values, left=left,
            label="U" if base == "T" else base, color=base_colors[base],
            height=0.8, edgecolor="white", linewidth=0.5)
    left = [l + v for l, v in zip(left, values)]

ax.set_yticks(range(len(data)))
# labels reflect the actual sample count/order passed in -- not a fixed
# 3-replicate assumption, so this works for any number of samples
tick_labels = [f"{s} ({'-' if condition_map[s] == 'minus' else '+'} DMS)" for s in sample_order]
ax.set_yticklabels(tick_labels)
plt.ylabel("Sample")
plt.xlabel("Proportion of stops")
ax.set_xlim(0, 1)
ax.invert_yaxis()  # top-to-bottom order matches input order
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

legend_patches = [mpatches.Patch(color=base_colors[b], label="U" if b == "T" else b) for b in bases]
ax.legend(handles=legend_patches, bbox_to_anchor=(1.01, 1), loc="upper left", fancybox=False)

plt.tight_layout()
plt.savefig(args.output, dpi=500, bbox_inches="tight")
plt.close()
print(f"Saved plot to {args.output}")
