import os
import pandas as pd
import matplotlib.pyplot as plt

chrom = os.environ["CHR"]
start = int(os.environ["START"])
end = int(os.environ["END"])
main_window = os.environ["MAIN_WINDOW"]

score_col = f"log2_insulation_score_{main_window}"
boundary_col = f"is_boundary_{main_window}"

samples = [
    ("sample1", "results/insulation/sample1_insulation.tsv"),
    ("sample2", "results/insulation/sample2_insulation.tsv"),
]

plt.figure(figsize=(11, 5))

for name, path in samples:
    df = pd.read_table(path)

    reg = df[
        (df["chrom"] == chrom) &
        (df["start"] >= start) &
        (df["end"] <= end)
    ].copy()

    reg["mid"] = (reg["start"] + reg["end"]) / 2

    plt.plot(reg["mid"], reg[score_col], label=name)

    boundaries = reg[reg[boundary_col].astype(str).isin(["True", "TRUE", "1"])]

    for x in boundaries["mid"]:
        plt.axvline(x, linestyle="--", alpha=0.25)

plt.xlabel(f"Position on {chrom}, bp")
plt.ylabel(f"log2 insulation score, window={main_window}")
plt.title(f"Insulation score and TAD boundaries: {chrom}:{start}-{end}")
plt.legend()
plt.tight_layout()
plt.savefig(f"figures/insulation_boundaries_w{main_window}.png", dpi=200)
