import os
import pandas as pd
import matplotlib.pyplot as plt

chrom = os.environ["CHR"]
start = int(os.environ["START"])
end = int(os.environ["END"])
windows = [int(x) for x in os.environ["WINDOWS"].split()]

samples = [
    ("sample1", "results/insulation/sample1_insulation.tsv"),
    ("sample2", "results/insulation/sample2_insulation.tsv"),
]

rows = []

for name, path in samples:
    df = pd.read_table(path)

    reg = df[
        (df["chrom"] == chrom) &
        (df["start"] >= start) &
        (df["end"] <= end)
    ]

    for w in windows:
        col = f"is_boundary_{w}"

        n_boundaries = int(
            reg[col].astype(str).isin(["True", "TRUE", "1"]).sum()
        )

        n_tads_approx = max(n_boundaries - 1, 0)

        rows.append({
            "sample": name,
            "window": w,
            "n_boundaries": n_boundaries,
            "n_tads_approx": n_tads_approx
        })

out = pd.DataFrame(rows)
out.to_csv("results/tables/tad_count_by_window.tsv", sep="\t", index=False)

plt.figure(figsize=(7, 5))

for name, sub in out.groupby("sample"):
    plt.plot(sub["window"], sub["n_boundaries"], marker="o", label=name)

plt.xlabel("Window size, bp")
plt.ylabel("Number of TAD boundaries")
plt.title(f"TAD boundary count vs window size: {chrom}:{start}-{end}")
plt.legend()
plt.tight_layout()
plt.savefig("figures/tad_boundaries_count_by_window.png", dpi=200)
