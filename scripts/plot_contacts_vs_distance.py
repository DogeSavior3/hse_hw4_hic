import os
import pandas as pd
import matplotlib.pyplot as plt

chrom = os.environ["CHR"]

samples = [
    ("sample1", "results/expected/sample1_expected_cis.tsv"),
    ("sample2", "results/expected/sample2_expected_cis.tsv"),
]

plt.figure(figsize=(7, 5))

for name, path in samples:
    df = pd.read_table(path)

    if "region1" in df.columns:
        df = df[df["region1"].astype(str) == chrom]

    df = df[df["dist_bp"] > 0]

    if "balanced.avg" in df.columns:
        ycol = "balanced.avg"
    else:
        ycol = "count.avg"

    plt.loglog(df["dist_bp"], df[ycol], marker=".", linestyle="-", label=name)

plt.xlabel("Genomic distance, bp")
plt.ylabel("Average contact frequency")
plt.title(f"Contacts vs distance, {chrom}")
plt.legend()
plt.tight_layout()
plt.savefig("figures/contacts_vs_distance.png", dpi=200)
