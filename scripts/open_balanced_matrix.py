import os
import cooler
import numpy as np

samples = [
    ("sample1", os.environ["URI1"]),
    ("sample2", os.environ["URI2"]),
]

chrom = os.environ["CHR"]
region = os.environ["REGION"]

for name, uri in samples:
    clr = cooler.Cooler(uri)

    mat_chr = clr.matrix(balance=True).fetch(chrom)
    mat_region = clr.matrix(balance=True).fetch(region)

    np.savetxt(
        f"results/tables/{name}_balanced_matrix_region.tsv",
        mat_region,
        delimiter="\t"
    )

    with open(f"results/info/{name}_balanced_matrix_info.txt", "w") as f:
        f.write(f"Sample: {name}\n")
        f.write(f"URI: {uri}\n")
        f.write(f"Chromosome: {chrom}\n")
        f.write(f"Region: {region}\n")
        f.write(f"Balanced chromosome matrix shape: {mat_chr.shape}\n")
        f.write(f"Balanced region matrix shape: {mat_region.shape}\n")
