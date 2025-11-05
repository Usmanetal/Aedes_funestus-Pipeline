import pandas as pd
import matplotlib.pyplot as plt

# --- Read header line manually ---
with open("GGAFUN_DELTA_2.genes.txt") as f:
    for line in f:
        if line.startswith("#GeneName"):
            columns = line.strip().lstrip("#").split("\t")
            break

# --- Read data ---
df = pd.read_csv("GGAFUN_DELTA_2.genes.txt", sep="\t", comment="#", header=None)
df.columns = columns

# --- Filter for high/moderate impact variants ---
high_mod = df[(df["variants_impact_HIGH"] > 0) | (df["variants_impact_MODERATE"] > 0)]

# --- Compute total impact and select top 20 genes ---
top_genes = (
    high_mod.assign(total_impact=lambda x: x["variants_impact_HIGH"] + x["variants_impact_MODERATE"])
    .sort_values("total_impact", ascending=False)
    .head(20)
)

# --- Plot ---
plt.figure(figsize=(10, 6))
plt.barh(top_genes["GeneName"], top_genes["total_impact"], color="steelblue")
plt.xlabel("Number of HIGH + MODERATE Variants")
plt.ylabel("Gene")
plt.title("Top 20 Genes with Highest Impact SNPs (GGAFUN_DELTA sample SnpEff Annotation)")
plt.gca().invert_yaxis()
plt.tight_layout()

# ✅ Save first, then show
plt.savefig("top_20_genes_high_moderate_variants11.png", dpi=300)
plt.show()

# --- Summary of counts ---
print("\n=== Total Variants by Impact ===")
print(high_mod[["variants_impact_HIGH", "variants_impact_MODERATE"]].sum())



import pandas as pd
import matplotlib.pyplot as plt

# --- 1️⃣ Extract header names directly from the file ---
header_line = None
with open("FANG.ann.genes.txt") as f:
    for line in f:
        if line.startswith("#GeneName"):
            header_line = line.strip().lstrip("#").split("\t")
            break

if not header_line:
    raise ValueError("Header line starting with #GeneName not found in FANG.ann.genes.txt")

# --- 2️⃣ Read the actual data (skip comment lines) ---
df = pd.read_csv("FANG.ann.genes.txt", sep="\t", comment="#", header=None)
df.columns = header_line

# --- 3️⃣ Filter for genes with HIGH or MODERATE impact variants ---
high_mod = df[
    (df["variants_impact_HIGH"] > 0) | (df["variants_impact_MODERATE"] > 0)
]

# --- 4️⃣ Compute total impact and extract top 20 genes ---
top_genes = (
    high_mod.assign(
        total_impact=lambda x: x["variants_impact_HIGH"] + x["variants_impact_MODERATE"]
    )
    .sort_values("total_impact", ascending=False)
    .head(20)
)

# --- 5️⃣ Plot horizontal bar chart ---
plt.figure(figsize=(10, 6))
plt.barh(top_genes["GeneName"], top_genes["total_impact"], color="darkcyan")
plt.xlabel("Number of HIGH + MODERATE Variants")
plt.ylabel("Gene")
plt.title("Top 20 Genes with Highest Impact SNPs (FANG sample, SnpEff annotation)")
plt.gca().invert_yaxis()
plt.tight_layout()

# --- 6️⃣ Save BEFORE showing the plot ---
plt.savefig("FANG_top20_high_moderate_variants.png", dpi=300)
plt.show()

# --- 7️⃣ Print summary counts for sanity check ---
print("\n=== Total Variants by Impact (FANG sample) ===")
print(
    high_mod[["variants_impact_HIGH", "variants_impact_MODERATE"]]
    .sum()
)


import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# --- Load both saved plot images ---
img1 = mpimg.imread("top_20_genes_high_moderate_variants11.png")   # GGAFUN_DELTA_2
img2 = mpimg.imread("FANG_top20_high_moderate_variants.png")       # FANG sample

# --- Create a figure with 2 subplots side by side ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

axes[0].imshow(img1)
axes[0].axis("off")
axes[0].set_title("GGAFUN_DELTA_2: Top 20 Genes (HIGH + MODERATE Variants)")

axes[1].imshow(img2)
axes[1].axis("off")
axes[1].set_title("FANG: Top 20 Genes (HIGH + MODERATE Variants)")

plt.tight_layout()

# --- Save combined comparison figure ---
plt.savefig("comparison_GGAFUN_vs_FANG.png", dpi=300)
plt.show()
