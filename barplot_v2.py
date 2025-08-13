#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import argparse
import sys
import numpy as np
import plotly.io as pio

# === PARSE ARGUMENTS ===
parser = argparse.ArgumentParser(description="Create an interactive relative abundance barplot from an OTU table (tab-delimited, requires 'Taxonomy' and 'Total' columns).")
parser.add_argument("-i", "--input", required=True, help="Input OTU table (tab-delimited)")
parser.add_argument("-o", "--output", required=True, help="Output HTML filename")
parser.add_argument("-id", "--id_column", required=True, choices=["OTU_ID", "OTU_assignment", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                    help="Taxon label: OTU_ID, OTU_assignment or taxonomy level")
parser.add_argument("-m", "--metadata", required=False, help="Optional sample metadata file (tab-delimited)")
parser.add_argument("-g", "--group_column", required=False, help="Column in metadata to use for sample grouping (facets/colors)")
parser.add_argument("-t", "--threshold", type=float, default=0.0,
                    help="Relative abundance threshold (%) below which taxa are grouped as 'Other'")
parser.add_argument("-n", "--top_n", type=int, default=None,
                    help="Keep only the top N most abundant taxa (others grouped as 'Other')")
parser.add_argument("-e", "--export", help="Optional static image export (e.g., .png, .pdf)")

args = parser.parse_args()

# === LOAD OTU TABLE ===
try:
    df = pd.read_csv(args.input, sep="\t")
except Exception as e:
    print(f"❌ Failed to read input file: {e}")
    sys.exit(1)

# === VERIFY TAXONOMY COLUMN ===
if "Taxonomy" not in df.columns:
    print("❌ 'Taxonomy' column is required.")
    sys.exit(1)

# === SPLIT TAXONOMY ===
taxonomy_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
df[taxonomy_levels] = df["Taxonomy"].str.split(",", expand=True)

# === IDENTIFY SAMPLE COLUMNS ===
if "Total" not in df.columns:
    print("❌ Missing 'Total' column.")
    sys.exit(1)
sample_cols = df.columns[df.columns.get_loc("Total") + 1:]

# === SET TAXON LABEL ===
df["Taxon"] = df[args.id_column]

# === COMPUTE RELATIVE ABUNDANCE ===
df_rel = df[sample_cols].div(df[sample_cols].sum(), axis=1) * 100  # percent
df_rel["Taxon"] = df["Taxon"]

# === MEAN ABUNDANCE PER TAXON ===
mean_abundance = df_rel.groupby("Taxon")[sample_cols].mean().mean(axis=1)

# === SELECT TAXA TO KEEP ===
keep_taxa = set(mean_abundance.index)
if args.top_n is not None:
    top_taxa = mean_abundance.nlargest(args.top_n).index
    keep_taxa = keep_taxa.intersection(top_taxa)

if args.threshold > 0.0:
    above_threshold_taxa = mean_abundance[mean_abundance >= args.threshold].index
    keep_taxa = keep_taxa.intersection(above_threshold_taxa)

# === LABEL LOW ABUNDANCE TAXA AS 'Other' ===
df_rel["Taxon"] = df_rel["Taxon"].apply(lambda x: x if x in keep_taxa else "Other")

# === AGGREGATE ===
df_grouped = df_rel.groupby("Taxon")[sample_cols].sum().reset_index()
df_long = df_grouped.melt(id_vars="Taxon", var_name="Sample", value_name="Relative Abundance")

# === OPTIONAL METADATA ===
if args.metadata and args.group_column:
    try:
        metadata = pd.read_csv(args.metadata, sep="\t")
        df_long = df_long.merge(metadata[[args.group_column, "Sample"]], on="Sample", how="left")
    except Exception as e:
        print(f"⚠️ Warning: Metadata merge failed: {e}")
        df_long[args.group_column] = "Unknown"
else:
    df_long[args.group_column] = "All Samples"

# === PLOT ===
fig = px.bar(
    df_long,
    x="Sample",
    y="Relative Abundance",
    color="Taxon",
    facet_col=args.group_column if args.group_column else None,
    title=f"Relative Abundance by {args.id_column} (Top-N: {args.top_n}, Threshold: {args.threshold}%)",
    labels={"Relative Abundance": "Relative Abundance (%)", "Sample": "Sample"},
    height=600
)

fig.update_layout(barmode="stack", xaxis_tickangle=-45)

# === SAVE HTML ===
fig.write_html(args.output)
print(f"✅ Interactive HTML plot saved to: {args.output}")

# === OPTIONAL STATIC EXPORT ===
if args.export:
    try:
        import plotly.io as pio
        fig.write_image(args.export)
        print(f"📦 Static image exported to: {args.export}")
    except Exception as e:
        print(f"⚠️ Failed to export static image: {e}")
        print("💡 Try installing kaleido: pip install -U kaleido")
