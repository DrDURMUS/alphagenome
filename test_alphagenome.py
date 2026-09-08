"""
AlphaGenome Test & Usage Examples
=================================
Tests API connectivity and demonstrates key AlphaGenome features.

Usage:
    conda activate alphagenome
    python test_alphagenome.py
"""

import os

import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.variant_scorers import (
    get_recommended_scorers,
    tidy_scores,
    GeneMaskLFCScorer,
)

# ─── API Key ─────────────────────────────────────────────────────────────────
API_KEY = os.environ.get("ALPHAGENOME_API_KEY")
if not API_KEY:
    print("⚠️  ALPHAGENOME_API_KEY ortam değişkeni tanımlı değil.")
    print("   Testi çalıştırmak için: export ALPHAGENOME_API_KEY=\"your_key\"")
    import sys
    sys.exit(0)

# ─── Create the model client ────────────────────────────────────────────────
print("=" * 60)
print("AlphaGenome - Connection Test")
print("=" * 60)

model = dna_client.create(API_KEY)
print("✅ Model client created successfully!")

# Common region: APOL4 gene on chr22 (~1Mb window)
interval = genome.Interval(
    chromosome="chr22",
    start=35677410,
    end=36725986,
)

# ─── Example 1: Forward Prediction (predict_interval) ───────────────────────
print("\n" + "=" * 60)
print("Example 1: Forward Prediction (predict_interval - RNA-seq)")
print("=" * 60)
print(f"  Region: {interval.chromosome}:{interval.start}-{interval.end}")
print(f"  Length: {interval.end - interval.start:,} bp")

outputs = model.predict_interval(
    interval=interval,
    ontology_terms=["UBERON:0001157"],  # colon tissue
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
)
rna_seq = outputs.rna_seq
print(f"  ✅ RNA-seq prediction successful!")
print(f"     Output interval: {rna_seq.interval}")
print(f"     Values shape: {rna_seq.values.shape}")
print(f"     Value range: [{rna_seq.values.min():.4f}, {rna_seq.values.max():.4f}]")

# ─── Example 2: Variant Effect Prediction ───────────────────────────────────
print("\n" + "=" * 60)
print("Example 2: Variant Effect Prediction (predict_variant)")
print("=" * 60)

variant = genome.Variant(
    chromosome="chr22",
    position=36201698,
    reference_bases="A",
    alternate_bases="C",
)
print(f"  Variant: {variant.chromosome}:{variant.position} "
      f"{variant.reference_bases}>{variant.alternate_bases}")

variant_outputs = model.predict_variant(
    interval=interval,
    variant=variant,
    ontology_terms=["UBERON:0001157"],
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
)
ref_rna = variant_outputs.reference.rna_seq
alt_rna = variant_outputs.alternate.rna_seq
diff = alt_rna.values - ref_rna.values
print(f"  ✅ Variant prediction successful!")
print(f"     REF range: [{ref_rna.values.min():.4f}, {ref_rna.values.max():.4f}]")
print(f"     ALT range: [{alt_rna.values.min():.4f}, {alt_rna.values.max():.4f}]")
print(f"     Max |effect|: {np.abs(diff).max():.6f}")
print(f"     Mean |effect|: {np.abs(diff).mean():.6f}")

# ─── Example 3: Variant Scoring (score_variant) ─────────────────────────────
print("\n" + "=" * 60)
print("Example 3: Variant Scoring with Recommended Scorers")
print("=" * 60)

from alphagenome.models.dna_model import Organism

recommended_scorers = get_recommended_scorers(Organism.HOMO_SAPIENS)
print(f"  Using {len(recommended_scorers)} recommended scorers")

scores = model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=recommended_scorers,
)
df = tidy_scores(scores)
if df is not None:
    print(f"  ✅ Variant scoring successful!")
    print(f"     Score table shape: {df.shape}")
    print(f"     Columns: {list(df.columns[:5])}...")
    print(f"\n  Top scores (by absolute value):")
    if "score" in df.columns:
        top = df.reindex(df["score"].abs().sort_values(ascending=False).index).head(5)
        print(top.to_string(index=False, max_colwidth=30))
else:
    print("  ⚠️  No scores returned (variant may not overlap scored features)")

# ─── Example 4: Multiple Modalities ─────────────────────────────────────────
print("\n" + "=" * 60)
print("Example 4: Multi-modality (RNA-seq + CAGE + ATAC)")
print("=" * 60)

multi_outputs = model.predict_interval(
    interval=interval,
    ontology_terms=["UBERON:0001157"],
    requested_outputs=[
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.CAGE,
        dna_client.OutputType.ATAC,
    ],
)
print(f"  ✅ Multi-modality prediction successful!")
if multi_outputs.rna_seq is not None:
    print(f"     RNA-seq shape: {multi_outputs.rna_seq.values.shape}")
if multi_outputs.cage is not None:
    print(f"     CAGE shape: {multi_outputs.cage.values.shape}")
if multi_outputs.atac is not None:
    print(f"     ATAC shape: {multi_outputs.atac.values.shape}")

# ─── Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("✅ All 4 examples completed successfully!")
print("=" * 60)
