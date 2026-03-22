#!/usr/bin/env python3
"""
DEK-AFF2 Fusion Neoantigen Binding Prediction
==============================================
Runs MHCflurry binding predictions for all junction peptides against
a panel of ~45 common HLA Class I alleles covering >95% of the global population.

Requires: pip install mhcflurry pandas
          mhcflurry-downloads fetch models_class1_presentation

Usage: python binding_prediction.py [--input-dir ./neoantigen_output] [--output-dir ./neoantigen_output]
"""

import argparse
import sys
import os
import warnings
import re
from pathlib import Path

# Suppress TensorFlow noise
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
warnings.filterwarnings("ignore")

import pandas as pd


# Common HLA Class I alleles covering >95% of global population
# Based on allele frequency data from AFND and TCGA
COMMON_HLA_ALLELES = [
    # HLA-A (17 alleles)
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06",
    "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02",
    "HLA-A*26:01", "HLA-A*29:02", "HLA-A*30:01", "HLA-A*30:02",
    "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01",
    "HLA-A*68:02",
    # HLA-B (23 alleles)
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*13:02", "HLA-B*14:02",
    "HLA-B*15:01", "HLA-B*15:02", "HLA-B*18:01", "HLA-B*27:05",
    "HLA-B*35:01", "HLA-B*38:01", "HLA-B*39:01", "HLA-B*40:01",
    "HLA-B*40:02", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*46:01",
    "HLA-B*48:01", "HLA-B*51:01", "HLA-B*52:01", "HLA-B*53:01",
    "HLA-B*54:01", "HLA-B*57:01", "HLA-B*58:01",
    # HLA-C (14 alleles)
    "HLA-C*01:02", "HLA-C*02:02", "HLA-C*03:03", "HLA-C*03:04",
    "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*07:01",
    "HLA-C*07:02", "HLA-C*08:02", "HLA-C*12:03", "HLA-C*14:02",
    "HLA-C*15:02", "HLA-C*16:01",
]


def parse_fasta(fasta_path):
    """Parse a FASTA file and return list of (header, sequence) tuples."""
    entries = []
    header = None
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def parse_junction_peptides(fasta_path):
    """Parse junction peptides FASTA, returning only class I peptides with metadata."""
    entries = parse_fasta(fasta_path)
    peptides = []
    for header, seq in entries:
        parts = header.split()
        name = parts[0]
        if not name.startswith("class_I_"):
            continue

        # Parse both the current DEK=7aa_AFF2=1aa format and a future
        # DEK=7aa AFF2=1aa format if the FASTA is regenerated.
        junction_match = re.search(r"junction_pos=(\d+)", header)
        residue_match = re.search(r"DEK=(\d+)aa(?:_|\s+)AFF2=(\d+)aa", header)

        peptides.append({
            "name": name,
            "sequence": seq,
            "length": len(seq),
            "junction_pos": int(junction_match.group(1)) if junction_match else 0,
            "dek_residues": residue_match.group(1) if residue_match else "",
            "aff2_residues": residue_match.group(2) if residue_match else "",
        })
    return peptides


def run_predictions(peptides, alleles):
    """Run MHCflurry predictions for all peptide x allele combinations."""
    from mhcflurry import Class1AffinityPredictor

    print(f"  Loading MHCflurry models...")
    predictor = Class1AffinityPredictor.load()

    # Filter to alleles actually supported by MHCflurry
    supported = set(predictor.supported_alleles)
    valid_alleles = [a for a in alleles if a in supported]
    skipped = [a for a in alleles if a not in supported]
    if skipped:
        print(f"  Skipped {len(skipped)} unsupported alleles: {', '.join(skipped[:5])}...")
    print(f"  Running {len(peptides)} peptides x {len(valid_alleles)} alleles = {len(peptides) * len(valid_alleles)} predictions...")

    # Build matched input lists: one prediction per peptide-allele pair
    all_peptide_seqs = []
    all_alleles = []
    all_meta = []

    for pep in peptides:
        for allele in valid_alleles:
            all_peptide_seqs.append(pep["sequence"])
            all_alleles.append(allele)
            all_meta.append(pep)

    # Run predictions in one batch
    df = predictor.predict_to_dataframe(
        peptides=all_peptide_seqs,
        alleles=all_alleles,
    )

    # Merge metadata
    results = []
    for i, row in df.iterrows():
        meta = all_meta[i]
        results.append({
            "peptide": row["peptide"],
            "length": meta["length"],
            "allele": row["allele"],
            "affinity_ic50_nM": round(row["prediction"], 2),
            "affinity_percentile": round(row["prediction_percentile"], 4),
            "junction_pos": meta["junction_pos"],
            "dek_residues": meta["dek_residues"],
            "aff2_residues": meta["aff2_residues"],
            "name": meta["name"],
        })

    return pd.DataFrame(results)


def classify_binder(row):
    """Classify binding strength based on IC50 and percentile."""
    if row["affinity_ic50_nM"] < 50 and row["affinity_percentile"] < 0.5:
        return "ELITE"
    elif row["affinity_ic50_nM"] < 500 and row["affinity_percentile"] < 2.0:
        return "STRONG"
    elif row["affinity_ic50_nM"] < 5000 and row["affinity_percentile"] < 5.0:
        return "WEAK"
    else:
        return "NO_BIND"


def generate_reports(results_df, outdir):
    """Generate all output files from prediction results."""

    # 1. Full results
    full_path = outdir / "binding_results_all_alleles.csv"
    results_df.to_csv(full_path, index=False)
    print(f"  Full results: {full_path} ({len(results_df)} rows)")

    # 2. Strong binders
    results_df["binding_class"] = results_df.apply(classify_binder, axis=1)
    binders = results_df[results_df["binding_class"].isin(["ELITE", "STRONG"])].copy()
    binders = binders.sort_values("affinity_ic50_nM")
    binders_path = outdir / "strong_binders_summary.csv"
    binders.to_csv(binders_path, index=False)
    print(f"  Strong binders: {binders_path} ({len(binders)} hits)")

    # 3. Promiscuous binders — peptides that bind many alleles
    if len(binders) > 0:
        promisc = binders.groupby("peptide").agg(
            n_alleles=("allele", "nunique"),
            best_ic50=("affinity_ic50_nM", "min"),
            best_percentile=("affinity_percentile", "min"),
            alleles=("allele", lambda x: ", ".join(sorted(set(x)))),
            n_elite=("binding_class", lambda x: (x == "ELITE").sum()),
            length=("length", "first"),
            junction_pos=("junction_pos", "first"),
        ).reset_index().sort_values("n_alleles", ascending=False)
    else:
        promisc = pd.DataFrame()

    promisc_path = outdir / "promiscuous_binders.csv"
    promisc.to_csv(promisc_path, index=False)
    print(f"  Promiscuous binders: {promisc_path} ({len(promisc)} peptides)")

    # 4. Human-readable report
    report_path = outdir / "binding_report.txt"
    with open(report_path, "w") as f:
        f.write("DEK-AFF2 Fusion Neoantigen — HLA Binding Prediction Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("WHAT THIS REPORT IS\n")
        f.write("-" * 40 + "\n")
        f.write("This report predicts which of the 34 DEK-AFF2 fusion junction\n")
        f.write("peptides can bind to common HLA Class I molecules. Binding is\n")
        f.write("necessary (but not sufficient) for immune recognition.\n\n")

        f.write("HOW TO READ THE RESULTS\n")
        f.write("-" * 40 + "\n")
        f.write("- IC50 (nM): Lower = stronger binding. <50 nM is elite, <500 nM is strong.\n")
        f.write("- Percentile: Lower = better. <0.5% is elite, <2% is strong.\n")
        f.write("- Promiscuous binder: A peptide that binds many HLA types — best vaccine candidate.\n\n")

        n_total = len(results_df)
        n_elite = len(results_df[results_df["binding_class"] == "ELITE"])
        n_strong = len(results_df[results_df["binding_class"] == "STRONG"])
        n_weak = len(results_df[results_df["binding_class"] == "WEAK"])
        n_alleles = results_df["allele"].nunique()
        n_peptides = results_df["peptide"].nunique()

        f.write("SUMMARY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Peptides tested:    {n_peptides}\n")
        f.write(f"HLA alleles tested: {n_alleles}\n")
        f.write(f"Total predictions:  {n_total}\n")
        f.write(f"Elite binders:      {n_elite} (IC50 < 50 nM, percentile < 0.5%)\n")
        f.write(f"Strong binders:     {n_strong} (IC50 < 500 nM, percentile < 2%)\n")
        f.write(f"Weak binders:       {n_weak} (IC50 < 5000 nM, percentile < 5%)\n\n")

        if len(promisc) > 0:
            f.write("TOP CANDIDATES (by number of HLA alleles bound)\n")
            f.write("-" * 40 + "\n")
            for _, row in promisc.head(10).iterrows():
                f.write(f"\n  Peptide: {row['peptide']} ({row['length']}mer)\n")
                f.write(f"  Junction position: {row['junction_pos']} "
                        f"(DEK|AFF2 split at position {row['junction_pos']})\n")
                f.write(f"  Binds {row['n_alleles']} alleles "
                        f"({row['n_elite']} elite)\n")
                f.write(f"  Best IC50: {row['best_ic50']:.1f} nM "
                        f"(percentile: {row['best_percentile']:.3f}%)\n")
                f.write(f"  Alleles: {row['alleles']}\n")

        f.write("\n\nNOTE ON CLASS II PEPTIDES\n")
        f.write("-" * 40 + "\n")
        f.write("The 14 class II peptides (15mers) are NOT included in this analysis.\n")
        f.write("MHCflurry only supports Class I (HLA-A, -B, -C).\n")
        f.write("For Class II predictions (HLA-DR, -DQ, -DP), use:\n")
        f.write("  https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/\n\n")

        f.write("NEXT STEPS\n")
        f.write("-" * 40 + "\n")
        f.write("1. Get your HLA type from clinical typing or Hartwig WGS\n")
        f.write("2. Run: python lookup_hla.py --alleles HLA-A*XX:XX,HLA-B*XX:XX,...\n")
        f.write("   This instantly filters these results to YOUR alleles\n")
        f.write("3. Discuss top candidates with your oncologist\n")

    print(f"  Report: {report_path}")

    return binders, promisc


def main():
    parser = argparse.ArgumentParser(
        description="DEK-AFF2 Fusion Neoantigen Binding Prediction")
    parser.add_argument("--input-dir", type=str, default="./neoantigen_output",
                       help="Directory with junction_peptides.fasta (default: ./neoantigen_output)")
    parser.add_argument("--output-dir", type=str, default="./neoantigen_output",
                       help="Output directory (default: ./neoantigen_output)")
    args = parser.parse_args()

    indir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("DEK-AFF2 Fusion — HLA Binding Prediction")
    print("=" * 60)

    # Step 1: Load peptides
    print(f"\n[1/3] Loading junction peptides from {indir}...")
    fasta_path = indir / "junction_peptides.fasta"
    if not fasta_path.exists():
        print(f"  ERROR: {fasta_path} not found. Run neoantigen_pipeline.py first.")
        sys.exit(1)

    peptides = parse_junction_peptides(fasta_path)
    print(f"  Loaded {len(peptides)} Class I peptides (8-11mers)")

    # Step 2: Run predictions
    print(f"\n[2/3] Running binding predictions against {len(COMMON_HLA_ALLELES)} common HLA alleles...")
    results_df = run_predictions(peptides, COMMON_HLA_ALLELES)

    # Step 3: Generate reports
    print(f"\n[3/3] Generating output files...")
    binders, promisc = generate_reports(results_df, outdir)

    # Print summary to terminal
    print(f"\n{'=' * 60}")
    print(f"RESULTS SUMMARY")
    print(f"{'=' * 60}")

    n_elite = len(results_df[results_df["binding_class"] == "ELITE"])
    n_strong = len(results_df[results_df["binding_class"] == "STRONG"])

    if len(binders) > 0:
        print(f"\n  {n_elite} elite + {n_strong} strong binding predictions found!")
        print(f"\n  Top peptides by number of HLA alleles bound:")
        for _, row in promisc.head(5).iterrows():
            print(f"    {row['peptide']:15s}  binds {row['n_alleles']:2d} alleles  "
                  f"(best IC50: {row['best_ic50']:7.1f} nM)")
    else:
        print(f"\n  No strong binders found across common HLA alleles.")
        print(f"  This doesn't necessarily mean no candidates exist —")
        print(f"  your specific HLA type may still yield binders.")

    print(f"\n  Output files in: {outdir.resolve()}")
    print(f"  Next: run lookup_hla.py with your HLA alleles when available")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
