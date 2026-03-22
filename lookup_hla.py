#!/usr/bin/env python3
"""
HLA Lookup — Filter pre-computed binding results to your specific alleles.
==========================================================================
Run this once you have your HLA typing results.

Usage:
  python lookup_hla.py --alleles HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02,HLA-C*07:02,HLA-C*05:01
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Filter binding predictions to your HLA alleles",
        epilog="Example: python lookup_hla.py --alleles HLA-A*02:01,HLA-B*07:02,HLA-C*07:02")
    parser.add_argument("--alleles", type=str, required=True,
                       help="Comma-separated HLA alleles (e.g., HLA-A*02:01,HLA-B*07:02)")
    parser.add_argument("--input-dir", type=str, default="./neoantigen_output",
                       help="Directory with binding_results_all_alleles.csv")
    parser.add_argument("--output-dir", type=str, default="./neoantigen_output",
                       help="Output directory")
    args = parser.parse_args()

    indir = Path(args.input_dir)
    outdir = Path(args.output_dir)

    # Parse alleles
    my_alleles = [a.strip() for a in args.alleles.split(",")]
    print(f"\nYour HLA alleles: {', '.join(my_alleles)}")

    # Load pre-computed results
    results_path = indir / "binding_results_all_alleles.csv"
    if not results_path.exists():
        print(f"ERROR: {results_path} not found. Run binding_prediction.py first.")
        sys.exit(1)

    df = pd.read_csv(results_path)
    print(f"Loaded {len(df)} pre-computed predictions")

    # Filter to user's alleles
    my_results = df[df["allele"].isin(my_alleles)].copy()

    if len(my_results) == 0:
        # Try without HLA- prefix or with different formatting
        alt_alleles = []
        for a in my_alleles:
            alt_alleles.append(a)
            alt_alleles.append(a.replace("HLA-", ""))
            alt_alleles.append("HLA-" + a if not a.startswith("HLA-") else a)
        my_results = df[df["allele"].isin(alt_alleles)].copy()

    if len(my_results) == 0:
        print(f"\nNo predictions found for your alleles.")
        print(f"Available alleles in the data: {', '.join(sorted(df['allele'].unique())[:10])}...")
        sys.exit(1)

    print(f"Found {len(my_results)} predictions for your alleles")

    # Classify binders
    def classify(row):
        if row["affinity_ic50_nM"] < 50 and row["affinity_percentile"] < 0.5:
            return "ELITE"
        elif row["affinity_ic50_nM"] < 500 and row["affinity_percentile"] < 2.0:
            return "STRONG"
        elif row["affinity_ic50_nM"] < 5000 and row["affinity_percentile"] < 5.0:
            return "WEAK"
        return "NO_BIND"

    my_results["binding_class"] = my_results.apply(classify, axis=1)
    my_results = my_results.sort_values("affinity_ic50_nM")

    # Save filtered results
    my_results_path = outdir / "my_hla_results.csv"
    my_results.to_csv(my_results_path, index=False)

    # Get binders
    binders = my_results[my_results["binding_class"].isin(["ELITE", "STRONG"])]

    # Print results
    print(f"\n{'=' * 60}")
    print(f"YOUR PERSONALIZED NEOANTIGEN CANDIDATES")
    print(f"{'=' * 60}")
    print(f"\nHLA alleles: {', '.join(my_alleles)}")

    n_elite = len(my_results[my_results["binding_class"] == "ELITE"])
    n_strong = len(my_results[my_results["binding_class"] == "STRONG"])
    n_weak = len(my_results[my_results["binding_class"] == "WEAK"])

    print(f"\n  Elite binders:  {n_elite}")
    print(f"  Strong binders: {n_strong}")
    print(f"  Weak binders:   {n_weak}")

    if len(binders) > 0:
        print(f"\n  STRONG/ELITE BINDERS:")
        print(f"  {'Peptide':<15} {'Length':<7} {'HLA Allele':<15} {'IC50 (nM)':<12} {'%ile':<8} {'Class':<8}")
        print(f"  {'-' * 65}")
        for _, row in binders.iterrows():
            print(f"  {row['peptide']:<15} {row['length']:<7} {row['allele']:<15} "
                  f"{row['affinity_ic50_nM']:<12.1f} {row['affinity_percentile']:<8.3f} {row['binding_class']:<8}")

        # Best candidates: unique peptides ranked by IC50
        best = binders.drop_duplicates(subset="peptide").head(5)
        print(f"\n  TOP VACCINE CANDIDATES (unique peptides):")
        for i, (_, row) in enumerate(best.iterrows(), 1):
            print(f"  {i}. {row['peptide']} ({row['length']}mer) — "
                  f"IC50: {row['affinity_ic50_nM']:.1f} nM on {row['allele']}")
    else:
        print(f"\n  No strong binders found for your specific HLA alleles.")
        print(f"  This may mean the fusion junction peptides don't fit your HLA molecules well.")
        print(f"  Consider:")
        print(f"  - Checking Class II binding (15mers on NetMHCIIpan)")
        print(f"  - Checking point mutation neoantigens (NF1 R2183Q, ATR D626N)")

    print(f"\n  Full results saved to: {my_results_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
