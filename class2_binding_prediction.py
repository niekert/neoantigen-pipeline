#!/usr/bin/env python3
"""
DEK-AFF2 Fusion — Class II HLA Binding Prediction
===================================================
Runs HLA Class II binding predictions for junction peptides (15mers) using
the IEDB REST API (which runs NetMHCIIpan server-side).

This tests against common HLA-DRB1 alleles covering the European/Dutch population.

Note: For higher accuracy, consider running NetMHCIIpan 4.3 locally once
academic license access is obtained (requires institutional email).

Requires: pip install requests pandas
Usage: python class2_binding_prediction.py [--input-dir ./neoantigen_output]
"""

import argparse
import sys
import time
import re
from pathlib import Path

import pandas as pd
import requests


IEDB_API_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# Common HLA-DRB1 alleles in European/Dutch populations
COMMON_CLASS2_ALLELES = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*04:04",
    "HLA-DRB1*04:05",
    "HLA-DRB1*07:01",
    "HLA-DRB1*08:01",
    "HLA-DRB1*09:01",
    "HLA-DRB1*10:01",
    "HLA-DRB1*11:01",
    "HLA-DRB1*12:01",
    "HLA-DRB1*13:01",
    "HLA-DRB1*13:02",
    "HLA-DRB1*14:01",
    "HLA-DRB1*15:01",
    "HLA-DRB1*16:02",
]


def parse_fasta(fasta_path):
    """Parse FASTA file, return list of (header, sequence) tuples."""
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


def parse_class2_peptides(fasta_path):
    """Parse Class II peptides from junction_peptides.fasta."""
    entries = parse_fasta(fasta_path)
    peptides = []
    for header, seq in entries:
        parts = header.split()
        name = parts[0]
        if not name.startswith("class_II_"):
            continue

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


def predict_iedb(peptide_seq, allele, method="netmhciipan_ba", retries=5):
    """Run a single prediction via IEDB API. Returns parsed result or None."""
    for attempt in range(retries):
        try:
            resp = requests.post(
                IEDB_API_URL,
                data={
                    "method": method,
                    "sequence_text": peptide_seq,
                    "allele": allele,
                    "length": str(len(peptide_seq)),
                },
                timeout=60,
            )
            if resp.status_code == 200 and not resp.text.startswith("Error"):
                return resp.text
            elif resp.status_code == 403:
                wait_s = 5 * (attempt + 1)
                print(f"    API rate limit / 403 for {allele}; retrying in {wait_s}s")
                if attempt < retries - 1:
                    time.sleep(wait_s)
            else:
                print(f"    API error for {allele}: {resp.text[:100]}")
                if attempt < retries - 1:
                    time.sleep(2)
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2)
            else:
                print(f"    Request failed for {allele}: {e}")
    return None


def run_predictions(peptides, alleles):
    """Run IEDB API predictions for all peptide x allele combinations."""
    total = len(peptides) * len(alleles)
    print(f"  Running {len(peptides)} peptides x {len(alleles)} alleles = {total} predictions via IEDB API...")
    print(f"  (This calls the IEDB server — may take 1-2 minutes)")

    results = []
    done = 0

    for allele in alleles:
        # Submit all peptides for this allele as one sequence
        # IEDB can handle multiple peptides submitted as separate sequences
        for pep in peptides:
            raw = predict_iedb(pep["sequence"], allele)
            done += 1

            if raw is None:
                continue

            # Parse tab-separated response
            lines = raw.strip().split("\n")
            if len(lines) < 2:
                continue

            header_cols = lines[0].split("\t")
            for line in lines[1:]:
                cols = line.split("\t")
                if len(cols) < len(header_cols):
                    continue

                row = dict(zip(header_cols, cols))
                try:
                    results.append({
                        "peptide": row.get("peptide", pep["sequence"]),
                        "length": pep["length"],
                        "allele": row.get("allele", allele),
                        "core_peptide": row.get("core_peptide", ""),
                        "affinity_ic50_nM": round(float(row.get("ic50", 99999)), 2),
                        "rank": round(float(row.get("rank", 100)), 4),
                        "junction_pos": pep["junction_pos"],
                        "dek_residues": pep["dek_residues"],
                        "aff2_residues": pep["aff2_residues"],
                        "name": pep["name"],
                    })
                except (ValueError, KeyError):
                    continue

            # Brief pause to be polite to the API and reduce 403 bursts
            time.sleep(0.2)
            if done % 5 == 0:
                sys.stdout.write(f"\r  Progress: {done}/{total} predictions...")
                sys.stdout.flush()
                time.sleep(0.8)

    print(f"\r  Completed {done}/{total} predictions.          ")
    return pd.DataFrame(results)


def classify_binder_class2(row):
    """Classify Class II binding strength. Thresholds differ from Class I."""
    if row["rank"] < 2.0:
        return "STRONG"
    elif row["rank"] < 10.0:
        return "WEAK"
    else:
        return "NO_BIND"


def generate_reports(results_df, outdir):
    """Generate output files from prediction results."""

    # 1. Full results
    full_path = outdir / "class2_binding_results.csv"
    results_df.to_csv(full_path, index=False)
    print(f"  Full results: {full_path} ({len(results_df)} rows)")

    # 2. Classify and filter
    results_df["binding_class"] = results_df.apply(classify_binder_class2, axis=1)
    binders = results_df[results_df["binding_class"].isin(["STRONG"])].copy()
    binders = binders.sort_values("affinity_ic50_nM")

    binders_path = outdir / "class2_strong_binders.csv"
    binders.to_csv(binders_path, index=False)
    print(f"  Strong binders: {binders_path} ({len(binders)} hits)")

    # 3. Human-readable report
    report_path = outdir / "class2_binding_report.txt"
    with open(report_path, "w") as f:
        f.write("DEK-AFF2 Fusion — Class II HLA Binding Prediction Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("PREDICTION METHOD\n")
        f.write("-" * 40 + "\n")
        f.write("Tool: IEDB API (NetMHCIIpan binding affinity)\n")
        f.write("Note: For higher accuracy, re-run with NetMHCIIpan 4.3 locally\n")
        f.write("      once academic license access is obtained.\n\n")

        f.write("WHY CLASS II MATTERS\n")
        f.write("-" * 40 + "\n")
        f.write("Class II HLA molecules (HLA-DR, -DQ, -DP) present peptides to\n")
        f.write("CD4+ helper T cells. These helper cells are critical for mounting\n")
        f.write("and sustaining a strong CD8+ killer T cell response. A vaccine\n")
        f.write("that activates both Class I (killer) and Class II (helper)\n")
        f.write("responses is more likely to be effective.\n\n")

        f.write("HOW TO READ THESE RESULTS\n")
        f.write("-" * 40 + "\n")
        f.write("- IC50 (nM): Lower = stronger binding.\n")
        f.write("- Rank (%): Percentile rank. Lower = better.\n")
        f.write("  < 2% = strong binder, < 10% = weak binder (Class II thresholds)\n")
        f.write("- Core peptide: The 9-residue binding core within the 15mer.\n\n")

        n_total = len(results_df)
        n_strong = len(results_df[results_df["binding_class"] == "STRONG"])
        n_weak = len(results_df[results_df["binding_class"] == "WEAK"])
        n_alleles = results_df["allele"].nunique()
        n_peptides = results_df["peptide"].nunique()

        f.write("SUMMARY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Peptides tested:    {n_peptides} (15mers)\n")
        f.write(f"HLA-DRB1 alleles:   {n_alleles}\n")
        f.write(f"Total predictions:  {n_total}\n")
        f.write(f"Strong binders:     {n_strong} (rank < 2%)\n")
        f.write(f"Weak binders:       {n_weak} (rank < 10%)\n\n")

        if len(binders) > 0:
            f.write("STRONG BINDERS\n")
            f.write("-" * 40 + "\n")
            for _, row in binders.iterrows():
                f.write(f"\n  Peptide: {row['peptide']} (15mer)\n")
                f.write(f"  HLA allele: {row['allele']}\n")
                f.write(f"  IC50: {row['affinity_ic50_nM']:.1f} nM | Rank: {row['rank']:.2f}%\n")
                f.write(f"  Core: {row['core_peptide']}\n")
                f.write(f"  Junction position: {row['junction_pos']}\n")

        f.write("\n\nNEXT STEPS\n")
        f.write("-" * 40 + "\n")
        f.write("1. Get HLA-DRB1 typing (included in standard HLA typing blood draw)\n")
        f.write("2. Cross-reference your DRB1 alleles with the strong binders above\n")
        f.write("3. Consider re-running with NetMHCIIpan 4.3 locally for higher accuracy\n")
        f.write("   (requires academic license from DTU: https://services.healthtech.dtu.dk)\n")

    print(f"  Report: {report_path}")
    return binders


def main():
    parser = argparse.ArgumentParser(
        description="DEK-AFF2 Fusion — Class II HLA Binding Prediction")
    parser.add_argument("--input-dir", type=str, default="./neoantigen_output")
    parser.add_argument("--output-dir", type=str, default="./neoantigen_output")
    args = parser.parse_args()

    indir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("DEK-AFF2 Fusion — Class II HLA Binding Prediction")
    print("=" * 60)

    # Step 1: Load peptides
    print(f"\n[1/3] Loading Class II peptides from {indir}...")
    fasta_path = indir / "junction_peptides.fasta"
    if not fasta_path.exists():
        print(f"  ERROR: {fasta_path} not found. Run neoantigen_pipeline.py first.")
        sys.exit(1)

    peptides = parse_class2_peptides(fasta_path)
    print(f"  Loaded {len(peptides)} Class II peptides (15mers)")

    if len(peptides) == 0:
        print("  No Class II peptides found in FASTA file.")
        sys.exit(1)

    # Step 2: Run predictions
    print(f"\n[2/3] Running binding predictions against {len(COMMON_CLASS2_ALLELES)} common HLA-DRB1 alleles...")
    results_df = run_predictions(peptides, COMMON_CLASS2_ALLELES)

    if len(results_df) == 0:
        print("  ERROR: No results returned from IEDB API.")
        sys.exit(1)

    # Step 3: Generate reports
    print(f"\n[3/3] Generating output files...")
    binders = generate_reports(results_df, outdir)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"CLASS II RESULTS SUMMARY")
    print(f"{'=' * 60}")

    n_strong = len(results_df[results_df["binding_class"] == "STRONG"])
    if len(binders) > 0:
        print(f"\n  {n_strong} strong Class II binding predictions found!")
        print(f"\n  Strong binders:")
        for _, row in binders.iterrows():
            print(f"    {row['peptide']:20s}  {row['allele']:<20s}  "
                  f"IC50: {row['affinity_ic50_nM']:8.1f} nM  rank: {row['rank']:.2f}%")
    else:
        print(f"\n  No strong Class II binders found.")
        print(f"  Weak binders may still be relevant for CD4+ helper response.")

    print(f"\n  Output in: {outdir.resolve()}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
