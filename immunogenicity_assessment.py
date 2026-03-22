#!/usr/bin/env python3
"""
DEK-AFF2 Fusion — Immunogenicity Assessment
=============================================
Adds exploratory differential-binding analysis on top of HLA binding predictions:

1. Agretopicity (DAI): Compares fusion vs wildtype binding as a heuristic
   prioritization signal
2. Combined ranking of all candidates

Requires: pip install mhcflurry pandas (same venv312 as binding_prediction.py)
Usage: python immunogenicity_assessment.py
"""

import argparse
import sys
import os
import warnings
from pathlib import Path

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
warnings.filterwarnings("ignore")

import pandas as pd


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


def parse_wt_peptides(fasta_path):
    """Parse wildtype peptides, return dict: name -> {DEK: seq, AFF2: seq}."""
    entries = parse_fasta(fasta_path)
    wt_map = {}
    for header, seq in entries:
        # Header format: class_I_8mer_1_WT_DEK or class_I_8mer_1_WT_AFF2
        if "_WT_DEK" in header:
            base_name = header.replace("_WT_DEK", "")
            wt_map.setdefault(base_name, {})["DEK"] = seq
        elif "_WT_AFF2" in header:
            base_name = header.replace("_WT_AFF2", "")
            wt_map.setdefault(base_name, {})["AFF2"] = seq
    return wt_map


def run_wt_predictions(wt_peptides, alleles):
    """Run MHCflurry on wildtype peptides for the given alleles."""
    from mhcflurry import Class1AffinityPredictor

    print(f"  Loading MHCflurry models...")
    predictor = Class1AffinityPredictor.load()
    supported = set(predictor.supported_alleles)

    all_peptides = []
    all_alleles = []
    all_meta = []

    for (name, side, seq) in wt_peptides:
        for allele in alleles:
            if allele not in supported:
                continue
            all_peptides.append(seq)
            all_alleles.append(allele)
            all_meta.append({"name": name, "side": side, "sequence": seq})

    if not all_peptides:
        return pd.DataFrame()

    print(f"  Running {len(all_peptides)} WT predictions...")
    df = predictor.predict_to_dataframe(peptides=all_peptides, alleles=all_alleles)

    results = []
    for i, row in df.iterrows():
        meta = all_meta[i]
        results.append({
            "wt_peptide": row["peptide"],
            "wt_name": meta["name"],
            "wt_side": meta["side"],
            "allele": row["allele"],
            "wt_ic50_nM": round(row["prediction"], 2),
            "wt_percentile": round(row["prediction_percentile"], 4),
        })

    return pd.DataFrame(results)


def compute_agretopicity(strong_binders_df, wt_results_df, wt_map):
    """Compute DAI for each strong binder vs its WT counterparts."""
    results = []

    for _, binder in strong_binders_df.iterrows():
        name = binder["name"]
        allele = binder["allele"]
        fusion_ic50 = binder["affinity_ic50_nM"]

        # Find WT counterparts
        wt_data = wt_map.get(name, {})

        for side in ["DEK", "AFF2"]:
            wt_seq = wt_data.get(side)
            if not wt_seq or len(wt_seq) != binder["length"]:
                continue

            # Find WT prediction for this sequence + allele
            wt_match = wt_results_df[
                (wt_results_df["wt_name"] == name) &
                (wt_results_df["wt_side"] == side) &
                (wt_results_df["allele"] == allele)
            ]

            if wt_match.empty:
                continue

            wt_ic50 = wt_match.iloc[0]["wt_ic50_nM"]

            # DAI = WT_IC50 / Fusion_IC50
            # Higher = fusion binds much better than WT = more immunogenic
            dai = wt_ic50 / fusion_ic50 if fusion_ic50 > 0 else 0

            results.append({
                "peptide": binder["peptide"],
                "name": name,
                "length": binder["length"],
                "allele": allele,
                "fusion_ic50_nM": fusion_ic50,
                "fusion_percentile": binder["affinity_percentile"],
                "binding_class": binder["binding_class"],
                "wt_side": side,
                "wt_peptide": wt_seq,
                "wt_ic50_nM": wt_ic50,
                "dai": round(dai, 2),
                "junction_pos": binder.get("junction_pos", ""),
            })

    return pd.DataFrame(results)


def summarize_agretopicity(agretopicity_df):
    """
    Collapse DEK/AFF2-side WT comparisons to one row per peptide-allele pair.

    We keep both side-specific DAI values and use the lower available DAI as the
    conservative ranking metric. This avoids overstating agretopicity by only
    keeping whichever WT side happens to give the most favorable comparison.
    """
    if len(agretopicity_df) == 0:
        return pd.DataFrame()

    group_cols = [
        "peptide",
        "name",
        "length",
        "allele",
        "fusion_ic50_nM",
        "fusion_percentile",
        "binding_class",
        "junction_pos",
    ]

    rows = []
    for keys, group in agretopicity_df.groupby(group_cols, dropna=False):
        row = dict(zip(group_cols, keys))

        dai_values = []
        for side in ["DEK", "AFF2"]:
            side_rows = group[group["wt_side"] == side]
            if side_rows.empty:
                row[f"wt_peptide_{side.lower()}"] = ""
                row[f"wt_ic50_nM_{side.lower()}"] = None
                row[f"dai_{side.lower()}"] = None
                continue

            side_row = side_rows.iloc[0]
            row[f"wt_peptide_{side.lower()}"] = side_row["wt_peptide"]
            row[f"wt_ic50_nM_{side.lower()}"] = side_row["wt_ic50_nM"]
            row[f"dai_{side.lower()}"] = side_row["dai"]
            dai_values.append(side_row["dai"])

        row["dai_conservative"] = round(min(dai_values), 2) if dai_values else None
        row["dai_optimistic"] = round(max(dai_values), 2) if dai_values else None
        row["dai_available_sides"] = ",".join(
            side for side in ["DEK", "AFF2"] if pd.notna(row.get(f"dai_{side.lower()}"))
        )
        rows.append(row)

    return pd.DataFrame(rows)


def generate_final_ranking(agretopicity_df, outdir):
    """Generate final candidate ranking and reports."""

    # Save agretopicity results
    agret_path = outdir / "agretopicity_results.csv"
    agretopicity_df.to_csv(agret_path, index=False)
    print(f"  Agretopicity results: {agret_path} ({len(agretopicity_df)} rows)")

    best_dai = summarize_agretopicity(agretopicity_df)

    # Create composite score: prioritize by DAI, then by binding strength
    if len(best_dai) > 0:
        # Normalize conservative DAI and IC50 to 0-1 scale for ranking
        max_dai = best_dai["dai_conservative"].max()
        max_ic50 = best_dai["fusion_ic50_nM"].max()

        best_dai["dai_normalized"] = (
            best_dai["dai_conservative"] / max_dai if max_dai > 0 else 0
        )
        best_dai["binding_normalized"] = 1 - (best_dai["fusion_ic50_nM"] / max_ic50)

        # Composite: 50% conservative DAI + 50% binding strength
        best_dai["composite_score"] = (
            0.5 * best_dai["dai_normalized"] +
            0.5 * best_dai["binding_normalized"]
        ).round(4)

        best_dai = best_dai.sort_values("composite_score", ascending=False)

    # Save final ranking
    final_path = outdir / "final_candidates.csv"
    best_dai.to_csv(final_path, index=False)
    print(f"  Final candidates: {final_path} ({len(best_dai)} entries)")

    # Generate human-readable report
    report_path = outdir / "immunogenicity_report.txt"
    with open(report_path, "w") as f:
        f.write("DEK-AFF2 Fusion — Immunogenicity Assessment Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("WHAT THIS REPORT IS\n")
        f.write("-" * 40 + "\n")
        f.write("This report extends the HLA binding output with an\n")
        f.write("exploratory agretopicity (DAI) comparison against WT\n")
        f.write("DEK- and AFF2-side peptides.\n\n")

        f.write("AGRETOPICITY (DAI) — WHAT IT MEANS\n")
        f.write("-" * 40 + "\n")
        f.write("DAI = IC50_wildtype / IC50_fusion\n\n")
        f.write("- DAI > 1: The fusion peptide binds HLA better than the\n")
        f.write("  chosen wildtype comparator for this allele. That can be\n")
        f.write("  useful for prioritization, but it does not prove\n")
        f.write("  presentation or T-cell recognition.\n")
        f.write("- DAI > 5: Larger differential vs the chosen comparator.\n")
        f.write("- DAI > 10: Very large differential vs the chosen comparator.\n")
        f.write("- DAI < 1: The chosen wildtype comparator binds similarly or\n")
        f.write("  better, which makes the fusion peptide less distinctive by\n")
        f.write("  this heuristic.\n\n")
        f.write("For fusion neoantigens, high DAI can occur because the\n")
        f.write("junction sequence differs substantially from each native side.\n")
        f.write("That still does not model processing, abundance, or TCR\n")
        f.write("recognition.\n\n")
        f.write("This report shows BOTH DEK-side and AFF2-side WT comparisons.\n")
        f.write("For ranking, it uses the LOWER available DAI as a conservative\n")
        f.write("summary rather than keeping only the most favorable value.\n\n")

        if len(best_dai) > 0:
            f.write("RANKED CANDIDATES\n")
            f.write("-" * 40 + "\n")
            f.write(
                f"{'Rank':<5} {'Peptide':<16} {'Allele':<15} {'IC50(nM)':<10} "
                f"{'Cons.DAI':<10} {'Max DAI':<10} {'Score':<8} {'Class':<8}\n"
            )
            f.write("-" * 88 + "\n")

            for rank, (_, row) in enumerate(best_dai.iterrows(), 1):
                dai_flag = ""
                if row["dai_conservative"] >= 10:
                    dai_flag = " ***"
                elif row["dai_conservative"] >= 5:
                    dai_flag = " **"
                elif row["dai_conservative"] >= 2:
                    dai_flag = " *"

                f.write(
                    f"{rank:<5} {row['peptide']:<16} {row['allele']:<15} "
                    f"{row['fusion_ic50_nM']:<10.1f} {row['dai_conservative']:<10.1f} "
                    f"{row['dai_optimistic']:<10.1f} {row['composite_score']:<8.4f} "
                    f"{row['binding_class']:<8}{dai_flag}\n"
                )

            f.write("\n  *** = conservative DAI >= 10 (high differential)\n")
            f.write("  **  = conservative DAI >= 5 (moderate-high differential)\n")
            f.write("  *   = conservative DAI >= 2 (modest differential)\n")

            # Highlight top 5 unique peptides
            top_peptides = best_dai.drop_duplicates(subset="peptide").head(5)
            f.write(f"\n\nTOP 5 WORKING CANDIDATES\n")
            f.write("-" * 40 + "\n")

            for rank, (_, row) in enumerate(top_peptides.iterrows(), 1):
                f.write(f"\n  #{rank}: {row['peptide']} ({row['length']}mer)\n")
                f.write(f"      Best allele: {row['allele']}\n")
                f.write(f"      Binding: {row['fusion_ic50_nM']:.1f} nM "
                        f"({row['binding_class']})\n")
                f.write(
                    f"      Conservative DAI: {row['dai_conservative']:.1f} "
                    f"(max observed DAI: {row['dai_optimistic']:.1f})"
                )
                if row["dai_conservative"] >= 10:
                    f.write(" — HIGH DIFFERENTIAL\n")
                elif row["dai_conservative"] >= 5:
                    f.write(" — MODERATE-HIGH DIFFERENTIAL\n")
                elif row["dai_conservative"] >= 2:
                    f.write(" — MODEST DIFFERENTIAL\n")
                else:
                    f.write(" — LOW DIFFERENTIAL / SELF-LIKE ON AT LEAST ONE SIDE\n")
                if pd.notna(row.get("dai_dek")):
                    f.write(
                        f"      DEK-side WT: {row['wt_peptide_dek']} | "
                        f"WT IC50: {row['wt_ic50_nM_dek']:.1f} nM | "
                        f"DAI: {row['dai_dek']:.1f}\n"
                    )
                if pd.notna(row.get("dai_aff2")):
                    f.write(
                        f"      AFF2-side WT: {row['wt_peptide_aff2']} | "
                        f"WT IC50: {row['wt_ic50_nM_aff2']:.1f} nM | "
                        f"DAI: {row['dai_aff2']:.1f}\n"
                    )
                f.write(f"      Composite score: {row['composite_score']:.4f}\n")

        f.write(f"\n\nADDITIONAL VALIDATION (MANUAL)\n")
        f.write("-" * 40 + "\n")
        f.write("For further confidence, submit top candidates to:\n\n")
        f.write("1. NetCTLpan 1.1 (proteasomal cleavage + TAP transport + binding):\n")
        f.write("   https://services.healthtech.dtu.dk/services/NetCTLpan-1.1/\n")
        f.write("   Input: peptide sequences + HLA alleles\n\n")
        f.write("2. IEDB Immunogenicity Predictor:\n")
        f.write("   http://tools.iedb.org/immunogenicity/\n")
        f.write("   Input: peptide sequences\n\n")
        f.write("3. DeepImmuno (web interface):\n")
        f.write("   https://deepimmuno.research.cchmc.org/\n")
        f.write("   Input: peptide + HLA allele (9-10mers only)\n")

    print(f"  Report: {report_path}")

    # Generate NetCTLpan input file
    netctlpan_path = outdir / "netctlpan_input.txt"
    with open(netctlpan_path, "w") as f:
        f.write("# Peptides for NetCTLpan 1.1 submission\n")
        f.write("# Submit at: https://services.healthtech.dtu.dk/services/NetCTLpan-1.1/\n")
        f.write("# Paste the sequences below into the input box\n")
        f.write("# Select relevant HLA alleles\n\n")
        if len(best_dai) > 0:
            seen = set()
            for _, row in best_dai.iterrows():
                if row["peptide"] not in seen:
                    f.write(f">{row['peptide']}_{row['length']}mer\n")
                    f.write(f"{row['peptide']}\n")
                    seen.add(row["peptide"])
    print(f"  NetCTLpan input: {netctlpan_path}")

    return best_dai


def main():
    parser = argparse.ArgumentParser(
        description="DEK-AFF2 Fusion — Immunogenicity Assessment")
    parser.add_argument("--input-dir", type=str, default="./neoantigen_output")
    parser.add_argument("--output-dir", type=str, default="./neoantigen_output")
    args = parser.parse_args()

    indir = Path(args.input_dir)
    outdir = Path(args.output_dir)

    print("=" * 60)
    print("DEK-AFF2 Fusion — Immunogenicity Assessment")
    print("=" * 60)

    # Step 1: Load strong binders
    print(f"\n[1/4] Loading strong binders...")
    binders_path = indir / "strong_binders_summary.csv"
    if not binders_path.exists():
        print(f"  ERROR: {binders_path} not found. Run binding_prediction.py first.")
        sys.exit(1)

    strong_binders = pd.read_csv(binders_path)
    print(f"  Loaded {len(strong_binders)} strong/elite binding predictions")

    # Step 2: Load wildtype peptides
    print(f"\n[2/4] Loading wildtype peptides...")
    wt_fasta = indir / "wildtype_peptides.fasta"
    if not wt_fasta.exists():
        print(f"  ERROR: {wt_fasta} not found.")
        sys.exit(1)

    wt_map = parse_wt_peptides(wt_fasta)
    print(f"  Loaded WT counterparts for {len(wt_map)} peptides")

    # Step 3: Run MHCflurry on WT peptides
    print(f"\n[3/4] Running MHCflurry on wildtype peptides...")

    # Collect unique (name, side, allele) combinations needed
    wt_to_predict = []
    alleles_needed = set()
    for _, binder in strong_binders.iterrows():
        name = binder["name"]
        allele = binder["allele"]
        alleles_needed.add(allele)
        wt_data = wt_map.get(name, {})
        for side in ["DEK", "AFF2"]:
            if side in wt_data and len(wt_data[side]) == binder["length"]:
                wt_to_predict.append((name, side, wt_data[side]))

    # Deduplicate
    wt_unique = list(set(wt_to_predict))
    print(f"  Need predictions for {len(wt_unique)} unique WT peptide-side pairs "
          f"across {len(alleles_needed)} alleles")

    wt_results = run_wt_predictions(wt_unique, list(alleles_needed))
    print(f"  Got {len(wt_results)} WT predictions")

    # Step 4: Compute agretopicity and generate ranking
    print(f"\n[4/4] Computing agretopicity and generating final ranking...")
    agretopicity = compute_agretopicity(strong_binders, wt_results, wt_map)

    if len(agretopicity) == 0:
        print("  WARNING: No agretopicity results computed.")
    else:
        print(f"  Computed DAI for {len(agretopicity)} peptide-allele-WT combinations")
        print(f"  DAI range: {agretopicity['dai'].min():.1f} - {agretopicity['dai'].max():.1f}")

    final = generate_final_ranking(agretopicity, outdir)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"IMMUNOGENICITY ASSESSMENT SUMMARY")
    print(f"{'=' * 60}")

    if len(final) > 0:
        top5 = final.drop_duplicates(subset="peptide").head(5)
        print(f"\n  Top candidates (by composite score):")
        for rank, (_, row) in enumerate(top5.iterrows(), 1):
            dai_label = "HIGH_DIFF" if row["dai_conservative"] >= 10 else \
                       "MID_DIFF" if row["dai_conservative"] >= 5 else \
                       "MODEST" if row["dai_conservative"] >= 2 else "LOW"
            print(
                f"    #{rank} {row['peptide']:16s} IC50: {row['fusion_ic50_nM']:7.1f} nM  "
                f"cons.DAI: {row['dai_conservative']:6.1f} ({dai_label})  on {row['allele']}"
            )

    print(f"\n  Output in: {outdir.resolve()}")
    print(f"  For additional validation, see netctlpan_input.txt")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
