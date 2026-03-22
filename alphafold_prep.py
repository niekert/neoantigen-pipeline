#!/usr/bin/env python3
"""
DEK-AFF2 Fusion — AlphaFold 3 Input Preparation
=================================================
Generates AlphaFold 3 JSON input files for peptide-HLA complex structure prediction.
Submitting these to alphafoldserver.com gives a 3D model showing how each neoantigen
peptide sits inside its HLA binding groove, and which residues face the T-cell receptor.

Requires: pip install requests
Usage:    python alphafold_prep.py [--output-dir ./neoantigen_output/alphafold_inputs]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import requests


# ---------------------------------------------------------------------------
# Peptide candidates (top 3 from immunogenicity_report.txt)
# ---------------------------------------------------------------------------

PEPTIDE_CANDIDATES = [
    {"peptide": "EAVEKAKPR",   "hla": "HLA-A*68:01", "ic50": 30.6,  "dai": 936},
    {"peptide": "SEEEAVEKA",   "hla": "HLA-B*40:02", "ic50": 65.3,  "dai": 345},
    {"peptide": "KESEEEAV",    "hla": "HLA-B*40:02", "ic50": 95.7,  "dai": 314},
]

# ---------------------------------------------------------------------------
# Beta-2 microglobulin (B2M) — invariant light chain of all HLA class I
# Mature protein (signal peptide removed), 99 aa. UniProt P61769.
# ---------------------------------------------------------------------------

B2M_SEQUENCE = (
    "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEH"
    "SDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
)

# ---------------------------------------------------------------------------
# HLA alpha chain sequences — extracellular domain (signal peptide removed)
# Source: IMGT/HLA database. Fetched via API below; hardcoded as fallback.
#
# HLA-A*68:01: 276 aa extracellular domain
# HLA-B*40:02: 276 aa extracellular domain
# ---------------------------------------------------------------------------

HLA_FALLBACK_SEQUENCES = {
    "HLA-A*68:01": (
        "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP"
    ),
    "HLA-B*40:02": (
        "GSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAARVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEP"
    ),
}


def fetch_hla_sequence(allele: str) -> str | None:
    """Try to fetch HLA alpha chain sequence from IMGT/HLA IPD API."""
    # Encode allele name for URL (e.g. HLA-A*68:01 -> HLA-A*68:01)
    url = f"https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/{allele}"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            # Navigate the IPD response structure
            seq = (
                data.get("sequence", {})
                    .get("protein", {})
                    .get("sequence", "")
            )
            if seq and len(seq) > 50:
                # Remove signal peptide (first 24 aa for HLA-A/B)
                return seq[24:] if len(seq) > 300 else seq
    except Exception:
        pass
    return None


def get_hla_sequence(allele: str) -> tuple[str, str]:
    """Return (sequence, source) for the given HLA allele."""
    print(f"  Fetching {allele} from IMGT/HLA API...", end=" ")
    seq = fetch_hla_sequence(allele)
    if seq:
        print(f"OK ({len(seq)} aa)")
        return seq, "IMGT/HLA API"

    print("failed — using hardcoded sequence")
    if allele in HLA_FALLBACK_SEQUENCES:
        seq = HLA_FALLBACK_SEQUENCES[allele]
        return seq, "hardcoded (IMGT/HLA)"

    print(f"  ERROR: No fallback sequence for {allele}.")
    print(f"  Please download it manually from: https://www.ebi.ac.uk/ipd/imgt/hla/")
    sys.exit(1)


def make_alphafold3_json(peptide: str, hla_allele: str, hla_seq: str) -> dict:
    """Build the AlphaFold 3 server JSON input for a peptide-HLA-B2M complex.

    Format per alphafoldserver.com spec:
    https://github.com/google-deepmind/alphafold/blob/main/server/README.md
    """
    job_name = f"{peptide}_{hla_allele.replace('*', '').replace(':', '')}"

    def protein_chain(seq: str, use_template: bool = True) -> dict:
        return {"proteinChain": {
            "sequence": seq,
            "count": 1,
            "glycans": [],
            "modifications": [],
            "useStructureTemplate": use_template,
        }}

    # Must be a list even for a single job — per alphafoldserver.com spec
    return [{
        "name": job_name,
        "modelSeeds": [],
        "sequences": [
            protein_chain(hla_seq,      use_template=True),
            protein_chain(B2M_SEQUENCE, use_template=True),
            protein_chain(peptide,      use_template=False),  # novel peptide, no template
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }]


def tcr_facing_analysis(peptide: str, hla_allele: str) -> list[dict]:
    """
    Annotate each peptide residue by its role in HLA binding vs TCR recognition.

    For HLA class I, standard 9mer binding conventions:
    - P2, P9: Primary HLA anchors (buried in the groove, allele-specific pockets)
    - P1:     N-terminal, partially exposed
    - P3-P8:  Solvent-exposed face — the part the T-cell receptor actually sees
              P5 is the central TCR contact (most immunogenic position)
    - P8:     Secondary anchor for some alleles (partially buried)

    For 8mers and 10mers, positions shift slightly but P2 and P-omega remain anchors.
    """
    n = len(peptide)
    residues = []

    for i, aa in enumerate(peptide):
        pos = i + 1  # 1-indexed (P1, P2, ...)
        p_omega = n   # last position

        if pos == 2:
            role = "HLA ANCHOR"
            note = "buried in B pocket — allele-specific"
        elif pos == p_omega:
            role = "HLA ANCHOR"
            note = "buried in F pocket — C-terminal anchor"
        elif pos == 1:
            role = "partial exposure"
            note = "N-terminal, partially visible"
        elif pos == n - 1:
            role = "partial anchor"
            note = "secondary anchor for some alleles"
        elif pos == (n // 2) + 1:
            role = "TCR-FACING ★"
            note = "central position — primary TCR contact"
        else:
            role = "TCR-FACING"
            note = "solvent-exposed, visible to T-cell receptor"

        residues.append({
            "position": f"P{pos}",
            "amino_acid": aa,
            "role": role,
            "note": note,
        })

    return residues


def write_tcr_report(candidates_data: list, outdir: Path) -> None:
    """Write human-readable TCR-facing residue analysis."""
    report_path = outdir / "tcr_facing_residues.txt"
    with open(report_path, "w") as f:
        f.write("DEK-AFF2 Fusion — TCR-Facing Residue Analysis\n")
        f.write("=" * 60 + "\n\n")
        f.write("WHAT THIS MEANS\n")
        f.write("-" * 40 + "\n")
        f.write("HLA class I molecules bind peptides in a groove. The peptide\n")
        f.write("is held by anchor residues buried at each end. The middle\n")
        f.write("residues stick up like a ridge — that's what the T-cell\n")
        f.write("receptor (TCR) scans for. If those residues look foreign,\n")
        f.write("the T cell fires.\n\n")
        f.write("Anchor positions are allele-specific (different HLA alleles\n")
        f.write("prefer different amino acids at P2 and P9). TCR-facing\n")
        f.write("positions determine immunogenicity.\n\n")

        for item in candidates_data:
            cand = item["candidate"]
            analysis = item["analysis"]
            f.write("=" * 60 + "\n")
            f.write(f"Peptide:  {cand['peptide']}\n")
            f.write(f"HLA:      {cand['hla']}\n")
            f.write(f"IC50:     {cand['ic50']} nM\n")
            f.write(f"DAI:      {cand['dai']}x\n")
            f.write("-" * 40 + "\n")
            f.write(f"{'Pos':<5} {'AA':<4} {'Role':<22} {'Note'}\n")
            f.write("-" * 60 + "\n")
            for r in analysis:
                f.write(f"{r['position']:<5} {r['amino_acid']:<4} {r['role']:<22} {r['note']}\n")
            f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("NEXT STEP\n")
        f.write("-" * 40 + "\n")
        f.write("1. Upload each JSON file to: https://alphafoldserver.com\n")
        f.write("   (free, requires Google account)\n")
        f.write("2. Download the resulting PDB file (~5-15 min per job)\n")
        f.write("3. Visualize at: https://molstar.org/viewer\n")
        f.write("   - Load PDB → color by chain\n")
        f.write("   - Chain A = HLA alpha (blue), Chain B = B2M (green),\n")
        f.write("     Chain C = peptide (red)\n")
        f.write("   - Zoom into the peptide — confirm TCR-facing residues\n")
        f.write("     are solvent-exposed (sticking up out of the groove)\n")
        f.write("4. Optional: search rcsb.org for 'HLA-A*68:01' to find\n")
        f.write("   existing crystal structures to compare against\n")

    print(f"  TCR analysis: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare AlphaFold 3 input files for peptide-HLA complex prediction")
    parser.add_argument("--output-dir", type=str,
                        default="./neoantigen_output/alphafold_inputs")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("DEK-AFF2 Fusion — AlphaFold 3 Input Preparation")
    print("=" * 60)
    print(f"\nGenerating inputs for {len(PEPTIDE_CANDIDATES)} peptide-HLA complexes")
    print(f"Output directory: {outdir.resolve()}\n")

    # Cache HLA sequences (avoid fetching same allele twice)
    hla_seq_cache: dict[str, tuple[str, str]] = {}

    candidates_data = []

    for cand in PEPTIDE_CANDIDATES:
        peptide = cand["peptide"]
        allele = cand["hla"]

        print(f"\n[{peptide} / {allele}]")

        # Get HLA sequence
        if allele not in hla_seq_cache:
            hla_seq_cache[allele] = get_hla_sequence(allele)
            time.sleep(0.5)  # polite pause between API calls
        hla_seq, source = hla_seq_cache[allele]
        print(f"  HLA sequence: {len(hla_seq)} aa ({source})")
        print(f"  B2M sequence: {len(B2M_SEQUENCE)} aa (hardcoded, UniProt P61769)")
        print(f"  Peptide:      {len(peptide)} aa")

        # Generate AlphaFold 3 JSON
        af3_json = make_alphafold3_json(peptide, allele, hla_seq)
        allele_clean = allele.replace("*", "").replace(":", "")
        json_path = outdir / f"{peptide}_{allele_clean}.json"
        with open(json_path, "w") as f:
            json.dump(af3_json, f, indent=2)
        print(f"  Saved: {json_path.name}")

        # TCR-facing analysis
        analysis = tcr_facing_analysis(peptide, allele)
        candidates_data.append({"candidate": cand, "analysis": analysis})

    # Write combined TCR report
    print(f"\n[TCR-facing residue analysis]")
    write_tcr_report(candidates_data, outdir)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"DONE — {len(PEPTIDE_CANDIDATES)} AlphaFold 3 input files generated")
    print(f"{'=' * 60}")
    print(f"\nFiles in {outdir.resolve()}/:")
    for cand in PEPTIDE_CANDIDATES:
        allele_clean = cand["hla"].replace("*", "").replace(":", "")
        print(f"  {cand['peptide']}_{allele_clean}.json")
    print(f"  tcr_facing_residues.txt")
    print(f"\nNext steps:")
    print(f"  1. Go to https://alphafoldserver.com (Google account required)")
    print(f"  2. Click 'New prediction' → upload each JSON file")
    print(f"  3. Download PDB when done (~5-15 min per job)")
    print(f"  4. Visualize at https://molstar.org/viewer")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
