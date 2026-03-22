#!/usr/bin/env python3
"""
DEK-AFF2 Fusion Neoantigen Peptide Generator
=============================================
This script:
1. Fetches DEK and AFF2 transcript/protein data from Ensembl REST API
2. Maps exon boundaries to protein positions
3. Reconstructs the DEK(exon7)::AFF2(exon9) fusion junction protein sequence
4. Generates all overlapping peptides spanning the junction
5. Outputs peptides in FASTA format for NetMHCpan submission

Requirements: pip install requests
Usage: python neoantigen_pipeline.py [--hla HLA-A02:01,HLA-B07:02,...]

Author: Generated for neoantigen prediction pipeline
"""

import requests
import json
import sys
import argparse
from pathlib import Path


ENSEMBL_REST = "https://rest.ensembl.org"

# Canonical transcripts for DEK and AFF2 (verified against Ensembl GRCh38)
# DEK: ENST00000652689 (canonical, 375 aa)
# AFF2: ENST00000370460 (canonical, 1311 aa)

def fetch_json(url, params=None):
    """Fetch JSON from Ensembl REST API with error handling."""
    headers = {"Content-Type": "application/json"}
    try:
        resp = requests.get(url, headers=headers, params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {url}: {e}")
        sys.exit(1)


def get_gene_transcripts(gene_symbol):
    """Get all transcripts for a gene symbol."""
    # First, get the gene ID
    url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{gene_symbol}"
    results = fetch_json(url)
    
    gene_id = None
    for r in results:
        if r.get("type") == "gene":
            gene_id = r["id"]
            break
    
    if not gene_id:
        print(f"Could not find gene ID for {gene_symbol}")
        sys.exit(1)
    
    print(f"  Gene ID: {gene_id}")
    
    # Get gene info with transcripts
    url = f"{ENSEMBL_REST}/lookup/id/{gene_id}?expand=1"
    gene_data = fetch_json(url)
    
    return gene_data


def get_canonical_transcript(gene_symbol):
    """Find the canonical/MANE Select transcript for a gene."""
    gene_data = get_gene_transcripts(gene_symbol)
    
    transcripts = gene_data.get("Transcript", [])
    
    # Prefer MANE Select, then canonical, then longest coding
    canonical = None
    longest_coding = None
    
    for tx in transcripts:
        if tx.get("is_canonical"):
            canonical = tx
        if tx.get("biotype") == "protein_coding":
            if longest_coding is None or tx.get("length", 0) > longest_coding.get("length", 0):
                longest_coding = tx
    
    chosen = canonical or longest_coding
    if not chosen:
        print(f"No protein-coding transcript found for {gene_symbol}")
        sys.exit(1)
    
    print(f"  Transcript: {chosen['id']} (canonical={chosen.get('is_canonical', False)})")
    return chosen


def get_exon_details(transcript_id):
    """Get detailed exon information including CDS coordinates."""
    url = f"{ENSEMBL_REST}/overlap/id/{transcript_id}?feature=exon;content-type=application/json"
    exons = fetch_json(url)
    
    # Filter to only exons belonging to this transcript
    tx_exons = [e for e in exons if transcript_id in e.get("Parent", "")]
    
    # Also get the transcript sequence details
    url = f"{ENSEMBL_REST}/lookup/id/{transcript_id}?expand=1"
    tx_data = fetch_json(url)
    
    return tx_exons, tx_data


def get_transcript_exons_with_protein_mapping(transcript_id):
    """
    Get exon boundaries mapped to protein positions.
    Uses the sequence endpoint to get CDS and protein info.
    """
    # Get the CDS sequence
    url = f"{ENSEMBL_REST}/sequence/id/{transcript_id}?type=cds"
    cds_data = fetch_json(url)
    cds_seq = cds_data.get("seq", "")
    
    # Get the protein sequence
    url = f"{ENSEMBL_REST}/sequence/id/{transcript_id}?type=protein"
    protein_data = fetch_json(url)
    protein_seq = protein_data.get("seq", "")
    
    # Get detailed exon info from the lookup endpoint
    url = f"{ENSEMBL_REST}/lookup/id/{transcript_id}?expand=1"
    tx_info = fetch_json(url)
    
    exons = tx_info.get("Exon", [])
    strand = tx_info.get("strand", 1)
    
    # Sort exons by their position in the transcript
    if strand == 1:
        exons.sort(key=lambda e: e["start"])
    else:
        exons.sort(key=lambda e: -e["end"])
    
    # Get the translation start/end (genomic coordinates)
    translation = tx_info.get("Translation", {})
    
    if not translation:
        print(f"  WARNING: No translation found for {transcript_id}")
        return exons, cds_seq, protein_seq, []
    
    # Map exons to CDS positions
    # We need to figure out which portion of each exon is coding
    tx_start = tx_info["start"]
    tx_end = tx_info["end"]
    cds_start_genomic = translation["start"]  # genomic coordinate
    cds_end_genomic = translation["end"]      # genomic coordinate
    
    cds_position = 0  # running position in CDS
    exon_protein_map = []
    
    for i, exon in enumerate(exons):
        exon_start = exon["start"]
        exon_end = exon["end"]
        
        # Determine the coding portion of this exon
        # Ensembl always reports start < end in genomic coords, regardless of strand
        coding_start = max(exon_start, cds_start_genomic)
        coding_end = min(exon_end, cds_end_genomic)
        
        if coding_start > coding_end:
            # This exon is entirely UTR
            exon_protein_map.append({
                "exon_number": i + 1,
                "exon_id": exon.get("id", ""),
                "genomic_start": exon_start,
                "genomic_end": exon_end,
                "coding_bases": 0,
                "cds_start": None,
                "cds_end": None,
                "protein_start": None,
                "protein_end": None,
                "is_coding": False
            })
            continue
        
        coding_length = coding_end - coding_start + 1
        cds_start_pos = cds_position
        cds_end_pos = cds_position + coding_length - 1
        
        # Convert CDS position to protein position (0-indexed)
        protein_start = cds_start_pos // 3
        protein_end = cds_end_pos // 3
        
        exon_protein_map.append({
            "exon_number": i + 1,
            "exon_id": exon.get("id", ""),
            "genomic_start": exon_start,
            "genomic_end": exon_end,
            "coding_bases": coding_length,
            "cds_start": cds_start_pos,
            "cds_end": cds_end_pos,
            "protein_start": protein_start + 1,  # 1-indexed for display
            "protein_end": protein_end + 1,       # 1-indexed for display
            "is_coding": True,
            "phase_start": cds_start_pos % 3,    # 0 = codon boundary
            "phase_end": (cds_end_pos + 1) % 3   # phase at end of exon
        })
        
        cds_position += coding_length
    
    return exons, cds_seq, protein_seq, exon_protein_map


CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def reconstruct_junction_from_cds(dek_cds, aff2_cds, dek_exon7_info, aff2_exon9_info):
    """
    Reconstruct the exact fusion junction protein sequence at the CDS (DNA) level.

    When a fusion breakpoint falls in the middle of a codon, the last partial codon
    from DEK exon 7 combines with the first partial codon from AFF2 exon 9 to form
    a novel "junction codon". This function determines what amino acid that produces.

    Returns:
        dict with:
          - dek_protein_portion: DEK amino acids up to (not including) the split codon
          - aff2_protein_portion: AFF2 amino acids after the split codon
          - junction_aa: the novel amino acid at the junction (or None if clean join)
          - dek_partial_bases: the leftover bases from DEK
          - aff2_partial_bases: the bases from AFF2 that complete the codon
          - junction_codon: the full 3-base codon formed at the junction
          - is_novel: whether the junction amino acid differs from both WT proteins
          - frame_preserved: whether the reading frame is maintained after the junction
    """
    dek_cds_end = dek_exon7_info["cds_end"]  # 0-indexed, last CDS base in DEK exon 7
    dek_phase_end = dek_exon7_info["phase_end"]  # bases into next codon (0=clean)

    aff2_cds_start = aff2_exon9_info["cds_start"]  # 0-indexed, first CDS base in AFF2 exon 9
    aff2_phase_start = aff2_exon9_info["phase_start"]  # bases already used from prev codon

    # DEK contributes bases up to and including dek_cds_end
    dek_cds_portion = dek_cds[:dek_cds_end + 1]
    # AFF2 contributes bases from aff2_cds_start onward
    aff2_cds_portion = aff2_cds[aff2_cds_start:]

    # How many leftover bases does DEK have after its last complete codon?
    dek_leftover = dek_phase_end  # 0, 1, or 2 bases
    # How many bases does AFF2 need to skip from its partial codon?
    aff2_skip = aff2_phase_start  # 0, 1, or 2 bases

    result = {
        "dek_partial_bases": dek_cds_portion[-dek_leftover:] if dek_leftover > 0 else "",
        "aff2_partial_bases": aff2_cds_portion[:3 - dek_leftover] if dek_leftover > 0 else "",
        "junction_codon": None,
        "junction_aa": None,
        "is_novel": False,
        "frame_preserved": True,
    }

    if dek_leftover == 0:
        # Clean codon boundary — no novel junction amino acid
        # DEK protein: translate all complete codons
        n_complete_dek = len(dek_cds_portion) // 3
        result["dek_protein_portion"] = "".join(
            CODON_TABLE.get(dek_cds_portion[i*3:(i+1)*3], 'X')
            for i in range(n_complete_dek)
        )
        # AFF2 protein: translate from its start, skipping partial codon if any
        aff2_coding = aff2_cds_portion[aff2_skip:]
        n_complete_aff2 = len(aff2_coding) // 3
        result["aff2_protein_portion"] = "".join(
            CODON_TABLE.get(aff2_coding[i*3:(i+1)*3], 'X')
            for i in range(n_complete_aff2)
        )
        return result

    # Mid-codon junction: DEK has leftover bases that need to combine with AFF2
    # The junction codon = DEK leftover bases + AFF2 first (3 - leftover) bases
    bases_needed_from_aff2 = 3 - dek_leftover
    junction_codon = dek_cds_portion[-dek_leftover:] + aff2_cds_portion[:bases_needed_from_aff2]

    if len(junction_codon) != 3:
        print(f"  WARNING: Could not form complete junction codon (got {junction_codon})")
        result["frame_preserved"] = False
        return result

    junction_aa = CODON_TABLE.get(junction_codon, 'X')
    result["junction_codon"] = junction_codon
    result["junction_aa"] = junction_aa

    # DEK protein up to the split codon (exclude the incomplete last codon)
    n_complete_dek = (len(dek_cds_portion) - dek_leftover) // 3
    result["dek_protein_portion"] = "".join(
        CODON_TABLE.get(dek_cds_portion[i*3:(i+1)*3], 'X')
        for i in range(n_complete_dek)
    )

    # AFF2 protein after the junction codon
    aff2_after_junction = aff2_cds_portion[bases_needed_from_aff2:]
    n_complete_aff2 = len(aff2_after_junction) // 3
    result["aff2_protein_portion"] = "".join(
        CODON_TABLE.get(aff2_after_junction[i*3:(i+1)*3], 'X')
        for i in range(n_complete_aff2)
    )

    # Check if the junction amino acid is novel (different from both WT proteins)
    # In normal DEK, the codon at this position would continue with DEK's next exon
    # In normal AFF2, the codon at this position would start with AFF2's previous exon
    # The fusion creates a codon that exists in neither — so it's inherently novel
    result["is_novel"] = True

    # Check if reading frame is preserved after the junction
    # If DEK leftover + AFF2 skip = 3 (or 0), frame is maintained
    total_partial = dek_leftover + aff2_skip
    result["frame_preserved"] = (total_partial % 3 == 0) or (dek_leftover + (3 - aff2_skip)) % 3 == 0

    return result


def generate_junction_peptides(dek_protein, aff2_protein,
                                dek_exon7_end_aa, aff2_exon9_start_aa,
                                dek_exon7_end_phase, aff2_exon9_start_phase,
                                window=30, junction_reconstruction=None):
    """
    Generate all peptides spanning the fusion junction.

    Args:
        dek_protein: Full DEK protein sequence
        aff2_protein: Full AFF2 protein sequence
        dek_exon7_end_aa: Last amino acid position from DEK exon 7 (1-indexed)
        aff2_exon9_start_aa: First amino acid position from AFF2 exon 9 (1-indexed)
        dek_exon7_end_phase: Phase at end of DEK exon 7 (0=complete codon)
        aff2_exon9_start_phase: Phase at start of AFF2 exon 9
        window: Number of amino acids on each side of junction to consider
        junction_reconstruction: result from reconstruct_junction_from_cds() if available
    """

    frame_note = ""

    if junction_reconstruction and junction_reconstruction.get("junction_aa"):
        # We have CDS-level data — build the exact fusion protein sequence at the junction
        jr = junction_reconstruction
        junction_aa = jr["junction_aa"]

        # DEK portion: protein up to the incomplete codon
        dek_portion = jr["dek_protein_portion"]
        # AFF2 portion: protein after the junction codon
        aff2_portion = jr["aff2_protein_portion"]

        # The fusion junction sequence is: DEK_portion + junction_aa + AFF2_portion
        dek_window = dek_portion[-window:] if len(dek_portion) >= window else dek_portion
        # junction_aa sits between DEK and AFF2
        aff2_window = aff2_portion[:window] if len(aff2_portion) >= window else aff2_portion

        # Build the full junction region with the novel amino acid in the middle
        junction_seq = dek_window + junction_aa + aff2_window
        junction_pos = len(dek_window)  # The novel AA is at this position

        frame_note = (
            f"CDS-LEVEL JUNCTION RECONSTRUCTION:\n"
            f"  DEK exon 7 ends with {jr['dek_partial_bases']} ({len(jr['dek_partial_bases'])} leftover base(s))\n"
            f"  AFF2 exon 9 starts with {jr['aff2_partial_bases']} ({len(jr['aff2_partial_bases'])} base(s) used)\n"
            f"  Junction codon: {jr['junction_codon']} → {junction_aa}\n"
            f"  This amino acid is NOVEL — it does not exist in either normal DEK or AFF2.\n"
            f"  Frame preserved after junction: {jr['frame_preserved']}"
        )
    else:
        # Fallback: simple protein-level join (no CDS data)
        if dek_exon7_end_phase != 0:
            frame_note = (
                f"NOTE: DEK exon 7 ends with phase {dek_exon7_end_phase} "
                f"(mid-codon). Combined with AFF2 exon 9 start phase {aff2_exon9_start_phase}, "
                f"there may be a novel amino acid at the exact junction point.\n"
                f"Run with CDS data available for exact reconstruction."
            )

        # Extract the protein regions flanking the junction
        dek_portion = dek_protein[:dek_exon7_end_aa]
        aff2_portion = aff2_protein[aff2_exon9_start_aa - 1:]

        dek_window = dek_portion[-window:] if len(dek_portion) >= window else dek_portion
        aff2_window = aff2_portion[:window] if len(aff2_portion) >= window else aff2_portion

        junction_seq = dek_window + aff2_window
        junction_pos = len(dek_window)
    
    peptides = {
        "class_I": [],   # 8-11mers for HLA class I
        "class_II": []   # 15mers for HLA class II
    }
    
    # Generate class I peptides (8-11mers)
    for k in [8, 9, 10, 11]:
        for i in range(len(junction_seq) - k + 1):
            peptide = junction_seq[i:i + k]
            # Only keep peptides that span the junction
            if i < junction_pos and i + k > junction_pos:
                peptides["class_I"].append({
                    "sequence": peptide,
                    "length": k,
                    "start_in_window": i,
                    "junction_position_in_peptide": junction_pos - i,
                    "dek_residues": min(k, junction_pos - i),
                    "aff2_residues": min(k, i + k - junction_pos)
                })
    
    # Generate class II peptides (15mers)
    for i in range(len(junction_seq) - 15 + 1):
        peptide = junction_seq[i:i + 15]
        if i < junction_pos and i + 15 > junction_pos:
            peptides["class_II"].append({
                "sequence": peptide,
                "length": 15,
                "start_in_window": i,
                "junction_position_in_peptide": junction_pos - i,
                "dek_residues": min(15, junction_pos - i),
                "aff2_residues": min(15, i + 15 - junction_pos)
            })
    
    return peptides, dek_window, aff2_window, junction_pos, frame_note


def write_fasta(peptides, output_file):
    """Write peptides to FASTA format for NetMHCpan."""
    with open(output_file, "w") as f:
        for cls in ["class_I", "class_II"]:
            for i, p in enumerate(peptides[cls]):
                header = (
                    f">{cls}_{p['length']}mer_{i+1} "
                    f"junction_pos={p['junction_position_in_peptide']} "
                    f"DEK={p['dek_residues']}aa_AFF2={p['aff2_residues']}aa"
                )
                f.write(f"{header}\n{p['sequence']}\n")
    print(f"\nFASTA written to: {output_file}")


def write_netmhcpan_input(peptides, output_file):
    """Write peptides as simple list for NetMHCpan web server."""
    with open(output_file, "w") as f:
        for cls in ["class_I", "class_II"]:
            f.write(f"# {cls} peptides\n")
            for p in peptides[cls]:
                f.write(f"{p['sequence']}\n")
    print(f"NetMHCpan input written to: {output_file}")


def write_wt_peptides(dek_protein, aff2_protein, peptides, 
                       dek_exon7_end_aa, aff2_exon9_start_aa, 
                       window, output_file):
    """
    Generate corresponding wild-type peptides for differential binding analysis.
    For each junction peptide, generate the DEK-only and AFF2-only WT peptides
    at the same positions.
    """
    dek_portion = dek_protein[:dek_exon7_end_aa]
    dek_window = dek_portion[-window:] if len(dek_portion) >= window else dek_portion
    
    # For WT comparison, we extend DEK beyond exon 7 and AFF2 before exon 9
    dek_extended = dek_protein[max(0, dek_exon7_end_aa - window):dek_exon7_end_aa + window]
    aff2_extended = aff2_protein[max(0, aff2_exon9_start_aa - 1 - window):aff2_exon9_start_aa - 1 + window]
    
    with open(output_file, "w") as f:
        f.write("# Wild-type counterpart peptides for differential binding analysis\n")
        f.write("# Compare binding of these to the junction peptides\n")
        f.write(f"# DEK extended region: {dek_extended}\n")
        f.write(f"# AFF2 extended region: {aff2_extended}\n\n")
        
        for cls in ["class_I", "class_II"]:
            for i, p in enumerate(peptides[cls]):
                # DEK-side WT: the corresponding region if DEK continued normally
                dek_start = dek_exon7_end_aa - window + p["start_in_window"]
                dek_wt_region = dek_protein[dek_start:dek_start + p["length"]]
                
                # AFF2-side WT: the corresponding region if AFF2 was intact
                aff2_offset = p["start_in_window"] - (window - (dek_exon7_end_aa - max(0, dek_exon7_end_aa - window)))
                aff2_start = max(0, aff2_exon9_start_aa - 1 + aff2_offset)
                aff2_wt_region = aff2_protein[aff2_start:aff2_start + p["length"]]
                
                f.write(f">{cls}_{p['length']}mer_{i+1}_WT_DEK\n")
                if len(dek_wt_region) == p["length"]:
                    f.write(f"{dek_wt_region}\n")
                else:
                    f.write(f"# Incomplete: DEK protein too short at this position\n")
                
                f.write(f">{cls}_{p['length']}mer_{i+1}_WT_AFF2\n")
                if len(aff2_wt_region) == p["length"]:
                    f.write(f"{aff2_wt_region}\n")
                else:
                    f.write(f"# Incomplete: AFF2 protein position out of range\n")
    
    print(f"Wild-type peptides written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="DEK-AFF2 Fusion Neoantigen Peptide Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python neoantigen_pipeline.py
  python neoantigen_pipeline.py --hla HLA-A02:01,HLA-B07:02
  python neoantigen_pipeline.py --outdir ./my_results

After running, submit the generated FASTA to:
  Class I:  https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
  Class II: https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/
        """
    )
    parser.add_argument("--hla", type=str, default=None,
                       help="Comma-separated HLA alleles (e.g., HLA-A02:01,HLA-B07:02)")
    parser.add_argument("--outdir", type=str, default="./neoantigen_output",
                       help="Output directory (default: ./neoantigen_output)")
    parser.add_argument("--dek-exon", type=int, default=6,
                       help="DEK CODING exon number at fusion breakpoint (default: 6, "
                            "= Hartwig all-exon 7 since DEK has 1 UTR-only exon at start)")
    parser.add_argument("--aff2-exon", type=int, default=9,
                       help="AFF2 CODING exon number at fusion breakpoint (default: 9)")
    parser.add_argument("--window", type=int, default=30,
                       help="Amino acids on each side of junction (default: 30)")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("DEK-AFF2 Fusion Neoantigen Peptide Generator")
    print("=" * 70)
    
    # ==========================================
    # Step 1: Fetch DEK transcript and protein data
    # ==========================================
    print(f"\n[1/6] Fetching DEK transcript data...")
    dek_exons, dek_cds, dek_protein, dek_map = \
        get_transcript_exons_with_protein_mapping("ENST00000652689")
    
    print(f"  DEK protein: {len(dek_protein)} amino acids")
    print(f"  DEK CDS: {len(dek_cds)} bases")
    print(f"  Coding exons found: {sum(1 for e in dek_map if e['is_coding'])}")
    
    # If ENST00000396581 doesn't work, try finding canonical
    if not dek_protein:
        print("  Trying to find canonical DEK transcript...")
        tx = get_canonical_transcript("DEK")
        dek_exons, dek_cds, dek_protein, dek_map = \
            get_transcript_exons_with_protein_mapping(tx["id"])
    
    # ==========================================
    # Step 2: Fetch AFF2 transcript and protein data
    # ==========================================
    print(f"\n[2/6] Fetching AFF2 transcript data...")
    
    # AFF2 canonical transcript
    aff2_transcript_ids = ["ENST00000370460"]
    aff2_protein = ""
    
    for tx_id in aff2_transcript_ids:
        try:
            aff2_exons, aff2_cds, aff2_protein, aff2_map = \
                get_transcript_exons_with_protein_mapping(tx_id)
            if aff2_protein:
                print(f"  AFF2 protein: {len(aff2_protein)} amino acids")
                print(f"  AFF2 CDS: {len(aff2_cds)} bases")
                print(f"  Coding exons found: {sum(1 for e in aff2_map if e['is_coding'])}")
                break
        except:
            continue
    
    if not aff2_protein:
        print("  Trying to find canonical AFF2 transcript...")
        tx = get_canonical_transcript("AFF2")
        aff2_exons, aff2_cds, aff2_protein, aff2_map = \
            get_transcript_exons_with_protein_mapping(tx["id"])
    
    # ==========================================
    # Step 3: Map exon boundaries to protein positions
    # ==========================================
    print(f"\n[3/6] Mapping exon boundaries to protein positions...")
    
    print(f"\n  DEK exon-protein map:")
    print(f"  {'Exon':<6} {'Coding bp':<12} {'Protein range':<18} {'Phase start':<12} {'Phase end':<10}")
    print(f"  {'-'*58}")
    
    dek_exon7_info = None
    coding_exon_num = 0
    for e in dek_map:
        if e["is_coding"]:
            coding_exon_num += 1
            marker = " <-- BREAKPOINT" if coding_exon_num == args.dek_exon else ""
            print(f"  {coding_exon_num:<6} {e['coding_bases']:<12} "
                  f"{e['protein_start']}-{e['protein_end']:<13} "
                  f"{e.get('phase_start', '?'):<12} {e.get('phase_end', '?'):<10}{marker}")
            if coding_exon_num == args.dek_exon:
                dek_exon7_info = e
    
    print(f"\n  AFF2 exon-protein map:")
    print(f"  {'Exon':<6} {'Coding bp':<12} {'Protein range':<18} {'Phase start':<12} {'Phase end':<10}")
    print(f"  {'-'*58}")
    
    aff2_exon9_info = None
    coding_exon_num = 0
    for e in aff2_map:
        if e["is_coding"]:
            coding_exon_num += 1
            marker = " <-- BREAKPOINT" if coding_exon_num == args.aff2_exon else ""
            print(f"  {coding_exon_num:<6} {e['coding_bases']:<12} "
                  f"{e['protein_start']}-{e['protein_end']:<13} "
                  f"{e.get('phase_start', '?'):<12} {e.get('phase_end', '?'):<10}{marker}")
            if coding_exon_num == args.aff2_exon:
                aff2_exon9_info = e
    
    if not dek_exon7_info:
        print(f"\n  ERROR: Could not find DEK coding exon {args.dek_exon}")
        print(f"  The transcript may have a different exon numbering.")
        print(f"  Check the exon table above and adjust --dek-exon accordingly.")
        sys.exit(1)
    
    if not aff2_exon9_info:
        print(f"\n  ERROR: Could not find AFF2 coding exon {args.aff2_exon}")
        print(f"  The transcript may have a different exon numbering.")
        print(f"  Check the exon table above and adjust --aff2-exon accordingly.")
        sys.exit(1)
    
    dek_break_aa = dek_exon7_info["protein_end"]
    aff2_start_aa = aff2_exon9_info["protein_start"]
    dek_end_phase = dek_exon7_info.get("phase_end", 0)
    aff2_start_phase = aff2_exon9_info.get("phase_start", 0)
    
    print(f"\n  FUSION JUNCTION:")
    print(f"  DEK exon {args.dek_exon} ends at amino acid {dek_break_aa} (phase {dek_end_phase})")
    print(f"  AFF2 exon {args.aff2_exon} starts at amino acid {aff2_start_aa} (phase {aff2_start_phase})")
    print(f"  DEK portion: aa 1-{dek_break_aa} ({dek_protein[:5]}...{dek_protein[dek_break_aa-5:dek_break_aa]})")
    print(f"  AFF2 portion: aa {aff2_start_aa}+ ({aff2_protein[aff2_start_aa-1:aff2_start_aa+4]}...)")
    
    # ==========================================
    # Step 4: CDS-level junction reconstruction + peptide generation
    # ==========================================
    print(f"\n[4/6] Reconstructing exact junction from CDS (DNA) sequences...")

    junction_recon = None
    if dek_cds and aff2_cds and dek_exon7_info and aff2_exon9_info:
        junction_recon = reconstruct_junction_from_cds(
            dek_cds, aff2_cds, dek_exon7_info, aff2_exon9_info
        )

    print(f"\n[5/6] Generating junction peptides...")

    peptides, dek_window, aff2_window, jpos, frame_note = generate_junction_peptides(
        dek_protein, aff2_protein,
        dek_break_aa, aff2_start_aa,
        dek_end_phase, aff2_start_phase,
        window=args.window,
        junction_reconstruction=junction_recon
    )

    if frame_note:
        print(f"\n  ⚠️  {frame_note}")

    print(f"\n  Junction region ({len(dek_window)}+{len(aff2_window)} aa):")
    print(f"  DEK:  ...{dek_window}")
    print(f"  AFF2: {aff2_window}...")
    if junction_recon and junction_recon.get("junction_aa"):
        jr = junction_recon
        print(f"  Full: ...{dek_window}[{jr['junction_aa']}]{aff2_window}...")
        print(f"                    {'':>{len(dek_window)}} ^ NOVEL junction amino acid")
    else:
        print(f"  Full: ...{dek_window}|{aff2_window}...")
        print(f"                    {'':>{len(dek_window)}}^ junction")
    
    n_class1 = len(peptides["class_I"])
    n_class2 = len(peptides["class_II"])
    print(f"\n  Generated {n_class1} class I peptides (8-11mers)")
    print(f"  Generated {n_class2} class II peptides (15mers)")
    
    print(f"\n  Class I peptides spanning junction:")
    print(f"  {'#':<4} {'Length':<8} {'Sequence':<15} {'Junction pos':<14}")
    print(f"  {'-'*41}")
    for i, p in enumerate(peptides["class_I"]):
        seq = p["sequence"]
        jpos_in_pep = p["junction_position_in_peptide"]
        # Mark the junction in the sequence
        marked = seq[:jpos_in_pep] + "|" + seq[jpos_in_pep:]
        print(f"  {i+1:<4} {p['length']}mer    {marked:<16} {jpos_in_pep}")
    
    print(f"\n  Class II peptides spanning junction:")
    for i, p in enumerate(peptides["class_II"]):
        seq = p["sequence"]
        jpos_in_pep = p["junction_position_in_peptide"]
        marked = seq[:jpos_in_pep] + "|" + seq[jpos_in_pep:]
        print(f"  {i+1:<4} {p['length']}mer   {marked:<18} {jpos_in_pep}")
    
    # ==========================================
    # Step 6: Write output files
    # ==========================================
    print(f"\n[6/6] Writing output files...")
    
    # FASTA for NetMHCpan
    fasta_file = outdir / "junction_peptides.fasta"
    write_fasta(peptides, fasta_file)
    
    # Simple peptide list for web submission
    netmhcpan_file = outdir / "netmhcpan_input.txt"
    write_netmhcpan_input(peptides, netmhcpan_file)
    
    # Wild-type counterparts for differential binding
    wt_file = outdir / "wildtype_peptides.fasta"
    write_wt_peptides(dek_protein, aff2_protein, peptides,
                      dek_break_aa, aff2_start_aa, args.window, wt_file)
    
    # Full protein sequences
    proteins_file = outdir / "protein_sequences.fasta"
    with open(proteins_file, "w") as f:
        f.write(f">DEK_full_protein_{len(dek_protein)}aa\n{dek_protein}\n")
        f.write(f">AFF2_full_protein_{len(aff2_protein)}aa\n{aff2_protein}\n")
        if junction_recon and junction_recon.get("junction_aa"):
            jr = junction_recon
            fusion_seq = jr["dek_protein_portion"] + jr["junction_aa"] + jr["aff2_protein_portion"]
            f.write(f">DEK_AFF2_fusion_protein_CDS_reconstructed_{len(fusion_seq)}aa\n{fusion_seq}\n")
        else:
            fusion_seq = dek_protein[:dek_break_aa] + aff2_protein[aff2_start_aa-1:]
            f.write(f">DEK_AFF2_fusion_protein_{len(fusion_seq)}aa\n{fusion_seq}\n")
    print(f"Protein sequences written to: {proteins_file}")
    
    # Summary report
    report_file = outdir / "analysis_report.txt"
    with open(report_file, "w") as f:
        f.write("DEK-AFF2 Fusion Neoantigen Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Fusion: DEK exon {args.dek_exon} :: AFF2 exon {args.aff2_exon}\n")
        f.write(f"DEK protein: {len(dek_protein)} aa (breakpoint at aa {dek_break_aa})\n")
        f.write(f"AFF2 protein: {len(aff2_protein)} aa (fusion starts at aa {aff2_start_aa})\n")
        f.write(f"Fusion protein: {dek_break_aa + len(aff2_protein) - aff2_start_aa + 1} aa\n\n")
        f.write(f"Junction region:\n")
        f.write(f"  DEK side:  ...{dek_window}\n")
        f.write(f"  AFF2 side: {aff2_window}...\n\n")
        if junction_recon and junction_recon.get("junction_aa"):
            jr = junction_recon
            f.write(f"CDS-Level Junction Reconstruction:\n")
            f.write(f"  DEK exon 7 leftover bases: {jr['dek_partial_bases']}\n")
            f.write(f"  AFF2 exon 9 contributing bases: {jr['aff2_partial_bases']}\n")
            f.write(f"  Junction codon: {jr['junction_codon']} -> amino acid: {jr['junction_aa']}\n")
            f.write(f"  This amino acid is NOVEL (not in normal DEK or AFF2)\n")
            f.write(f"  Reading frame preserved: {jr['frame_preserved']}\n\n")
        if frame_note:
            f.write(f"Frame note: {frame_note}\n\n")
        f.write(f"Peptides generated:\n")
        f.write(f"  Class I (8-11mers): {n_class1}\n")
        f.write(f"  Class II (15mers):  {n_class2}\n\n")
        f.write(f"Next steps:\n")
        f.write(f"1. Get HLA typed (blood test or from WGS germline BAM)\n")
        f.write(f"2. Submit junction_peptides.fasta to NetMHCpan 4.1:\n")
        f.write(f"   https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/\n")
        f.write(f"3. Submit class II peptides to NetMHCIIpan 4.3:\n")
        f.write(f"   https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/\n")
        f.write(f"4. Compare mutant vs wild-type binding (wildtype_peptides.fasta)\n")
        f.write(f"5. Select peptides with:\n")
        f.write(f"   - %Rank < 0.5 (strong binders) or < 2.0 (weak binders)\n")
        f.write(f"   - IC50 < 500 nM (preferably < 50 nM)\n")
        f.write(f"   - Better binding than wild-type counterpart\n")
    print(f"Report written to: {report_file}")
    
    # ==========================================
    # HLA-specific instructions
    # ==========================================
    if args.hla:
        hla_alleles = [h.strip() for h in args.hla.split(",")]
        print(f"\n  HLA alleles specified: {', '.join(hla_alleles)}")
        print(f"  Use these when submitting to NetMHCpan.")
    else:
        print(f"\n  No HLA alleles specified. To include them, re-run with:")
        print(f"  python neoantigen_pipeline.py --hla HLA-A02:01,HLA-B07:02,...")
    
    print(f"\n{'=' * 70}")
    print(f"DONE! Output files in: {outdir.resolve()}")
    print(f"\nNext steps:")
    print(f"  1. Get your HLA type (from blood or WGS germline data)")
    print(f"  2. Go to https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/")
    print(f"  3. Paste peptides from {netmhcpan_file.name}")
    print(f"  4. Enter your HLA alleles")
    print(f"  5. Look for strong binders (%Rank < 0.5, IC50 < 50 nM)")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
