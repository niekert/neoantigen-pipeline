# DEK-AFF2 Fusion Neoantigen Discovery

Computational pipeline to identify candidate neoantigen peptides from a DEK-AFF2 gene fusion for personalized cancer immunotherapy.

## Clinical Context

- **Diagnosis:** Stage IV NSCLC with neuroendocrine differentiation
- **Driver:** DEK-AFF2 gene fusion (DEK exon 7 → AFF2 exon 9, in-frame)
- **Source:** Whole genome sequencing by Hartwig Medical Foundation (COREDB010694T2, hg38)
- **TMB:** 1.4 mut/Mb | **PD-L1:** <1% | **MSI:** Stable

## What This Pipeline Does

1. **Reconstructs the fusion junction** — fetches DEK and AFF2 protein sequences from Ensembl, maps exon boundaries, and determines the exact amino acid sequence at the fusion point
2. **Generates candidate peptides** — creates all 8-15mer peptides that span the junction (48 candidates)
3. **Predicts HLA binding** — runs MHCflurry against 54 common HLA Class I alleles to identify which peptides could be presented to the immune system
4. **Provides instant HLA lookup** — once HLA typing is available, results are filtered to the patient's specific alleles in seconds

## Key Results

**Fusion junction:** `...ESEEE|AVEKA...` (DEK aa 254 | AFF2 aa 454, 1112 aa fusion protein)

**Binding predictions across 54 common HLA alleles:**

| Candidate   | Type  | Best IC50 | HLA Alleles with Strong Binding        |
| ----------- | ----- | --------: | -------------------------------------- |
| EAVEKAKPR   | 9mer  |   30.6 nM | A\*68:01, A\*33:01, A\*31:01           |
| SEEEAVEKA   | 9mer  |   65.3 nM | B\*40:01, B\*40:02, B\*44:02, B\*44:03 |
| KESEEEAVEKA | 11mer |   70.7 nM | B\*40:01, B\*40:02, B\*44:02, B\*44:03 |
| ESEEEAVEK   | 9mer  |  136.8 nM | A\*68:01, C\*02:02, C\*03:03, C\*06:02 |

2 elite binders (IC50 < 50 nM), 23 strong binders (IC50 < 500 nM) identified across the allele panel.

## What Is Needed Next

**HLA Class I typing (HLA-A, -B, -C)** — required to determine which of these candidates match the patient's immune system. Once available:

```bash
source venv312/bin/activate
python lookup_hla.py --alleles HLA-A*XX:XX,HLA-A*XX:XX,HLA-B*XX:XX,HLA-B*XX:XX,HLA-C*XX:XX,HLA-C*XX:XX
```

## Output Files

All results are in `neoantigen_output/`:

| File                              | Description                                         |
| --------------------------------- | --------------------------------------------------- |
| `full_project_summary.txt`        | Complete human-readable project summary             |
| `binding_report.txt`              | Binding prediction results and interpretation guide |
| `strong_binders_summary.csv`      | All strong/elite binding predictions                |
| `promiscuous_binders.csv`         | Peptides ranked by number of HLA alleles bound      |
| `binding_results_all_alleles.csv` | Full results matrix (2592 predictions)              |
| `junction_peptides.fasta`         | All 48 candidate peptides                           |
| `wildtype_peptides.fasta`         | Normal counterpart peptides for comparison          |
| `protein_sequences.fasta`         | Full DEK, AFF2, and fusion protein sequences        |
| `analysis_report.txt`             | Fusion junction reconstruction details              |

## Technical Notes

- **Exon numbering:** Hartwig/LINX reports "DEK exon 7" counting all exons including UTR. DEK has 1 UTR-only exon at the 5' end, so this maps to coding exon 6 (ENST00000652689). This interpretation is confirmed by frame analysis — only coding exon 6 (phase 0) produces an in-frame fusion with AFF2 exon 9 (phase 0).
- **Transcripts used:** DEK `ENST00000652689` (canonical, 375 aa, RefSeq NM_003472.4), AFF2 `ENST00000370460` (canonical, 1311 aa)
- **Binding predictions:** MHCflurry 2.1.5, Class I only. Class II (HLA-DR/DQ/DP) predictions pending.
- **Not yet done:** LINX breakpoint verification, Class II binding, point mutation neoantigens (NF1 R2183Q, ATR D626N), immunogenicity prediction.

## How to Run

```bash
# Phase 1-2: Fusion reconstruction + peptide generation
source venv/bin/activate
python neoantigen_pipeline.py

# Phase 4: Binding prediction
source venv312/bin/activate
python binding_prediction.py

# When HLA type is known:
python lookup_hla.py --alleles HLA-A*XX:XX,...
```

Requires Python 3.14 (venv/) and Python 3.12 (venv312/) due to MHCflurry compatibility.

## Tools & References

- [Ensembl REST API](https://rest.ensembl.org) — gene/protein/exon data (GRCh38)
- [MHCflurry 2.0](https://github.com/openvax/mhcflurry) — O'Donnell et al., Cell Systems 2020
- [NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) — Reynisson et al., Nucleic Acids Res 2020
- Wells et al., "Key Parameters of Tumor Epitope Immunogenicity", Genome Medicine 2020
- Kuo et al., "DEK-AFF2 fusion in sinonasal carcinoma", Mod Pathol 2021

## Disclaimer

This is a computational research project, not a clinical recommendation. All findings are predictions based on published algorithms and public databases. Results should be reviewed and validated by qualified immunologists and oncologists before any clinical application.
