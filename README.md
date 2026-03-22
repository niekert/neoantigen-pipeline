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
3. **Predicts HLA binding** — runs [MHCflurry](https://github.com/openvax/mhcflurry) against 54 common HLA Class I alleles to identify which peptides could be presented to the immune system
4. **Provides instant HLA lookup** — once HLA typing is available, results are filtered to the patient's specific alleles in seconds

## Key Results

**Fusion junction:** `...ESEEE|AVEKA...` (DEK aa 254 | AFF2 aa 454, 1112 aa fusion protein)

The fusion protein sequence was reconstructed by fetching DEK (`ENST00000652689`, 375 aa) and AFF2 (`ENST00000370460`, 1311 aa) protein and exon data from the Ensembl REST API (GRCh38). Exon-to-protein boundary mapping confirmed that DEK coding exon 6 (Hartwig exon 7, which includes a 5' UTR exon) joins AFF2 exon 9 in-frame — both exons are phase 0, producing a clean junction. All 48 overlapping peptides (8-15mers) spanning this junction were generated as candidate neoantigens.

**Binding predictions across 54 common HLA alleles:**

Binding affinity was predicted using [MHCflurry 2.1.5](https://github.com/openvax/mhcflurry) (Class1AffinityPredictor), testing all 48 junction peptides against a panel of 54 common HLA Class I alleles (17 HLA-A, 23 HLA-B, 14 HLA-C) covering >95% of the global population. This produced 2,592 peptide-allele predictions. Binders were classified using a dual threshold: **elite** (IC50 < 50 nM *and* percentile < 0.5%), **strong** (IC50 < 500 nM *and* percentile < 2%).

| Candidate   | Type  | Best IC50 | HLA Alleles with Strong Binding        |
| ----------- | ----- | --------: | -------------------------------------- |
| EAVEKAKPR   | 9mer  |   30.6 nM | A\*68:01, A\*33:01, A\*31:01           |
| SEEEAVEKA   | 9mer  |   65.3 nM | B\*40:01, B\*40:02, B\*44:02, B\*44:03 |
| KESEEEAVEKA | 11mer |   70.7 nM | B\*40:01, B\*40:02, B\*44:02, B\*44:03 |
| ESEEEAVEK   | 9mer  |  136.8 nM | A\*68:01, C\*02:02, C\*03:03, C\*06:02 |

2 elite binders (IC50 < 50 nM), 23 strong binders (IC50 < 500 nM) identified across the allele panel. Peptides are ranked by the number of HLA alleles they bind (promiscuous binders), as these make the strongest candidates regardless of the patient's specific HLA type.

## What Is Needed Next

1. **HLA Class I typing (HLA-A, -B, -C)** — required to determine which of these candidates match the patient's immune system. Standard blood draw. Once available:

```bash
source venv312/bin/activate
python lookup_hla.py --alleles HLA-A*XX:XX,HLA-A*XX:XX,HLA-B*XX:XX,HLA-B*XX:XX,HLA-C*XX:XX,HLA-C*XX:XX
```

2. **HLA Class II binding predictions — DONE (no strong binders found).** All 14 junction 15mers were tested against 16 common HLA-DRB1 alleles via IEDB API (NetMHCIIpan). No strong (rank < 2%) or weak (rank < 10%) binders were found. This means the fusion junction peptides are unlikely to activate CD4+ helper T cells through Class II presentation. The Class I candidates remain the primary focus. For higher-accuracy validation, re-run with [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) locally (requires academic license from DTU — needs an institutional email e.g. from NKI-AVL or LUMC).

3. **Hartwig LINX breakpoint coordinates** — to definitively confirm the exon numbering interpretation (see Technical Notes).

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

---

## Full Project Summary

**Date:** 2026-03-22 | **Sample:** COREDB010694T2 (Hartwig Medical Foundation) | **Diagnosis:** Stage IV NSCLC with neuroendocrine differentiation

### Background

The tumor is driven by a DEK-AFF2 gene fusion (DEK exon 7 → AFF2 exon 9, as reported by Hartwig/LINX). This fusion creates a chimeric protein that does not exist in normal human cells. The junction between DEK and AFF2 produces novel peptide sequences that could be recognized by the immune system — making them ideal targets for a personalized mRNA cancer vaccine.

### Phase 1-2: Fusion Junction Reconstruction & Peptide Generation

**Tool:** Custom Python script using Ensembl REST API (`neoantigen_pipeline.py`)

1. **Fetched protein sequences** for DEK (375 aa) and AFF2 (1311 aa) from Ensembl using their canonical transcripts:
   - DEK: `ENST00000652689` (RefSeq NM_003472.4)
   - AFF2: `ENST00000370460`

2. **Resolved exon numbering discrepancy:**
   - Hartwig reports "DEK exon 7" counting ALL exons (including UTR)
   - DEK has 11 total exons, 1 of which is a 5' UTR-only exon
   - Hartwig's "exon 7" = coding exon 6 in Ensembl numbering
   - Coding exon 6 ends at a clean codon boundary (phase 0), confirming the fusion is **in-frame**
   - Coding exon 7 would have been out-of-frame (phase 1), producing a premature stop codon after just 2 amino acids

3. **Confirmed the fusion is in-frame:**
   - DEK coding exon 6: phase_end = 0 (clean codon boundary)
   - AFF2 coding exon 9: phase_start = 0 (clean codon boundary)
   - No novel junction amino acid, no frameshift
   - Fusion protein: 1112 amino acids (functional chimeric protein)

4. **Reconstructed the fusion junction:**
   ```
   DEK side:  ...ILSDESSSDEDEKKNKEESSDDEDKESEEE
   AFF2 side: AVEKAKPRNNPVNPPLATPQPPPAVQASGG...
   Junction:  ...ESEEE|AVEKA...
              (DEK aa 254 | AFF2 aa 454)
   ```

5. **Generated 48 candidate peptides** spanning the junction:
   - 34 Class I peptides (8-11mers) for HLA-A, -B, -C presentation
   - 14 Class II peptides (15mers) for HLA-DR, -DQ, -DP presentation
   - All peptides contain at least 1 amino acid from each side of the junction, ensuring they are truly tumor-specific

### Phase 4: HLA Class I Binding Prediction

**Tool:** MHCflurry 2.1.5 (Class1AffinityPredictor) via `binding_prediction.py`

Since the patient's HLA type is not yet known, binding predictions were run against a panel of 54 common HLA Class I alleles covering >95% of the global population. This allows instant lookup once HLA typing arrives.

**Results:** 2,592 predictions (48 peptides x 54 alleles)
- 2 elite binders (IC50 < 50 nM, percentile < 0.5%)
- 23 strong binders (IC50 < 500 nM, percentile < 2%)
- 138 weak binders (IC50 < 5000 nM, percentile < 5%)

### Top Neoantigen Candidates

Ranked by number of HLA alleles showing strong binding:

| Rank | Peptide | Type | Junction | Alleles Bound | Best IC50 | Best HLA | Notes |
|------|---------|------|----------|---------------|-----------|----------|-------|
| #1 | EAVEKAKPR | 9mer | pos 1 (1 DEK, 8 AFF2) | 3 (2 elite) | 30.6 nM | A\*68:01 | Strongest binder overall. Also: A\*33:01 (43.5 nM, elite), A\*31:01 (350 nM) |
| #2 | ESEEEAVEK | 9mer | pos 5 (5 DEK, 4 AFF2) | 4 (0 elite) | 136.8 nM | A\*68:01 | Also: C\*02:02, C\*03:03, C\*06:02 |
| #3 | SEEEAVEKA | 9mer | pos 4 (4 DEK, 5 AFF2) | 4 (0 elite) | 65.3 nM | B\*40:02 | Also: B\*40:01, B\*44:02, B\*44:03 |
| #4 | KESEEEAVEKA | 11mer | pos 6 (6 DEK, 5 AFF2) | 4 (0 elite) | 70.7 nM | B\*40:02 | Also: B\*40:01, B\*44:02, B\*44:03 |
| #5 | KESEEEAVEK | 10mer | pos 6 (6 DEK, 4 AFF2) | 3 (0 elite) | 253.1 nM | B\*27:05 | Also: B\*40:02, C\*03:03 |

### HLA Class II Binding Prediction

**Tool:** IEDB REST API (NetMHCIIpan server-side) via `class2_binding_prediction.py`

All 14 junction 15mers were tested against 16 common HLA-DRB1 alleles. **No strong (rank < 2%) or weak (rank < 10%) binders were found.** This means the fusion junction peptides are unlikely to activate CD4+ helper T cells through Class II presentation. The Class I candidates remain the primary focus.

### How to Read These Results

- **IC50 (nM)** — Binding affinity. The concentration of peptide needed to occupy 50% of HLA molecules. Lower = stronger. <50 nM elite, <500 nM strong, <5000 nM weak.
- **Percentile (%)** — How this peptide ranks vs. random peptides for the same HLA. Lower = better. <0.5% elite, <2% strong.
- **Junction position** — Where the DEK|AFF2 split falls within the peptide. Peptides with the junction near the middle are most novel (more residues from both sides).
- **Promiscuous binder** — A peptide that binds strongly to many different HLA alleles. Best vaccine candidates because they work regardless of HLA type.

### What Is Still Needed

1. **HLA typing (critical, blocking)** — The patient's HLA-A, -B, -C alleles determine which candidates will actually work. Standard blood draw. Once available, run `lookup_hla.py`.
2. **Hartwig LINX data (important, not blocking)** — Exact genomic breakpoint coordinates to confirm exon numbering interpretation. Look for `*.linx.fusion.tsv` file.
3. **Point mutation neoantigens (secondary)** — NF1 R2183Q and ATR D626N could also produce neoantigens, though lower priority than the fusion (subclonal, lower VAF).
4. **Immunogenicity prediction (stretch goal)** — Not all HLA-binding peptides trigger T-cell responses. Tools like IEDB immunogenicity predictor or DeepImmuno can estimate this.

### Important Caveats

- These are computational predictions, not experimental results. Binding prediction accuracy is ~80-90% for strong binders.
- Binding to HLA is necessary but not sufficient for immune recognition. The peptide must also be processed, presented, and recognized by T cells.
- The exon numbering interpretation (Hartwig exon 7 = coding exon 6) is based on frame analysis and is highly likely correct, but can only be 100% confirmed with the LINX breakpoint coordinates.
