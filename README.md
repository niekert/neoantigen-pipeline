# DEK-AFF2 Fusion Neoantigen Discovery

Exploratory computational pipeline to generate and rank candidate neoantigen peptides from a reported DEK-AFF2 gene fusion.

## Clinical Context

- **Diagnosis:** Stage IV NSCLC with neuroendocrine differentiation
- **Driver:** Reported DEK-AFF2 gene fusion (DEK exon 7 → AFF2 exon 9)
- **Source:** Whole genome sequencing by Hartwig Medical Foundation (hg38)
- **TMB:** 1.4 mut/Mb | **PD-L1:** <1% | **MSI:** Stable

> **Status note**
> This repository is hypothesis-generating, not clinically validated. The exact fusion breakpoint still needs confirmation from LINX and/or RNA evidence. The current output files reflect the present code and should still be treated as exploratory.

## What This Pipeline Does

1. **Reconstructs the fusion junction** — fetches DEK and AFF2 protein sequences from Ensembl, maps exon boundaries, and determines the exact amino acid sequence at the fusion point
2. **Generates candidate peptides** — creates all 8-15mer peptides that span the junction (48 candidates)
3. **Predicts HLA binding** — runs [MHCflurry](https://github.com/openvax/mhcflurry) against 54 common HLA Class I alleles to identify which peptides could be presented to the immune system
4. **Assesses differential binding heuristics** — computes agretopicity (DAI) by comparing fusion vs wildtype binding as a rough ranking signal, not a direct measure of clinical immunogenicity
5. **Provides instant HLA lookup** — once HLA typing is available, results are filtered to the patient's specific alleles in seconds

## Current Working Model

**Proposed fusion junction under the current transcript/exon assumptions:** `...ESEEE|AVEKA...` (DEK aa 254 | AFF2 aa 454, 1112 aa fusion protein)

The fusion protein sequence was reconstructed by fetching DEK (`ENST00000652689`, 375 aa) and AFF2 (`ENST00000370460`, 1311 aa) protein and exon data from the Ensembl REST API (GRCh38). Under this model, DEK coding exon 6 (Hartwig exon 7, if Hartwig counted the 5' UTR exon) joins AFF2 exon 9 in-frame, with both exons at phase 0. This is a reasonable working assumption, but it is not yet definitive without exact LINX breakpoint coordinates or RNA fusion support. All 48 overlapping peptides (34 Class I 8-11mers and 14 Class II 15mers) spanning this assumed junction were generated as candidate neoantigens.

**Class I binding signal across 54 common HLA alleles:**

Binding affinity was predicted using [MHCflurry 2.1.5](https://github.com/openvax/mhcflurry) (Class1AffinityPredictor) against the 34 Class I junction peptides and 54 common HLA Class I alleles (17 HLA-A, 23 HLA-B, 14 HLA-C).

| Candidate   | Type  | Best IC50 | HLA Alleles with Strong Binding        | Interpretation |
| ----------- | ----- | --------: | -------------------------------------- | -------------- |
| EAVEKAKPR   | 9mer  |   30.6 nM | A\*68:01, A\*33:01, A\*31:01           | Strongest current Class I signal in this modeled junction |
| SEEEAVEKA   | 9mer  |   65.3 nM | B\*40:01, B\*40:02, B\*44:02, B\*44:03 | Plausible candidate, pending patient HLA match |
| KESEEEAVEKA | 11mer |   70.7 nM | B\*40:01, B\*40:02, B\*44:02, B\*44:03 | Plausible candidate, pending patient HLA match |
| ESEEEAVEK   | 9mer  |  136.8 nM | A\*68:01, C\*02:02, C\*03:03, C\*06:02 | Plausible candidate, pending patient HLA match |

Current results show 2 elite binders and 21 additional strong binders (23 strong/elite hits total) across the allele panel. Agretopicity remains a rough prioritization aid rather than evidence that the peptides will be processed, presented, and recognized by T cells in vivo.

## What Is Needed Next

1. **Exact fusion breakpoint / transcript confirmation** — the most important unresolved issue. The current peptide set is only as good as the assumed junction. Confirm with Hartwig LINX coordinates and, ideally, RNA fusion evidence.

2. **HLA Class I typing (HLA-A, -B, -C)** — required to determine which of these candidates could match the patient's immune system. Standard blood draw. Once available:

```bash
source venv312/bin/activate
python lookup_hla.py --alleles HLA-A*XX:XX,HLA-A*XX:XX,HLA-B*XX:XX,HLA-B*XX:XX,HLA-C*XX:XX,HLA-C*XX:XX
```

3. **HLA Class II binding predictions — initial DRB1 screen only.** All 14 junction 15mers were tested against 16 common HLA-DRB1 alleles via IEDB API (NetMHCIIpan). No obvious binders were seen in the returned panel, although one peptide/allele call is still missing because the IEDB API intermittently returned `403` responses. That is useful negative evidence, but it is not definitive: the run is DRB1-only, and negative binding predictions do not rule out all possible CD4 help. For higher-accuracy validation, re-run with [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) locally and consider DQ/DP as well.

4. **Immunogenicity assessment — exploratory heuristic only.** The agretopicity/DAI output is useful for rough ranking, but it is not a validated immunogenicity model and should not be read as proof that the peptides are clinically meaningful neoantigens. The current implementation reports both DEK-side and AFF2-side comparisons and ranks candidates using the lower available DAI as a conservative summary. For additional validation, submit candidates to [NetCTLpan 1.1](https://services.healthtech.dtu.dk/services/NetCTLpan-1.1/) (proteasomal cleavage + TAP transport) and [DeepImmuno](https://deepimmuno.research.cchmc.org/) (T-cell response prediction).

5. **Class II re-run with NetMHCIIpan 4.3 locally, including HLA-DQ and HLA-DP alleles** — the current Class II results used the IEDB REST API (NetMHCIIpan via web) against DRB1 only. A local run would cover DQ/DP and give higher accuracy.

   > **BLOCKER:** NetMHCIIpan 4.3 requires a download from DTU (healthtech.dtu.dk) which is restricted to academic/whitelisted institutional email addresses. No `.edu` or whitelisted address available. Resolution: request access via NKI-AVL or LUMC collaborator, or ask oncologist to facilitate.

6. **Hartwig LINX breakpoint coordinates** — still needed to definitively confirm the exon numbering interpretation (see Technical Notes).

## Output Files

Current exploratory outputs are in `neoantigen_output/`:

| File                              | Description                                         |
| --------------------------------- | --------------------------------------------------- |
| `full_project_summary.txt`        | Human-readable project summary                    |
| `binding_report.txt`              | Class I binding summary                             |
| `strong_binders_summary.csv`      | Strong/elite predictions                            |
| `promiscuous_binders.csv`         | Peptides ranked by number of HLA alleles bound      |
| `binding_results_all_alleles.csv` | Class I binding matrix (34 peptides x 54 alleles)   |
| `junction_peptides.fasta`         | All 48 candidate peptides                           |
| `wildtype_peptides.fasta`         | Normal counterpart peptides for comparison          |
| `protein_sequences.fasta`         | Full DEK, AFF2, and fusion protein sequences        |
| `analysis_report.txt`             | Fusion junction reconstruction details              |
| `agretopicity_results.csv`        | Fusion vs wildtype binding comparison (DAI scores)  |
| `final_candidates.csv`            | Ranked candidates by composite score                |
| `immunogenicity_report.txt`       | Immunogenicity assessment with top 5 candidates     |
| `netctlpan_input.txt`             | Ready-to-submit file for NetCTLpan web server       |
| `alphafold_outputs/EAVEKAKPR_HLA-A6801.png` | AlphaFold 3 structure render — peptide (yellow) in HLA groove |
| `alphafold_outputs/EAVEKAKPR_HLA-A6801.zip` | Full AlphaFold 3 output (CIF models, confidence scores) |
| `alphafold_inputs/*.json`         | AlphaFold 3 input files for top 3 candidates        |

## Technical Notes

- **Exon numbering:** Hartwig/LINX reports "DEK exon 7" counting all exons including UTR. DEK has 1 UTR-only exon at the 5' end, so this may map to coding exon 6 (ENST00000652689). That interpretation is plausible from frame analysis, but it is still an inference until the exact breakpoint is checked.
- **Transcripts used:** DEK `ENST00000652689` (canonical, 375 aa, RefSeq NM_003472.4), AFF2 `ENST00000370460` (canonical, 1311 aa)
- **Binding predictions:** MHCflurry 2.1.5 was used for Class I affinity prediction across the 34 generated Class I junction peptides and 54 HLA-I alleles.
- **Immunogenicity:** Agretopicity (DAI) was computed as an exploratory differential-binding heuristic. The current implementation reports both DEK-side and AFF2-side values and ranks candidates using the lower available DAI as a conservative summary.
- **AlphaFold 3 structure prediction:** Top candidates were modeled as peptide + HLA + B2M complexes using AlphaFold Server. The peptide-HLA chain-pair ipTM is the relevant confidence metric (overall ipTM is dominated by the well-known HLA-B2M interface). Results in `neoantigen_output/alphafold_outputs/`.

## AlphaFold 3 Structure Predictions

Peptide-HLA-B2M complexes modeled via [AlphaFold Server](https://alphafoldserver.com). Each model requires HLA typing confirmation before clinical interpretation — these are only meaningful if the patient carries the listed allele.

| Peptide | HLA allele | IC50 | Cons. DAI | peptide-HLA ipTM | Structure | Notes |
|---------|-----------|-----:|----------:|:----------------:|:---------:|-------|
| EAVEKAKPR | A\*68:01 | 30.6 nM | 14.1x | 0.64 | [view](neoantigen_output/alphafold_outputs/EAVEKAKPR_HLA-A6801.png) | Flat in groove. Best binder overall. **Requires A\*68:01.** |
| SEEEAVEKA | B\*40:02 | 65.3 nM | 1.4x | 0.86 | [view](neoantigen_output/alphafold_outputs/SEEEAVEKA_HLA-B4002.png) | Flat in groove. Best structural confidence. Low DAI. **Requires B\*40:02.** |
| KESEEEAV | B\*40:02 | 95.7 nM | 103.4x | 0.64 | [view](neoantigen_output/alphafold_outputs/KESEEEAV_HLA-B4002.png) | Flat in groove. Best conservative DAI. **Requires B\*40:02.** |
| ESEEEAVEKAK | A\*68:01 | 176.7 nM | 54.5x | 0.78 | [view](neoantigen_output/alphafold_outputs/ESEEEAVEKAK_HLA-A6801.png) | 11mer — bulging loop above groove. Weaker fit. **Requires A\*68:01.** |

**ipTM guidance:** >0.8 = high confidence, 0.6–0.8 = moderate (plausible positioning), <0.6 = low confidence. The headline ipTM shown on the AlphaFold server (~0.96 for all models) is dominated by the HLA-B2M interface; the peptide-HLA value above is what matters for binding assessment.
- **Not yet done:** definitive LINX breakpoint verification, RNA-level fusion confirmation, point mutation neoantigens (NF1 R2183Q, ATR D626N), NetCTLpan pathway validation, deeper immunogenicity modeling.

## How to Run

### Setup (first time after cloning)

MHCflurry requires Python 3.12 (it does not support 3.13+). All other scripts work with 3.12+ as well, so a single Python 3.12 venv is the simplest setup.

```bash
# Create virtual environment (Python 3.12 recommended)
python3.12 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
pip install -r requirements-mhcflurry.txt

# Download MHCflurry models (~150 MB, one-time)
mhcflurry-downloads fetch
```

### Running the pipeline

```bash
source venv/bin/activate

# Phase 1-2: Fusion reconstruction + peptide generation
python neoantigen_pipeline.py

# Phase 3: HLA Class I binding prediction
python binding_prediction.py

# Phase 4: HLA Class II binding prediction (uses IEDB API, no local models needed)
python class2_binding_prediction.py

# Phase 5: Immunogenicity assessment (agretopicity / DAI)
python immunogenicity_assessment.py

# When HLA type is known:
python lookup_hla.py --alleles HLA-A*XX:XX,HLA-A*XX:XX,HLA-B*XX:XX,HLA-B*XX:XX,HLA-C*XX:XX,HLA-C*XX:XX
```

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

**Date:** 2026-03-22 | **Source:** Hartwig Medical Foundation WGS | **Diagnosis:** Stage IV NSCLC with neuroendocrine differentiation

### Background

The tumor is reported to be driven by a DEK-AFF2 gene fusion (DEK exon 7 → AFF2 exon 9, as reported by Hartwig/LINX). If the assumed breakpoint and transcript model are correct, this fusion would create a chimeric protein absent from normal cells and could generate novel junction peptides worth screening as candidate neoantigens.

### Phase 1-2: Fusion Junction Reconstruction & Peptide Generation

**Tool:** Custom Python script using Ensembl REST API (`neoantigen_pipeline.py`)

1. **Fetched protein sequences** for DEK (375 aa) and AFF2 (1311 aa) from Ensembl using their canonical transcripts:
   - DEK: `ENST00000652689` (RefSeq NM_003472.4)
   - AFF2: `ENST00000370460`

2. **Resolved exon numbering discrepancy:**
   - Hartwig reports "DEK exon 7" counting ALL exons (including UTR)
   - DEK has 11 total exons, 1 of which is a 5' UTR-only exon
   - Hartwig's "exon 7" = coding exon 6 in Ensembl numbering
   - Coding exon 6 ends at a clean codon boundary (phase 0), supporting an **in-frame working model**
   - Coding exon 7 would have been out-of-frame (phase 1), producing a premature stop codon after just 2 amino acids

3. **Built an in-frame working model of the fusion:**
   - DEK coding exon 6: phase_end = 0 (clean codon boundary)
   - AFF2 coding exon 9: phase_start = 0 (clean codon boundary)
   - No novel junction amino acid is predicted under this model
   - Fusion protein: 1112 amino acids under this model

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

Since the patient's HLA type is not yet known, binding predictions were run against a panel of 54 common HLA Class I alleles covering >95% of the global population. This is a reasonable broad pre-screen, but it does not replace patient-specific HLA matching.

**Current results:** 1,836 Class I predictions are reported in the current output files (34 peptides x 54 alleles). The shortlist below reflects the present modeled junction and current code, but it remains exploratory until the breakpoint/transcript model is confirmed.

### Top Neoantigen Candidates

Ranked by number of HLA alleles showing strong binding in the current output:

| Rank | Peptide | Type | Junction | Alleles Bound | Best IC50 | Best HLA | Notes |
|------|---------|------|----------|---------------|-----------|----------|-------|
| #1 | EAVEKAKPR | 9mer | pos 1 (1 DEK, 8 AFF2) | 3 (2 elite) | 30.6 nM | A\*68:01 | Strongest binder overall. Also: A\*33:01 (43.5 nM, elite), A\*31:01 (350 nM) |
| #2 | ESEEEAVEK | 9mer | pos 5 (5 DEK, 4 AFF2) | 4 (0 elite) | 136.8 nM | A\*68:01 | Also: C\*02:02, C\*03:03, C\*06:02 |
| #3 | SEEEAVEKA | 9mer | pos 4 (4 DEK, 5 AFF2) | 4 (0 elite) | 65.3 nM | B\*40:02 | Also: B\*40:01, B\*44:02, B\*44:03 |
| #4 | KESEEEAVEKA | 11mer | pos 6 (6 DEK, 5 AFF2) | 4 (0 elite) | 70.7 nM | B\*40:02 | Also: B\*40:01, B\*44:02, B\*44:03 |
| #5 | KESEEEAVEK | 10mer | pos 6 (6 DEK, 4 AFF2) | 3 (0 elite) | 253.1 nM | B\*27:05 | Also: B\*40:02, C\*03:03 |

### HLA Class II Binding Prediction

**Tool:** IEDB REST API (NetMHCIIpan server-side) via `class2_binding_prediction.py`

All 14 junction 15mers were tested against 16 common HLA-DRB1 alleles in an initial screen. **No strong (rank < 2%) or weak (rank < 10%) binders were seen in the returned DRB1 results.** This argues against obvious DRB1-mediated Class II presentation for these exact 15mers, but it is not definitive because DQ/DP were not tested and one expected API row did not return.

### Immunogenicity Assessment (Agretopicity)

**Tool:** MHCflurry 2.1.5 via `immunogenicity_assessment.py`

For each strong-binding fusion peptide, the corresponding wildtype peptides (from normal DEK and AFF2) were also run through MHCflurry. The Differential Agretopicity Index (DAI = IC50_wildtype / IC50_fusion) is used here only as a rough ranking heuristic: higher DAI means the fusion peptide binds better than the chosen wildtype comparator. It does **not** by itself establish immunogenicity.

**Current results:** the script reports both DEK-side and AFF2-side wildtype comparisons and uses the lower available DAI as the conservative summary for ranking. These values are still heuristic and should be interpreted cautiously.

| Rank | Peptide     | Best HLA    | IC50 (nM) | Conservative DAI | DEK/AFF2 WT IC50 (nM) |
|------|-------------|-------------|----------:|-----------------:|----------------------:|
| #1   | KESEEEAV    | B\*40:02    |      95.7 |            103.4 |        9,895 / 30,068 |
| #2   | ESEEEAVEKAK | A\*68:01    |     176.7 |             54.5 |       26,091 / 9,635  |
| #3   | EAVEKAKPR   | A\*68:01    |      30.6 |             14.1 |          431 / 28,625 |
| #4   | KESEEEAVEKA | B\*40:02    |      70.7 |              3.6 |         255 / 10,789  |
| #5   | SEEEAVEKA   | B\*40:02    |      65.3 |              1.4 |          95 / 22,505  |

Large DAI values can occur for fusion peptides because the junction creates sequences absent from the native proteins. That can be a useful prioritization signal, but it should not be read as proof of presentation or clinical immune recognition.

### How to Read These Results

- **IC50 (nM)** — Binding affinity. The concentration of peptide needed to occupy 50% of HLA molecules. Lower = stronger. <50 nM elite, <500 nM strong, <5000 nM weak.
- **Percentile (%)** — How this peptide ranks vs. random peptides for the same HLA. Lower = better. <0.5% elite, <2% strong.
- **Junction position** — Where the DEK|AFF2 split falls within the peptide. Peptides with the junction near the middle are most novel (more residues from both sides).
- **Promiscuous binder** — A peptide that binds strongly to many different HLA alleles. Best vaccine candidates because they work regardless of HLA type.

### What Is Still Needed

1. **Breakpoint confirmation (critical)** — Exact LINX coordinates, and ideally RNA fusion support, are needed to confirm that the assumed junction is the right one.
2. **HLA typing (critical)** — The patient's HLA-A, -B, -C alleles determine which candidates could plausibly matter. Once available, run `lookup_hla.py`.
3. **Point mutation neoantigens (secondary)** — NF1 R2183Q and ATR D626N could also produce neoantigens, though lower priority than the fusion (subclonal, lower VAF).
4. **Additional validation (optional)** — Submit top candidates to [NetCTLpan 1.1](https://services.healthtech.dtu.dk/services/NetCTLpan-1.1/) for proteasomal cleavage + TAP transport pathway prediction, and [DeepImmuno](https://deepimmuno.research.cchmc.org/) for T-cell response prediction (9-10mers only).

### Important Caveats

- These are computational predictions, not experimental results. Even strong predicted binders often fail at the stages of processing, presentation, or T-cell recognition.
- Binding to HLA is necessary but not sufficient for immune recognition. This repository currently uses MHC binding affinity as the main screen, not a full presentation or TCR-recognition model.
- The exon numbering interpretation (Hartwig exon 7 = coding exon 6) is a plausible working model based on frame analysis, but it can only be confirmed with the LINX breakpoint coordinates or RNA fusion evidence.
