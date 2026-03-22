# DEK-AFF2 Fusion Neoantigen Discovery Pipeline

## Background

Stage IV NSCLC with neuroendocrine differentiation, driven by a **DEK-AFF2 gene fusion** (DEK exon 7 → AFF2 exon 9). Detected by whole genome sequencing (Hartwig Medical Foundation, sample COREDB010694T2).

The fusion creates a novel protein junction that doesn't exist in normal human biology — an ideal target for a personalized mRNA cancer vaccine. It's clonal, functionally essential, and the tumor can't easily lose it.

The goal is to computationally identify candidate neoantigen peptides from this fusion that could be presented by HLA molecules and recognized by the immune system.

---

## Key Data

- **Fusion:** DEK exon 7 → AFF2 exon 9 (in-frame, high driver likelihood)
- **Other somatic variants:**
  - NF1 c.6548G>A p.Arg2183Gln (tVAF 14%)
  - ATR c.1876G>A p.Asp626Asn (tVAF 46%)
- **TMB:** 1.4 mut/Mb (low)
- **PD-L1:** <1%
- **MSI:** Stable
- **HLA type:** NOT YET KNOWN — pending clinical typing
- **Hartwig sample ID:** COREDB010694T2
- **Reference genome:** hg38 (GRCh38)

---

## Pipeline Phases & Status

### Phase 1: Fusion Junction Reconstruction — DONE

**Goal:** Determine the exact amino acid sequence at the DEK-AFF2 fusion junction.

**What was done:**
- Used Ensembl REST API to fetch DEK and AFF2 protein/CDS/exon data
- Canonical transcripts: DEK `ENST00000652689` (375 aa), AFF2 `ENST00000370460` (1311 aa)
- Resolved exon numbering: Hartwig's "DEK exon 7" = coding exon 6 (DEK has 1 UTR-only exon at the 5' end that Hartwig counts)
- Confirmed in-frame fusion: both sides at phase 0 (clean codon boundary)
- CDS-level junction reconstruction confirms no novel junction amino acid

**Key finding:**
```
Junction: ...ESEEE|AVEKA...
           DEK aa 254 | AFF2 aa 454
Fusion protein: 1112 amino acids
```

**Output:** `neoantigen_output/protein_sequences.fasta`

---

### Phase 2: Junction Peptide Generation — DONE

**Goal:** Generate all candidate peptides (8-11mers for HLA class I, 15mers for class II) spanning the fusion junction.

**What was done:**
- Sliding window over 60 aa junction region (30 aa each side)
- Generated all peptides containing at least 1 residue from each side
- Also generated wild-type counterpart peptides for differential binding comparison

**Result:** 48 candidate peptides
- 34 Class I peptides (8-11mers)
- 14 Class II peptides (15mers)

**Output:**
- `neoantigen_output/junction_peptides.fasta` — all candidate neoepitopes
- `neoantigen_output/wildtype_peptides.fasta` — normal counterparts
- `neoantigen_output/netmhcpan_input.txt` — simple list for web submission

---

### Phase 3: HLA Typing — BLOCKED

**Goal:** Determine HLA-A, HLA-B, HLA-C (class I) and HLA-DRB1 (class II) alleles.

**Status:** Waiting on clinical HLA typing from blood draw.

**Alternative path:** Computational extraction from Hartwig germline WGS BAM using OptiType (class I) or HLA-HD (class II) — requires Hartwig data access.

---

### Phase 4: Binding Prediction — DONE (pre-computed)

**Goal:** Predict which junction peptides bind HLA molecules.

**What was done:**
- Ran MHCflurry 2.1.5 (Class1AffinityPredictor) locally
- Tested all 34 class I peptides against 54 common HLA alleles (>95% population coverage)
- 2592 total predictions

**Results:**
- **2 elite binders** (IC50 < 50 nM, percentile < 0.5%)
- **23 strong binders** (IC50 < 500 nM, percentile < 2%)
- **138 weak binders** (IC50 < 5000 nM, percentile < 5%)

**Top candidates:**

| Peptide | # Alleles | Best IC50 | Key HLA alleles |
|---------|:---------:|----------:|-----------------|
| EAVEKAKPR | 3 | 30.6 nM | A*68:01, A*33:01, A*31:01 |
| SEEEAVEKA | 4 | 65.3 nM | B*40:01, B*40:02, B*44:02, B*44:03 |
| KESEEEAVEKA | 4 | 70.7 nM | B*40:01, B*40:02, B*44:02, B*44:03 |
| ESEEEAVEK | 4 | 136.8 nM | A*68:01, C*02:02, C*03:03, C*06:02 |
| KESEEEAVEK | 3 | 253.1 nM | B*27:05, B*40:02, C*03:03 |

**When HLA typing arrives, run:**
```bash
source venv312/bin/activate
python lookup_hla.py --alleles HLA-A*XX:XX,HLA-A*XX:XX,HLA-B*XX:XX,HLA-B*XX:XX,HLA-C*XX:XX,HLA-C*XX:XX
```

**Output:**
- `neoantigen_output/binding_results_all_alleles.csv` — all 2592 predictions
- `neoantigen_output/strong_binders_summary.csv` — 25 strong/elite hits
- `neoantigen_output/promiscuous_binders.csv` — peptides ranked by # alleles
- `neoantigen_output/binding_report.txt` — human-readable summary

**Class II binding predictions — DONE (no strong binders):**
Tested all 14 junction 15mers against 16 common HLA-DRB1 alleles via IEDB API (NetMHCIIpan). No strong or weak binders found (best rank: 74%). The fusion junction peptides do not appear to bind Class II HLA molecules well. Tool: `class2_binding_prediction.py` (uses IEDB REST API). For validation, consider re-running with NetMHCIIpan 4.3 locally (requires academic license from DTU — needs institutional email e.g. NKI-AVL or LUMC).

---

### Phase 5: Point Mutation Neoantigens — TODO (secondary)

**Goal:** Check NF1 R2183Q and ATR D626N for neoantigen potential.

Lower priority than the fusion (subclonal, lower VAF), but worth checking once HLA type is known.

---

### Phase 6: Immunogenicity Assessment — TODO (stretch)

**Goal:** Further prioritize candidates by predicted immunogenicity.

Not all peptides that bind HLA trigger a T-cell response. Tools like IEDB immunogenicity predictor, Seq2Neo, or DeepImmuno can estimate this.

---

## Codebase

| Script | Python | Purpose |
|--------|--------|---------|
| `neoantigen_pipeline.py` | 3.14 (venv/) | Phase 1-2: fusion reconstruction + peptide generation |
| `binding_prediction.py` | 3.12 (venv312/) | Phase 4: HLA binding prediction against all common alleles |
| `lookup_hla.py` | 3.12 (venv312/) | Filter pre-computed results to specific HLA alleles |
| `class2_binding_prediction.py` | 3.12 (venv312/) | Class II binding prediction via IEDB API |

Two venvs because MHCflurry requires Python 3.12 (incompatible with 3.14).

---

## Still Needed

- [ ] HLA type (pending clinical typing or Hartwig BAM access)
- [ ] Hartwig LINX breakpoint coordinates (to confirm exon numbering interpretation)
- [x] Class II binding predictions — done, no strong binders found
- [ ] Re-validate Class II with NetMHCIIpan 4.3 locally (requires academic license)
- [ ] Point mutation neoantigen analysis (NF1, ATR)
- [ ] Immunogenicity prediction (stretch)
- [ ] RNA-seq data (would confirm fusion expression — may not exist)

---

## References

- MHCflurry 2.0: O'Donnell et al., Cell Systems 2020
- NetMHCpan 4.1: Reynisson et al., Nucleic Acids Res 2020
- Best practices for neoantigen prediction: Wells et al., Genome Medicine 2020
- pVACtools: Hundal et al., Cancer Immunol Res 2020
- DEK-AFF2 in sinonasal carcinoma: Kuo et al., Mod Pathol 2021
