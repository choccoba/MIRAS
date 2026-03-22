# MIRAS
### miRNA Integrated Regulatory Analysis across disease Spectrum of chronic HBV infection

**Principal Investigator**: YoungMin (choccoba)  
**Companion project**: [ITLAS](https://github.com/choccoba/ITLAS) — HBV-IT immunopathogenesis (CMH manuscript)

---

## Project Overview

MIRAS is a multi-pipeline bioinformatics project investigating **miRNA regulation** across the full disease spectrum of chronic HBV infection:

> **IT** (Immune Tolerant) → **IA** (Immune Active) → **CR** (Chronic Resolved) → **HCC**

Companion to **ITLAS/v18**: explores the miRNA regulatory layer underlying the 6-layer effector suppression framework in the HBV-IT phase. The ultimate goal is to build a **complete molecular disease map of HBV pathogenesis** and find the path from IT phase to functional cure.

---

## Analysis Pipeline

| Step | Directory | Tools | Purpose |
|------|-----------|-------|---------|
| 1 | 01_Step1_GEO_DEG | R (GEOquery, limma) | GEO datasets + IT-specific DEG |
| 2 | 02_Step2_miRNA_Mapping | R (multiMiR, HMDD v4.0) | IT-gene seed → miRNA mapping |
| 3 | 03_Step3_MetaAnalysis | Python / Colab | Literature meta-analysis |
| 4 | 04_Step4_Integration | R + Python | Integrated IT-miRNA network |

---

## Key GEO Datasets

| Dataset | Content | Phase annotation |
|---------|---------|------------------|
| GSE65359 | Liver transcriptome | IT / IA / IC (primary) |
| GSE83148 | HBV liver tissue n=128 | HBV vs Normal |
| GSE84044 | CHB fibrosis cohort n=124 | Fibrosis grade |
| GSE69580 | HBV-HCC miRNA microarray | HCC vs non-tumor |

---

## v18 Connection

MIRAS targets miRNA regulators of ITLAS v18 key genes:
AIM2, LGALS9, TGFB1, DNMT1, SOCS1, FOXP3, IL10, PDCD1, CD274, CTLA4

Hypothesized miRNA axes:
- miR-152 / miR-148a → DNMT1 : epigenetic regulation in IT phase
- miR-155 → SOCS1 : JAK-STAT immune suppression
- miR-21 → TGFB1 : TGF-beta tolerogenic pathway
- miR-146a → SOCS1/TRAF6 : innate immune modulation

---

## Project Philosophy

> HBV pathogenesis를 완벽하게 이해하고 IT 단계부터 functional cure로 이르는 길을 찾는 것
> Papers are milestones, not destinations.

---

## Local Setup

git clone https://github.com/choccoba/MIRAS "G:/내 드라이브/MIRAS"

---

## Related Projects

- ITLAS https://github.com/choccoba/ITLAS — HBV-IT immunopathogenesis (CMH submission)
- hbv-tahoe-analysis https://github.com/choccoba/hbv-tahoe-analysis
- alphafold2_hbv_core https://github.com/choccoba/alphafold2_hbv_core
