# HGSOC Virtual Screening Pipeline
### AI-driven drug repurposing for RelB-dependent targets in High-Grade Serous Ovarian Cancer

> **TARVIA Lab** · Omar Lujano Olazaba, PhD · June 2026

---

## Overview

This pipeline screens 3,311 FDA-approved small molecules against a classifier trained on 26 druggable targets derived from a 352-gene RelB-dependent transcriptomic signature in HGSOC. It bridges the gap between the upstream biomarker-discovery-pipeline (which identifies RelB-regulated genes) and experimental validation by ranking repurposable compounds by predicted activity.

---

## Scientific Context

**RelB** is an NF-κB transcription factor subunit that drives platinum-resistance and spheroid survival in HGSOC. Using RNA-seq differential expression data from RelB-knockdown vs. control experiments, we identified 352 RelB-dependent genes (logFC-ranked, Benjamini–Hochberg corrected). Of these:

- **73 high-priority** (|RelB logFC| ≥ 2.0, p < 1×10⁻⁶)
- **102 medium-priority** (|RelB logFC| ≥ 1.0, p < 1×10⁻⁴)
- **177 low-priority** (remaining)

The 175 high+medium genes were queried against ChEMBL to identify druggable targets with sufficient training data (≥ 100 known ligands).

---

## Pipeline Stages

```
RelB 352-gene signature (Excel)
         │
         ▼
[01] Target Acquisition
     175 high/medium genes → ChEMBL target mapping → 26 druggable targets
         │
         ▼
[02] Ligand Retrieval
     ChEMBL bioactivity data (IC50, Ki, binding assays)
     Binary labels: active (pChEMBL ≥ 6) / inactive (pChEMBL < 5)
         │
         ▼
[03] Feature Engineering
     Lipinski Ro5 filter → ECFP4 Morgan fingerprints (radius=2, 2048 bits)
         │
         ▼
[04] Classifier Training
     Random Forest (200 trees, class_weight=balanced)
     Stratified 5-fold cross-validation
         │
         ▼
[05] Virtual Screening
     3,311 FDA-approved compounds (ChEMBL max_phase=4)
     Ranked by p(active)
         │
         ▼
    Top 50 ranked hits
```

---

## Druggable RelB Targets (26 identified)

| Gene | ChEMBL ID | Priority | RelB logFC | Known Ligands | Biological Role |
|------|-----------|----------|-----------|--------------|-----------------|
| KLKB1 | CHEMBL2000 | medium | -1.75 | 3,592 | Plasma kallikrein / coagulation |
| CYP19A1 | CHEMBL1978 | medium | -2.44 | 2,991 | Aromatase / estrogen synthesis |
| MMP13 | CHEMBL280 | **high** | -2.60 | 2,829 | Matrix metalloproteinase / ECM remodeling |
| PRL | CHEMBL2014 | medium | -2.67 | 2,612 | Prolactin / cytokine signaling |
| GNAL | CHEMBL4026 | medium | -1.50 | 1,163 | G-protein signaling |
| ABCG2 | CHEMBL5393 | **high** | -2.37 | 1,045 | Multidrug efflux transporter |
| AKR1C3 | CHEMBL4681 | medium | -1.78 | 1,009 | Steroid metabolism |
| TRAC | CHEMBL1825 | **high** | -2.14 | 905 | T-cell receptor / immune |
| PPARA | CHEMBL239 | medium | -1.52 | 889 | Peroxisome proliferator receptor |
| CYP17A1 | CHEMBL3522 | medium | -3.21 | 729 | Steroidogenesis |
| MMP12 | CHEMBL4393 | **high** | -3.12 | 712 | Macrophage elastase / invasion |
| HSD17B2 | CHEMBL2789 | medium | -1.58 | 671 | 17β-HSD / steroid inactivation |
| ALB | CHEMBL2083 | **high** | -4.79 | 603 | Serum albumin |
| CYP1A1 | CHEMBL2231 | **high** | -2.40 | 569 | Xenobiotic metabolism |
| ANPEP | CHEMBL1907 | **high** | -3.29 | 446 | Aminopeptidase N / angiogenesis |
| PLA2G7 | CHEMBL3514 | medium | -1.55 | 437 | Phospholipase / lipid signaling |
| GIF | CHEMBL2085 | medium | -1.87 | 399 | Gastric intrinsic factor |
| LIPC | CHEMBL2127 | **high** | -3.68 | 378 | Hepatic lipase |
| CES1 | CHEMBL2265 | **high** | -4.46 | 350 | Carboxylesterase / prodrug activation |
| MGAM | CHEMBL2074 | **high** | -3.76 | 313 | Maltase-glucoamylase |
| P2RX4 | CHEMBL2104 | medium | -1.51 | 291 | Purinergic receptor / inflammation |
| TNNI3 | CHEMBL5260 | **high** | -2.41 | 212 | Cardiac troponin I |
| ALPI | CHEMBL5573 | **high** | -4.68 | 163 | Intestinal alkaline phosphatase |
| FABP1 | CHEMBL3344 | **high** | -5.85 | 128 | Fatty acid binding protein |
| HP | CHEMBL1861 | **high** | -2.96 | 113 | Haptoglobin |
| CMKLR1 | CHEMBL3540 | medium | -1.95 | 103 | Chemerin receptor / inflammation |

---

## Model Performance

| Metric | 5-Fold CV | Held-out Test |
|--------|-----------|---------------|
| **ROC-AUC** | 0.960 ± 0.005 | **0.958** |
| **PR-AUC** | 0.992 ± 0.001 | **0.990** |

Training set: compound–activity pairs from 19 RelB-pathway targets (after Lipinski filter and SMILES validation). Random Forest, 200 estimators, class_weight=balanced to correct active/inactive imbalance.

---

## Top 10 Screening Hits (FDA-Approved)

| Rank | Drug | ChEMBL ID | p(active) | Drug Class | HGSOC Relevance |
|------|------|-----------|-----------|------------|-----------------|
| **1** | **Buprenorphine** | CHEMBL511142 | **0.9950** | Opioid partial agonist | Reported anti-proliferative activity in ovarian cancer cells |
| 2 | Flurbiprofen | CHEMBL563 | 0.9800 | NSAID / COX inhibitor | Anti-inflammatory; COX-2 overexpressed in HGSOC |
| 3 | Estrone | CHEMBL1405 | 0.9650 | Estrogen | CYP19A1 (aromatase) substrate — validates aromatase axis |
| 4 | **Anastrozole** | CHEMBL1399 | 0.9550 | Aromatase inhibitor | FDA-approved; breast cancer → HGSOC repurposing candidate |
| 5 | **Letrozole** | CHEMBL1444 | 0.9550 | Aromatase inhibitor | FDA-approved; being evaluated in ovarian cancer trials |
| 6 | **Exemestane** | CHEMBL1200374 | 0.9550 | Aromatase inhibitor | FDA-approved; steroidal; irreversible CYP19A1 inhibitor |
| 7 | Naproxen | CHEMBL154 | 0.9500 | NSAID | NF-κB inhibitory activity reported |
| 8 | Nalmefene | CHEMBL982 | 0.9500 | Opioid antagonist | Low-dose naltrexone studied in ovarian cancer |
| 9 | Naltrexone | CHEMBL19019 | 0.9350 | Opioid antagonist | Low-dose naltrexone — immune modulation in cancer |
| 10 | Nalmefene HCl | CHEMBL1201152 | 0.9350 | Opioid antagonist | Salt form of rank 8 |

---

## Key Biological Insights

### 1. Aromatase axis emerges as top repurposing opportunity
Three FDA-approved aromatase inhibitors (anastrozole, letrozole, exemestane) rank 4–6 with p(active) = 0.955. CYP19A1 (aromatase) was the second most ligand-rich RelB-dependent target (2,991 compounds, logFC = -2.44). This is biologically coherent: estrogen signaling and NF-κB/RelB share bidirectional crosstalk, and aromatase inhibitors are already evaluated in hormone receptor-positive ovarian cancer. **Priority experimental target.**

### 2. Opioid receptor ligands cluster in the top hits
Buprenorphine (#1), nalmefene (#8, #10), naltrexone (#9) suggest opioid receptor involvement. Buprenorphine has published anti-proliferative activity in OVCAR cell lines (δ-opioid receptor pathway). Delta-opioid receptors modulate NF-κB signaling in tumor cells.

### 3. COX/NF-κB anti-inflammatory axis (flurbiprofen, naproxen)
NSAIDs at ranks 2 and 7 are consistent with the known NF-κB-suppressing effects of COX inhibition. Flurbiprofen also inhibits microsomal PGE2 synthesis — relevant to HGSOC ascites-driven inflammation.

### 4. MMP inhibition as a secondary therapeutic strategy
MMP13 and MMP12 were among the highest-priority druggable targets (logFC = -2.60 and -3.12). Both are RelB-suppressed, suggesting that RelB loss enables matrix invasion. MMP inhibitors may synergize with platinum re-sensitization.

---

## Repository Structure

```
HGSOC-virtual-screening/
├── config.yaml                        # pipeline parameters
├── requirements.txt                   # Python dependencies
├── run_pipeline.py                    # single-command end-to-end runner
├── create_notebooks.py                # generates .ipynb files
├── data/
│   ├── inputs/
│   │   ├── relb_gene_signature_352.csv    # 352-gene RelB signature (source)
│   │   └── example_biomarker_targets.csv  # placeholder
│   └── processed/
│       ├── target_manifest.csv        # 26 druggable targets
│       ├── combined_bioactivity.csv   # all ChEMBL training data
│       ├── screening_hits.csv         # top 50 FDA-approved hits
│       └── models/                    # trained classifier (git-ignored)
├── src/
│   ├── acquisition.py                 # ChEMBL data fetching
│   ├── features.py                    # ECFP4 fingerprints + Lipinski
│   ├── models.py                      # Random Forest classifier
│   └── screening.py                   # screening pipeline
└── notebooks/
    ├── 01_target_acquisition.ipynb
    ├── 02_ligand_retrieval.ipynb
    ├── 03_feature_engineering.ipynb
    ├── 04_classifier_training.ipynb
    └── 05_virtual_screening.ipynb
```

---

## Quickstart

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run full pipeline (uses cached ChEMBL data on repeat runs)
~/miniconda3/bin/python3 run_pipeline.py

# 3. Or explore interactively
jupyter lab
```

---

## Dependencies

- `chembl_webresource_client` — ChEMBL REST API client
- `rdkit` — cheminformatics (ECFP4 fingerprints, Lipinski filter)
- `scikit-learn` — Random Forest classifier
- `umap-learn` — chemical space visualization
- `pandas`, `numpy`, `matplotlib` — data processing and plotting

---

## Citation

If you use this pipeline, please cite:

> Lujano Olazaba O. *HGSOC Virtual Screening Pipeline: AI-driven drug repurposing for RelB-dependent targets in high-grade serous ovarian cancer.* TARVIA Lab, 2026. https://github.com/TARVIA-lab/HGSOC-virtual-screening

---

*Built with the methodology from: Manning — Build AI Drug Discovery Pipelines (MEAP, 2026)*
