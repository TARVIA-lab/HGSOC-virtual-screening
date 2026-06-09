"""
Run once to generate the five starter notebooks:
    python create_notebooks.py
"""

import json, uuid
from pathlib import Path

NB_DIR = Path("notebooks")
NB_DIR.mkdir(exist_ok=True)

_cid = 0

def cid():
    global _cid
    _cid += 1
    return f"{_cid:08x}-{uuid.uuid4().hex[:4]}-{uuid.uuid4().hex[:4]}-{uuid.uuid4().hex[:4]}-{uuid.uuid4().hex[:12]}"

def md(source):
    return {"cell_type": "markdown", "id": cid(), "metadata": {}, "source": source}

def code(source):
    return {"cell_type": "code", "execution_count": None, "id": cid(), "metadata": {}, "outputs": [], "source": source}

def nb(cells):
    return {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.10.0"},
        },
        "cells": cells,
    }


# ─── Notebook 1: Target Acquisition ──────────────────────────────────────────

nb1 = nb([
    md("""# 01 — Target Acquisition

Map RelB-dependent HGSOC biomarker genes to ChEMBL targets and rank by druggability.

**Input:** `../data/inputs/example_biomarker_targets.csv`
**Output:** `../data/processed/target_manifest.csv`"""),

    code("""\
import sys, yaml
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, '..')
from src.acquisition import get_chembl_id_for_gene

with open('../config.yaml') as f:
    config = yaml.safe_load(f)

Path('../data/processed').mkdir(parents=True, exist_ok=True)
print("Ready.")"""),

    md("""## Load biomarker pipeline output

Replace `example_biomarker_targets.csv` with the ranked output from your
`biomarker-discovery-pipeline`. Required column: `gene_symbol`.
Optional columns: `priority`, `survival_hr`, `survival_p`."""),

    code("""\
targets_df = pd.read_csv('../data/inputs/example_biomarker_targets.csv')
print(f"{len(targets_df)} candidate targets loaded")
targets_df"""),

    md("## Map gene symbols → ChEMBL target IDs"),

    code("""\
if 'chembl_id' not in targets_df.columns:
    targets_df['chembl_id'] = None

unmapped = targets_df['chembl_id'].isna()
print(f"Mapping {unmapped.sum()} genes to ChEMBL IDs…")

for idx, row in targets_df[unmapped].iterrows():
    targets_df.at[idx, 'chembl_id'] = get_chembl_id_for_gene(row['gene_symbol'])

targets_df = targets_df.dropna(subset=['chembl_id'])
print(f"\\n{len(targets_df)} targets successfully mapped")
targets_df[['gene_symbol', 'chembl_id', 'priority']]"""),

    md("""## Druggability proxy: count known ChEMBL ligands

Targets with ≥ 100 known ligands have enough data to train a reliable classifier."""),

    code("""\
from chembl_webresource_client.new_client import new_client
activity_api = new_client.activity

def count_ligands(chembl_id):
    acts = activity_api.filter(
        target_chembl_id=chembl_id,
        standard_type__in=['IC50', 'Ki'],
        pchembl_value__isnull=False,
    ).only(['molecule_chembl_id'])
    return len({a['molecule_chembl_id'] for a in acts})

targets_df['n_ligands'] = targets_df['chembl_id'].apply(count_ligands)
targets_df = targets_df.sort_values('n_ligands', ascending=False)

fig, ax = plt.subplots(figsize=(9, 4))
colors = ['#2ecc71' if n >= 100 else '#e74c3c' for n in targets_df['n_ligands']]
ax.barh(targets_df['gene_symbol'], targets_df['n_ligands'], color=colors)
ax.axvline(100, color='black', ls='--', label='Min threshold (100)')
ax.set_xlabel('Known ChEMBL ligands')
ax.set_title('Druggability proxy: HGSOC target library')
ax.legend()
plt.tight_layout()
plt.savefig('../data/processed/druggability_bar.png', dpi=150)
plt.show()"""),

    code("""\
MIN_LIGANDS = 100
manifest = targets_df[targets_df['n_ligands'] >= MIN_LIGANDS].copy()
manifest.to_csv('../data/processed/target_manifest.csv', index=False)
print(f"Saved {len(manifest)} screeable targets → data/processed/target_manifest.csv")
manifest"""),
])


# ─── Notebook 2: Ligand Retrieval ────────────────────────────────────────────

nb2 = nb([
    md("""# 02 — Ligand Retrieval

Fetch compound–activity data from ChEMBL for each prioritized target and assign binary labels.

**Input:** `../data/processed/target_manifest.csv`
**Output:** `../data/processed/combined_bioactivity.csv`"""),

    code("""\
import sys, yaml
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '..')
from src.acquisition import fetch_bioactivity_data, label_activity

with open('../config.yaml') as f:
    config = yaml.safe_load(f)

manifest = pd.read_csv('../data/processed/target_manifest.csv')
print(f"Processing {len(manifest)} targets")
manifest[['gene_symbol', 'chembl_id', 'n_ligands']]"""),

    md("## Fetch bioactivity data from ChEMBL binding assays"),

    code("""\
all_dfs = []
for _, row in manifest.iterrows():
    df = fetch_bioactivity_data(
        chembl_id=row['chembl_id'],
        activity_types=config['acquisition']['activity_types'],
        assay_type=config['acquisition']['assay_type'],
        max_compounds=config['acquisition']['max_compounds_per_target'],
    )
    if not df.empty:
        df['target_gene'] = row['gene_symbol']
        all_dfs.append(df)
        print(f"  {row['gene_symbol']:10s}: {len(df):>5d} compounds")

combined = pd.concat(all_dfs, ignore_index=True)
print(f"\\nTotal: {len(combined)} compound–target pairs across {len(all_dfs)} targets")"""),

    md("""## Binary activity labels

pChEMBL ≥ 6 (IC₅₀ ≤ 1 µM) → **active = 1**
pChEMBL < 5 (IC₅₀ > 10 µM) → **inactive = 0**
Gray zone (5–6) is dropped to keep class boundaries sharp."""),

    code("""\
labeled = label_activity(
    combined,
    active_threshold=config['acquisition']['pchembl_cutoff_active'],
    inactive_threshold=config['acquisition']['pchembl_cutoff_inactive'],
)

n_active   = labeled['active'].sum()
n_inactive = (labeled['active'] == 0).sum()
print(f"Active: {n_active} | Inactive: {n_inactive} | Ratio: {n_inactive/n_active:.1f}:1")

labeled.to_csv('../data/processed/combined_bioactivity.csv', index=False)"""),

    md("## Exploratory data analysis"),

    code("""\
fig, axes = plt.subplots(1, 2, figsize=(13, 4))

# pChEMBL distribution by target
for gene, grp in labeled.groupby('target_gene'):
    axes[0].hist(grp['pchembl_value'], bins=25, alpha=0.5, label=gene)
axes[0].axvline(config['acquisition']['pchembl_cutoff_active'],   color='green', ls='--', label='Active threshold')
axes[0].axvline(config['acquisition']['pchembl_cutoff_inactive'], color='red',   ls='--', label='Inactive threshold')
axes[0].set_xlabel('pChEMBL value')
axes[0].set_title('Activity distributions by target')
axes[0].legend(fontsize=7)

# Class balance per target
balance = labeled.groupby(['target_gene', 'active']).size().unstack(fill_value=0)
balance.columns = ['Inactive', 'Active']
balance.plot(kind='bar', ax=axes[1], color=['#e74c3c', '#2ecc71'])
axes[1].set_title('Class balance per target')
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('../data/processed/eda_activity.png', dpi=150)
plt.show()"""),
])


# ─── Notebook 3: Feature Engineering ─────────────────────────────────────────

nb3 = nb([
    md("""# 03 — Feature Engineering

Compute Morgan (ECFP4) fingerprints, apply Lipinski's Rule of Five, and visualize
the chemical space with UMAP.

**Input:** `../data/processed/combined_bioactivity.csv`
**Output:** `features_X.npy`, `features_y.npy`, `features_meta.csv`"""),

    code("""\
import sys, yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '..')
from src.features import lipinski_filter, compute_fingerprints, compute_rdkit_descriptors

with open('../config.yaml') as f:
    config = yaml.safe_load(f)

df = pd.read_csv('../data/processed/combined_bioactivity.csv')
print(f"Loaded {len(df)} labeled compounds")"""),

    md("""## Lipinski's Rule of Five

Filters out compounds unlikely to be orally bioavailable leads.
Rule: MW ≤ 500, LogP ≤ 5, H-bond donors ≤ 5, H-bond acceptors ≤ 10."""),

    code("""\
df_ro5 = lipinski_filter(df)
print(f"After Ro5 filter: {len(df_ro5)} / {len(df)} compounds retained ({len(df_ro5)/len(df):.1%})")"""),

    md("""## Compute Morgan fingerprints (ECFP4)

Radius = 2, 2048 bits. Each set bit encodes a circular substructure
centred on an atom out to radius 2 bonds. These are the features the
ML classifier will train on."""),

    code("""\
fp_cfg = config['features']
X, valid_df = compute_fingerprints(
    df_ro5,
    radius=fp_cfg['morgan_radius'],
    n_bits=fp_cfg['morgan_n_bits'],
)
y = valid_df['active'].values

print(f"Feature matrix: {X.shape}  |  Active: {y.sum()}  |  Inactive: {(y==0).sum()}")

np.save('../data/processed/features_X.npy', X)
np.save('../data/processed/features_y.npy', y)
valid_df.to_csv('../data/processed/features_meta.csv', index=False)
print("Saved features_X.npy, features_y.npy, features_meta.csv")"""),

    md("""## UMAP: visualise chemical space

Project 2048-D fingerprints to 2D. Well-separated active/inactive clusters
suggest the fingerprints carry discriminative signal."""),

    code("""\
from umap import UMAP

reducer = UMAP(n_components=2, random_state=42, n_jobs=1)
embedding = reducer.fit_transform(X)

fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(embedding[:, 0], embedding[:, 1], c=y, cmap='RdYlGn', alpha=0.4, s=8)
plt.colorbar(sc, label='Active (1) / Inactive (0)')
ax.set_title('Chemical space: UMAP of ECFP4 fingerprints')
ax.set_xlabel('UMAP-1')
ax.set_ylabel('UMAP-2')
plt.tight_layout()
plt.savefig('../data/processed/umap_chemical_space.png', dpi=150)
plt.show()"""),

    md("## Physicochemical descriptor distributions"),

    code("""\
descs = compute_rdkit_descriptors(valid_df)
descs['active'] = y

fig, axes = plt.subplots(2, 3, figsize=(14, 7))
for ax, col in zip(axes.flatten(), ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'RotBonds']):
    for label, color in [(0, '#e74c3c'), (1, '#2ecc71')]:
        subset = descs[descs['active'] == label][col].dropna()
        ax.hist(subset, bins=30, alpha=0.5, color=color, label='Active' if label else 'Inactive')
    ax.set_title(col)
    ax.legend(fontsize=8)

plt.suptitle('Physicochemical properties: active vs inactive', y=1.01)
plt.tight_layout()
plt.savefig('../data/processed/descriptor_distributions.png', dpi=150)
plt.show()"""),
])


# ─── Notebook 4: Classifier Training ─────────────────────────────────────────

nb4 = nb([
    md("""# 04 — Classifier Training

Train a Random Forest on ECFP4 fingerprints. Evaluate with stratified 5-fold
cross-validation and a held-out test set.

**Input:** `features_X.npy`, `features_y.npy`
**Output:** `../data/processed/models/combined_model.joblib`"""),

    code("""\
import sys, yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from pathlib import Path

sys.path.insert(0, '..')
from src.models import (
    build_classifier, cross_validate, train_and_evaluate,
    save_model, get_feature_importance,
)

with open('../config.yaml') as f:
    config = yaml.safe_load(f)

Path('../data/processed/models').mkdir(parents=True, exist_ok=True)

X = np.load('../data/processed/features_X.npy')
y = np.load('../data/processed/features_y.npy')

print(f"Feature matrix : {X.shape}")
print(f"Active         : {y.sum()} ({y.mean():.1%})")
print(f"Inactive       : {(y==0).sum()}")"""),

    md("## Train / test split"),

    code("""\
cfg = config['models']
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=cfg['test_size'], stratify=y, random_state=cfg['random_state'],
)
print(f"Train: {len(X_train)} | Test: {len(X_test)}")"""),

    md("""## Stratified cross-validation

Estimate generalisation before touching the test set.
`class_weight='balanced'` corrects for the active/inactive imbalance."""),

    code("""\
model = build_classifier(model_type=cfg['type'], n_estimators=cfg['n_estimators'])
cv = cross_validate(model, X_train, y_train, cv_folds=cfg['cv_folds'])

print(f"CV  ROC-AUC : {cv['roc_auc_mean']:.3f} ± {cv['roc_auc_std']:.3f}")
print(f"CV  PR-AUC  : {cv['pr_auc_mean']:.3f} ± {cv['pr_auc_std']:.3f}")"""),

    md("## Final model training and held-out test evaluation"),

    code("""\
metrics = train_and_evaluate(model, X_train, y_train, X_test, y_test)

print(f"Test ROC-AUC : {metrics['roc_auc']:.3f}")
print(f"Test PR-AUC  : {metrics['pr_auc']:.3f}")
print()
print(metrics['classification_report'])

model_path = '../data/processed/models/combined_model.joblib'
save_model(model, model_path)"""),

    md("## ROC and Precision–Recall curves"),

    code("""\
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

fpr, tpr, _ = metrics['roc_curve']
axes[0].plot(fpr, tpr, lw=2, color='steelblue', label=f"AUC = {metrics['roc_auc']:.3f}")
axes[0].plot([0, 1], [0, 1], 'k--')
axes[0].set_xlabel('False Positive Rate')
axes[0].set_ylabel('True Positive Rate')
axes[0].set_title('ROC Curve')
axes[0].legend()

prec, rec, _ = metrics['pr_curve']
axes[1].plot(rec, prec, lw=2, color='darkorange', label=f"AUC = {metrics['pr_auc']:.3f}")
axes[1].set_xlabel('Recall')
axes[1].set_ylabel('Precision')
axes[1].set_title('Precision–Recall Curve')
axes[1].legend()

plt.tight_layout()
plt.savefig('../data/processed/model_evaluation_curves.png', dpi=150)
plt.show()"""),

    md("""## Feature importance (top 20 Morgan bits)

Each bar is a fingerprint bit. High importance = that circular substructure
pattern strongly separates active from inactive compounds."""),

    code("""\
importances = get_feature_importance(model, top_n=20)

fig, ax = plt.subplots(figsize=(8, 5))
importances.sort_values().plot(kind='barh', ax=ax, color='steelblue')
ax.set_xlabel('Feature importance')
ax.set_ylabel('ECFP4 bit index')
ax.set_title('Top 20 most discriminative substructure bits')
plt.tight_layout()
plt.savefig('../data/processed/feature_importance.png', dpi=150)
plt.show()"""),
])


# ─── Notebook 5: Virtual Screening ───────────────────────────────────────────

nb5 = nb([
    md("""# 05 — Virtual Screening

Score the FDA-approved compound library with the trained classifier.
Rank hits by p(active) and inspect their 2-D structures.

**Input:** trained model + ChEMBL FDA-approved library
**Output:** `../data/processed/screening_hits.csv`, `top_hits_structures.png`"""),

    code("""\
import sys, yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, '..')
from src.screening import fetch_fda_approved_library, screen_library
from src.features import compute_fingerprints, lipinski_filter
from src.models import load_model

with open('../config.yaml') as f:
    config = yaml.safe_load(f)"""),

    md("""## Download FDA-approved compound library

Pulls `max_phase = 4` small molecules from ChEMBL (~2,500 compounds).
The CSV is cached after the first run."""),

    code("""\
library_path = '../data/processed/fda_library.csv'

if Path(library_path).exists():
    library_df = pd.read_csv(library_path)
    print(f"Loaded cached library: {len(library_df)} compounds")
else:
    library_df = fetch_fda_approved_library(library_path)
    print(f"Downloaded and cached: {len(library_df)} compounds")

library_df.head()"""),

    md("""## Screen against the trained model

Pipeline: Lipinski filter → ECFP4 fingerprints → predict p(active) → rank."""),

    code("""\
model_path = '../data/processed/models/combined_model.joblib'

hits = screen_library(
    library_df=library_df,
    model_path=model_path,
    apply_lipinski=True,
    top_n=config['screening']['top_n_hits'],
)

hits.to_csv('../data/processed/screening_hits.csv', index=False)
print(f"Top {len(hits)} hits saved to data/processed/screening_hits.csv")
hits[['rank', 'pref_name', 'molecule_chembl_id', 'p_active']].head(20)"""),

    md("## Score distribution across the full library"),

    code("""\
model = load_model(model_path)
lib_filt = lipinski_filter(library_df)
X_lib, _ = compute_fingerprints(lib_filt)
all_probs = model.predict_proba(X_lib)[:, 1]

cutoff = hits['p_active'].min()

fig, ax = plt.subplots(figsize=(9, 4))
ax.hist(all_probs, bins=60, color='steelblue', alpha=0.75)
ax.axvline(cutoff, color='red', ls='--', label=f'Top-{len(hits)} cutoff ({cutoff:.3f})')
ax.set_xlabel('p(active)')
ax.set_ylabel('Compound count')
ax.set_title('Predicted activity scores: FDA-approved library')
ax.legend()
plt.tight_layout()
plt.savefig('../data/processed/screening_score_dist.png', dpi=150)
plt.show()"""),

    md("## Visualise top hit structures"),

    code("""\
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display

top12 = hits.head(12)
mols   = [Chem.MolFromSmiles(s) for s in top12['canonical_smiles']]
labels = [
    f"{row['pref_name'] or row['molecule_chembl_id']}\\np={row['p_active']:.3f}"
    for _, row in top12.iterrows()
]

img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(320, 260), legends=labels)
img.save('../data/processed/top_hits_structures.png')
display(img)"""),

    md("""## Summary

| Notebook | Output |
|---|---|
| 01 Target Acquisition | `target_manifest.csv` |
| 02 Ligand Retrieval | `combined_bioactivity.csv` |
| 03 Feature Engineering | `features_X.npy`, UMAP plot |
| 04 Classifier Training | `models/combined_model.joblib`, ROC/PR curves |
| 05 Virtual Screening | `screening_hits.csv`, `top_hits_structures.png` |

**Next step →** share the top 10–20 hits with your HGSOC collaborators.
Prioritise compounds already approved for other indications (drug repurposing)."""),
])


# ─── Write all notebooks ──────────────────────────────────────────────────────

notebooks = {
    "01_target_acquisition.ipynb": nb1,
    "02_ligand_retrieval.ipynb":   nb2,
    "03_feature_engineering.ipynb": nb3,
    "04_classifier_training.ipynb": nb4,
    "05_virtual_screening.ipynb":  nb5,
}

for name, notebook in notebooks.items():
    path = NB_DIR / name
    with open(path, "w", encoding="utf-8") as f:
        json.dump(notebook, f, indent=1, ensure_ascii=False)
    print(f"Created: {path}")
