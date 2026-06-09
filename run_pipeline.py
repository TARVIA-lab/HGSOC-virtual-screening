"""
End-to-end HGSOC virtual screening pipeline.
Run with: ~/miniconda3/bin/python3 run_pipeline.py
"""

import sys, time, logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.insert(0, str(Path(__file__).parent))
from src.acquisition import get_chembl_id_for_gene, fetch_bioactivity_data, label_activity
from src.features import lipinski_filter, compute_fingerprints
from src.models import build_classifier, cross_validate, train_and_evaluate, save_model
from src.screening import fetch_fda_approved_library, screen_library

logging.basicConfig(level=logging.WARNING)  # suppress verbose API logs
Path("data/processed/models").mkdir(parents=True, exist_ok=True)

t0 = time.time()


# ── STAGE 1: Load RelB gene signature ────────────────────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 1  Load RelB gene signature")
print("═══════════════════════════════════════════════")

genes_df = pd.read_csv("data/inputs/relb_gene_signature_352.csv")
high_med = genes_df[genes_df["priority"].isin(["high", "medium"])].copy()
print(f"  {len(genes_df)} total RelB genes  →  {len(high_med)} high/medium priority selected")


# ── STAGE 2: Map gene symbols → ChEMBL IDs ───────────────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 2  ChEMBL target mapping")
print("═══════════════════════════════════════════════")

cache_path = Path("data/processed/chembl_id_cache.csv")
if cache_path.exists():
    cache = pd.read_csv(cache_path).set_index("gene_symbol")["chembl_id"].to_dict()
    print(f"  Loaded cached mapping for {len(cache)} genes")
else:
    cache = {}

mapped, failed = [], []
for _, row in high_med.iterrows():
    gene = row["gene_symbol"]
    if gene in cache:
        cid = cache[gene]
    else:
        cid = get_chembl_id_for_gene(gene)
        cache[gene] = cid
    if cid:
        mapped.append({"gene_symbol": gene, "chembl_id": cid, "priority": row["priority"],
                       "relb_logfc": row["relb_logfc"]})
    else:
        failed.append(gene)
    sys.stdout.write(f"\r  Mapped {len(mapped)+len(failed)}/{len(high_med)} genes...")
    sys.stdout.flush()

# Save cache
pd.DataFrame([{"gene_symbol": k, "chembl_id": v} for k, v in cache.items()]).to_csv(cache_path, index=False)

mapped_df = pd.DataFrame(mapped)
print(f"\n  {len(mapped_df)} mapped  |  {len(failed)} unmapped (non-protein/lncRNA)")


# ── STAGE 3: Filter by druggability (ligand count) ───────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 3  Druggability filter (≥100 known ligands)")
print("═══════════════════════════════════════════════")

from chembl_webresource_client.new_client import new_client
activity_api = new_client.activity

ligand_cache_path = Path("data/processed/ligand_count_cache.csv")
if ligand_cache_path.exists():
    lc = pd.read_csv(ligand_cache_path).set_index("chembl_id")["n_ligands"].to_dict()
    print(f"  Loaded cached ligand counts for {len(lc)} targets")
else:
    lc = {}

counts = []
for i, row in mapped_df.iterrows():
    cid = row["chembl_id"]
    if cid not in lc:
        acts = activity_api.filter(
            target_chembl_id=cid,
            standard_type__in=["IC50", "Ki"],
            pchembl_value__isnull=False,
        ).only(["molecule_chembl_id"])
        lc[cid] = len({a["molecule_chembl_id"] for a in acts})
    counts.append(lc[cid])
    sys.stdout.write(f"\r  Checked {len(counts)}/{len(mapped_df)}... {row['gene_symbol']} → {lc[cid]} ligands")
    sys.stdout.flush()

pd.DataFrame([{"chembl_id": k, "n_ligands": v} for k, v in lc.items()]).to_csv(ligand_cache_path, index=False)

mapped_df["n_ligands"] = counts
druggable = mapped_df[mapped_df["n_ligands"] >= 100].sort_values("n_ligands", ascending=False)
druggable.to_csv("data/processed/target_manifest.csv", index=False)

print(f"\n\n  {'Gene':<14} {'ChEMBL ID':<18} {'Priority':<10} {'N Ligands':>10}")
print(f"  {'-'*54}")
for _, r in druggable.iterrows():
    print(f"  {r['gene_symbol']:<14} {r['chembl_id']:<18} {r['priority']:<10} {r['n_ligands']:>10,}")


# ── STAGE 4: Fetch bioactivity data ──────────────────────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 4  Fetch compound–activity data from ChEMBL")
print("═══════════════════════════════════════════════")

all_dfs = []
for _, row in druggable.iterrows():
    gene, cid = row["gene_symbol"], row["chembl_id"]
    df = fetch_bioactivity_data(cid, activity_types=["IC50", "Ki"], assay_type="B", max_compounds=3000)
    if df.empty:
        continue
    df = label_activity(df, active_threshold=6.0, inactive_threshold=5.0)
    df["target_gene"] = gene
    all_dfs.append(df)
    n_act = df["active"].sum()
    n_ina = (df["active"] == 0).sum()
    print(f"  {gene:<14} {n_act:>5} active  {n_ina:>5} inactive")

combined = pd.concat(all_dfs, ignore_index=True)
combined.to_csv("data/processed/combined_bioactivity.csv", index=False)
print(f"\n  Total: {len(combined):,} labeled compound–target pairs across {len(all_dfs)} targets")


# ── STAGE 5: Feature engineering ─────────────────────────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 5  ECFP4 fingerprints + Lipinski filter")
print("═══════════════════════════════════════════════")

df_ro5 = lipinski_filter(combined)
print(f"  Lipinski filter: {len(df_ro5):,} / {len(combined):,} compounds retained")

X, valid_df = compute_fingerprints(df_ro5, radius=2, n_bits=2048)
y = valid_df["active"].values

np.save("data/processed/features_X.npy", X)
np.save("data/processed/features_y.npy", y)
valid_df.to_csv("data/processed/features_meta.csv", index=False)
print(f"  Feature matrix: {X.shape}  |  Active: {y.sum():,}  |  Inactive: {(y==0).sum():,}")


# ── STAGE 6: Train classifier ─────────────────────────────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 6  Random Forest classifier (5-fold CV)")
print("═══════════════════════════════════════════════")

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

model = build_classifier("random_forest", n_estimators=200)
cv = cross_validate(model, X_train, y_train, cv_folds=5)
print(f"  CV  ROC-AUC : {cv['roc_auc_mean']:.3f} ± {cv['roc_auc_std']:.3f}")
print(f"  CV  PR-AUC  : {cv['pr_auc_mean']:.3f} ± {cv['pr_auc_std']:.3f}")

metrics = train_and_evaluate(model, X_train, y_train, X_test, y_test)
print(f"  Test ROC-AUC: {metrics['roc_auc']:.3f}")
print(f"  Test PR-AUC : {metrics['pr_auc']:.3f}")

model_path = "data/processed/models/combined_model.joblib"
save_model(model, model_path)


# ── STAGE 7: Virtual screening — FDA-approved library ─────────────────────────
print("\n═══════════════════════════════════════════════")
print(" STAGE 7  Virtual screening: FDA-approved library")
print("═══════════════════════════════════════════════")

library_path = "data/processed/fda_library.csv"
if Path(library_path).exists():
    library_df = pd.read_csv(library_path)
    print(f"  Loaded cached FDA library: {len(library_df):,} compounds")
else:
    print("  Downloading FDA-approved compounds from ChEMBL...")
    library_df = fetch_fda_approved_library(library_path)

hits = screen_library(library_df, model_path, apply_lipinski=True, top_n=50)
hits.to_csv("data/processed/screening_hits.csv", index=False)


# ── RESULTS ───────────────────────────────────────────────────────────────────
elapsed = (time.time() - t0) / 60
print("\n═══════════════════════════════════════════════")
print(f" DONE  ({elapsed:.1f} min)")
print("═══════════════════════════════════════════════")

print("\n TOP 10 PREDICTED HITS (FDA-approved compounds)\n")
print(f"  {'Rank':<5} {'Name':<35} {'ChEMBL ID':<18} {'p(active)':>10}")
print(f"  {'-'*70}")
for _, row in hits.head(10).iterrows():
    name = (row.get("pref_name") or row["molecule_chembl_id"])[:33]
    print(f"  {int(row['rank']):<5} {name:<35} {row['molecule_chembl_id']:<18} {row['p_active']:>10.4f}")

top = hits.iloc[0]
print(f"\n ★  TOP HIT: {top.get('pref_name') or top['molecule_chembl_id']}")
print(f"    ChEMBL ID : {top['molecule_chembl_id']}")
print(f"    p(active) : {top['p_active']:.4f}")
print(f"    SMILES    : {top['canonical_smiles'][:80]}")
print()
