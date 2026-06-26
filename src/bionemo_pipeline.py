#!/usr/bin/env python3
"""
RelB Drug Discovery Pipeline
GenMol → DiffDock → Boltz-2

Target: Human RelB (NF-κB subunit, UniProt Q01201)
Goal: Generate novel small molecules, dock them, and predict binding affinity.
"""

import requests
import json
import os
import time
import csv
from pathlib import Path

NGC_API_KEY = os.environ.get("NGC_API_KEY") or os.environ["NVIDIA_API_KEY"]

HEADERS = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {NGC_API_KEY}",
}

# Human RelB Rel Homology Domain (RHD) — the druggable DNA-binding region (residues 120-400)
RELB_SEQUENCE = (
    "MDSLESSFMS QDFSRHQARQ GPYENEYSEF QDSPSGSEQS PHNVDALQGR"
    "EECLRAGKPS ISDIIQDQFE RALTEPHHEQ AMRKFAQHQQ PTQHRSEKRN"
    "RLAPQWQLHH VQAIFPPNWH KFAADRRFRF RQLVEAQAGK LPECLHISAA"
    "YRGTKMTEQE LQERYDQLLR ESGFNLQWKQ PRGFHLLTKN SPDQYPACEA"
    "HMRNLYEGEC VKGMKCGQEE LNITIVKLQS ATQNFKFTKP AAVLKRTCQE"
    "GVDAQAAYPF LSTPQVALFG SGDGHGHHDG SLDSLGTLSA PLVTAPPVES"
    "FCQHAFIRNR FKESGEFSYA LSAVKQPTVQ TIKPRVDKQG SPVYTSASSD"
    "VTHVFKRRGS TSTKGFPQRE GMEQKLISEE DLSSLAGTPV SSLGIISSAS"
    "PTDIDSPFGH SAGLSRTTAQ SESLVHQEPS PPTITMSELS VSKLHSSAQK"
    "DTAPYSGQSS AQVSISGHGT SLPVASSGQT LHQHDQQHKV HIQRPQQQGL"
    "SPFASSTFSS SSLSFPQHAD YLLHTHAPAP PAPCPHDIYS DLLDTAPPSP"
).replace(" ", "").replace("\n", "")

# Clean sequence — remove any non-amino-acid characters
RELB_SEQUENCE = "".join(c for c in RELB_SEQUENCE if c.isalpha())

print(f"RelB sequence length: {len(RELB_SEQUENCE)} residues")
print(f"First 50 residues: {RELB_SEQUENCE[:50]}...")

# ═══════════════════════════════════════════════════════════════
# STEP 1: Generate molecules with GenMol
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 1: Generating novel molecules with GenMol")
print("=" * 70)

genmol_url = "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate"

payload = {
    "smiles": "[*{20-30}]",
    "num_molecules": 30,
    "scoring": "QED",
    "unique": True,
    "temperature": "1.0",
    "noise": "1.0",
}

t0 = time.time()
r = requests.post(genmol_url, headers=HEADERS, json=payload, timeout=180)
r.raise_for_status()
genmol_result = r.json()
t1 = time.time()

molecules = sorted(genmol_result["molecules"], key=lambda x: x["score"], reverse=True)
top_20 = molecules[:20]

print(f"Generated {len(molecules)} molecules in {t1 - t0:.1f}s")
print(f"\nTop 10 by QED score:")
print(f"{'Rank':>4}  {'QED':>6}  {'SMILES'}")
print("-" * 70)
for i, mol in enumerate(top_20[:10], 1):
    print(f"{i:4d}  {mol['score']:6.4f}  {mol['smiles'][:55]}")

# ═══════════════════════════════════════════════════════════════
# STEP 2: Dock molecules with DiffDock
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 2: Docking top 20 molecules against RelB with DiffDock")
print("=" * 70)

# Fetch RelB PDB structure (5OG6 — RelB RHD crystal structure)
print("Fetching RelB crystal structure (PDB: 5OG6)...")
pdb_r = requests.get("https://files.rcsb.org/download/5OG6.pdb", timeout=30)
pdb_r.raise_for_status()
receptor_pdb = "\n".join(
    line for line in pdb_r.text.splitlines() if line.startswith("ATOM")
)
print(f"Loaded receptor: {len(receptor_pdb.splitlines())} ATOM lines")

diffdock_url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"

docking_results = []
for i, mol in enumerate(top_20):
    payload = {
        "protein": receptor_pdb,
        "ligand": mol["smiles"],
        "ligand_file_type": "txt",
        "num_poses": 5,
        "time_divisions": 20,
        "steps": 18,
        "save_trajectory": False,
    }

    try:
        r = requests.post(diffdock_url, headers=HEADERS, json=payload, timeout=300)
        r.raise_for_status()
        result = r.json()

        best_conf = result["position_confidence"][0]
        best_pose = result["ligand_positions"][0]

        docking_results.append({
            "smiles": mol["smiles"],
            "qed_score": mol["score"],
            "docking_confidence": best_conf,
            "best_pose_sdf": best_pose,
        })
        print(f"  Mol {i + 1:2d}/{len(top_20)}: QED={mol['score']:.3f}  dock_conf={best_conf:.4f}  ✓")
    except Exception as e:
        print(f"  Mol {i + 1:2d}/{len(top_20)}: QED={mol['score']:.3f}  FAILED: {e}")

docking_results.sort(key=lambda x: x["docking_confidence"], reverse=True)

print(f"\nSuccessfully docked: {len(docking_results)}/{len(top_20)}")
print(f"\nTop 5 by docking confidence:")
print(f"{'Rank':>4}  {'Dock Conf':>10}  {'QED':>6}  {'SMILES'}")
print("-" * 70)
for i, d in enumerate(docking_results[:5], 1):
    print(f"{i:4d}  {d['docking_confidence']:10.4f}  {d['qed_score']:6.4f}  {d['smiles'][:45]}")

# ═══════════════════════════════════════════════════════════════
# STEP 3: Predict binding affinity with Boltz-2
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 3: Predicting binding affinity with Boltz-2 (top 5 docking hits)")
print("=" * 70)

boltz_url = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

# Use a trimmed RHD sequence for Boltz-2 (full-length is too long for fast inference)
relb_rhd = RELB_SEQUENCE[120:400]
print(f"Using RelB RHD region: {len(relb_rhd)} residues")

affinity_results = []
for i, d in enumerate(docking_results[:5]):
    payload = {
        "polymers": [
            {"id": "A", "molecule_type": "protein", "sequence": relb_rhd}
        ],
        "ligands": [
            {"id": "L1", "smiles": d["smiles"], "predict_affinity": True}
        ],
        "recycling_steps": 3,
        "sampling_steps": 50,
        "diffusion_samples": 1,
        "sampling_steps_affinity": 200,
        "diffusion_samples_affinity": 5,
        "output_format": "mmcif",
    }

    try:
        r = requests.post(boltz_url, headers=HEADERS, json=payload, timeout=600)
        r.raise_for_status()
        result = r.json()

        aff = result["affinities"]["L1"]
        pic50 = aff["affinity_pic50"][0]
        prob_binding = aff["affinity_probability_binary"][0]

        affinity_results.append({
            **d,
            "pic50": pic50,
            "probability_binding": prob_binding,
        })
        print(f"  Hit {i + 1}: pIC50={pic50:.2f}  P(bind)={prob_binding:.3f}  dock_conf={d['docking_confidence']:.4f}  ✓")
    except Exception as e:
        print(f"  Hit {i + 1}: FAILED — {e}")

affinity_results.sort(key=lambda x: x["pic50"], reverse=True)

# ═══════════════════════════════════════════════════════════════
# FINAL RESULTS
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL RESULTS — RelB Drug Candidates Ranked by Predicted Affinity")
print("=" * 70)
print(f"\n{'Rank':>4}  {'pIC50':>6}  {'P(bind)':>8}  {'Dock Conf':>10}  {'QED':>6}  {'SMILES'}")
print("-" * 80)
for i, ar in enumerate(affinity_results, 1):
    print(
        f"{i:4d}  {ar['pic50']:6.2f}  {ar['probability_binding']:8.3f}  "
        f"{ar['docking_confidence']:10.4f}  {ar['qed_score']:6.4f}  {ar['smiles'][:40]}"
    )

print("\n\nInterpretation:")
print("  pIC50 > 6  → sub-micromolar binding (potent)")
print("  pIC50 > 8  → low-nanomolar binding (very potent)")
print("  P(bind) > 0.7 → likely binder")
print("  QED > 0.5 → drug-like properties")

# Save results
output_path = Path("/Users/omarlujanoolazaba/Desktop/Amazon/relb_pipeline_results.csv")
with open(output_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["rank", "smiles", "pic50", "probability_binding", "docking_confidence", "qed_score"],
    )
    writer.writeheader()
    for i, ar in enumerate(affinity_results, 1):
        writer.writerow({
            "rank": i,
            "smiles": ar["smiles"],
            "pic50": ar["pic50"],
            "probability_binding": ar["probability_binding"],
            "docking_confidence": ar["docking_confidence"],
            "qed_score": ar["qed_score"],
        })

print(f"\nResults saved to: {output_path}")

# Save top SDF poses
sdf_dir = Path("/Users/omarlujanoolazaba/Desktop/Amazon/relb_docking_poses")
sdf_dir.mkdir(exist_ok=True)
for i, ar in enumerate(affinity_results, 1):
    sdf_path = sdf_dir / f"relb_hit_{i}.sdf"
    sdf_path.write_text(ar["best_pose_sdf"])
print(f"Docking poses saved to: {sdf_dir}/")
