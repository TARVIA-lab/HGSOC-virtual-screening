#!/usr/bin/env python3
"""
Lead Optimization: RelB Hit #1 (pyridine-piperazine amide)
Multiple GenMol strategies → DiffDock → Boltz-2

Parent: Cc1cc(C(=O)N2CCN(c3ccccc3)CC2)ccn1  (pIC50=6.93, QED=0.848)
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

PARENT_SMILES = "Cc1cc(C(=O)N2CCN(c3ccccc3)CC2)ccn1"
PARENT_PIC50 = 6.93
PARENT_QED = 0.848

RELB_SEQUENCE = (
    "MDSLESSFMSQDFSRHQARQGPYENEYSEFQDSPSGSEQSPHNVDALQGR"
    "EECLRAGKPSISDIIQDQFERALTEPHHEQAMRKFAQHQQPTQHRSEKRN"
    "RLAPQWQLHHVQAIFPPNWHKFAADRRFRFRQLVEAQAGKLPECLHISAA"
    "YRGTKMTEQELQERYDQLLRESGFNLQWKQPRGFHLLTKNSPDQYPACEA"
    "HMRNLYEGECVKGMKCGQEELNITIVKLQSATQNFKFTKPAAVLKRTCQE"
    "GVDAQAAYPFLSTPQVALFGSGDGHGHHDGSLDSLGTLSAPLVTAPPVES"
    "FCQHAFIRNRFKESGEFSYALSAVKQPTVQTIKPRVDKQGSPVYTSASSD"
    "VTHVFKRRGSTSTKGFPQREGMEQKLISEEDLSSLAGTPVSSLGIISSAS"
    "PTDIDSPFGHSAGLSRTTAQSESLVHQEPSPPTITMLSELVSKLHSSAQK"
    "DTAPYSGQSSAQVSISGHGTSLPVASSGQTLHQHDQQHKVHIQRPQQQGL"
    "SPFASSTFSSSSLSFPQHADYLLHTHAPAPPAPCPHDIYSDLLDTAPPSP"
)
RELB_RHD = RELB_SEQUENCE[120:400]

genmol_url = "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate"

print("=" * 70)
print("LEAD OPTIMIZATION — Hit #1: Pyridine-piperazine amide")
print(f"Parent: {PARENT_SMILES}")
print(f"Parent pIC50: {PARENT_PIC50}  QED: {PARENT_QED}")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════
# STEP 1: Multi-strategy analog generation
# ═══════════════════════════════════════════════════════════════
print("\nSTEP 1: Multi-strategy analog generation with GenMol")

# Strategy A: Use key substructure fragments as seeds
# The parent has: pyridine ring + amide + piperazine + phenyl
# We'll use partial scaffolds with masks

strategies = [
    # Keep pyridine-amide core, vary piperazine substituent
    ("Pyridine-amide + varied amine", "C(=O)N1CCN(CC1).[*{10-20}]"),
    # Keep piperazine-phenyl, vary acyl group
    ("Piperazine-phenyl + varied acyl", "N1CCN(c2ccccc2)CC1.[*{10-18}]"),
    # De novo around same size (drug-like 20-30 atoms)
    ("De novo drug-like (batch 1)", "[*{20-30}]"),
    # De novo slightly larger
    ("De novo drug-like (batch 2)", "[*{22-32}]"),
    # De novo with tighter QED
    ("De novo compact (batch 3)", "[*{18-25}]"),
    # Heterocyclic-biased size range
    ("De novo (batch 4)", "[*{20-28}]"),
]

all_molecules = []
for name, mask in strategies:
    print(f"\n  Strategy: {name}")
    print(f"  SAFE input: {mask}")

    for temp, noise in [("0.8", "0.8"), ("1.2", "1.0")]:
        payload = {
            "smiles": mask,
            "num_molecules": 100,
            "scoring": "QED",
            "unique": True,
            "temperature": temp,
            "noise": noise,
        }

        try:
            r = requests.post(genmol_url, headers=HEADERS, json=payload, timeout=180)
            r.raise_for_status()
            result = r.json()
            mols = result.get("molecules", [])
            for m in mols:
                m["strategy"] = name
            all_molecules.extend(mols)
            print(f"    temp={temp} noise={noise}: {len(mols)} molecules")
        except Exception as e:
            print(f"    temp={temp} noise={noise}: FAILED — {e}")

# Deduplicate, filter parent, and rank
seen = {PARENT_SMILES}
unique_mols = []
for m in sorted(all_molecules, key=lambda x: x["score"], reverse=True):
    if m["smiles"] not in seen:
        seen.add(m["smiles"])
        unique_mols.append(m)

print(f"\nTotal unique novel analogs: {len(unique_mols)}")
top_30 = unique_mols[:30]

print(f"\nTop 15 candidates by QED:")
print(f"{'Rank':>4}  {'QED':>6}  {'Strategy':<35}  {'SMILES'}")
print("-" * 100)
for i, m in enumerate(top_30[:15], 1):
    print(f"{i:4d}  {m['score']:6.4f}  {m['strategy']:<35}  {m['smiles'][:50]}")

# ═══════════════════════════════════════════════════════════════
# STEP 2: Dock with DiffDock
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 2: Docking top 30 candidates against RelB (PDB: 5OG6)")
print("=" * 70)

print("Fetching RelB crystal structure...")
pdb_r = requests.get("https://files.rcsb.org/download/5OG6.pdb", timeout=30)
pdb_r.raise_for_status()
receptor_pdb = "\n".join(
    line for line in pdb_r.text.splitlines() if line.startswith("ATOM")
)
print(f"Loaded receptor: {len(receptor_pdb.splitlines())} ATOM lines")

diffdock_url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
docking_results = []

for i, mol in enumerate(top_30):
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
            "strategy": mol["strategy"],
            "docking_confidence": best_conf,
            "best_pose_sdf": best_pose,
        })
        print(f"  {i+1:2d}/{len(top_30)}: QED={mol['score']:.3f}  dock_conf={best_conf:.4f}  ✓")
    except Exception as e:
        print(f"  {i+1:2d}/{len(top_30)}: FAILED — {str(e)[:60]}")

docking_results.sort(key=lambda x: x["docking_confidence"], reverse=True)
print(f"\nDocked: {len(docking_results)}/{len(top_30)}")

# ═══════════════════════════════════════════════════════════════
# STEP 3: Boltz-2 affinity (top 10 docking hits)
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 3: Boltz-2 affinity prediction for top 10 docking hits")
print("=" * 70)

boltz_url = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
affinity_results = []

for i, d in enumerate(docking_results[:10]):
    payload = {
        "polymers": [
            {"id": "A", "molecule_type": "protein", "sequence": RELB_RHD}
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

        improvement = pic50 - PARENT_PIC50
        tag = "★ IMPROVED" if improvement > 0 else "~ similar" if improvement > -0.3 else "  worse"

        affinity_results.append({
            **d,
            "pic50": pic50,
            "probability_binding": prob_binding,
            "improvement": improvement,
        })
        print(f"  {i+1:2d}: pIC50={pic50:.2f} ({improvement:+.2f})  P(bind)={prob_binding:.3f}  {tag}")
    except Exception as e:
        print(f"  {i+1:2d}: FAILED — {str(e)[:60]}")

affinity_results.sort(key=lambda x: x["pic50"], reverse=True)

# ═══════════════════════════════════════════════════════════════
# FINAL RESULTS
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL OPTIMIZATION RESULTS")
print(f"Parent: {PARENT_SMILES}  pIC50={PARENT_PIC50}  QED={PARENT_QED}")
print("=" * 70)

improved = [a for a in affinity_results if a["improvement"] > 0]
print(f"\nImproved over parent: {len(improved)}/{len(affinity_results)}")

print(f"\n{'#':>2}  {'pIC50':>6}  {'Δ':>6}  {'P(bind)':>8}  {'Dock':>7}  {'QED':>6}  {'SMILES'}")
print("-" * 95)
for i, ar in enumerate(affinity_results, 1):
    delta = ar["improvement"]
    marker = "★" if delta > 0 else " "
    print(
        f"{i:2d}  {ar['pic50']:6.2f}  {delta:+5.2f}  {ar['probability_binding']:8.3f}  "
        f"{ar['docking_confidence']:7.3f}  {ar['qed_score']:6.3f}  {marker}{ar['smiles'][:38]}"
    )

# Save
output_path = Path("/Users/omarlujanoolazaba/Desktop/Amazon/relb_hit1_optimization_results.csv")
with open(output_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["rank", "smiles", "pic50", "delta_pic50", "probability_binding",
                     "docking_confidence", "qed_score", "strategy"],
    )
    writer.writeheader()
    for i, ar in enumerate(affinity_results, 1):
        writer.writerow({
            "rank": i,
            "smiles": ar["smiles"],
            "pic50": round(ar["pic50"], 3),
            "delta_pic50": round(ar["improvement"], 3),
            "probability_binding": round(ar["probability_binding"], 4),
            "docking_confidence": round(ar["docking_confidence"], 4),
            "qed_score": round(ar["qed_score"], 4),
            "strategy": ar["strategy"],
        })
print(f"\nResults saved to: {output_path}")

sdf_dir = Path("/Users/omarlujanoolazaba/Desktop/Amazon/relb_hit1_optimized_poses")
sdf_dir.mkdir(exist_ok=True)
for i, ar in enumerate(affinity_results, 1):
    sdf_path = sdf_dir / f"hit1_analog_{i}.sdf"
    sdf_path.write_text(ar["best_pose_sdf"])
print(f"Docking poses saved to: {sdf_dir}/")
