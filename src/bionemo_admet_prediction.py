#!/usr/bin/env python3
"""
ADMET Property Prediction — Top 2 Improved RelB Analogs + Parent
Uses RDKit descriptors for comprehensive drugability profiling.
"""

import csv
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

compounds = [
    {
        "name": "Analog #1 — Dimethylphenyl oxadiazole-THF",
        "smiles": "Cc1ccc(C(=O)Nc2noc([C@@H]3CCCO3)n2)c(C)c1",
        "pic50": 7.08,
        "delta": +0.15,
    },
    {
        "name": "Analog #2 — Methoxyquinoline aminopiperidine",
        "smiles": "COc1ccc2cc(NC[C@H]3CCCCO3)cnc2c1",
        "pic50": 6.99,
        "delta": +0.06,
    },
    {
        "name": "Parent — Pyridine-piperazine amide",
        "smiles": "Cc1cc(C(=O)N2CCN(c3ccccc3)CC2)ccn1",
        "pic50": 6.93,
        "delta": 0.0,
    },
]

print("=" * 80)
print("ADMET PROPERTY PREDICTION — RelB Drug Candidates")
print("=" * 80)

# PAINS filter
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog.FilterCatalog(params)

# Brenk filter (unwanted substructures)
params_brenk = FilterCatalogParams()
params_brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
brenk_catalog = FilterCatalog.FilterCatalog(params_brenk)

all_results = []

for compound in compounds:
    mol = Chem.MolFromSmiles(compound["smiles"])
    if mol is None:
        print(f"\n  ERROR: Could not parse {compound['smiles']}")
        continue

    # ── Physicochemical properties ──
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hba = Descriptors.NumHAcceptors(mol)
    hbd = Descriptors.NumHDonors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Descriptors.NumRotatableBonds(mol)
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    rings = Descriptors.RingCount(mol)
    qed_score = QED.qed(mol)
    mr = Descriptors.MolMR(mol)  # molar refractivity
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)

    # ── Lipinski Rule of 5 ──
    lipinski_violations = 0
    if mw > 500: lipinski_violations += 1
    if logp > 5: lipinski_violations += 1
    if hba > 10: lipinski_violations += 1
    if hbd > 5: lipinski_violations += 1
    lipinski_pass = lipinski_violations <= 1

    # ── Veber rules (oral bioavailability) ──
    veber_pass = tpsa <= 140 and rotatable <= 10

    # ── Ghose filter ──
    ghose_pass = (160 <= mw <= 480 and -0.4 <= logp <= 5.6
                  and 40 <= mr <= 130 and 20 <= heavy_atoms <= 70)

    # ── Egan filter (oral absorption) ──
    egan_pass = tpsa <= 131.6 and logp <= 5.88

    # ── Muegge filter ──
    muegge_pass = (200 <= mw <= 600 and -2 <= logp <= 5
                   and tpsa <= 150 and rings <= 7
                   and hba <= 10 and hbd <= 5 and rotatable <= 15)

    # ── PAINS alerts ──
    pains_matches = pains_catalog.GetMatches(mol)
    pains_alerts = [m.GetDescription() for m in pains_matches]
    pains_pass = len(pains_alerts) == 0

    # ── Brenk alerts ──
    brenk_matches = brenk_catalog.GetMatches(mol)
    brenk_alerts = [m.GetDescription() for m in brenk_matches]

    # ── Absorption estimates ──
    # GI absorption (Egan egg model approximation)
    gi_absorption = "High" if (tpsa <= 131.6 and logp <= 5.88) else "Low"

    # BBB permeation (simplified Egan model)
    bbb_permeant = "Yes" if (tpsa <= 79 and logp >= 1) else "No"

    # Pgp substrate likelihood (rough heuristic)
    pgp_substrate = "Likely" if (mw > 400 or hbd > 3) else "Unlikely"

    # ── Solubility estimate (ESOL model - Delaney) ──
    log_s = (0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotatable
             - 0.74 * aromatic_rings)
    solubility_class = ("Highly soluble" if log_s >= 0
                        else "Soluble" if log_s >= -2
                        else "Moderately soluble" if log_s >= -4
                        else "Poorly soluble" if log_s >= -6
                        else "Insoluble")

    # ── Synthetic accessibility (SA Score) ──
    from rdkit.Chem import RDConfig
    import os, sys
    sa_path = os.path.join(RDConfig.RDContribDir, 'SA_Score')
    if sa_path not in sys.path:
        sys.path.insert(0, sa_path)
    try:
        import sascorer
        sa_score = sascorer.calculateScore(mol)
    except Exception:
        sa_score = None

    # ── Lead-likeness (Teague) ──
    leadlike = mw <= 350 and logp <= 3.5 and rotatable <= 7

    result = {
        "name": compound["name"],
        "smiles": compound["smiles"],
        "pic50": compound["pic50"],
        "delta_pic50": compound["delta"],
        "mw": round(mw, 1),
        "logp": round(logp, 2),
        "hba": hba,
        "hbd": hbd,
        "tpsa": round(tpsa, 1),
        "rotatable_bonds": rotatable,
        "aromatic_rings": aromatic_rings,
        "heavy_atoms": heavy_atoms,
        "fsp3": round(fsp3, 3),
        "qed": round(qed_score, 3),
        "lipinski_violations": lipinski_violations,
        "lipinski_pass": lipinski_pass,
        "veber_pass": veber_pass,
        "ghose_pass": ghose_pass,
        "egan_pass": egan_pass,
        "muegge_pass": muegge_pass,
        "gi_absorption": gi_absorption,
        "bbb_permeant": bbb_permeant,
        "pgp_substrate": pgp_substrate,
        "log_s_esol": round(log_s, 2),
        "solubility_class": solubility_class,
        "pains_pass": pains_pass,
        "pains_alerts": pains_alerts,
        "brenk_alerts": brenk_alerts,
        "sa_score": round(sa_score, 2) if sa_score else None,
        "leadlike": leadlike,
    }
    all_results.append(result)

    # Print report
    print(f"\n{'─' * 80}")
    print(f"  {compound['name']}")
    print(f"  {compound['smiles']}")
    print(f"  pIC50: {compound['pic50']}  (Δ {compound['delta']:+.2f})")
    print(f"{'─' * 80}")

    print(f"\n  PHYSICOCHEMICAL PROPERTIES")
    print(f"    MW:              {mw:.1f} Da    {'✓' if mw <= 500 else '✗ >500'}")
    print(f"    LogP:            {logp:.2f}       {'✓' if logp <= 5 else '✗ >5'}")
    print(f"    HB Acceptors:    {hba}           {'✓' if hba <= 10 else '✗ >10'}")
    print(f"    HB Donors:       {hbd}           {'✓' if hbd <= 5 else '✗ >5'}")
    print(f"    TPSA:            {tpsa:.1f} Å²   {'✓' if tpsa <= 140 else '✗ >140'}")
    print(f"    Rotatable bonds: {rotatable}           {'✓' if rotatable <= 10 else '✗ >10'}")
    print(f"    Aromatic rings:  {aromatic_rings}")
    print(f"    Fsp3:            {fsp3:.3f}       {'✓ 3D-rich' if fsp3 >= 0.25 else '✗ flat'}")
    print(f"    QED:             {qed_score:.3f}")

    print(f"\n  DRUG-LIKENESS FILTERS")
    print(f"    Lipinski Ro5:    {'PASS' if lipinski_pass else 'FAIL'} ({lipinski_violations} violations)")
    print(f"    Veber:           {'PASS' if veber_pass else 'FAIL'}")
    print(f"    Ghose:           {'PASS' if ghose_pass else 'FAIL'}")
    print(f"    Egan:            {'PASS' if egan_pass else 'FAIL'}")
    print(f"    Muegge:          {'PASS' if muegge_pass else 'FAIL'}")
    print(f"    Lead-like:       {'YES' if leadlike else 'NO'}")

    print(f"\n  ABSORPTION & DISTRIBUTION")
    print(f"    GI absorption:   {gi_absorption}")
    print(f"    BBB permeant:    {bbb_permeant}")
    print(f"    Pgp substrate:   {pgp_substrate}")

    print(f"\n  SOLUBILITY (ESOL)")
    print(f"    Log S:           {log_s:.2f}")
    print(f"    Class:           {solubility_class}")

    print(f"\n  MEDICINAL CHEMISTRY")
    print(f"    PAINS alerts:    {'CLEAN ✓' if pains_pass else f'FLAGGED: {pains_alerts}'}")
    print(f"    Brenk alerts:    {'CLEAN ✓' if not brenk_alerts else f'{brenk_alerts}'}")
    if sa_score:
        sa_label = "Easy" if sa_score <= 3 else "Moderate" if sa_score <= 6 else "Difficult"
        print(f"    SA Score:        {sa_score:.2f} ({sa_label})")

# ═══════════════════════════════════════════════════════════════
# COMPARATIVE SUMMARY
# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 80}")
print("COMPARATIVE SUMMARY")
print(f"{'=' * 80}")
print(f"\n{'Property':<22} {'Analog #1':>16} {'Analog #2':>16} {'Parent':>16}")
print(f"{'─' * 72}")

props = [
    ("pIC50", "pic50", ""),
    ("Δ pIC50", "delta_pic50", ""),
    ("MW (Da)", "mw", ""),
    ("LogP", "logp", ""),
    ("TPSA (Å²)", "tpsa", ""),
    ("HB Acceptors", "hba", ""),
    ("HB Donors", "hbd", ""),
    ("Rot. Bonds", "rotatable_bonds", ""),
    ("Fsp3", "fsp3", ""),
    ("QED", "qed", ""),
    ("Lipinski", "lipinski_pass", "bool"),
    ("Veber", "veber_pass", "bool"),
    ("GI Absorption", "gi_absorption", ""),
    ("BBB Permeant", "bbb_permeant", ""),
    ("Solubility", "solubility_class", ""),
    ("PAINS Clean", "pains_pass", "bool"),
    ("SA Score", "sa_score", ""),
]

for label, key, fmt in props:
    vals = []
    for r in all_results:
        v = r[key]
        if fmt == "bool":
            vals.append("✓ PASS" if v else "✗ FAIL")
        elif isinstance(v, float):
            vals.append(f"{v:+.2f}" if key == "delta_pic50" else f"{v:.2f}" if v != int(v) else f"{v:.1f}")
        else:
            vals.append(str(v))
    print(f"  {label:<20} {vals[0]:>16} {vals[1]:>16} {vals[2]:>16}")

# Save results
output_path = Path("/Users/omarlujanoolazaba/Desktop/Amazon/relb_admet_results.json")
with open(output_path, "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\nFull results saved to: {output_path}")

output_csv = Path("/Users/omarlujanoolazaba/Desktop/Amazon/relb_admet_results.csv")
with open(output_csv, "w", newline="") as f:
    flat_keys = [k for k in all_results[0] if k not in ("pains_alerts", "brenk_alerts")]
    writer = csv.DictWriter(f, fieldnames=flat_keys + ["pains_alerts", "brenk_alerts"])
    writer.writeheader()
    for r in all_results:
        row = {k: r[k] for k in flat_keys}
        row["pains_alerts"] = "; ".join(r["pains_alerts"]) if r["pains_alerts"] else ""
        row["brenk_alerts"] = "; ".join(r["brenk_alerts"]) if r["brenk_alerts"] else ""
        writer.writerow(row)
print(f"CSV saved to: {output_csv}")
