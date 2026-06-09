"""End-to-end virtual screening pipeline."""

import logging
from pathlib import Path

import pandas as pd

from .features import compute_fingerprints, lipinski_filter
from .models import load_model

logger = logging.getLogger(__name__)


def fetch_fda_approved_library(output_path: str = "data/processed/fda_library.csv") -> pd.DataFrame:
    """
    Download all FDA-approved (max_phase=4) small molecules from ChEMBL.
    Results are cached at output_path after the first run.
    """
    from chembl_webresource_client.new_client import new_client

    approved = new_client.molecule.filter(
        max_phase=4, molecule_type="Small molecule"
    ).only(["molecule_chembl_id", "molecule_structures", "pref_name"])

    records = []
    for mol in approved:
        structs = mol.get("molecule_structures") or {}
        smi = structs.get("canonical_smiles")
        if smi:
            records.append({
                "molecule_chembl_id": mol["molecule_chembl_id"],
                "pref_name": mol.get("pref_name", ""),
                "canonical_smiles": smi,
            })

    df = pd.DataFrame(records)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.info("Downloaded %d FDA-approved compounds", len(df))
    return df


def screen_library(
    library_df: pd.DataFrame,
    model_path: str,
    smiles_col: str = "canonical_smiles",
    apply_lipinski: bool = True,
    top_n: int = 50,
    radius: int = 2,
    n_bits: int = 2048,
) -> pd.DataFrame:
    """
    Score a compound library with a trained classifier.

    Returns the top_n hits ranked by p(active).
    """
    model = load_model(model_path)

    if apply_lipinski:
        library_df = lipinski_filter(library_df, smiles_col=smiles_col)

    X, valid_df = compute_fingerprints(library_df, smiles_col=smiles_col, radius=radius, n_bits=n_bits)
    probs = model.predict_proba(X)[:, 1]

    valid_df = valid_df.copy()
    valid_df["p_active"] = probs

    hits = valid_df.nlargest(top_n, "p_active").reset_index(drop=True)
    hits.insert(0, "rank", hits.index + 1)

    logger.info(
        "Screening done. Top hit: %s (p=%.3f)",
        hits["molecule_chembl_id"].iloc[0],
        hits["p_active"].iloc[0],
    )
    return hits
