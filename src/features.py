"""Molecular feature engineering: fingerprints, physicochemical descriptors, filters."""

import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def lipinski_filter(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """Drop compounds that fail Lipinski's Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    def passes(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return False
        return (
            Descriptors.ExactMolWt(mol) <= 500
            and Descriptors.MolLogP(mol) <= 5
            and Descriptors.NumHDonors(mol) <= 5
            and Descriptors.NumHAcceptors(mol) <= 10
        )

    mask = df[smiles_col].apply(passes)
    logger.info("Lipinski filter: kept %d / %d", mask.sum(), len(df))
    return df[mask].copy()


def smiles_to_morgan(
    smiles: str,
    radius: int = 2,
    n_bits: int = 2048,
    use_chirality: bool = True,
) -> Optional[np.ndarray]:
    """Convert a SMILES string to an ECFP Morgan fingerprint as a float32 array."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=radius, nBits=n_bits, useChirality=use_chirality
    )
    return np.array(fp, dtype=np.float32)


def compute_fingerprints(
    df: pd.DataFrame,
    smiles_col: str = "canonical_smiles",
    radius: int = 2,
    n_bits: int = 2048,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Compute Morgan fingerprints for all compounds.

    Returns (X [n_valid, n_bits], valid_df) — rows with invalid SMILES are dropped.
    """
    fps, valid_idx = [], []
    for idx, row in df.iterrows():
        fp = smiles_to_morgan(row[smiles_col], radius=radius, n_bits=n_bits)
        if fp is not None:
            fps.append(fp)
            valid_idx.append(idx)

    X = np.vstack(fps)
    valid_df = df.loc[valid_idx].reset_index(drop=True)
    logger.info("Fingerprints: %d valid / %d total", len(valid_df), len(df))
    return X, valid_df


def compute_rdkit_descriptors(
    df: pd.DataFrame, smiles_col: str = "canonical_smiles"
) -> pd.DataFrame:
    """Return a DataFrame of 7 standard 2-D physicochemical descriptors per compound."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    fns = {
        "MW": Descriptors.ExactMolWt,
        "LogP": Descriptors.MolLogP,
        "HBD": Descriptors.NumHDonors,
        "HBA": Descriptors.NumHAcceptors,
        "TPSA": Descriptors.TPSA,
        "RotBonds": Descriptors.NumRotatableBonds,
        "AromaticRings": Descriptors.NumAromaticRings,
    }
    rows = []
    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        rows.append({k: fn(mol) for k, fn in fns.items()} if mol else {k: np.nan for k in fns})
    return pd.DataFrame(rows)
