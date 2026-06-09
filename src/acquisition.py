"""ChEMBL data acquisition for the HGSOC virtual screening pipeline."""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def get_chembl_id_for_gene(gene_symbol: str) -> Optional[str]:
    """Return the primary ChEMBL target ID for a human gene symbol."""
    from chembl_webresource_client.new_client import new_client

    results = list(
        new_client.target.filter(
            target_synonym__icontains=gene_symbol,
            target_type="SINGLE PROTEIN",
            organism="Homo sapiens",
        ).only(["target_chembl_id", "pref_name"])
    )
    if not results:
        logger.warning("No ChEMBL target found for %s", gene_symbol)
        return None
    logger.info("Mapped %s → %s", gene_symbol, results[0]["target_chembl_id"])
    return results[0]["target_chembl_id"]


def fetch_bioactivity_data(
    chembl_id: str,
    activity_types: list = None,
    assay_type: str = "B",
    max_compounds: int = 5000,
) -> pd.DataFrame:
    """
    Pull compound–activity pairs from ChEMBL for a given target.

    Returns columns: molecule_chembl_id, canonical_smiles, standard_type,
    standard_value, standard_units, pchembl_value.
    """
    from chembl_webresource_client.new_client import new_client

    if activity_types is None:
        activity_types = ["IC50", "Ki", "Kd"]

    records = list(
        new_client.activity.filter(
            target_chembl_id=chembl_id,
            standard_type__in=activity_types,
            standard_relation__in=["=", "<"],
            assay_type=assay_type,
            pchembl_value__isnull=False,
        ).only(
            [
                "molecule_chembl_id",
                "canonical_smiles",
                "standard_type",
                "standard_value",
                "standard_units",
                "pchembl_value",
            ]
        )[:max_compounds]
    )

    if not records:
        logger.warning("No activity data for %s", chembl_id)
        return pd.DataFrame()

    df = pd.DataFrame.from_records(records)
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    df = df.dropna(subset=["canonical_smiles", "pchembl_value"])
    df = df.drop_duplicates(subset=["molecule_chembl_id"])
    logger.info("Fetched %d compounds for %s", len(df), chembl_id)
    return df


def label_activity(
    df: pd.DataFrame,
    active_threshold: float = 6.0,
    inactive_threshold: float = 5.0,
) -> pd.DataFrame:
    """
    Assign binary activity labels.

    pChEMBL >= active_threshold  → 1 (active)
    pChEMBL <  inactive_threshold → 0 (inactive)
    Compounds between thresholds are dropped (ambiguous gray zone).
    """
    df = df.copy()
    df["active"] = np.nan
    df.loc[df["pchembl_value"] >= active_threshold, "active"] = 1
    df.loc[df["pchembl_value"] < inactive_threshold, "active"] = 0
    df = df.dropna(subset=["active"])
    df["active"] = df["active"].astype(int)
    logger.info(
        "Labeled: %d active, %d inactive", df["active"].sum(), (df["active"] == 0).sum()
    )
    return df


def build_training_dataset(
    targets: list,
    config: dict,
    output_dir: str = "data/processed",
) -> pd.DataFrame:
    """
    Acquire and label compounds for all targets, save per-target and combined CSVs.

    targets: list of dicts with keys 'gene_symbol' and optionally 'chembl_id'.
    config:  loaded config.yaml dict.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    all_dfs = []

    for target in targets:
        gene = target.get("gene_symbol", "unknown")
        chembl_id = target.get("chembl_id") or get_chembl_id_for_gene(gene)
        if not chembl_id:
            continue

        df = fetch_bioactivity_data(
            chembl_id,
            activity_types=config["acquisition"]["activity_types"],
            assay_type=config["acquisition"]["assay_type"],
            max_compounds=config["acquisition"]["max_compounds_per_target"],
        )
        if df.empty:
            continue

        df = label_activity(
            df,
            active_threshold=config["acquisition"]["pchembl_cutoff_active"],
            inactive_threshold=config["acquisition"]["pchembl_cutoff_inactive"],
        )
        df["target_gene"] = gene
        df["target_chembl_id"] = chembl_id

        df.to_csv(Path(output_dir) / f"{gene}_bioactivity.csv", index=False)
        all_dfs.append(df)

    if not all_dfs:
        return pd.DataFrame()

    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(Path(output_dir) / "combined_bioactivity.csv", index=False)
    logger.info("Total: %d compounds across %d targets", len(combined), len(all_dfs))
    return combined
