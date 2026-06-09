"""ML classifier training and evaluation for virtual screening."""

import logging
from pathlib import Path
from typing import Union

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    average_precision_score,
    classification_report,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold, cross_val_score

logger = logging.getLogger(__name__)

Model = Union[RandomForestClassifier, GradientBoostingClassifier]


def build_classifier(model_type: str = "random_forest", **kwargs) -> Model:
    defaults = {
        "random_forest": dict(n_estimators=200, n_jobs=-1, class_weight="balanced", random_state=42),
        "gradient_boosting": dict(n_estimators=200, learning_rate=0.05, max_depth=4, random_state=42),
    }
    if model_type not in defaults:
        raise ValueError(f"Unknown model type: {model_type!r}. Choose 'random_forest' or 'gradient_boosting'.")
    params = {**defaults[model_type], **kwargs}
    cls = RandomForestClassifier if model_type == "random_forest" else GradientBoostingClassifier
    return cls(**params)


def cross_validate(
    model: Model, X: np.ndarray, y: np.ndarray, cv_folds: int = 5
) -> dict:
    """Stratified k-fold CV. Returns mean ± std for ROC-AUC and PR-AUC."""
    skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
    roc = cross_val_score(model, X, y, cv=skf, scoring="roc_auc", n_jobs=-1)
    pr = cross_val_score(model, X, y, cv=skf, scoring="average_precision", n_jobs=-1)
    results = {
        "roc_auc_mean": roc.mean(), "roc_auc_std": roc.std(),
        "pr_auc_mean": pr.mean(),   "pr_auc_std": pr.std(),
    }
    logger.info(
        "CV ROC-AUC %.3f ± %.3f | PR-AUC %.3f ± %.3f",
        results["roc_auc_mean"], results["roc_auc_std"],
        results["pr_auc_mean"],  results["pr_auc_std"],
    )
    return results


def train_and_evaluate(
    model: Model,
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray,  y_test: np.ndarray,
) -> dict:
    """Fit model, evaluate on held-out test set, return metrics dict."""
    model.fit(X_train, y_train)
    y_prob = model.predict_proba(X_test)[:, 1]
    y_pred = model.predict(X_test)
    metrics = {
        "roc_auc": roc_auc_score(y_test, y_prob),
        "pr_auc": average_precision_score(y_test, y_prob),
        "classification_report": classification_report(y_test, y_pred),
        "roc_curve": roc_curve(y_test, y_prob),
        "pr_curve": precision_recall_curve(y_test, y_prob),
        "y_prob": y_prob,
    }
    logger.info("Test ROC-AUC %.3f | PR-AUC %.3f", metrics["roc_auc"], metrics["pr_auc"])
    return metrics


def get_feature_importance(model: RandomForestClassifier, top_n: int = 20) -> pd.Series:
    """Return top-N feature importances. Only valid for tree-based models."""
    return pd.Series(model.feature_importances_).nlargest(top_n)


def save_model(model: Model, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, path)
    logger.info("Saved model → %s", path)


def load_model(path: str) -> Model:
    return joblib.load(path)
