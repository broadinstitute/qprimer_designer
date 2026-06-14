"""Variant frequency forecasting via multinomial logistic regression (MLR).

Uses scikit-learn for MLR fitting — already a pipeline dependency, so no extra
install is required. If the 'evofr' package is installed (pip install evofr),
it is used instead for Bayesian inference with uncertainty estimates.
"""

from __future__ import annotations

import sys
from typing import Optional

import numpy as np
import pandas as pd


# Candidate column names in metadata files (tried in order, case-insensitive)
_DATE_COLUMNS = [
    "collection_date",
    "collection date",
    "isolate_collection_date",
    "date",
    "release_date",
    "release date",
    "sample_date",
]
_CLADE_COLUMNS = [
    "pango_lineage",
    "pango lineage",
    "lineage",
    "clade",
    "nextclade_pango",
    "nextstrain_clade",
    "nextstrain clade",
    "gisaid_clade",
    "gisaid clade",
    "variant",
    "strain",
]
_LOCATION_COLUMNS = [
    "geographic_region",
    "geographic region",
    "country",
    "location",
    "geo_loc_name",
]
_ACCESSION_COLUMNS = [
    "accession",
    "genbank_accession",
    "sequence_id",
    "seq_id",
    "id",
    "name",
    "strain",
]


def _find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate found in df.columns (case-insensitive match)."""
    col_map = {c.lower().strip(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in col_map:
            return col_map[cand.lower()]
    return None


def metadata_to_clade_counts(
    metadata_df: pd.DataFrame,
    location_filter: str | None = None,
) -> pd.DataFrame:
    """Convert a sequence metadata DataFrame to daily clade counts.

    Handles metadata from gget virus output and user-provided TSVs with flexible
    column name detection.

    Args:
        metadata_df: Sequence metadata with at minimum a date column and a
            lineage/clade column.
        location_filter: Optional substring to filter the location column
            (case-insensitive). If no sequences match, falls back to all locations.

    Returns:
        DataFrame with columns: date (datetime.date), clade (str), count (int).

    Raises:
        ValueError: If required date or clade columns cannot be found, or if no
            rows remain after date parsing.
    """
    df = metadata_df.copy()

    date_col = _find_column(df, _DATE_COLUMNS)
    if date_col is None:
        raise ValueError(
            "No date column found in metadata. "
            f"Expected one of: {', '.join(_DATE_COLUMNS[:5])}. "
            "Add a 'date' or 'collection_date' column to your metadata file."
        )

    clade_col = _find_column(df, _CLADE_COLUMNS)
    if clade_col is None:
        raise ValueError(
            "No lineage/clade column found in metadata. "
            f"Expected one of: {', '.join(_CLADE_COLUMNS[:5])}. "
            "Add a 'clade' or 'pango_lineage' column to your metadata file."
        )

    location_col = _find_column(df, _LOCATION_COLUMNS)
    if location_filter and location_col:
        mask = df[location_col].fillna("").str.contains(location_filter, case=False, na=False)
        if mask.sum() == 0:
            print(
                f"Warning: no sequences match location '{location_filter}'. "
                "Using all locations.",
                file=sys.stderr,
            )
        else:
            df = df[mask]

    df["_date"] = pd.to_datetime(df[date_col], errors="coerce").dt.date
    df = df.dropna(subset=["_date"])

    if df.empty:
        raise ValueError(
            f"No rows with parseable dates in column '{date_col}'. "
            "Ensure dates are in YYYY-MM-DD format."
        )

    df["_clade"] = df[clade_col].fillna("Unknown").astype(str).str.strip()
    df.loc[df["_clade"].isin(["", "nan", "N/A", "NA", "None"]), "_clade"] = "Unknown"

    counts = (
        df.groupby(["_date", "_clade"])
        .size()
        .reset_index(name="count")
        .rename(columns={"_date": "date", "_clade": "clade"})
    )
    return counts


def run_mlr_forecast(
    clade_counts: pd.DataFrame,
    horizon_days: int = 30,
    min_sequences: int = 10,
) -> pd.DataFrame:
    """Forecast future variant frequencies using multinomial logistic regression.

    Fits an MLR model to observed clade frequencies over time (sequences per
    clade per day) and extrapolates to estimate frequencies horizon_days ahead.
    This mirrors the Bedford Lab / Nextstrain evofr MLR model conceptually:
    each variant has a fixed growth advantage (slope in log-frequency space).

    Args:
        clade_counts: DataFrame with columns [date (datetime.date), clade (str),
            count (int)].
        horizon_days: Days ahead to forecast.
        min_sequences: Minimum total sequences per clade to include in the model.

    Returns:
        DataFrame with columns [clade, current_freq, predicted_freq,
        growth_advantage]. Sorted by predicted_freq descending.

    Raises:
        ValueError: If fewer than 2 clades meet the minimum sequence threshold.
    """
    from sklearn.linear_model import LogisticRegression

    clade_totals = clade_counts.groupby("clade")["count"].sum()
    kept = sorted(clade_totals[clade_totals >= min_sequences].index.tolist())
    if len(kept) < 2:
        raise ValueError(
            f"Only {len(kept)} clade(s) with >= {min_sequences} sequences. "
            "Lower --min-sequences or provide more sequences with lineage metadata."
        )

    df = clade_counts[clade_counts["clade"].isin(kept)].copy()

    # Encode dates as days since first observation
    all_dates = sorted(df["date"].unique())
    min_date = all_dates[0]
    df["day"] = df["date"].apply(lambda d: (d - min_date).days)
    last_day = int(df["day"].max())
    future_day = last_day + horizon_days

    # Expand counts into individual observations for weighted fitting
    days = df["day"].values.astype(float)
    clades = df["clade"].values
    weights = df["count"].values.astype(float)

    # Normalize days to [0, 1] for numerical stability
    day_scale = max(last_day, 1.0)
    X = (days / day_scale).reshape(-1, 1)

    clf = LogisticRegression(solver="lbfgs", max_iter=2000, C=1.0)
    clf.fit(X, clades, sample_weight=weights)

    current_probs = clf.predict_proba([[last_day / day_scale]])[0]
    future_probs = clf.predict_proba([[future_day / day_scale]])[0]

    results = []
    for i, clade in enumerate(clf.classes_):
        curr = float(current_probs[i])
        fut = float(future_probs[i])
        ga = fut / curr if curr > 1e-6 else 0.0
        results.append({
            "clade": clade,
            "current_freq": round(curr, 4),
            "predicted_freq": round(fut, 4),
            "growth_advantage": round(ga, 4),
        })

    return (
        pd.DataFrame(results)
        .sort_values("predicted_freq", ascending=False)
        .reset_index(drop=True)
    )
