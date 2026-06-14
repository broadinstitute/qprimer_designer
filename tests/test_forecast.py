"""Tests for the variant frequency forecasting module."""

import random
import tempfile
from datetime import date, timedelta
from pathlib import Path

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Helpers for building synthetic test data
# ---------------------------------------------------------------------------

def _make_clade_counts(
    clades: dict[str, tuple[int, float]],
    n_days: int = 60,
    start: date | None = None,
) -> pd.DataFrame:
    """Build synthetic clade counts with specified growth rates.

    Args:
        clades: {clade_name: (initial_count, daily_growth_factor)}
        n_days: Number of days of data to generate.
        start: Start date (defaults to 60 days before today).

    Returns:
        DataFrame with [date, clade, count] columns.
    """
    if start is None:
        start = date.today() - timedelta(days=n_days)

    rows = []
    for d in range(n_days):
        current_date = start + timedelta(days=d)
        for clade, (init, factor) in clades.items():
            count = max(1, int(init * (factor ** d)))
            rows.append({"date": current_date, "clade": clade, "count": count})
    return pd.DataFrame(rows)


def _make_metadata_df(
    clades: dict[str, tuple[int, float]],
    n_days: int = 60,
) -> pd.DataFrame:
    """Build a metadata DataFrame mimicking gget virus output."""
    start = date.today() - timedelta(days=n_days)
    rows = []
    seq_id = 0
    for d in range(n_days):
        current_date = start + timedelta(days=d)
        for clade, (init, factor) in clades.items():
            count = max(1, int(init * (factor ** d)))
            for _ in range(count):
                rows.append({
                    "accession": f"SEQ{seq_id:05d}",
                    "collection_date": current_date.isoformat(),
                    "pango_lineage": clade,
                    "country": "USA",
                })
                seq_id += 1
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Tests for metadata_to_clade_counts
# ---------------------------------------------------------------------------

class TestMetadataToCladeCountss:
    def test_basic_conversion(self):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        meta = _make_metadata_df({"BA.2": (10, 1.0), "XBB.1.5": (5, 1.02)}, n_days=30)
        counts = metadata_to_clade_counts(meta)

        assert set(counts.columns) >= {"date", "clade", "count"}
        assert "BA.2" in counts["clade"].values
        assert "XBB.1.5" in counts["clade"].values
        assert (counts["count"] > 0).all()

    def test_flexible_column_names(self):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        # Uses alternate column names
        meta = pd.DataFrame({
            "id": ["S1", "S2", "S3"],
            "date": ["2024-01-01", "2024-01-02", "2024-01-01"],
            "lineage": ["BA.2", "BA.2", "XBB"],
        })
        counts = metadata_to_clade_counts(meta)
        assert "BA.2" in counts["clade"].values

    def test_location_filter(self):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        meta = _make_metadata_df({"BA.2": (10, 1.0)}, n_days=30)
        counts_usa = metadata_to_clade_counts(meta, location_filter="USA")
        counts_all = metadata_to_clade_counts(meta)
        # All sequences are from USA, so counts should be equal
        assert counts_usa["count"].sum() == counts_all["count"].sum()

    def test_location_filter_no_match_falls_back(self, capsys):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        meta = _make_metadata_df({"BA.2": (5, 1.0)}, n_days=10)
        # Filter for a location that doesn't exist — should warn and fall back
        counts = metadata_to_clade_counts(meta, location_filter="Antarctica")
        assert len(counts) > 0  # fell back to all locations
        captured = capsys.readouterr()
        assert "Warning" in captured.err

    def test_missing_date_column_raises(self):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        meta = pd.DataFrame({"accession": ["S1"], "lineage": ["BA.2"]})
        with pytest.raises(ValueError, match="date"):
            metadata_to_clade_counts(meta)

    def test_missing_clade_column_raises(self):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        meta = pd.DataFrame({"accession": ["S1"], "date": ["2024-01-01"]})
        with pytest.raises(ValueError, match="clade"):
            metadata_to_clade_counts(meta)

    def test_invalid_dates_dropped(self):
        from qprimer_designer.external.evofr import metadata_to_clade_counts

        meta = pd.DataFrame({
            "accession": ["S1", "S2", "S3"],
            "date": ["2024-01-01", "not-a-date", "2024-01-02"],
            "clade": ["BA.2", "BA.2", "XBB"],
        })
        counts = metadata_to_clade_counts(meta)
        # Row with bad date is dropped; two valid rows remain
        assert counts["count"].sum() == 2


# ---------------------------------------------------------------------------
# Tests for run_mlr_forecast
# ---------------------------------------------------------------------------

class TestRunMlrForecast:
    def test_basic_forecast_shape(self):
        from qprimer_designer.external.evofr import run_mlr_forecast

        counts = _make_clade_counts(
            {"Alpha": (50, 0.97), "Delta": (30, 1.0), "Omicron": (20, 1.05)},
            n_days=60,
        )
        result = run_mlr_forecast(counts, horizon_days=30, min_sequences=10)

        assert set(result.columns) == {"clade", "current_freq", "predicted_freq", "growth_advantage"}
        assert len(result) == 3
        assert abs(result["current_freq"].sum() - 1.0) < 0.01
        assert abs(result["predicted_freq"].sum() - 1.0) < 0.01

    def test_growing_clade_has_higher_future_freq(self):
        from qprimer_designer.external.evofr import run_mlr_forecast

        counts = _make_clade_counts(
            {"slow": (50, 0.99), "fast": (10, 1.05)},
            n_days=90,
        )
        result = run_mlr_forecast(counts, horizon_days=30, min_sequences=5)
        result = result.set_index("clade")

        # Fast-growing clade should have higher predicted frequency relative to current
        assert result.loc["fast", "growth_advantage"] > result.loc["slow", "growth_advantage"]

    def test_sorted_by_predicted_freq_descending(self):
        from qprimer_designer.external.evofr import run_mlr_forecast

        counts = _make_clade_counts(
            {"A": (40, 1.0), "B": (30, 1.02), "C": (20, 0.98)},
            n_days=60,
        )
        result = run_mlr_forecast(counts, horizon_days=30)
        freqs = result["predicted_freq"].tolist()
        assert freqs == sorted(freqs, reverse=True)

    def test_too_few_clades_raises(self):
        from qprimer_designer.external.evofr import run_mlr_forecast

        counts = _make_clade_counts({"OnlyClade": (100, 1.0)}, n_days=60)
        with pytest.raises(ValueError, match="clade"):
            run_mlr_forecast(counts, min_sequences=10)

    def test_min_sequences_filter(self):
        from qprimer_designer.external.evofr import run_mlr_forecast

        counts = _make_clade_counts(
            {"common": (50, 1.0), "rare": (2, 1.0), "alsorare": (3, 1.0)},
            n_days=10,
        )
        # rare and alsorare have < 50 total sequences
        with pytest.raises(ValueError):
            run_mlr_forecast(counts, min_sequences=50)

        # Should work with lower threshold
        result = run_mlr_forecast(counts, min_sequences=5)
        assert "common" in result["clade"].values


# ---------------------------------------------------------------------------
# Tests for forecast_variants command
# ---------------------------------------------------------------------------

class TestForecastVariantsCommand:
    def _make_fasta(self, accessions: list[str], tmp_path: Path) -> Path:
        fasta = tmp_path / "seqs.fa"
        with open(fasta, "w") as f:
            for acc in accessions:
                f.write(f">{acc}\nATCGATCGATCGATCGATCG\n")
        return fasta

    def _make_metadata(self, accessions: list[str], clade_map: dict[str, str], tmp_path: Path) -> Path:
        meta_file = tmp_path / "metadata.csv"
        rows = []
        for i, acc in enumerate(accessions):
            days_ago = len(accessions) - i
            rows.append({
                "accession": acc,
                "collection_date": (date.today() - timedelta(days=days_ago)).isoformat(),
                "pango_lineage": clade_map.get(acc, "Unknown"),
                "country": "USA",
            })
        pd.DataFrame(rows).to_csv(meta_file, index=False)
        return meta_file

    def test_basic_run_produces_outputs(self, tmp_path):
        from qprimer_designer.commands.forecast_variants import run

        # 50 sequences split across 3 clades with enough days of data
        accessions = [f"ACC{i:03d}" for i in range(150)]
        clades = ["BA.2"] * 50 + ["XBB.1.5"] * 60 + ["EG.5"] * 40
        clade_map = dict(zip(accessions, clades))

        # Spread dates over 60 days
        rows = []
        for i, acc in enumerate(accessions):
            days_ago = 60 - (i % 60)
            rows.append({
                "accession": acc,
                "collection_date": (date.today() - timedelta(days=days_ago)).isoformat(),
                "pango_lineage": clade_map[acc],
                "country": "USA",
            })
        meta_df = pd.DataFrame(rows)
        meta_file = tmp_path / "metadata.csv"
        meta_df.to_csv(meta_file, index=False)

        fasta_file = self._make_fasta(accessions, tmp_path)
        out_fa = tmp_path / "future.fa"
        out_tsv = tmp_path / "forecast.tsv"

        class Args:
            sequences = str(fasta_file)
            metadata = str(meta_file)
            horizon = 30
            min_sequences = 10
            top_n = 3
            top_n_per_clade = 20
            location = None
            out_sequences = str(out_fa)
            out_forecast = str(out_tsv)
            param_file = "params.txt"

        run(Args())

        assert out_fa.exists() and out_fa.stat().st_size > 0
        assert out_tsv.exists()

        forecast = pd.read_csv(out_tsv, sep="\t")
        assert set(forecast.columns) >= {"clade", "current_freq", "predicted_freq", "growth_advantage"}
        assert len(forecast) >= 3
