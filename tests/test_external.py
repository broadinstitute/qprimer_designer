"""Tests for external tool wrappers (mocked - no external tools needed)."""

import subprocess
import warnings
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

from qprimer_designer.external.vienna import (
    _validate_sequence,
    compute_dimer_dg,
    compute_batch_dimer_dg,
    compute_self_dimer_dg,
    find_rnaduplex,
    VALID_DNA_CHARS,
    _IUPAC_AMBIGUITY_CHARS,
)
from qprimer_designer.external.bowtie import (
    find_bowtie2,
    find_bowtie2_build,
    build_index,
    align_primers,
)
from qprimer_designer.external.mafft import (
    find_mafft,
    align_sequences,
)


# --- Vienna RNA tests ---


class TestValidateSequence:
    """Tests for _validate_sequence."""

    def test_valid_sequence(self):
        assert _validate_sequence("ATCGATCG") == "ATCGATCG"

    def test_strips_whitespace(self):
        assert _validate_sequence("  ATCG  ") == "ATCG"

    def test_none_raises_type_error(self):
        with pytest.raises(TypeError, match="cannot be None"):
            _validate_sequence(None)

    def test_not_string_raises_type_error(self):
        with pytest.raises(TypeError, match="must be a string"):
            _validate_sequence(123)

    def test_empty_raises_value_error(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            _validate_sequence("")

    def test_whitespace_only_raises_value_error(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            _validate_sequence("   ")

    def test_internal_whitespace_raises_value_error(self):
        with pytest.raises(ValueError, match="contains internal whitespace"):
            _validate_sequence("AT CG")

    def test_invalid_chars_raises_value_error(self):
        with pytest.raises(ValueError, match="contains invalid characters"):
            _validate_sequence("ATCGX")

    def test_iupac_codes_replaced_with_n(self):
        """IUPAC ambiguity codes should be replaced with N and warn."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _validate_sequence("ARCG")
            assert result == "ANCG"
            assert len(w) == 1
            assert "IUPAC ambiguity" in str(w[0].message)

    def test_valid_dna_chars(self):
        """All valid DNA characters should pass."""
        result = _validate_sequence("ATCGatcgNn-")
        assert result == "ATCGatcgNn-"

    def test_gap_char(self):
        assert _validate_sequence("-") == "-"

    def test_n_char(self):
        assert _validate_sequence("N") == "N"

    def test_custom_name_in_error(self):
        """Custom name should appear in error messages."""
        with pytest.raises(TypeError, match="primer1"):
            _validate_sequence(None, name="primer1")


class TestComputeDimerDg:
    """Tests for compute_dimer_dg (mocked)."""

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    @patch("qprimer_designer.external.vienna.subprocess.run")
    def test_parses_output(self, mock_run, mock_find):
        mock_run.return_value = MagicMock(
            stdout="..((.....))&..((.....)). : 1,12 : 1,12 ( -5.30)\n",
            returncode=0,
        )
        result = compute_dimer_dg("ATCGATCGATCG", "CGATCGATCGAT")
        assert result == -5.3

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    @patch("qprimer_designer.external.vienna.subprocess.run")
    def test_positive_dg(self, mock_run, mock_find):
        mock_run.return_value = MagicMock(
            stdout=".&. : 0,0 : 0,0 ( 0.00)\n",
            returncode=0,
        )
        result = compute_dimer_dg("AAAA", "AAAA")
        assert result == 0.0

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    @patch("qprimer_designer.external.vienna.subprocess.run")
    def test_unparseable_output_raises(self, mock_run, mock_find):
        mock_run.return_value = MagicMock(
            stdout="no match\n",
            returncode=0,
        )
        with pytest.raises(ValueError, match="Could not parse"):
            compute_dimer_dg("ATCG", "GCTA")

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    def test_none_seq_raises(self, mock_find):
        with pytest.raises(TypeError):
            compute_dimer_dg(None, "ATCG")

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    def test_empty_seq_raises(self, mock_find):
        with pytest.raises(ValueError):
            compute_dimer_dg("", "ATCG")


class TestComputeBatchDimerDg:
    """Tests for compute_batch_dimer_dg (mocked)."""

    def test_empty_list(self):
        assert compute_batch_dimer_dg([]) == []

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    @patch("qprimer_designer.external.vienna.subprocess.run")
    def test_multiple_pairs(self, mock_run, mock_find):
        mock_run.return_value = MagicMock(
            stdout=(
                "..((.....))&..((.....)). : 1,12 : 1,12 ( -5.30)\n"
                "..((.....))&..((.....)). : 1,12 : 1,12 ( -3.20)\n"
            ),
            returncode=0,
        )
        pairs = [("ATCGATCGATCG", "CGATCGATCGAT"), ("AAAA", "TTTT")]
        result = compute_batch_dimer_dg(pairs)
        assert len(result) == 2
        assert result[0] == -5.3
        assert result[1] == -3.2

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    @patch("qprimer_designer.external.vienna.subprocess.run")
    def test_count_mismatch_raises(self, mock_run, mock_find):
        mock_run.return_value = MagicMock(
            stdout="..((.....))&..((.....)). : 1,12 : 1,12 ( -5.30)\n",
            returncode=0,
        )
        pairs = [("ATCG", "GCTA"), ("AAAA", "TTTT")]
        with pytest.raises(ValueError, match="Expected 2"):
            compute_batch_dimer_dg(pairs)


class TestComputeSelfDimerDg:
    """Tests for compute_self_dimer_dg (mocked)."""

    @patch("qprimer_designer.external.vienna.find_rnaduplex", return_value="/usr/bin/RNAduplex")
    @patch("qprimer_designer.external.vienna.subprocess.run")
    def test_calls_compute_dimer_dg_with_same_seq(self, mock_run, mock_find):
        mock_run.return_value = MagicMock(
            stdout="..((.....))&..((.....)). : 1,12 : 1,12 ( -2.50)\n",
            returncode=0,
        )
        result = compute_self_dimer_dg("ATCGATCGATCG")
        assert result == -2.5
        # Verify it was called with the same sequence twice
        call_args = mock_run.call_args
        input_str = call_args.kwargs.get("input") or call_args[1].get("input")
        lines = input_str.strip().split("\n")
        assert lines[0] == lines[1]


class TestFindRnaduplex:
    """Tests for find_rnaduplex."""

    def setup_method(self):
        """Clear LRU cache before each test."""
        find_rnaduplex.cache_clear()

    @patch("qprimer_designer.external.vienna.shutil.which", return_value="/usr/bin/RNAduplex")
    def test_found(self, mock_which):
        assert find_rnaduplex() == "/usr/bin/RNAduplex"

    @patch("qprimer_designer.external.vienna.shutil.which", return_value=None)
    def test_not_found(self, mock_which):
        with pytest.raises(FileNotFoundError, match="RNAduplex not found"):
            find_rnaduplex()


# --- Bowtie2 tests ---


class TestFindBowtie2:
    """Tests for find_bowtie2."""

    @patch("qprimer_designer.external.bowtie.shutil.which", return_value="/usr/bin/bowtie2")
    def test_found(self, mock_which):
        assert find_bowtie2() == "/usr/bin/bowtie2"

    @patch("qprimer_designer.external.bowtie.shutil.which", return_value=None)
    def test_not_found(self, mock_which):
        with pytest.raises(FileNotFoundError, match="bowtie2 not found"):
            find_bowtie2()


class TestFindBowtie2Build:
    """Tests for find_bowtie2_build."""

    @patch("qprimer_designer.external.bowtie.shutil.which", return_value="/usr/bin/bowtie2-build")
    def test_found(self, mock_which):
        assert find_bowtie2_build() == "/usr/bin/bowtie2-build"

    @patch("qprimer_designer.external.bowtie.shutil.which", return_value=None)
    def test_not_found(self, mock_which):
        with pytest.raises(FileNotFoundError, match="bowtie2-build not found"):
            find_bowtie2_build()


class TestBuildIndex:
    """Tests for build_index (mocked)."""

    @patch("qprimer_designer.external.bowtie.find_bowtie2_build", return_value="/usr/bin/bowtie2-build")
    @patch("qprimer_designer.external.bowtie.subprocess.run")
    def test_success(self, mock_run, mock_find, tmp_path):
        fasta = tmp_path / "test.fa"
        fasta.write_text(">seq1\nATCG\n")
        mock_run.return_value = MagicMock(returncode=0)

        build_index(fasta, tmp_path / "idx")
        mock_run.assert_called_once()

    @patch("qprimer_designer.external.bowtie.find_bowtie2_build", return_value="/usr/bin/bowtie2-build")
    def test_missing_input(self, mock_find, tmp_path):
        with pytest.raises(FileNotFoundError, match="Input FASTA file not found"):
            build_index(tmp_path / "nonexistent.fa", tmp_path / "idx")

    @patch("qprimer_designer.external.bowtie.find_bowtie2_build", return_value="/usr/bin/bowtie2-build")
    @patch("qprimer_designer.external.bowtie.subprocess.run")
    def test_subprocess_failure(self, mock_run, mock_find, tmp_path):
        fasta = tmp_path / "test.fa"
        fasta.write_text(">seq1\nATCG\n")
        mock_run.side_effect = subprocess.CalledProcessError(1, "bowtie2-build")

        with pytest.raises(subprocess.CalledProcessError):
            build_index(fasta, tmp_path / "idx")


class TestAlignPrimers:
    """Tests for align_primers (mocked)."""

    @patch("qprimer_designer.external.bowtie.find_bowtie2", return_value="/usr/bin/bowtie2")
    @patch("qprimer_designer.external.bowtie.subprocess.run")
    def test_success(self, mock_run, mock_find, tmp_path):
        query = tmp_path / "primers.fa"
        query.write_text(">p1\nATCG\n")
        mock_run.return_value = MagicMock(returncode=0)

        align_primers(tmp_path / "idx", query, tmp_path / "out.sam")
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "--local" in cmd
        assert "--very-sensitive-local" in cmd

    @patch("qprimer_designer.external.bowtie.find_bowtie2", return_value="/usr/bin/bowtie2")
    @patch("qprimer_designer.external.bowtie.subprocess.run")
    def test_no_local(self, mock_run, mock_find, tmp_path):
        query = tmp_path / "primers.fa"
        query.write_text(">p1\nATCG\n")
        mock_run.return_value = MagicMock(returncode=0)

        align_primers(tmp_path / "idx", query, tmp_path / "out.sam", local=False, very_sensitive=True)
        cmd = mock_run.call_args[0][0]
        assert "--local" not in cmd
        assert "--very-sensitive" in cmd

    @patch("qprimer_designer.external.bowtie.find_bowtie2", return_value="/usr/bin/bowtie2")
    def test_missing_query(self, mock_find, tmp_path):
        with pytest.raises(FileNotFoundError, match="Query FASTA file not found"):
            align_primers(tmp_path / "idx", tmp_path / "nonexistent.fa", tmp_path / "out.sam")


# --- MAFFT tests ---


class TestFindMafft:
    """Tests for find_mafft."""

    @patch("qprimer_designer.external.mafft.shutil.which", return_value="/usr/bin/mafft")
    def test_found(self, mock_which):
        assert find_mafft() == "/usr/bin/mafft"

    @patch("qprimer_designer.external.mafft.shutil.which", return_value=None)
    def test_not_found(self, mock_which):
        with pytest.raises(FileNotFoundError, match="mafft not found"):
            find_mafft()


class TestAlignSequences:
    """Tests for align_sequences (mocked)."""

    @patch("qprimer_designer.external.mafft.find_mafft", return_value="/usr/bin/mafft")
    def test_missing_input(self, mock_find, tmp_path):
        with pytest.raises(FileNotFoundError, match="Input FASTA file not found"):
            align_sequences(tmp_path / "nonexistent.fa", tmp_path / "out.fa")

    @patch("qprimer_designer.external.mafft.find_mafft", return_value="/usr/bin/mafft")
    def test_empty_input(self, mock_find, tmp_path):
        fasta = tmp_path / "empty.fa"
        fasta.write_text("")
        with pytest.raises(ValueError, match="empty"):
            align_sequences(fasta, tmp_path / "out.fa")

    @patch("qprimer_designer.external.mafft.find_mafft", return_value="/usr/bin/mafft")
    def test_invalid_format(self, mock_find, tmp_path):
        fasta = tmp_path / "bad.fa"
        fasta.write_text("not a fasta file\n")
        with pytest.raises(ValueError, match="does not appear to be in FASTA format"):
            align_sequences(fasta, tmp_path / "out.fa")

    @patch("qprimer_designer.external.mafft.find_mafft", return_value="/usr/bin/mafft")
    @patch("qprimer_designer.external.mafft.subprocess.run")
    def test_success(self, mock_run, mock_find, tmp_path):
        fasta = tmp_path / "input.fa"
        fasta.write_text(">seq1\nATCG\n>seq2\nGCTA\n")
        output = tmp_path / "output.fa"

        mock_run.return_value = MagicMock(returncode=0)

        align_sequences(fasta, output)
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "mafft" in cmd[0]
        assert "--auto" in cmd
        assert "--quiet" in cmd

    @patch("qprimer_designer.external.mafft.find_mafft", return_value="/usr/bin/mafft")
    @patch("qprimer_designer.external.mafft.subprocess.run")
    def test_no_auto_no_quiet(self, mock_run, mock_find, tmp_path):
        fasta = tmp_path / "input.fa"
        fasta.write_text(">seq1\nATCG\n>seq2\nGCTA\n")
        output = tmp_path / "output.fa"

        mock_run.return_value = MagicMock(returncode=0)

        align_sequences(fasta, output, auto=False, quiet=False)
        cmd = mock_run.call_args[0][0]
        assert "--auto" not in cmd
        assert "--quiet" not in cmd


class TestViennaConstants:
    """Tests for vienna module constants."""

    def test_valid_dna_chars(self):
        """Valid chars include standard bases, N, and gap."""
        assert 'A' in VALID_DNA_CHARS
        assert 'T' in VALID_DNA_CHARS
        assert 'C' in VALID_DNA_CHARS
        assert 'G' in VALID_DNA_CHARS
        assert 'N' in VALID_DNA_CHARS
        assert 'n' in VALID_DNA_CHARS
        assert '-' in VALID_DNA_CHARS
        assert 'X' not in VALID_DNA_CHARS

    def test_iupac_ambiguity_chars(self):
        """IUPAC codes should include R, Y, W, S, M, K, B, D, H, V."""
        for c in "RYWSMKBDHVryswmkbdhv":
            assert c in _IUPAC_AMBIGUITY_CHARS
        assert 'A' not in _IUPAC_AMBIGUITY_CHARS
        assert 'N' not in _IUPAC_AMBIGUITY_CHARS
