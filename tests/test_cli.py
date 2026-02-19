"""Tests for CLI entry point."""

import sys
from unittest.mock import patch, MagicMock

import pytest

from qprimer_designer.cli import main


class TestCLI:
    """Tests for the main CLI entry point."""

    def test_no_args_exits(self):
        """Running with no args should print help and exit with code 1."""
        with patch("sys.argv", ["qprimer"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1

    def test_version(self):
        """--version should print version and exit."""
        with patch("sys.argv", ["qprimer", "--version"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0

    def test_help(self):
        """--help should exit with code 0."""
        with patch("sys.argv", ["qprimer", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0

    def test_generate_help(self):
        """generate --help should exit with code 0."""
        with patch("sys.argv", ["qprimer", "generate", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0

    def test_evaluate_help(self):
        """evaluate --help should exit with code 0."""
        with patch("sys.argv", ["qprimer", "evaluate", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0

    def test_filter_help(self):
        """filter --help should exit with code 0."""
        with patch("sys.argv", ["qprimer", "filter", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0

    def test_unknown_command_exits(self):
        """Unknown subcommand should cause argparse error."""
        with patch("sys.argv", ["qprimer", "nonexistent"]):
            with pytest.raises(SystemExit):
                main()

    def test_subcommands_registered(self):
        """All expected subcommands should be parseable."""
        subcommands = [
            "generate", "generate-probe", "parse-probe-mapping",
            "prepare-features", "pick-representatives", "prepare-input",
            "evaluate", "filter", "build-output", "select-multiplex",
            "export-report",
        ]
        for cmd in subcommands:
            with patch("sys.argv", ["qprimer", cmd, "--help"]):
                with pytest.raises(SystemExit) as exc_info:
                    main()
                assert exc_info.value.code == 0, f"Failed for subcommand: {cmd}"

    def test_generate_missing_required_args(self):
        """generate without required args should exit with error."""
        with patch("sys.argv", ["qprimer", "generate"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 2
