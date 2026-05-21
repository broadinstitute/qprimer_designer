"""NCBI BLAST+ wrapper for remote primer specificity checking."""

import functools
import os
import shutil
import subprocess
import tempfile

BLAST_OUTFMT_FIELDS = [
    "qseqid", "sseqid", "stitle", "pident", "length",
    "mismatch", "evalue", "bitscore", "staxids",
]

_NUMERIC_FIELDS = {
    "pident": float,
    "length": int,
    "mismatch": int,
    "evalue": float,
    "bitscore": float,
}


@functools.lru_cache(maxsize=1)
def find_blastn() -> str:
    path = shutil.which("blastn")
    if path is None:
        raise FileNotFoundError(
            "blastn not found in PATH. "
            "Install BLAST+: conda install -c bioconda blast"
        )
    return path


def parse_blast_results(raw_output: str, fields: list[str] | None = None) -> list[dict]:
    if fields is None:
        fields = BLAST_OUTFMT_FIELDS

    results = []
    for line in raw_output.strip().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != len(fields):
            continue
        row = {}
        for field, value in zip(fields, parts):
            if field in _NUMERIC_FIELDS:
                try:
                    row[field] = _NUMERIC_FIELDS[field](value)
                except (ValueError, TypeError):
                    row[field] = value
            else:
                row[field] = value
        results.append(row)
    return results


def run_blastn_remote(
    sequences: dict[str, str],
    negative_taxids: list[int] | None = None,
    evalue: float = 10,
    word_size: int = 11,
    max_target_seqs: int = 50,
    timeout: int = 600,
) -> list[dict]:
    blastn = find_blastn()

    if not sequences:
        raise ValueError("sequences must be a non-empty dict")
    for name, seq in sequences.items():
        if not seq or not seq.strip():
            raise ValueError(f"Empty sequence for '{name}'")

    if negative_taxids is not None:
        for tid in negative_taxids:
            if not isinstance(tid, int) or tid <= 0:
                raise ValueError(f"Invalid TaxID: {tid} (must be a positive integer)")

    query_file = None
    try:
        query_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, prefix="blast_query_",
        )
        for name, seq in sequences.items():
            query_file.write(f">{name}\n{seq.strip()}\n")
        query_file.flush()
        query_file.close()

        outfmt_str = "6 " + " ".join(BLAST_OUTFMT_FIELDS)

        cmd = [
            blastn, "-remote", "-db", "nt",
            "-query", query_file.name,
            "-task", "blastn-short",
            "-evalue", str(evalue),
            "-word_size", str(word_size),
            "-dust", "no",
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", outfmt_str,
        ]

        # -negative_taxidlist is incompatible with -remote; use -entrez_query
        if negative_taxids:
            exclusion = " AND ".join(f"NOT txid{tid}[ORGN]" for tid in negative_taxids)
            cmd.extend(["-entrez_query", exclusion])

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"blastn exited with code {result.returncode}: {result.stderr.strip()}"
            )

        return parse_blast_results(result.stdout, BLAST_OUTFMT_FIELDS)

    finally:
        if query_file is not None:
            try:
                os.unlink(query_file.name)
            except OSError:
                pass
