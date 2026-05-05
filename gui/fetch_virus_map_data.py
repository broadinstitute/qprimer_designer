#!/usr/bin/env python3
"""Fetch virus sequence geographic distribution from NCBI for the Virus Map feature.

Queries NCBI datasets summary API per continent, parses country + date per record,
and saves results as JSON for each virus. The GUI can then filter by date client-side.

Usage:
    python fetch_virus_map_data.py [--output-dir DIR] [--viruses NAME1 NAME2 ...]
"""
import argparse
import json
import subprocess
import sys
from fnmatch import fnmatch
from pathlib import Path

CONTINENTS = [
    "Africa",
    "Asia",
    "Europe",
    "North America",
    "South America",
    "Oceania",
]

# All supported viruses
# Tuple format: (name, taxid, subtype_filter or "")
DEFAULT_VIRUSES = [
    # Coronaviruses (SARS-CoV-2 excluded — too large for summary API)
    ("SARS_CoV", 694009, ""),
    ("MERS_CoV", 1335626, ""),
    ("HCoV_229E", 11137, ""),
    ("HCoV_OC43", 31631, ""),
    ("HCoV_NL63", 277944, ""),
    ("HCoV_HKU1", 290028, ""),
    # Influenza (Influenza A excluded — too large for summary API)
    ("Influenza_B", 11520, ""),
    # Respiratory
    ("RSV_A", 208893, ""),
    ("RSV_B", 208895, ""),
    ("Human_metapneumovirus", 162145, ""),
    ("Parainfluenza_virus_1", 12730, ""),
    ("Parainfluenza_virus_3", 11216, ""),
    ("Human_bocavirus_1", 329641, ""),
    ("Human_rhinovirus_A", 147711, ""),
    ("Adenovirus_HAdV_C", 129951, ""),
    ("Measles_virus", 11234, ""),
    ("Mumps_virus", 2560602, ""),
    ("Rubella_virus", 11041, ""),
    # Hemorrhagic fever
    ("Ebola_virus", 186538, ""),
    ("Marburg_virus", 3052505, ""),
    ("Lassa_virus", 3052310, ""),
    ("CCHF_virus", 3052518, ""),
    ("Rift_Valley_Fever_virus", 11588, ""),
    ("Hantaan_virus", 3052480, ""),
    ("SFTS_virus", 1003835, ""),
    # Mosquito/tick-borne
    ("Dengue_virus_1", 11053, ""),
    ("Dengue_virus_2", 11060, ""),
    ("Dengue_virus_3", 11069, ""),
    ("Dengue_virus_4", 11070, ""),
    ("Zika_virus", 64320, ""),
    ("Chikungunya_virus", 37124, ""),
    ("Oropouche_virus", 118655, ""),
    ("West_Nile_virus", 11082, ""),
    ("Yellow_Fever_virus", 11089, ""),
    ("Japanese_Encephalitis_virus", 11072, ""),
    ("Tick_borne_encephalitis_virus", 11084, ""),
    # Blood-borne / hepatitis (HIV-1 excluded — too large for summary API)
    #("HIV_1", 11676, ""),
    ("HIV_2", 11709, ""),
    ("Hepatitis_A_virus", 12092, ""),
    ("Hepatitis_B_virus", 10407, ""),
    ("Hepatitis_C_virus", 3052230, ""),
    ("Hepatitis_E_virus", 1678143, ""),
    # Other
    ("MPOX", 10244, ""),
    ("Nipah_virus", 3052225, ""),
    ("Rabies_virus", 11292, ""),
    ("Enterovirus_A71", 39054, ""),
    ("Enterovirus_D68", 42789, ""),
    ("Norovirus_GII", 142786, ""),
    ("Rotavirus_A", 28875, ""),
]


def fetch_continent(taxid: int, continent: str) -> list[dict]:
    """Run datasets summary for one taxid + continent, return parsed records."""
    cmd = [
        "datasets", "summary", "virus", "genome",
        "taxon", str(taxid),
        "--geo-location", continent,
        "--as-json-lines",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=1800,
        )
    except subprocess.TimeoutExpired:
        print(f"  [WARN] Timeout for {continent}", file=sys.stderr)
        return []

    records = []
    for line in result.stdout.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            continue

    return records


def _matches_subtype(record: dict, subtype_filter: str) -> bool:
    """Check if a record matches the subtype filter pattern."""
    if not subtype_filter:
        return True
    organism = record.get("virus", {}).get("organism_name", "")
    isolate_name = record.get("isolate", {}).get("name", "")
    header = f"{organism} {isolate_name}"
    pattern = subtype_filter.strip()
    if not pattern.startswith("*"):
        pattern = "*" + pattern
    if not pattern.endswith("*"):
        pattern = pattern + "*"
    return fnmatch(header.lower(), pattern.lower())


def extract_record(record: dict) -> dict | None:
    """Extract country, collection_date, and release_date from a record."""
    loc = record.get("location", {})
    geo = loc.get("geographic_location", "")
    if not geo:
        return None
    country = geo.split(":")[0].strip()
    if not country:
        return None

    # Collection date: "2024-08", "2024-08-15", or missing
    collection_date = record.get("isolate", {}).get("collection_date", "")
    # Release date: "2025-04-08T00:00:00Z"
    release_date = record.get("release_date", "")
    if release_date:
        release_date = release_date[:10]  # keep YYYY-MM-DD only

    return {
        "country": country,
        "collection_date": collection_date or "",
        "release_date": release_date or "",
        "continent": "",  # filled by fetch_virus()
    }


def _fetch_raw_all_continents(taxid: int, label: str) -> list[tuple[dict, str]]:
    """Fetch raw records for all continents. Returns [(record, continent), ...]."""
    all_raw = []
    for continent in CONTINENTS:
        raw = fetch_continent(taxid, continent)
        print(f"  {continent}: {len(raw)} records")
        for rec in raw:
            all_raw.append((rec, continent))
    print(f"  {label} raw total: {len(all_raw)} records")
    return all_raw


def _build_result(name: str, taxid: int, subtype_filter: str,
                  raw_with_continent: list[tuple[dict, str]]) -> dict:
    """Build result dict from raw records, optionally filtering by subtype."""
    records = []
    for rec, continent in raw_with_continent:
        if subtype_filter and not _matches_subtype(rec, subtype_filter):
            continue
        parsed = extract_record(rec)
        if parsed:
            parsed["continent"] = continent
            records.append(parsed)

    countries = set(r["country"] for r in records)
    print(f"  {name}: {len(records)} sequences, {len(countries)} countries")
    return {
        "name": name,
        "taxid": taxid,
        "subtype_filter": subtype_filter,
        "total_sequences": len(records),
        "num_countries": len(countries),
        "records": records,
    }


def fetch_virus(name: str, taxid: int, subtype_filter: str = "") -> dict:
    """Fetch per-record country + date data for a single virus."""
    label = f"{name} (taxid {taxid})"
    if subtype_filter:
        label += f" [subtype: {subtype_filter}]"
    print(f"Fetching {label}...")

    all_raw = _fetch_raw_all_continents(taxid, label)
    return _build_result(name, taxid, subtype_filter, all_raw)


def fetch_virus_group(entries: list[tuple[str, int, str]]) -> list[dict]:
    """Fetch once for a shared taxid and split into multiple results by subtype.

    All entries must share the same taxid.
    """
    taxid = entries[0][1]
    subtypes = [s for _, _, s in entries]
    print(f"Fetching taxid {taxid} (subtypes: {', '.join(subtypes)})...")

    all_raw = _fetch_raw_all_continents(taxid, f"taxid {taxid}")

    results = []
    for name, _, subtype_filter in entries:
        result = _build_result(name, taxid, subtype_filter, all_raw)
        results.append(result)
    return results


def main():
    parser = argparse.ArgumentParser(description="Fetch virus geographic distribution")
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).parent / "virus_map_data",
        help="Output directory for JSON files",
    )
    parser.add_argument(
        "--viruses", nargs="+",
        help="Virus names to fetch (default: all)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download even if file already exists",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    viruses = DEFAULT_VIRUSES
    if args.viruses:
        name_set = set(args.viruses)
        viruses = [(n, t, s) for n, t, s in DEFAULT_VIRUSES if n in name_set]

    # Group entries by taxid so we fetch once for shared taxids (e.g. Influenza A subtypes)
    from collections import OrderedDict
    taxid_groups: OrderedDict[int, list[tuple[str, int, str]]] = OrderedDict()
    for name, taxid, subtype in viruses:
        taxid_groups.setdefault(taxid, []).append((name, taxid, subtype))

    for taxid, entries in taxid_groups.items():
        # Check which entries still need downloading
        needed = []
        for name, tid, sub in entries:
            out_path = args.output_dir / f"{name}.json"
            if out_path.exists() and not args.force:
                print(f"Skipping {name} (already exists, use --force to re-download)")
            else:
                needed.append((name, tid, sub))

        if not needed:
            continue

        # If multiple entries share a taxid, fetch once and split
        if len(entries) > 1:
            results = fetch_virus_group(needed)
        else:
            name, tid, sub = needed[0]
            results = [fetch_virus(name, tid, sub)]

        for data in results:
            out_path = args.output_dir / f"{data['name']}.json"
            with open(out_path, "w") as f:
                json.dump(data, f)
            print(f"  Saved to {out_path}")

    # Also save a combined index
    index_path = args.output_dir / "index.json"
    index = []
    for p in sorted(args.output_dir.glob("*.json")):
        if p.name == "index.json":
            continue
        with open(p) as f:
            d = json.load(f)
        index.append({
            "name": d["name"],
            "taxid": d["taxid"],
            "subtype_filter": d.get("subtype_filter", ""),
            "total_sequences": d["total_sequences"],
            "num_countries": d["num_countries"],
        })
    with open(index_path, "w") as f:
        json.dump(index, f, indent=2)
    print(f"\nIndex saved to {index_path}")


if __name__ == "__main__":
    main()
