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
from pathlib import Path

CONTINENTS = [
    "Africa",
    "Asia",
    "Europe",
    "North America",
    "South America",
    "Oceania",
]

# All supported viruses (excluding SARS-CoV-2 and Influenza A — too large for summary API)
DEFAULT_VIRUSES = [
    # Coronaviruses
    ("SARS_CoV", 694009),
    ("MERS_CoV", 1335626),
    ("HCoV_229E", 11137),
    ("HCoV_OC43", 31631),
    ("HCoV_NL63", 277944),
    ("HCoV_HKU1", 290028),
    # Influenza
    ("Influenza_B", 11520),
    # Respiratory
    ("RSV_A", 208893),
    ("RSV_B", 208895),
    ("Human_metapneumovirus", 162145),
    ("Parainfluenza_virus_1", 12730),
    ("Parainfluenza_virus_3", 11216),
    ("Human_bocavirus_1", 329641),
    ("Human_rhinovirus_A", 147711),
    ("Adenovirus_HAdV_C", 129951),
    ("Measles_virus", 11234),
    ("Mumps_virus", 2560602),
    ("Rubella_virus", 11041),
    # Hemorrhagic fever
    ("Ebola_virus", 186538),
    ("Marburg_virus", 3052505),
    ("Lassa_virus", 3052310),
    ("CCHF_virus", 3052518),
    ("Rift_Valley_Fever_virus", 11588),
    ("Hantaan_virus", 3052480),
    ("SFTS_virus", 1003835),
    # Mosquito/tick-borne
    ("Dengue_virus_1", 11053),
    ("Dengue_virus_2", 11060),
    ("Dengue_virus_3", 11069),
    ("Dengue_virus_4", 11070),
    ("Zika_virus", 64320),
    ("Chikungunya_virus", 37124),
    ("Oropouche_virus", 118655),
    ("West_Nile_virus", 11082),
    ("Yellow_Fever_virus", 11089),
    ("Japanese_Encephalitis_virus", 11072),
    ("Tick_borne_encephalitis_virus", 11084),
    # Blood-borne / hepatitis
    ("HIV_1", 11676),
    ("HIV_2", 11709),
    ("Hepatitis_A_virus", 12092),
    ("Hepatitis_B_virus", 10407),
    ("Hepatitis_C_virus", 3052230),
    ("Hepatitis_E_virus", 1678143),
    # Other
    ("MPOX", 10244),
    ("Nipah_virus", 3052225),
    ("Rabies_virus", 11292),
    ("Enterovirus_A71", 39054),
    ("Enterovirus_D68", 42789),
    ("Norovirus_GII", 142786),
    ("Rotavirus_A", 28875),
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
            cmd, capture_output=True, text=True, timeout=120,
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
    }


def fetch_virus(name: str, taxid: int) -> dict:
    """Fetch per-record country + date data for a single virus."""
    print(f"Fetching {name} (taxid {taxid})...")
    records = []
    total_raw = 0

    for continent in CONTINENTS:
        raw = fetch_continent(taxid, continent)
        print(f"  {continent}: {len(raw)} records")
        total_raw += len(raw)
        for rec in raw:
            parsed = extract_record(rec)
            if parsed:
                records.append(parsed)

    countries = set(r["country"] for r in records)
    result = {
        "name": name,
        "taxid": taxid,
        "total_sequences": len(records),
        "num_countries": len(countries),
        "records": records,
    }
    print(f"  Total: {len(records)} sequences, {len(countries)} countries")
    return result


def main():
    parser = argparse.ArgumentParser(description="Fetch virus geographic distribution")
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).parent / "virus_map_data",
        help="Output directory for JSON files",
    )
    parser.add_argument(
        "--viruses", nargs="+",
        help="Virus names to fetch (default: 5 PoC viruses)",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    viruses = DEFAULT_VIRUSES
    if args.viruses:
        name_set = set(args.viruses)
        viruses = [(n, t) for n, t in DEFAULT_VIRUSES if n in name_set]

    for name, taxid in viruses:
        out_path = args.output_dir / f"{name}.json"
        if out_path.exists():
            print(f"Skipping {name} (already exists)")
            continue
        data = fetch_virus(name, taxid)
        out_path = args.output_dir / f"{name}.json"
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
            "total_sequences": d["total_sequences"],
            "num_countries": d["num_countries"],
        })
    with open(index_path, "w") as f:
        json.dump(index, f, indent=2)
    print(f"\nIndex saved to {index_path}")


if __name__ == "__main__":
    main()
