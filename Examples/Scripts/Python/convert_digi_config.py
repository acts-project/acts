#!/usr/bin/env python3
"""Convert a digitization config JSON from one geometry to another.

Uses a geometry ID map CSV (produced by ``generate_geoid_map.py``) to remap
the volume IDs in a digitization config JSON file.  The digi config entries
are keyed by volume, so the conversion builds a volume→volume mapping from
the full surface-level geo ID map and applies it.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/convert_digi_config.py \\
    --input  Examples/Configs/odd-digi-smearing-config-notime.json \\
    --geoid-map geoid_map.csv \\
    --source-prefix gen1 --target-prefix gen3 \\
    --output odd-digi-smearing-config-notime-gen3.json
"""

import argparse
import csv
import json
import tempfile
from pathlib import Path


def convert_digi_config(
    digi_config_path,
    geoid_map_path,
    source_prefix="gen1",
    target_prefix="gen3",
    output_path=None,
):
    """Remap volume IDs in a digi config JSON using a geo ID map CSV.

    Parameters
    ----------
    digi_config_path : Path
        Input digitization config JSON file.
    geoid_map_path : Path
        Geo ID map CSV produced by ``generate_geoid_map.py``.
    source_prefix, target_prefix : str
        Column name prefixes in the CSV for source and target geometries.
    output_path : Path or None
        Output JSON path.  If None, writes to a temporary file and returns
        the path.

    Returns
    -------
    Path
        Path to the output JSON file.
    """
    volume_map = _build_volume_map(geoid_map_path, source_prefix, target_prefix)

    with open(digi_config_path) as f:
        config = json.load(f)

    new_entries = []
    for entry in config["entries"]:
        src_vol = entry.get("volume")
        if src_vol is None:
            new_entries.append(entry)
            continue

        if src_vol not in volume_map:
            raise RuntimeError(
                f"Volume {src_vol} in digi config has no mapping in geo ID "
                f"map {geoid_map_path}"
            )

        for tgt_vol in sorted(volume_map[src_vol]):
            new_entry = dict(entry)
            new_entry["volume"] = tgt_vol
            new_entries.append(new_entry)

    config["entries"] = new_entries

    if output_path is None:
        fd, output_path = tempfile.mkstemp(
            suffix=".json", prefix="digi_config_converted_"
        )
        with open(fd, "w") as f:
            json.dump(config, f, indent=4)
        output_path = Path(output_path)
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(config, f, indent=4)

    return output_path


def _build_volume_map(geoid_map_path, source_prefix, target_prefix):
    """Build a source_volume → set of target_volumes mapping from the geo ID CSV."""
    volume_map = {}
    with open(geoid_map_path, newline="") as f:
        reader = csv.DictReader(f)
        src_col = f"{source_prefix}_volume"
        tgt_col = f"{target_prefix}_volume"
        for row in reader:
            src_vol = int(row[src_col])
            tgt_vol = int(row[tgt_col])
            volume_map.setdefault(src_vol, set()).add(tgt_vol)
    return volume_map


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Input digitization config JSON",
    )
    parser.add_argument(
        "--geoid-map",
        type=Path,
        required=True,
        help="Geo ID map CSV (from generate_geoid_map.py)",
    )
    parser.add_argument(
        "--source-prefix",
        default="gen1",
        help="Source geometry column prefix (default: gen1)",
    )
    parser.add_argument(
        "--target-prefix",
        default="gen3",
        help="Target geometry column prefix (default: gen3)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output JSON path (default: auto-named temp file)",
    )
    args = parser.parse_args()

    out = convert_digi_config(
        args.input,
        args.geoid_map,
        source_prefix=args.source_prefix,
        target_prefix=args.target_prefix,
        output_path=args.output,
    )
    print(f"Written converted digi config to {out}")


if __name__ == "__main__":
    main()
