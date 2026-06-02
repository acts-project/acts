#!/usr/bin/env python3
"""Generate a ColliderML → ACTS geometry ID mapping CSV file.

For each unique (detector, volume, layer, surface) tuple found in a ColliderML
tracker-hits parquet file, this script computes the centroid of all hits with
that tuple, then uses a KDTree over ACTS sensitive surface centres to find the
closest surface.  Each match is validated by calling globalToLocal on the
surface; unmatched tuples are skipped with a warning.

Output: a CSV file with columns
    detector, volume, layer, surface, acts_geo_id
where acts_geo_id is the hex-encoded Acts::GeometryIdentifier value.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/generate_colliderml_geo_map.py \\
    --input  colliderml-sample/CERN__ColliderML-Release-1 \\
    --output colliderml_geo_map.csv
"""

import argparse
import csv
import pathlib
import sys

import numpy as np
import pyarrow.parquet as pq

try:
    from scipy.spatial import KDTree
except ImportError:
    print(
        "ERROR: scipy is required.  Install with:  pip install scipy", file=sys.stderr
    )
    sys.exit(1)

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


def collect_surfaces(tracking_geometry, gctx):
    """Return list of (GeometryIdentifier, centre_xyz) for all sensitive surfaces."""
    surfaces = []

    def visitor(surface):
        gid = surface.geometryId()
        if gid.sensitive() == 0:
            return
        centre = surface.center(gctx)
        surfaces.append((gid, centre))

    tracking_geometry.visitSurfaces(visitor)
    return surfaces


def load_hit_centroids(parquet_path: pathlib.Path):
    """Return dict (det, vol, layer, surf) -> mean_xyz from the first shard."""
    table = pq.read_table(parquet_path)

    from collections import defaultdict

    sums = defaultdict(lambda: np.zeros(3))
    counts = defaultdict(int)

    for row in range(table.num_rows):
        dets = table.column("detector")[row].as_py()
        vols = table.column("volume_id")[row].as_py()
        layers = table.column("layer_id")[row].as_py()
        surfs = table.column("surface_id")[row].as_py()
        xs = table.column("true_x")[row].as_py()
        ys = table.column("true_y")[row].as_py()
        zs = table.column("true_z")[row].as_py()

        for d, v, l, s, x, y, z in zip(dets, vols, layers, surfs, xs, ys, zs):
            key = (int(d), int(v), int(l), int(s))
            sums[key] += np.array([x, y, z])
            counts[key] += 1

    return {k: sums[k] / counts[k] for k in sums}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        "-i",
        type=pathlib.Path,
        required=True,
        help="Root of the ColliderML sample " "(contains ttbar_pu200_tracker_hits/)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=pathlib.Path,
        default=pathlib.Path("colliderml_geo_map.csv"),
        help="Output CSV path (default: colliderml_geo_map.csv)",
    )
    parser.add_argument(
        "--odd-dir",
        type=pathlib.Path,
        default=None,
        help="ODD XML directory (default: auto-detect)",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=5.0,
        help="Maximum KDTree match distance in mm (default: 5)",
    )
    args = parser.parse_args()

    hits_dir = (
        args.input / "ttbar_pu200_tracker_hits" / "data" / "ttbar_pu200_tracker_hits"
    )
    parquet_files = sorted(hits_dir.glob("*.parquet"))
    if not parquet_files:
        print(f"ERROR: no parquet files found under {hits_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading hits from {parquet_files[0]} ...")
    centroids = load_hit_centroids(parquet_files[0])
    print(f"  {len(centroids)} unique (det, vol, layer, surface) tuples")

    gctx = acts.GeometryContext()
    odd_dir = args.odd_dir or getOpenDataDetectorDirectory()
    print("Building ODD geometry ...")
    odd = getOpenDataDetector(odd_dir=odd_dir)
    tgeo = odd.trackingGeometry()

    print("Collecting sensitive surfaces ...")
    surf_list = collect_surfaces(tgeo, gctx)
    print(f"  {len(surf_list)} sensitive surfaces")

    centres = np.array([c for _, c in surf_list])
    tree = KDTree(centres)

    rows = []
    n_matched = 0
    n_unmatched = 0

    for (det, vol, layer, surf), centroid in sorted(centroids.items()):
        dist, idx = tree.query(centroid)
        if dist > args.max_distance:
            print(
                f"  WARN: ({det},{vol},{layer},{surf}) nearest surface is "
                f"{dist:.1f} mm away (>{args.max_distance} mm), skipping"
            )
            n_unmatched += 1
            continue

        acts_geo_id = surf_list[idx][0]
        acts_surface = tgeo.findSurface(acts_geo_id)
        if acts_surface is None:
            n_unmatched += 1
            continue

        local_result = acts_surface.globalToLocal(gctx, centroid, acts.Vector3(0, 0, 0))
        if not local_result.ok():
            print(
                f"  WARN: ({det},{vol},{layer},{surf}) globalToLocal failed "
                f"on nearest surface {acts_geo_id}, skipping"
            )
            n_unmatched += 1
            continue

        rows.append((det, vol, layer, surf, acts_geo_id.value()))
        n_matched += 1

    print(f"Matched {n_matched} / {n_matched + n_unmatched} tuples")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["detector", "volume", "layer", "surface", "acts_geo_id"])
        for det, vol, layer, surf, geo_id_val in rows:
            writer.writerow([det, vol, layer, surf, hex(geo_id_val)])

    print(f"Written to {args.output}")


if __name__ == "__main__":
    main()
