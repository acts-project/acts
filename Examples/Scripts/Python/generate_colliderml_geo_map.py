#!/usr/bin/env python3
"""Generate a ColliderML → ACTS geometry ID mapping Parquet file.

For each unique (detector, volume, layer, surface) tuple found in a ColliderML
tracker-hits parquet file, this script finds the matching ACTS sensitive surface
using a KDTree over surface centres combined with a normal-vector perpendicular
distance check.  The perpendicular distance |(hit - centre) · normal| is nearly
zero for any hit generated on a surface, regardless of how far the hit is from
the surface centre along the surface (important for long strip sensors).

Output: a Parquet file with columns
    detector  (uint8)   ColliderML detector index
    volume    (uint8)   ColliderML volume_id
    layer     (uint16)  ColliderML layer_id
    surface   (uint32)  ColliderML surface_id
    acts_geo_id (uint64) ACTS GeometryIdentifier raw value

Read in C++ via loadColliderMLGeoIdMap() in ColliderMLInputConverter.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/generate_colliderml_geo_map.py \\
    --input  colliderml-sample/CERN__ColliderML-Release-1 \\
    --output colliderml_geo_map.parquet
"""

import argparse
import pathlib
import sys


import numpy as np
import pyarrow as pa
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
    """Return arrays of (GeometryIdentifier, centre_xyz, normal_xyz) for sensitive surfaces."""
    gids = []
    centres = []
    normals = []

    def visitor(surface):
        gid = surface.geometryId
        if gid.sensitive == 0:
            return
        centre = surface.center(gctx)
        # normal(gctx, global_position, direction) — direction is a dummy for orientation
        normal = surface.normal(gctx, centre, centre)
        gids.append(gid)
        centres.append(np.array(centre))
        normals.append(np.array(normal))

    tracking_geometry.visitSurfaces(visitor)
    return gids, np.array(centres), np.array(normals)


def load_hit_representatives(parquet_path: pathlib.Path):
    """Return dict (det, vol, layer, surf) -> first_hit_xyz from the first shard.

    Uses the FIRST hit position per surface rather than the centroid.  A single
    hit lies exactly on the surface plane; the centroid of many hits on a long
    strip sensor can be far from the surface centre yet still in the plane.
    """
    table = pq.read_table(parquet_path)

    representatives = {}

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
            if key not in representatives:
                representatives[key] = np.array([x, y, z])

    return representatives


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
        default=pathlib.Path("colliderml_geo_map.parquet"),
        help="Output Parquet path (default: colliderml_geo_map.parquet)",
    )
    parser.add_argument(
        "--odd-dir",
        type=pathlib.Path,
        default=None,
        help="ODD XML directory (default: auto-detect)",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=10,
        help="Number of KDTree neighbours to probe (default: 10)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0,
        help="Normal-vector perpendicular distance tolerance in mm (default: 1.0)",
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
    representatives = load_hit_representatives(parquet_files[0])
    print(f"  {len(representatives)} unique (det, vol, layer, surface) tuples")

    gctx = acts.GeometryContext()
    odd_dir = args.odd_dir or getOpenDataDetectorDirectory()
    print("Building ODD geometry ...")
    odd = getOpenDataDetector(odd_dir=odd_dir)
    tgeo = odd.trackingGeometry()

    print("Collecting sensitive surfaces ...")
    surf_gids, centres, normals = collect_surfaces(tgeo, gctx)
    print(f"  {len(surf_gids)} sensitive surfaces")

    tree = KDTree(centres)

    rows = []
    n_matched = 0
    n_unmatched = 0

    for (det, vol, layer, surf), hit_pos in sorted(representatives.items()):
        # Query k nearest neighbours, accept the first whose surface plane contains hit_pos
        dists, indices = tree.query(hit_pos, k=min(args.k, len(surf_gids)))
        matched_idx = None
        for idx in indices:
            delta = hit_pos - centres[idx]
            perp_dist = abs(float(np.dot(delta, normals[idx])))
            if perp_dist < args.tolerance:
                matched_idx = idx
                break

        if matched_idx is None:
            best_perp = min(
                abs(float(np.dot(hit_pos - centres[i], normals[i]))) for i in indices
            )
            print(
                f"  WARN: ({det},{vol},{layer},{surf}) no surface within "
                f"{args.tolerance:.1f} mm normal-distance "
                f"(best perp dist among {args.k} neighbours: {best_perp:.2f} mm)"
            )
            n_unmatched += 1
            continue

        acts_geo_id = surf_gids[matched_idx]
        rows.append((det, vol, layer, surf, acts_geo_id.value))
        n_matched += 1

    print(f"Matched {n_matched} / {n_matched + n_unmatched} tuples")
    if n_unmatched > 0:
        print(
            f"  {n_unmatched} unmatched tuples — these hits will be skipped "
            f"in ColliderMLInputConverter. Consider increasing --tolerance or --k."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Write as Parquet using the same type widths that C++ reads back:
    # detector/volume → uint8, layer → uint16, surface → uint32, acts_geo_id → uint64
    schema = pa.schema(
        [
            pa.field("detector", pa.uint8()),
            pa.field("volume", pa.uint8()),
            pa.field("layer", pa.uint16()),
            pa.field("surface", pa.uint32()),
            pa.field("acts_geo_id", pa.uint64()),
        ]
    )
    table = pa.table(
        {
            "detector": pa.array([r[0] for r in rows], type=pa.uint8()),
            "volume": pa.array([r[1] for r in rows], type=pa.uint8()),
            "layer": pa.array([r[2] for r in rows], type=pa.uint16()),
            "surface": pa.array([r[3] for r in rows], type=pa.uint32()),
            "acts_geo_id": pa.array([r[4] for r in rows], type=pa.uint64()),
        },
        schema=schema,
    )
    pq.write_table(table, args.output, compression="snappy")
    print(f"Written to {args.output}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
