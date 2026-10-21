#!/usr/bin/env python3
"""Generate a geometry ID mapping between two ACTS tracking geometries.

For each sensitive surface in geometry A, the script finds the matching
surface in geometry B by comparing surface centre positions.  The match
is deterministic and does not require any data files.

Output: a CSV file with 10 columns (5 per geometry, column names use the
user-supplied prefixes):

    {a}_packed   {a}_volume   {a}_layer   {a}_sensitive   {a}_extra
    {b}_packed   {b}_volume   {b}_layer   {b}_sensitive   {b}_extra

where ``_packed`` is the raw uint64 ``GeometryIdentifier.value``.

Usage
-----
./run_in_env.sh python3 Examples/Scripts/Python/generate_geoid_map.py \\
    --output geoid_map.csv --prefix-a gen1 --prefix-b gen3
"""

import argparse
import csv
from pathlib import Path

import numpy as np

import acts
from acts.examples.odd import getOpenDataDetector


def _collect_surfaces(tracking_geometry, gctx):
    """Return {raw_geoId: (GeometryIdentifier, centre_array)} for sensitive surfaces."""
    surfaces = {}

    def visitor(surface):
        gid = surface.geometryId
        if gid.sensitive == 0:
            return
        centre = np.array(surface.center(gctx))
        surfaces[gid.value] = (gid, centre)

    tracking_geometry.visitSurfaces(visitor)
    return surfaces


def generate_geoid_map(
    geo_a, geo_b, output_path, prefix_a="a", prefix_b="b", tolerance=0.1
):
    """Match sensitive surfaces between two geometries by centre position.

    Parameters
    ----------
    geo_a, geo_b : acts.TrackingGeometry
        Source and target tracking geometries.
    output_path : Path
        Output CSV file path.
    prefix_a, prefix_b : str
        Column name prefixes for the two geometries.
    tolerance : float
        Maximum allowed Euclidean distance (mm) between matched centres.

    Raises
    ------
    RuntimeError
        If any surface in geo_a has no unique match in geo_b.
    """
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    surfs_a = _collect_surfaces(geo_a, gctx)
    surfs_b = _collect_surfaces(geo_b, gctx)

    centres_b = np.array([c for _, c in surfs_b.values()])
    gids_b = [gid for gid, _ in surfs_b.values()]

    rows = []
    unmatched = []

    for raw_a, (gid_a, centre_a) in sorted(surfs_a.items()):
        dists = np.linalg.norm(centres_b - centre_a, axis=1)
        best_idx = np.argmin(dists)
        best_dist = dists[best_idx]

        if best_dist > tolerance:
            unmatched.append((gid_a, best_dist))
            continue

        n_close = np.sum(dists <= tolerance)
        if n_close > 1:
            raise RuntimeError(
                f"Ambiguous match for surface {gid_a.value} "
                f"(volume={gid_a.volume}, layer={gid_a.layer}, "
                f"sensitive={gid_a.sensitive}): "
                f"{n_close} surfaces within tolerance {tolerance} mm"
            )

        gid_b = gids_b[best_idx]
        rows.append((gid_a, gid_b))

    if unmatched:
        msg_parts = [
            f"  geoId={gid.value} (vol={gid.volume}, lay={gid.layer}, "
            f"sen={gid.sensitive}) — nearest={dist:.2f} mm"
            for gid, dist in unmatched[:10]
        ]
        raise RuntimeError(
            f"{len(unmatched)} surfaces in geometry A have no match in "
            f"geometry B (tolerance={tolerance} mm). First entries:\n"
            + "\n".join(msg_parts)
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        f"{prefix_a}_packed",
        f"{prefix_a}_volume",
        f"{prefix_a}_layer",
        f"{prefix_a}_sensitive",
        f"{prefix_a}_extra",
        f"{prefix_b}_packed",
        f"{prefix_b}_volume",
        f"{prefix_b}_layer",
        f"{prefix_b}_sensitive",
        f"{prefix_b}_extra",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for gid_a, gid_b in rows:
            writer.writerow(
                {
                    f"{prefix_a}_packed": gid_a.value,
                    f"{prefix_a}_volume": gid_a.volume,
                    f"{prefix_a}_layer": gid_a.layer,
                    f"{prefix_a}_sensitive": gid_a.sensitive,
                    f"{prefix_a}_extra": gid_a.extra,
                    f"{prefix_b}_packed": gid_b.value,
                    f"{prefix_b}_volume": gid_b.volume,
                    f"{prefix_b}_layer": gid_b.layer,
                    f"{prefix_b}_sensitive": gid_b.sensitive,
                    f"{prefix_b}_extra": gid_b.extra,
                }
            )

    print(f"Written {len(rows)} surface mappings to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("geoid_map.csv"),
        help="Output CSV path (default: geoid_map.csv)",
    )
    parser.add_argument(
        "--prefix-a",
        default="gen1",
        help="Column name prefix for geometry A (default: gen1)",
    )
    parser.add_argument(
        "--prefix-b",
        default="gen3",
        help="Column name prefix for geometry B (default: gen3)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.1,
        help="Match tolerance in mm (default: 0.1)",
    )
    args = parser.parse_args()

    print("Building ODD Gen1 geometry ...")
    detector_a = getOpenDataDetector()
    geo_a = detector_a.trackingGeometry()

    print("Building ODD Gen3 geometry ...")
    detector_b = getOpenDataDetector(gen3=True)
    geo_b = detector_b.trackingGeometry()

    generate_geoid_map(
        geo_a,
        geo_b,
        output_path=args.output,
        prefix_a=args.prefix_a,
        prefix_b=args.prefix_b,
        tolerance=args.tolerance,
    )


if __name__ == "__main__":
    main()
