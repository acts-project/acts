#!/usr/bin/env python3
"""2D display of tracking-geometry surfaces for ODD or Generic detectors."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
import math
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, FrozenSet, Tuple

import acts
import acts.examples
import numpy as np

from acts.utilities import (
    auto_axis_limits,
    collect_material_surfaces,
    collect_sensitive_surfaces,
    filter_surfaces_by_layer_ids,
    filter_surfaces_by_volume_ids,
    group_surfaces_by_volume_layer,
    parse_layer_selection_spec,
    parse_subsystem_region_spec,
    sort_surfaces_by_geometry_id,
    surface_to_2d_patches,
    volume_ids_for_selection,
)

# Volume ids for legacy ODD (Gen1 DD4hep)
_ODD_SUBSYSTEM_VOLUME_IDS: Dict[str, FrozenSet[int]] = {
    "all": frozenset({16, 17, 18, 23, 24, 25, 28, 29, 30}),
    "pixels": frozenset({16, 17, 18}),
    "short_strips": frozenset({23, 24, 25}),
    "long_strips": frozenset({28, 29, 30}),
}
_ODD_REGION_VOLUME_IDS: Dict[str, FrozenSet[int]] = {
    "n_endcap": frozenset({16, 23, 28}),
    "p_endcap": frozenset({18, 25, 30}),
    "barrel": frozenset({17, 24, 29}),
}
_ODD_IDENTIFICATION_MAP: Dict[str, Any] = {
    "subsystems": _ODD_SUBSYSTEM_VOLUME_IDS,
    "regions": _ODD_REGION_VOLUME_IDS,
}

# GenericDetector uses the same subsystem naming convention:
# pixels, short_strips, long_strips with n_endcap/barrel/p_endcap partitioning.
_GENERIC_SUBSYSTEM_VOLUME_IDS: Dict[str, FrozenSet[int]] = {
    "all": frozenset({7, 8, 9, 12, 13, 14, 16, 17, 18}),
    "pixels": frozenset({7, 8, 9}),
    "short_strips": frozenset({12, 13, 14}),
    "long_strips": frozenset({16, 17, 18}),
}
_GENERIC_REGION_VOLUME_IDS: Dict[str, FrozenSet[int]] = {
    "n_endcap": frozenset({7, 12, 16}),
    "p_endcap": frozenset({9, 14, 18}),
    "barrel": frozenset({8, 13, 17}),
}
_GENERIC_IDENTIFICATION_MAP: Dict[str, Any] = {
    "subsystems": _GENERIC_SUBSYSTEM_VOLUME_IDS,
    "regions": _GENERIC_REGION_VOLUME_IDS,
}

_IDENTIFICATION_MAPS: Dict[str, Dict[str, Any]] = {
    "odd": _ODD_IDENTIFICATION_MAP,
    "generic": _GENERIC_IDENTIFICATION_MAP,
}


def _identification_map(detector_name: str) -> Dict[str, Any]:
    return _IDENTIFICATION_MAPS[detector_name]


def _parse_sub_system_selector(
    spec: str, identification_map: Dict[str, Any]
) -> Tuple[str, str]:
    """Parse CLI ``--sub-system`` selector against chosen detector map."""

    subsystem, region = parse_subsystem_region_spec(spec)
    aliases = {"short_strip": "short_strips", "long_strip": "long_strips"}
    subsystem = aliases.get(subsystem, subsystem)

    subsystems = identification_map["subsystems"]
    if subsystem not in subsystems:
        raise ValueError(
            f"Unknown sub-system {subsystem!r}; expected one of "
            f"{sorted(subsystems.keys())}"
        )

    allowed_regions = ("all", "barrel", "endcaps", "n_endcap", "p_endcap")
    if region not in allowed_regions:
        raise ValueError(
            f"Unknown region {region!r}; expected one of {allowed_regions}"
        )
    return subsystem, region


def _volume_ids(detector_name: str, subsystem: str, region: str) -> FrozenSet[int]:
    return frozenset(
        volume_ids_for_selection(subsystem, region, _identification_map(detector_name))
    )


def _build_detector(detector_name: str, *, gen3: bool):
    if detector_name == "odd":
        from acts.examples.odd import getOpenDataDetector

        return getOpenDataDetector(gen3=gen3)
    if detector_name == "generic":
        return acts.examples.GenericDetector()
    raise ValueError(f"Unsupported detector {detector_name!r}")


def _parse_range(
    view: str, raw: list[float] | None
) -> tuple[float, float, float, float] | None:
    if raw is None:
        return None
    if len(raw) != 4:
        raise SystemExit(f"--range expects 4 values for view {view}; got {len(raw)}")
    return float(raw[0]), float(raw[1]), float(raw[2]), float(raw[3])


def _gather_patches(
    view,
    surfaces,
    gctx,
    quarter_segments: int,
    *,
    triangular_mesh: bool,
):
    polys = []
    failed = 0
    for srf in surfaces:
        try:
            ph = srf.polyhedronRepresentation(gctx, quarter_segments)
            polys.extend(
                surface_to_2d_patches(view, ph, prefer_triangles=triangular_mesh)
            )
        except Exception:
            failed += 1
    if failed:
        warnings.warn(f"Skipped {failed} surfaces (polyhedron failed)")
    return polys


def _clip_polygon_horizontal(poly, y_cut: float, keep_below: bool):
    """Clip polygon against y <= y_cut or y >= y_cut (Sutherland-Hodgman)."""

    if len(poly) < 3:
        return np.zeros((0, 2))

    def inside(p):
        return p[1] <= y_cut if keep_below else p[1] >= y_cut

    def intersect(p1, p2):
        dy = p2[1] - p1[1]
        if abs(dy) < 1e-15:
            return np.array([p2[0], y_cut], dtype=float)
        t = (y_cut - p1[1]) / dy
        return np.array([p1[0] + t * (p2[0] - p1[0]), y_cut], dtype=float)

    output = []
    prev = poly[-1]
    prev_in = inside(prev)
    for curr in poly:
        curr_in = inside(curr)
        if curr_in:
            if not prev_in:
                output.append(intersect(prev, curr))
            output.append(curr)
        elif prev_in:
            output.append(intersect(prev, curr))
        prev = curr
        prev_in = curr_in

    if len(output) < 3:
        return np.zeros((0, 2))
    return np.asarray(output, dtype=float)


def _split_zphi_wrapped_polygon(poly):
    """Split one zphi polygon across ±pi wrap if needed."""

    if len(poly) < 3:
        return [poly], False, None

    phi = poly[:, 1]
    if (np.max(phi) - np.min(phi)) <= math.pi:
        return [poly], False, None

    unwrapped = np.array(poly, copy=True)
    unwrapped[:, 1] = np.where(
        unwrapped[:, 1] < 0.0,
        unwrapped[:, 1] + 2.0 * math.pi,
        unwrapped[:, 1],
    )

    lower = _clip_polygon_horizontal(unwrapped, math.pi, keep_below=True)
    upper = _clip_polygon_horizontal(unwrapped, math.pi, keep_below=False)

    split = []
    if len(lower) >= 3:
        split.append(lower)
    if len(upper) >= 3:
        upper[:, 1] -= 2.0 * math.pi
        split.append(upper)

    if len(split) == 0:
        return [poly], False, None

    z_center = float(np.mean(poly[:, 0]))
    z_span = float(np.max(poly[:, 0]) - np.min(poly[:, 0]))
    return split, True, (z_center, z_span)


def _split_zphi_wrapped_polygons(polygons):
    """Split all wrap-crossing polygons and return additional join indicators."""

    split_polys = []
    join_markers = []
    for poly in polygons:
        pieces, wrapped, marker_info = _split_zphi_wrapped_polygon(poly)
        split_polys.extend(pieces)
        if wrapped and marker_info is not None:
            z_center, z_span = marker_info
            z_half = max(2.0, 0.06 * z_span)
            phi_delta = 0.08
            join_markers.append(
                np.array(
                    [
                        [z_center - z_half, math.pi],
                        [z_center + z_half, math.pi],
                        [z_center, math.pi - phi_delta],
                    ],
                    dtype=float,
                )
            )
            join_markers.append(
                np.array(
                    [
                        [z_center - z_half, -math.pi],
                        [z_center + z_half, -math.pi],
                        [z_center, -math.pi + phi_delta],
                    ],
                    dtype=float,
                )
            )
    return split_polys, join_markers


def _dedupe_surfaces_by_geometry_id(surfaces):
    unique = {}
    for srf in surfaces:
        unique[int(srf.geometryId.value)] = srf
    return list(unique.values())


def _select_surface_category(tracking_geometry, category: str):
    sensitive = collect_sensitive_surfaces(tracking_geometry)
    material = collect_material_surfaces(tracking_geometry)
    if category == "sensitive":
        return sensitive
    if category == "material":
        return material
    if category == "all":
        return _dedupe_surfaces_by_geometry_id([*sensitive, *material])
    raise ValueError(f"Unknown surface category: {category}")


def _layers_in_volumes(grouped, volume_ids):
    lids = set()
    for vid in volume_ids:
        layers = grouped.get(vid, {})
        lids.update(layers.keys())
    return sorted(lids)


def _format_layer_list(layer_ids):
    if not layer_ids:
        return "-"
    return ":".join(str(v) for v in layer_ids)


def _print_detector_inspect(grouped_all, select_category: str, detector_name: str):
    identification_map = _identification_map(detector_name)
    subsystems = identification_map["subsystems"]

    print(
        f"Detector inspect [{detector_name}] " f"(surface category: {select_category})"
    )
    print("=" * 72)

    region_order = ("all", "barrel", "n_endcap", "p_endcap", "endcaps")
    subsystem_names = [name for name in subsystems.keys() if name != "all"]

    for subsystem in subsystem_names:
        sub_vols = sorted(subsystems[subsystem])
        sub_layers = _layers_in_volumes(grouped_all, sub_vols)
        print(f"{subsystem}")
        print(f"  volumes : {','.join(str(v) for v in sub_vols)}")
        print(f"  layers  : {_format_layer_list(sub_layers)}")
        print("  parts:")
        for region in region_order:
            vol_ids = sorted(_volume_ids(detector_name, subsystem, region))
            lay_ids = _layers_in_volumes(grouped_all, vol_ids)
            print(
                f"    - {region:<8} volumes={','.join(str(v) for v in vol_ids) or '-':<12} "
                f"layers={_format_layer_list(lay_ids)}"
            )
        print("")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--detector",
        choices=("odd", "generic"),
        default="odd",
        help="Detector model to build",
    )
    parser.add_argument(
        "--sub-system",
        default="pixels::all",
        help=(
            "subsystem::region with subsystem in "
            "{all, pixels, short_strips, long_strips} and region in "
            "{all, barrel, endcaps, n_endcap, p_endcap}"
        ),
    )
    parser.add_argument(
        "--view",
        choices=("xy", "zr", "zphi"),
        default="xy",
        help="2D projection plane",
    )
    parser.add_argument(
        "--select",
        choices=("all", "material", "sensitive"),
        default="all",
        help=(
            "Surface category: sensitive only, material only, or union of both "
            "(portal/boundary surfaces not handled yet)"
        ),
    )
    parser.add_argument(
        "--layer",
        default="all",
        help="Layer selection: 'all' or colon/comma separated list, e.g. 2:4:6",
    )
    parser.add_argument(
        "--range",
        nargs=4,
        type=float,
        metavar=("LOW_0", "LOW_1", "HIGH_0", "HIGH_1"),
        help=(
            "Axis limits: for xy → x_min y_min x_max y_max; "
            "for zr → z_min r_min z_max r_max; "
            "for zphi → z_min phi_min z_max phi_max (radians)"
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show interactive window (matplotlib)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Save figure to this path (png/pdf/svg…)",
    )
    parser.add_argument(
        "--segments",
        type=int,
        default=8,
        help="Surface.polyhedronRepresentation quarterSegments",
    )
    parser.add_argument(
        "--triangular-mesh",
        action="store_true",
        help=(
            "Use polyhedron triangularMesh (tessellated triangles); "
            "default uses native polygon faces"
        ),
    )
    parser.add_argument(
        "--gen3",
        action="store_true",
        help="Use OpenDataDetector gen3 constructor (only meaningful for --detector odd)",
    )
    parser.add_argument(
        "--detector-inspect",
        action="store_true",
        dest="detector_inspect",
        help=(
            "Print detector/sub-part/layer overview for the selected surface "
            "category and exit"
        ),
    )
    parser.add_argument(
        "--list-hierarchy",
        action="store_true",
        help="Print volume/layer surface counts and exit",
    )
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
    except ImportError as e:
        raise SystemExit(
            "matplotlib is required for display_tracking_geometry.py "
            "(install e.g. pip install matplotlib)"
        ) from e

    if args.detector != "odd" and args.gen3:
        warnings.warn("--gen3 ignored for detector != odd")

    try:
        subsystem, region = _parse_sub_system_selector(
            args.sub_system, _identification_map(args.detector)
        )
    except ValueError as e:
        raise SystemExit(str(e)) from e
    try:
        layer_ids = parse_layer_selection_spec(args.layer)
    except ValueError as e:
        raise SystemExit(str(e)) from e

    limits = _parse_range(args.view, args.range)
    vol_ids = _volume_ids(args.detector, subsystem, region)

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    detector_obj = _build_detector(args.detector, gen3=args.gen3)
    detector_ctx = (
        detector_obj
        if hasattr(detector_obj, "__enter__")
        else nullcontext(detector_obj)
    )

    with detector_ctx as detector:
        tracking_geometry = detector.trackingGeometry()

        surfaces = _select_surface_category(tracking_geometry, args.select)
        grouped_all = group_surfaces_by_volume_layer(surfaces)
        if args.detector_inspect:
            _print_detector_inspect(grouped_all, args.select, args.detector)
            return

        selected = filter_surfaces_by_volume_ids(surfaces, vol_ids)
        selected = filter_surfaces_by_layer_ids(selected, layer_ids)
        selected = sort_surfaces_by_geometry_id(selected)

        if args.list_hierarchy:
            grouped = group_surfaces_by_volume_layer(selected)
            print(
                f"Geometry hierarchy for {args.sub_system!r} "
                f"({len(selected)} surfaces):"
            )
            for vol, layers in grouped.items():
                nlay = len(layers)
                ns = sum(len(v) for v in layers.values())
                print(f"  volume {vol}: {nlay} layers, {ns} surfaces")
            return

        if not selected:
            print(
                f"No surfaces for {args.sub_system!r} "
                f"(check subsystem and region).",
                file=sys.stderr,
            )
            return

        polys = _gather_patches(
            args.view,
            selected,
            gctx,
            max(1, args.segments),
            triangular_mesh=args.triangular_mesh,
        )
        if not polys:
            print("No drawable patches produced.", file=sys.stderr)
            return

        join_markers = []
        if args.view == "zphi":
            polys, join_markers = _split_zphi_wrapped_polygons(polys)

        if limits is None:
            xmin, ymin, xmax, ymax = auto_axis_limits(polys)
            if args.view == "zphi":
                ymin, ymax = -math.pi, math.pi
        else:
            xmin, ymin, xmax, ymax = limits

        fig, ax = plt.subplots(figsize=(9, 8))
        collection = PolyCollection(
            polys,
            facecolors=(0.2, 0.5, 0.95, 0.35),
            edgecolors=(0.1, 0.2, 0.4, 0.8),
            linewidths=0.2,
        )
        ax.add_collection(collection)
        if join_markers:
            marker_collection = PolyCollection(
                join_markers,
                facecolors="none",
                edgecolors=(0.1, 0.2, 0.4, 0.9),
                linewidths=0.8,
            )
            ax.add_collection(marker_collection)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        if args.view == "zphi":
            ax.set_aspect("auto")
            if limits is None:
                ax.set_yticks([-math.pi, 0.0, math.pi])
                ax.set_yticklabels(["-π", "0", "+π"])
        else:
            ax.set_aspect("equal", adjustable="box")
        labels = {
            "xy": ("x [mm]", "y [mm]"),
            "zr": ("z [mm]", "r [mm]"),
            "zphi": ("z [mm]", "φ [rad]"),
        }
        ax.set_xlabel(labels[args.view][0])
        ax.set_ylabel(labels[args.view][1])
        ax.set_title(f"{args.detector.upper()} {args.sub_system} — {args.view}")
        fig.tight_layout()

        if args.output is not None:
            fig.savefig(args.output, dpi=150)
        if args.show:
            plt.show()
        elif args.output is None:
            plt.close(fig)
            print(
                "Figure not shown (--show not set) and not saved (--output). "
                "Use --show and/or --output.",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
