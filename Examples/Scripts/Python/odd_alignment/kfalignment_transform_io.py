"""Alignment transform I/O: six-column delta text, index maps, GeoID maps."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

import acts
import acts.examples.alignment

# Local DoF names / size (Acts::eAlignmentSize).
ALIGNMENT_DOF_NAMES = ("dx", "dy", "dz", "rx", "ry", "rz")
ALIGNMENT_SIZE = 6


def parse_geoid(id_str: str) -> dict:
    """Parse ``vol=X|lay=Y|sen=Z[|ext=W]`` into a dict."""
    parts = {}
    for part in id_str.split("|"):
        key, value = part.split("=")
        parts[key] = int(value)
    return parts


def geoid_to_id_string(gid: acts.GeometryIdentifier) -> str:
    """Format GeometryIdentifier as ``vol=X|lay=Y|sen=Z[|ext=W]`` (omit ext=0)."""
    parts = [f"vol={gid.volume}", f"lay={gid.layer}", f"sen={gid.sensitive}"]
    if gid.extra != 0:
        parts.append(f"ext={gid.extra}")
    return "|".join(parts)


def copy_transform3(trf: acts.Transform3) -> acts.Transform3:
    """Return a mutable copy of an Acts Transform3."""
    t = trf.translation
    r = trf.rotation
    return acts.Transform3(
        acts.Vector3(float(t[0]), float(t[1]), float(t[2])),
        acts.RotationMatrix3(
            acts.Vector3(r[0, 0], r[1, 0], r[2, 0]),
            acts.Vector3(r[0, 1], r[1, 1], r[2, 1]),
            acts.Vector3(r[0, 2], r[1, 2], r[2, 2]),
        ),
    )


def sort_surfaces_alignment_order(matches: list) -> list:
    """Sort ``(surface, placement)`` by layer, then volume / sensitive / extra."""
    return sorted(
        matches,
        key=lambda sp: (
            sp[0].geometryId.layer,
            sp[0].geometryId.volume,
            sp[0].geometryId.sensitive,
            sp[0].geometryId.extra,
        ),
    )


def build_geoid_surface_index_map(matches: list) -> list:
    """Build dense surfaceIndex → GeoID records from sorted matches."""
    index_map = []
    for surface_index, (surface, _placement) in enumerate(matches):
        gid = surface.geometryId
        index_map.append(
            {
                "surface_index": surface_index,
                "geo_id": geoid_to_id_string(gid),
                "gid": gid,
            }
        )
    return index_map


def write_misalignment_index_map(path: Path, index_map: list) -> None:
    """Write surfaceIndex ↔ GeoID sidecar text (label, index, dof, name, geo_id)."""
    with open(path, "w") as f:
        f.write("# label  surface_index  dof  dof_name  geo_id\n")
        for entry in index_map:
            sidx = entry["surface_index"]
            geo_id = entry["geo_id"]
            for dof in range(ALIGNMENT_SIZE):
                label = sidx * ALIGNMENT_SIZE + dof + 1
                f.write(
                    f"{label:<8d} {sidx:<14d} {dof:<4d} "
                    f"{ALIGNMENT_DOF_NAMES[dof]:<8s} {geo_id}\n"
                )


def write_alignment_delta_res(
    path: Path,
    surface_deltas: list,
) -> None:
    """Write per-surface DoF deltas as alignment result text.

    Columns: ``label  value  0.  difference  error  99`` with
    ``label = surfaceIndex * 6 + dof + 1`` (1-based). For truth injection,
    ``value == difference`` and ``error`` is ``0.``.
    """
    with open(path, "w") as f:
        for surface_index, deltas in enumerate(surface_deltas):
            for dof, value in enumerate(deltas):
                label = surface_index * ALIGNMENT_SIZE + dof + 1
                f.write(
                    f"{label:8d}   {value:12g}   {0.0:4g}   "
                    f"{value:12g}   {0.0:12g}    99 \n"
                )


def geometry_identifier_from_id_string(id_str: str) -> acts.GeometryIdentifier:
    """Build GeometryIdentifier from ``vol=X|lay=Y|sen=Z[|ext=W]``."""
    parts = parse_geoid(id_str)
    gid_kwargs = {
        "volume": parts.get("vol", 0),
        "layer": parts.get("lay", 0),
        "sensitive": parts.get("sen", 0),
    }
    if "ext" in parts:
        gid_kwargs["extra"] = parts["ext"]
    return acts.GeometryIdentifier(**gid_kwargs)


def read_misalignment_index_map(path: Path) -> dict:
    """Read ``{surface_index: geo_id_string}`` from an index-map text file."""
    surface_to_geoid = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            surface_index = int(parts[1])
            geo_id = parts[4]
            surface_to_geoid[surface_index] = geo_id
    return surface_to_geoid


def read_alignment_delta_res(path: Path) -> dict:
    """Read ``{surface_index: [dx, dy, dz, rx, ry, rz]}`` from result text (value col)."""
    surface_deltas = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            label = int(parts[0])
            value = float(parts[1])
            surface_index = (label - 1) // ALIGNMENT_SIZE
            dof = (label - 1) % ALIGNMENT_SIZE
            if surface_index not in surface_deltas:
                surface_deltas[surface_index] = [0.0] * ALIGNMENT_SIZE
            surface_deltas[surface_index][dof] = value
    return surface_deltas


def apply_local_misalignment_deltas(
    trf: acts.Transform3, deltas: list
) -> acts.Transform3:
    """Apply ``[dx, dy, dz, rx, ry, rz]`` in place (local frame; ``T *= Rz*Ry*Rx``)."""
    dx, dy, dz, rotx, roty, rotz = deltas

    if dx != 0.0:
        gen = acts.examples.alignment.AlignmentGeneratorLocalShift()
        gen.axisDirection = acts.AxisDirection.AxisX
        gen.shift = dx
        gen(trf)
    if dy != 0.0:
        gen = acts.examples.alignment.AlignmentGeneratorLocalShift()
        gen.axisDirection = acts.AxisDirection.AxisY
        gen.shift = dy
        gen(trf)
    if dz != 0.0:
        gen = acts.examples.alignment.AlignmentGeneratorLocalShift()
        gen.axisDirection = acts.AxisDirection.AxisZ
        gen.shift = dz
        gen(trf)

    if rotz != 0.0:
        gen = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        gen.axis = acts.Vector3(0.0, 0.0, 1.0)
        gen.angle = rotz
        gen(trf)
    if roty != 0.0:
        gen = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        gen.axis = acts.Vector3(0.0, 1.0, 0.0)
        gen.angle = roty
        gen(trf)
    if rotx != 0.0:
        gen = acts.examples.alignment.AlignmentGeneratorLocalRotation()
        gen.axis = acts.Vector3(1.0, 0.0, 0.0)
        gen.angle = rotx
        gen(trf)

    return trf


def sensitive_surfaces_by_id_string(
    trackingGeometry: acts.TrackingGeometry,
) -> dict:
    """Map ``geoid_to_id_string(gid)`` -> sensitive ``Surface``."""
    id_to_surface = {}

    def visit(surface: acts.Surface) -> bool:
        if surface.isSensitive:
            id_to_surface[geoid_to_id_string(surface.geometryId)] = surface
        return True

    trackingGeometry.visitSurfaces(visit)
    return id_to_surface


def geo_id_map_from_delta_maps(
    trackingGeometry: acts.TrackingGeometry,
    surface_to_geoid: dict,
    surface_deltas: dict,
    id_to_surface: Optional[dict] = None,
) -> dict:
    """Build GeoID -> Transform3 from already-read alignment index/delta maps."""
    if id_to_surface is None:
        id_to_surface = sensitive_surfaces_by_id_string(trackingGeometry)

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    geoIdMap = {}

    for surface_index, geo_id_str in sorted(surface_to_geoid.items()):
        if surface_index not in surface_deltas:
            continue
        surface = id_to_surface.get(geo_id_str)
        if surface is None:
            print(
                f"Warning: surface for {geo_id_str} (index={surface_index}) "
                f"not found in trackingGeometry; skipping"
            )
            continue
        gid = surface.geometryId
        trf = copy_transform3(surface.localToGlobalTransform(gctx))
        apply_local_misalignment_deltas(trf, surface_deltas[surface_index])
        geoIdMap[gid] = trf

    return geoIdMap


def geo_id_map_from_misalignment(
    trackingGeometry: acts.TrackingGeometry,
    applied_path: Path,
    index_path: Path,
) -> dict:
    """Build GeoID → Transform3 from misalignment applied + index files."""
    surface_to_geoid = read_misalignment_index_map(index_path)
    surface_deltas = read_alignment_delta_res(applied_path)
    return geo_id_map_from_delta_maps(
        trackingGeometry, surface_to_geoid, surface_deltas
    )


def local_deltas_between_transforms(
    trf_from: acts.Transform3, trf_to: acts.Transform3
) -> list:
    """Local ``[dx, dy, dz, rx, ry, rz]`` taking ``trf_from`` to ``trf_to``."""
    R0 = np.array(
        [[float(trf_from.rotation[i, j]) for j in range(3)] for i in range(3)]
    )
    t0 = np.array([float(trf_from.translation[i]) for i in range(3)])
    R1 = np.array([[float(trf_to.rotation[i, j]) for j in range(3)] for i in range(3)])
    t1 = np.array([float(trf_to.translation[i]) for i in range(3)])

    delta_center_local = R0.T @ (t1 - t0)
    R_delta = R0.T @ R1  # Rz(rz) * Ry(ry) * Rx(rx)

    ry = float(np.arcsin(np.clip(-R_delta[2, 0], -1.0, 1.0)))
    cy = float(np.cos(ry))
    if abs(cy) > 1e-8:
        rx = float(np.arctan2(R_delta[2, 1], R_delta[2, 2]))
        rz = float(np.arctan2(R_delta[1, 0], R_delta[0, 0]))
    else:
        rx = float(np.arctan2(-R_delta[1, 2], R_delta[1, 1]))
        rz = 0.0

    return [
        float(delta_center_local[0]),
        float(delta_center_local[1]),
        float(delta_center_local[2]),
        rx,
        ry,
        rz,
    ]


def load_misalignment_for_alignment(
    trackingGeometry: acts.TrackingGeometry,
    applied_path: Path,
    index_path: Path,
) -> tuple:
    """Load misalignment for ``runAlignment``.

    Returns ``(geoIdMap, alignment_placements, surface_to_geoid)``.
    """
    surface_to_geoid = read_misalignment_index_map(index_path)
    surface_deltas = read_alignment_delta_res(applied_path)
    id_to_surface = sensitive_surfaces_by_id_string(trackingGeometry)
    geoIdMap = geo_id_map_from_delta_maps(
        trackingGeometry,
        surface_to_geoid,
        surface_deltas,
        id_to_surface=id_to_surface,
    )

    alignment_placements = []
    for surface_index in sorted(surface_to_geoid.keys()):
        geo_id_str = surface_to_geoid[surface_index]
        surface = id_to_surface.get(geo_id_str)
        if surface is None:
            print(
                f"Warning: no SurfacePlacement for {geo_id_str} "
                f"(index={surface_index}); skipping in alignedDetElements"
            )
            continue
        placement = acts.examples.alignment.surfacePlacement(surface)
        if placement is None:
            print(
                f"Warning: no SurfacePlacement for {geo_id_str} "
                f"(index={surface_index}); skipping in alignedDetElements"
            )
            continue
        alignment_placements.append(placement)

    return geoIdMap, alignment_placements, surface_to_geoid


def geo_id_map_keyed_by_id_string(geoIdMap: dict) -> dict:
    """Re-index a GeoID→Transform map by ``geoid_to_id_string`` (stable across pybind)."""
    return {geoid_to_id_string(gid): trf for gid, trf in geoIdMap.items()}


def rebase_geo_id_map_to_tracking_geometry(
    trackingGeometry: acts.TrackingGeometry, geoIdMap: dict
) -> dict:
    """Re-key transforms with live ``surface.geometryId`` from ``trackingGeometry``.

    Needed when map keys come from another construction (pybind GeoIDs may not
    compare equal). Maps built via ``geo_id_map_from_*`` already use live keys.
    """
    by_id = geo_id_map_keyed_by_id_string(geoIdMap)
    out = {}

    def visit(surface: acts.Surface) -> bool:
        if not surface.isSensitive:
            return True
        id_str = geoid_to_id_string(surface.geometryId)
        if id_str in by_id:
            out[surface.geometryId] = copy_transform3(by_id[id_str])
        return True

    trackingGeometry.visitSurfaces(visit)
    return out


def apply_delta_maps_to_geo_id_map(
    trackingGeometry: acts.TrackingGeometry,
    geoIdMap: dict,
    surface_to_geoid: dict,
    surface_deltas: dict,
    id_to_surface: Optional[dict] = None,
) -> dict:
    """Stack delta maps onto ``geoIdMap`` (missing surfaces start from nominal)."""
    if id_to_surface is None:
        id_to_surface = sensitive_surfaces_by_id_string(trackingGeometry)
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

    existing_by_id = geo_id_map_keyed_by_id_string(geoIdMap)
    out = {}
    for surface_index, geo_id_str in sorted(surface_to_geoid.items()):
        if surface_index not in surface_deltas:
            continue
        surface = id_to_surface.get(geo_id_str)
        if surface is None:
            continue
        if geo_id_str in existing_by_id:
            trf = copy_transform3(existing_by_id[geo_id_str])
        else:
            trf = copy_transform3(surface.localToGlobalTransform(gctx))
        apply_local_misalignment_deltas(trf, surface_deltas[surface_index])
        out[surface.geometryId] = trf
    return out


def apply_deltas_to_geo_id_map(
    trackingGeometry: acts.TrackingGeometry,
    geoIdMap: dict,
    applied_path: Path,
    index_path: Path,
) -> dict:
    """Apply local deltas from applied + index files onto an existing GeoID map."""
    return apply_delta_maps_to_geo_id_map(
        trackingGeometry,
        geoIdMap,
        read_misalignment_index_map(index_path),
        read_alignment_delta_res(applied_path),
    )


def geo_id_map_from_aligned(
    trackingGeometry: acts.TrackingGeometry,
    misalignment_applied_path: Path,
    misalignment_index_path: Path,
    aligned_result_path: Path,
    aligned_index_path: Path = None,
) -> dict:
    """Build GeoID map: nominal → misalignment deltas → alignment deltas (one visit)."""
    mis_to_geoid = read_misalignment_index_map(misalignment_index_path)
    mis_deltas = read_alignment_delta_res(misalignment_applied_path)

    ali_index_path = (
        aligned_index_path
        if aligned_index_path is not None
        else misalignment_index_path
    )
    if Path(ali_index_path).resolve() == Path(misalignment_index_path).resolve():
        ali_to_geoid = mis_to_geoid
    else:
        ali_to_geoid = read_misalignment_index_map(ali_index_path)
    ali_deltas = read_alignment_delta_res(aligned_result_path)

    id_to_surface = sensitive_surfaces_by_id_string(trackingGeometry)
    geoIdMap = geo_id_map_from_delta_maps(
        trackingGeometry,
        mis_to_geoid,
        mis_deltas,
        id_to_surface=id_to_surface,
    )
    return apply_delta_maps_to_geo_id_map(
        trackingGeometry,
        geoIdMap,
        ali_to_geoid,
        ali_deltas,
        id_to_surface=id_to_surface,
    )
