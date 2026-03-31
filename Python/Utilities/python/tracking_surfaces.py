"""Detector-agnostic helpers for surface selection and ``GeometryIdentifier`` grouping."""

from __future__ import annotations

from collections import defaultdict
from typing import Collection, DefaultDict, Dict, Iterable, List, Mapping, Set, Tuple


def parse_subsystem_region_spec(spec: str) -> Tuple[str, str]:
    """Parse ``subsystem::region`` (case-normalised to lower).

    Validates only syntax, not names — subsystem / region meaning comes from the
    identification map supplied to :func:`volume_ids_for_selection`.
    """

    parts = spec.split("::", 1)
    if len(parts) != 2:
        raise ValueError(f"Selector must be 'subsystem::region', got {spec!r}")
    subsystem, region = parts[0].strip(), parts[1].strip()
    if not subsystem or not region:
        raise ValueError(f"Empty subsystem or region in {spec!r}")
    return subsystem.lower(), region.lower()


def volume_ids_for_selection(
    subsystem: str,
    region: str,
    identification_map: Mapping[str, object],
) -> Set[int]:
    """Resolve allowed tracking volume ids for a subsystem/region.

    Parameters
    ----------
    subsystem
        Key into ``identification_map["subsystems"]`` (e.g. ``pixels``).
    region
        One of ``all``, ``endcaps``, or a key of ``identification_map["regions"]``.
    identification_map
        Mapping with:

        - ``"subsystems"``: ``dict[str, Collection[int]]`` — volume id sets per
          logical detector (plus e.g. ``all`` for a union you define).
        - ``"regions"``: ``dict[str, Collection[int]]`` — volume id sets for
          geometric regions (``barrel``, ``n_endcap``, ``p_endcap``, …).

    Returns
    -------
    set[int]
        Volume ids to keep; intersection of subsystem volumes with region
        volumes. For ``region == "all"``, returns the subsystem set only.
        For ``region == "endcaps"``, uses
        ``regions["n_endcap"] | regions["p_endcap"]`` intersected with the
        subsystem set (both keys must be present).
    """

    try:
        subsystems = identification_map["subsystems"]
        regions = identification_map["regions"]
    except KeyError as e:
        raise KeyError(
            'identification_map must have "subsystems" and "regions" mappings'
        ) from e

    if subsystem not in subsystems:
        raise KeyError(
            f"Unknown subsystem {subsystem!r}; known: {sorted(subsystems.keys())}"
        )

    base = set(subsystems[subsystem])
    region = region.lower()

    if region == "all":
        return set(base)

    if region == "endcaps":
        try:
            cap = set(regions["n_endcap"]) | set(regions["p_endcap"])
        except KeyError as e:
            raise KeyError(
                'region "endcaps" requires regions["n_endcap"] and '
                'regions["p_endcap"] in identification_map'
            ) from e
        return set(v for v in base if v in cap)

    if region not in regions:
        raise KeyError(
            f"Unknown region {region!r}; known: all, endcaps, "
            f'{", ".join(sorted(regions.keys()))}'
        )
    reg = set(regions[region])
    return set(v for v in base if v in reg)


def filter_surfaces_by_volume_ids(
    surfaces: Iterable, volume_ids: Collection[int]
) -> List:
    """Keep surfaces whose ``geometryId.volume`` is in ``volume_ids``."""

    vidset = set(volume_ids)
    return [s for s in surfaces if s.geometryId.volume in vidset]


def parse_layer_selection_spec(spec: str) -> Set[int] | None:
    """Parse layer selection string.

    ``"all"`` returns ``None``. Otherwise parse colon/comma-separated integers,
    e.g. ``"1:2:3"`` or ``"1,2,3"``.
    """

    normalized = spec.strip().lower()
    if normalized == "all":
        return None
    tokens = [
        t.strip() for t in normalized.replace(",", ":").split(":") if t.strip() != ""
    ]
    if not tokens:
        raise ValueError("Layer selection is empty; use 'all' or e.g. '1:2:3'")
    out: Set[int] = set()
    for token in tokens:
        try:
            value = int(token)
        except ValueError as e:
            raise ValueError(
                f"Invalid layer token {token!r}; use integers or 'all'"
            ) from e
        if value < 0:
            raise ValueError(f"Layer id must be non-negative, got {value}")
        out.add(value)
    return out


def filter_surfaces_by_layer_ids(
    surfaces: Iterable, layer_ids: Collection[int] | None
) -> List:
    """Keep surfaces whose ``geometryId.layer`` is in ``layer_ids``.

    If ``layer_ids`` is ``None``, input is returned as a list unchanged.
    """

    if layer_ids is None:
        return list(surfaces)
    lidset = set(layer_ids)
    return [s for s in surfaces if s.geometryId.layer in lidset]


def collect_sensitive_surfaces(tracking_geometry) -> List:
    """Collect sensitive surfaces (same as default ``TrackingGeometry.visitSurfaces``)."""

    surfs: List = []

    def visit(s):
        surfs.append(s)

    tracking_geometry.visitSurfaces(visit)
    return surfs


def collect_material_surfaces(tracking_geometry) -> List:
    """Collect all surfaces that carry material information."""

    return list(tracking_geometry.extractMaterialSurfaces())


def sort_surfaces_by_geometry_id(surfaces: Iterable) -> List:
    """Sort by encoded ``GeometryIdentifier`` value."""

    return sorted(surfaces, key=lambda s: s.geometryId.value)


def group_surfaces_by_volume_layer(
    surfaces: Iterable,
) -> Dict[int, Dict[int, List]]:
    """``result[volume][layer] -> [surfaces…]`` from each surface's ``geometryId``."""

    nested: DefaultDict[int, DefaultDict[int, List]] = defaultdict(
        lambda: defaultdict(list)
    )
    for srf in surfaces:
        gid = srf.geometryId
        nested[gid.volume][gid.layer].append(srf)
    return {
        v: {ly: surfs for ly, surfs in sorted(layers.items())}
        for v, layers in sorted(nested.items())
    }
