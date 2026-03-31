"""Pure-Python helpers shipped with Acts (geometry views, utilities)."""

from .geometry_identifier import GeometryIdentifierMasks
from .projection2d import (
    auto_axis_limits,
    project_polyhedron_vertices,
    surface_to_2d_patches,
)
from .tracking_surfaces import (
    collect_material_surfaces,
    collect_sensitive_surfaces,
    filter_surfaces_by_layer_ids,
    filter_surfaces_by_volume_ids,
    group_surfaces_by_volume_layer,
    parse_layer_selection_spec,
    parse_subsystem_region_spec,
    sort_surfaces_by_geometry_id,
    volume_ids_for_selection,
)

__all__ = [
    "GeometryIdentifierMasks",
    "auto_axis_limits",
    "collect_material_surfaces",
    "collect_sensitive_surfaces",
    "filter_surfaces_by_layer_ids",
    "filter_surfaces_by_volume_ids",
    "group_surfaces_by_volume_layer",
    "parse_layer_selection_spec",
    "parse_subsystem_region_spec",
    "project_polyhedron_vertices",
    "sort_surfaces_by_geometry_id",
    "surface_to_2d_patches",
    "volume_ids_for_selection",
]
