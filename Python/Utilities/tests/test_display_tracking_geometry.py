from __future__ import annotations

import importlib.util
from pathlib import Path

import acts
import acts.examples
import numpy as np
import pytest


def _load_display_module():
    repo_root = Path(__file__).resolve().parents[3]
    script = (
        repo_root / "Examples" / "Scripts" / "Python" / "display_tracking_geometry.py"
    )
    spec = importlib.util.spec_from_file_location("display_tracking_geometry", script)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def display_module():
    return _load_display_module()


@pytest.fixture(scope="module")
def generic_detector():
    return acts.examples.GenericDetector()


@pytest.fixture(scope="module")
def generic_tracking_geometry(generic_detector):
    return generic_detector.trackingGeometry()


def test_generic_subsystem_selector_and_volume_map(display_module):
    subsystem, region = display_module._parse_sub_system_selector(
        "pixels::barrel", display_module._GENERIC_IDENTIFICATION_MAP
    )
    assert subsystem == "pixels"
    assert region == "barrel"

    subsystem2, region2 = display_module._parse_sub_system_selector(
        "short_strip::endcaps", display_module._GENERIC_IDENTIFICATION_MAP
    )
    assert subsystem2 == "short_strips"
    assert region2 == "endcaps"

    assert display_module._volume_ids("generic", "pixels", "barrel") == {8}
    assert display_module._volume_ids("generic", "pixels", "n_endcap") == {7}
    assert display_module._volume_ids("generic", "pixels", "p_endcap") == {9}


def test_generic_sensitive_surface_selection(display_module, generic_tracking_geometry):
    surfaces = display_module._select_surface_category(
        generic_tracking_geometry, "sensitive"
    )
    assert len(surfaces) > 0

    pixel_barrel = display_module.filter_surfaces_by_volume_ids(surfaces, {8})
    assert len(pixel_barrel) > 0
    assert all(s.geometryId.volume == 8 for s in pixel_barrel)

    layers = display_module.parse_layer_selection_spec("2:4")
    layer_filtered = display_module.filter_surfaces_by_layer_ids(pixel_barrel, layers)
    assert len(layer_filtered) > 0
    assert all(s.geometryId.layer in {2, 4} for s in layer_filtered)


@pytest.mark.parametrize("view", ["xy", "zr", "zphi"])
def test_gather_patches_for_generic_views(
    view, display_module, generic_tracking_geometry
):
    surfaces = display_module._select_surface_category(
        generic_tracking_geometry, "sensitive"
    )
    selected = display_module.filter_surfaces_by_volume_ids(surfaces, {8})
    selected = display_module.sort_surfaces_by_geometry_id(selected)[:30]

    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    patches = display_module._gather_patches(
        view, selected, gctx, 4, triangular_mesh=False
    )
    assert len(patches) > 0
    assert all(len(poly) >= 3 for poly in patches)


def test_zphi_wrap_polygon_splitting(display_module):
    poly = np.array(
        [
            [0.0, -3.12],
            [10.0, -3.02],
            [10.0, 3.03],
            [0.0, 3.11],
        ],
        dtype=float,
    )
    pieces, wrapped, marker_info = display_module._split_zphi_wrapped_polygon(poly)
    assert wrapped is True
    assert marker_info is not None
    assert len(pieces) == 2

    split, markers = display_module._split_zphi_wrapped_polygons([poly])
    assert len(split) == 2
    assert len(markers) == 2


def test_detector_inspect_output_for_generic(
    capsys, display_module, generic_tracking_geometry
):
    surfaces = display_module._select_surface_category(
        generic_tracking_geometry, "sensitive"
    )
    grouped = display_module.group_surfaces_by_volume_layer(surfaces)
    display_module._print_detector_inspect(grouped, "sensitive", "generic")
    out = capsys.readouterr().out
    assert "Detector inspect [generic]" in out
    assert "pixels" in out
    assert "short_strips" in out
    assert "long_strips" in out
