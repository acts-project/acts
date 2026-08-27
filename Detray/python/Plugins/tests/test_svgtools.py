# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import gc
import os
import weakref
from xml.etree import ElementTree

import pytest

import detray.core
import detray.io
import detray.svgtools

# Every draw method is bound for every view, so covering all of the views in a
# single draw method is enough.
VIEWS = [
    (detray.svgtools.ViewXY, "xy"),
    (detray.svgtools.ViewZR, "zr"),
    (detray.svgtools.ViewZPhi, "zphi"),
    (detray.svgtools.ViewZRPhi, "zrphi"),
]


@pytest.fixture
def illustrator(detector):
    det, names = detector
    return detray.svgtools.Illustrator(
        det, names, detray.svgtools.tableauColorblindStyle()
    )


def test_write_svg(tmp_path):
    axes = detray.svgtools.drawAxes("axes", [-1.0, 1.0], [-1.0, 1.0], "x", "y", 12)
    assert axes.id == "axes"

    path = os.path.join(str(tmp_path), axes.id)
    detray.svgtools.writeSvg(path, [axes])

    ElementTree.parse(path + ".svg")


def test_style_properties():
    style = detray.svgtools.tableauColorblindStyle()
    style.fontSize = 20
    style.materialGradientPos = [800.0, -200.0]

    assert style.fontSize == 20
    assert style.materialGradientPos == pytest.approx([800.0, -200.0], rel=1e-6)


@pytest.mark.parametrize("view, suffix", VIEWS)
def test_draw_detector(illustrator, view, suffix):
    svg = illustrator.drawDetector(view(), detray.core.GeometryContext())

    assert svg.id == f"toy_detector_{suffix}"


def test_draw_volumes_by_index(illustrator):
    svg, grids = illustrator.drawVolumes(
        [0, 1], detray.svgtools.ViewXY(), detray.core.GeometryContext()
    )

    assert svg.id == "toy_detector_volumes_xy"
    # One entry per volume, empty for a volume that has no grid in this view.
    assert len(grids) == 2


def test_draw_volumes_by_name(illustrator, detector):
    _, names = detector
    svg, _ = illustrator.drawVolumes(
        [names[0], names[1]], detray.svgtools.ViewXY(), detray.core.GeometryContext()
    )

    assert svg.id == "toy_detector_volumes_xy"


def test_draw_surfaces(illustrator):
    svg, material = illustrator.drawSurfaces(
        [0, 1], detray.svgtools.ViewXY(), detray.core.GeometryContext()
    )

    assert svg.id == "toy_detector_surfaces_xy"
    assert material.id == "toy_detector_material_xy"


def test_illustrator_keeps_detector_and_names_alive(toy_detector_files):
    reader_config = detray.io.DetectorReaderConfig().addFile(
        toy_detector_files["geometry"]
    )
    memory_resource = detray.core.HostMemoryResource()
    det, names = detray.io.readDetector(memory_resource, reader_config)

    mr_ref = weakref.ref(memory_resource)
    det_ref = weakref.ref(det)
    names_ref = weakref.ref(names)

    il = detray.svgtools.Illustrator(
        det, names, detray.svgtools.tableauColorblindStyle()
    )

    del memory_resource, det, names
    gc.collect()

    assert mr_ref() is not None
    assert det_ref() is not None
    assert names_ref() is not None
    assert il.detName == "toy_detector"

    del il
    gc.collect()

    assert mr_ref() is None
    assert det_ref() is None
    assert names_ref() is None
