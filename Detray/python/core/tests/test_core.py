# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import gc
import weakref

import pytest

import detray.core
import detray.tests


def test_detector_keeps_memory_resource_alive():
    memory_resource = detray.core.HostMemoryResource()
    mr_ref = weakref.ref(memory_resource)

    det, names = detray.tests.build_toy_detector(
        memory_resource, detray.tests.ToyDetectorConfig()
    )

    del memory_resource
    gc.collect()

    assert mr_ref() is not None

    del det, names
    gc.collect()

    assert mr_ref() is None


def test_detector_outlives_memory_resource(toy_detector_files):
    reader_config = detray.io.DetectorReaderConfig()
    reader_config.add_files(
        toy_detector_files["geometry"], toy_detector_files["material"]
    )
    memory_resource = detray.core.HostMemoryResource()
    det, names = detray.io.read_detector_json(memory_resource, reader_config)

    del memory_resource
    gc.collect()

    assert len(det.volumes) == 22
    assert det.name(names) == "toy_detector"


def test_detector_volume_container_bounds(toy_detector):
    det, _ = toy_detector
    n = len(det.volumes)
    assert n == 22

    det.volumes[0]
    det.volumes[-1]  # Wraps to last element.

    with pytest.raises(IndexError):
        det.volumes[n]
    with pytest.raises(IndexError):
        det.volumes[-(n + 1)]


def test_name_map():
    names = detray.core.NameMap()
    names[0] = "volume_0"

    assert names[0] == "volume_0"
    assert names["volume_0"] == 0
    assert 0 in names
    assert "volume_0" in names
    with pytest.raises(KeyError):
        names[99]
    with pytest.raises(KeyError):
        names["does_not_exist"]
