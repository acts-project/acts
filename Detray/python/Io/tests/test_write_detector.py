# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os

import detray.core
import detray.io
import detray.tests


def test_detector_writer_config():
    config = detray.io.DetectorWriterConfig()
    config.path = "./json_writer_test/"
    config.replaceFiles = True
    config.compactifyJson = True
    config.writeMaterial = False
    config.writeGrids = False

    assert config.path == "./json_writer_test/"
    assert config.replaceFiles
    assert config.compactifyJson
    assert not config.writeMaterial
    assert not config.writeGrids


def test_write_geometry_only(tmp_path):
    """Both writer flags reach the C++ writer and suppress their output."""
    out_dir = str(tmp_path) + os.sep

    det, names = detray.tests.buildToyDetector(
        detray.core.HostMemoryResource(), detray.tests.ToyDetectorConfig()
    )

    config = detray.io.DetectorWriterConfig()
    config.path = out_dir

    config.writeMaterial = False
    config.writeGrids = False

    detray.io.writeDetector(det, names, config)

    assert os.listdir(out_dir) == ["toy_detector_geometry.json"]
