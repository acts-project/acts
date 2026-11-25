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
    config.replace_files = True
    config.compactify_json = True
    config.write_material = False
    config.write_grids = False

    assert config.path == "./json_writer_test/"
    assert config.replace_files
    assert config.compactify_json
    assert not config.write_material
    assert not config.write_grids


def test_write_geometry_only(tmp_path):
    """Both writer flags reach the C++ writer and suppress their output."""
    out_dir = str(tmp_path) + os.sep

    det, names = detray.tests.build_toy_detector(
        detray.core.HostMemoryResource(), detray.tests.ToyDetectorConfig()
    )

    config = detray.io.DetectorWriterConfig()
    config.path = out_dir

    config.write_material = False
    config.write_grids = False

    detray.io.write_detector_json(det, names, config)

    assert os.listdir(out_dir) == ["toy_detector_geometry.json"]
