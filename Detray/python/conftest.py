# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os

import pytest

import detray.core
import detray.examples
import detray.io


@pytest.fixture(scope="session")
def toy_detector_files(tmp_path_factory):
    """Write a toy detector to disk."""
    out_dir = str(tmp_path_factory.mktemp("toy_detector")) + os.sep
    detray.examples.generateToyDetector(outputDir=out_dir)

    files = {
        "geometry": os.path.join(out_dir, "toy_detector_geometry.json"),
        "material": os.path.join(out_dir, "toy_detector_homogeneous_material.json"),
    }
    for path in files.values():
        assert os.path.exists(path)
    return files


@pytest.fixture(scope="session")
def detector(toy_detector_files):
    """Read the toy detector."""
    reader_config = (
        detray.io.DetectorReaderConfig()
        .addFile(toy_detector_files["geometry"])
        .addFile(toy_detector_files["material"])
    )
    return detray.io.readDetector(detray.core.HostMemoryResource(), reader_config)
