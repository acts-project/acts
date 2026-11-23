# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os
import tempfile

import detray.core
import detray.io
import detray.tests


def test_generate_and_read_detector():
    # Use the toy generator to produce detector files to read back.
    out_dir = tempfile.mkdtemp() + os.sep

    det, det_names = detray.tests.build_toy_detector(
        detray.core.HostMemoryResource(), detray.tests.ToyDetectorConfig()
    )

    writer_config = detray.io.DetectorWriterConfig()
    writer_config.path = out_dir
    detray.io.write_detector_json(det, det_names, writer_config)

    geometry_file = os.path.join(out_dir, "toy_detector_geometry.json")
    assert os.path.exists(geometry_file)

    reader_config = detray.io.DetectorReaderConfig().add_file(geometry_file)

    detector, names = detray.io.read_detector_json(
        detray.core.HostMemoryResource(), reader_config
    )
    assert isinstance(names, detray.core.NameMap)
    assert len(detector.volumes) == 22
    assert len(detector.surfaces) == 3230
