# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import detray.io


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
