# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pytest

import detray.core
import detray.tests


def test_material_scan_config():
    cfg = detray.tests.MaterialScanConfig()
    cfg.trackGenerator.etaSteps = 9
    assert cfg.trackGenerator.etaSteps == 9


def test_material_validation_config():
    cfg = detray.tests.MaterialValidationConfig()
    cfg.materialFile = "navigation_output"
    cfg.nTracks = 42
    cfg.propagation.stepping.pathLimit = 1234.0

    assert cfg.materialFile == "navigation_output"
    assert cfg.nTracks == 42
    assert cfg.propagation.stepping.pathLimit == pytest.approx(1234.0, rel=1e-6)


def test_toy_detector_config():
    cfg = detray.tests.ToyDetectorConfig()
    cfg.nBarrelLayers = 2
    cfg.nEndcapLayers = 0
    cfg.useMaterialMaps = True

    assert cfg.nBarrelLayers == 2
    assert cfg.nEndcapLayers == 0
    assert cfg.useMaterialMaps


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
