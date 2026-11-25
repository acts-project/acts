# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pytest

import detray.tests


def test_toy_detector_config():
    cfg = detray.tests.ToyDetectorConfig()
    cfg.n_brl_layers = 2
    cfg.n_edc_layers = 0
    cfg.use_material_maps = True

    assert cfg.n_brl_layers == 2
    assert cfg.n_edc_layers == 0
    assert cfg.use_material_maps
