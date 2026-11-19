# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pytest
import detray.core


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


def test_propagation_configuration():
    config = detray.core.PropagationConfig()
    config.navigation.searchWindow = [2, 3]
    config.stepping.pathLimit = 42
    assert list(config.navigation.searchWindow) == [2, 3]
    assert config.stepping.pathLimit == pytest.approx(42)
