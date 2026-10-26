# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pytest

import acts


def test_proto_volume_binning():
    default = acts.ProtoVolumeMaterial()
    assert len(default.binning) == 3
    assert default.binning.deferred
    assert [default.binning[i].nBins for i in range(3)] == [1, 1, 1]

    binning = acts.MultiAxisSpec(
        [
            acts.AxisSpec.DeferredEquidistant(10, acts.AxisDirection.AxisR),
            acts.AxisSpec.DeferredEquidistant(20, acts.AxisDirection.AxisPhi),
            acts.AxisSpec.DeferredVariable(
                [0.0, 0.1, 1.0], direction=acts.AxisDirection.AxisZ
            ),
        ]
    )
    proxy = acts.ProtoVolumeMaterial(binning)
    assert proxy.binning == binning
    assert proxy.binning[2].nBins == 2

    with pytest.raises(ValueError, match="three axes"):
        acts.ProtoVolumeMaterial(
            acts.MultiAxisSpec([acts.AxisSpec.DeferredEquidistant(1)])
        )
