# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import gc
import weakref

import detray.core
import detray.tests


def test_detector_keeps_memory_resource_alive():
    memory_resource = detray.core.HostMemoryResource()
    mr_ref = weakref.ref(memory_resource)

    det, names = detray.tests.buildToyDetector(
        memory_resource, detray.tests.ToyDetectorConfig()
    )

    del memory_resource
    gc.collect()

    assert mr_ref() is not None

    del det, names
    gc.collect()

    assert mr_ref() is None
