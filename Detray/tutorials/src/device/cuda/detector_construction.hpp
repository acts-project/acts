// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tutorial/types.hpp"

namespace detray::tutorial {

// Detector
using metadata_t = detray::tutorial::toy_metadata;
using host_detector_t = host::detector<metadata_t>;
using device_detector_t = device::detector<metadata_t>;

using mask_id = typename host_detector_t::masks::id;
using acc_id = typename host_detector_t::accel::id;

/// Detector construction tutorial function (prints some detector statistics)
void print(typename host_detector_t::view_type det_data);

}  // namespace detray::tutorial
