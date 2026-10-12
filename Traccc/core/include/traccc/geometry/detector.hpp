// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// Detray include(s).

// Turn off excessive detector building debug logs (only works with gcc!)
// clang-format off
#include <detray/utils/quiet_log_start.hpp>
#include <detray/core/detector.hpp>
#include <detray/utils/quiet_log_end.hpp>
// clang-format on

#include <detray/detectors/default_metadata.hpp>
#include <detray/detectors/itk_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/detectors/telescope_metadata.hpp>
#include <detray/detectors/toy_metadata.hpp>
#include <detray/detectors/wire_chamber_metadata.hpp>

namespace traccc {

/// Default detector (can contain (almost) any detector data)
using default_detector =
    detray::host::detector<detray::default_metadata<traccc::default_algebra>>;

using default_detector_traits = detray::detector_traits<default_detector>;

/// ATLAS Inner Tracker (ITk) detector
using itk_detector =
    detray::host::detector<detray::itk_metadata<traccc::default_algebra>>;

/// Open Data Detector (ODD) detector
using odd_detector =
    detray::host::detector<detray::odd_metadata<traccc::default_algebra>>;

/// Detray telescope detector (test detector)
using telescope_detector = detray::host::detector<
    detray::telescope_metadata<traccc::default_algebra, detray::rectangle2D>>;

/// Detray toy detector (test detector)
using toy_detector =
    detray::host::detector<detray::toy_metadata<traccc::default_algebra>>;

/// Detray wire chamber detector (test detector)
using wire_chamber = detray::host::detector<
    detray::wire_chamber_metadata<traccc::default_algebra>>;

}  // namespace traccc
