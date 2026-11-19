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

#include <detray/detectors/metadata_registry.hpp>

namespace traccc {

/// Concrete aliases for every available metadata family. Backend instantiation
/// remains controlled by TRACCC_SUPPORTED_DETECTORS.
// clang-format off
#define TRACCC_DEFINE_DETECTOR(ID) \
  using ID##_detector = detray::host::detector< \
      detray::registry::ID##_metadata<traccc::default_algebra>>;
DETRAY_FOR_EACH_METADATA(TRACCC_DEFINE_DETECTOR)
#undef TRACCC_DEFINE_DETECTOR
// clang-format on

using available_detector_type_list =
    detray::registry::detector_types<traccc::default_algebra>;
using available_metadata_list =
    detray::registry::metadata_types<traccc::default_algebra>;

using default_detector_traits = detray::detector_traits<default_detector>;
using wire_chamber = wire_chamber_detector;

}  // namespace traccc
