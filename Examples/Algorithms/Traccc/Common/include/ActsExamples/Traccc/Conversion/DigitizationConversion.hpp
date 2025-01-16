// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Plugins/Traccc/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Creates a traccc digitalization config from an Acts geometry hierarchy map
/// that contains the digitization configuration.
/// @param config the Acts geometry hierarchy map that contains the digitization configuration.
/// @return a traccc digitization config.
Acts::TracccPlugin::DigitizationConfig tracccConfig(
    const Acts::GeometryHierarchyMap<DigiComponentsConfig>& config);

}  // namespace ActsExamples::Traccc::Common::Conversion
