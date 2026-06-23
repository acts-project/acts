// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace ActsExamples {

using StripModulePairMap =
    std::unordered_map<Acts::GeometryIdentifier, Acts::GeometryIdentifier>;

StripModulePairMap pairStripModules(
    const Acts::TrackingGeometry &trackingGeometry,
    const std::vector<Acts::GeometryIdentifier> &stripGeometrySelection,
    const Acts::Logger &logger);

}  // namespace ActsExamples
