// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Acts Examples include(s)
#include "ActsExamples/Traccc/Common/Util/MapUtil.hpp"

// Traccc include(s)
#include "traccc/edm/cell.hpp"
#include "traccc/geometry/geometry.hpp"
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <map>
#include <vector>

namespace ActsExamples::Traccc::Common::Conversion {

auto createModules(std::vector<Acts::GeometryIdentifier::Value> geomentryIds, std::map<Acts::GeometryIdentifier::Value, detray::geometry::barcode> barcodeMap){
    return 0;
}

}  // namespace ActsExamples::Traccc::Common::Conversion
