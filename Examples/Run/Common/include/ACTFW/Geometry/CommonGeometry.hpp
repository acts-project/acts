// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "ACTFW/Utilities/OptionsFwd.hpp"

namespace Acts {
class TrackingGeometry;
}

namespace FW {
class IBaseDetector;
class IContextDecorator;
namespace Geometry {

/// @brief helper method to setup the geometry
///
/// @tparam options_map_t Type of the options to be read
/// @tparam geometry_setupt_t Type of the callable geometry setup
///
/// @param vm the parsed options map
/// @param geometrySetup the callable geometry setup
///
/// @return a pair of TrackingGeometry and context decorators
std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
          std::vector<std::shared_ptr<FW::IContextDecorator>>>
build(const boost::program_options::variables_map& vm, IBaseDetector& detector);

}  // namespace Geometry
}  // namespace FW
