// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace Acts::detail::GeoModelExentHelper {

/// @brief Helper function to find out which ones are constraint from internal / external
///
/// @param boundsEntry the bounds entry from the database
/// @param boundsPos the position of the bound in the entry
/// @param ctype the type of the constraint as string from the database
///
/// @return a vector
std::vector<BinningValue> readConstaints(
    const std::vector<std::string>& boundsEntry, std::size_t boundsPos,
    const std::string& ctype = "i");

/// @brief Helper function to create the extent from database volume entry
///
/// @param boundsEntry the bounds entry from the database
/// @param boundsPos the position of the bound in the entry
/// @param externalExtent the extend from external constraints (marked "e" in the database)
/// @param internalExtent the extend of the internal objects (marked "i" in the database)
///
/// @return a tuple of the bounds type, the extent, and a list of binning values to be determined from the internals
std::tuple<VolumeBounds::BoundsType, Extent> extentFromTable(
    const std::vector<std::string>& boundsEntry, std::size_t boundsPos,
    const Extent& externalExtent = Extent(),
    const Extent& internalExtent = Extent());

}  // namespace Acts::detail::GeoModelExentHelper
