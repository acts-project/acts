// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp"

#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <DDRec/DetectorData.h>

namespace ActsPlugins {
/// @addtogroup dd4hep_plugin
/// @{

/// Allowed binning directions for DD4hep conversion
static const std::vector<std::tuple<std::string, Acts::AxisDirection>>
    allowedBinnings = {{"x", Acts::AxisDirection::AxisX},
                       {"y", Acts::AxisDirection::AxisY},
                       {"z", Acts::AxisDirection::AxisZ},
                       {"phi", Acts::AxisDirection::AxisPhi},
                       {"r", Acts::AxisDirection::AxisR}};

/// Helper method to convert the string to binning value
///
/// @param binningString
///
/// @return a binningValue
inline Acts::AxisDirection stringToAxisDirection(
    const std::string &binningString) {
  using enum Acts::AxisDirection;
  if (binningString == "x") {
    return AxisX;
  } else if (binningString == "y") {
    return AxisY;
  } else if (binningString == "z") {
    return AxisZ;
  } else if (binningString == "phi") {
    return AxisPhi;
  } else if (binningString == "r") {
    return AxisR;
  } else {
    throw std::invalid_argument("DD4hepBinningHelpers: Binning value " +
                                binningString + " not allowed.");
  }
}

/// Helper method to cenvert a binning list string to a vector of binning values
/// e.g. "r,z" -> {AxisR, AxisZ}
///
/// @param binningString
/// @param del the delimiter for the splitting
///
/// @return a vector of binninng values
inline std::vector<Acts::AxisDirection> stringToAxisDirections(
    const std::string &binningString, const char &del = ',') {
  if (binningString.empty()) {
    return {};
  }
  std::vector<Acts::AxisDirection> bBinning;
  std::stringstream s(binningString);
  std::string b = "";
  while (getline(s, b, del)) {
    bBinning.push_back(stringToAxisDirection(b));
  }
  return bBinning;
}

namespace DD4hepBinningHelpers {
/// @ingroup dd4hep_plugin

/// @brief This method converts the DD4hep binning into the Acts ProtoAxis
///
/// @param dd4hepElement the element which has a binning description attached
/// @param bname the binning base name, e.g. surface_binning, material_binning
///
/// @return a vector of proto binning descriptions
std::vector<std::tuple<Acts::DirectedProtoAxis, std::size_t>> convertBinning(
    const dd4hep::DetElement &dd4hepElement, const std::string &bname);

}  // namespace DD4hepBinningHelpers
/// @}
}  // namespace ActsPlugins
