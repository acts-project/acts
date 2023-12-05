// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Objects.h>
#include <DDRec/DetectorData.h>
#include <XML/Utilities.h>

namespace Acts {

static std::vector<std::tuple<std::string, BinningValue>> allowedBinnings = {
    {"x", binX}, {"y", binY}, {"z", binZ}, {"phi", binPhi}, {"r", binR}};

/// Helper method to convert the string to binning value
///
/// @param binningString
///
/// @return a binningValue
inline BinningValue stringToBinningValue(const std::string &binningString) {
  if (binningString == "x") {
    return binX;
  } else if (binningString == "y") {
    return binY;
  } else if (binningString == "z") {
    return binZ;
  } else if (binningString == "phi") {
    return binPhi;
  } else if (binningString == "r") {
    return binR;
  } else {
    throw std::invalid_argument("DD4hepBinningHelpers: Binning value " +
                                binningString + " not allowed.");
  }
}

/// Helper method to cenvert a binning list string to a vector of binning values
/// e.g. "r,z" -> {binR, binZ}
///
/// @param binningString
/// @param del the delimiter for the splitting
///
/// @return a vector of binninng values
inline std::vector<BinningValue> stringToBinningValues(
    const std::string &binningString, const char &del = ',') {
  if (binningString.empty()) {
    return {};
  }
  std::vector<BinningValue> bBinning;
  std::stringstream s(binningString);
  std::string b = "";
  while (getline(s, b, del)) {
    bBinning.push_back(stringToBinningValue(b));
  }
  return bBinning;
}

namespace DD4hepBinningHelpers {

/// Helper method to decode the binning from what would appear in the
/// xml into variant parameters, such that it can be understood in the
/// downstream processing.
///
/// This parses the dediced \< surface_binning \> tag
/// - allowed/understood binnings are x,y,z,phi,r
/// - allowed/unserstood types are equidistant/variable (those are
/// auto-detected)
///
/// Example for e.g. bname = \"surface_binning\":
///
/// - Equidistant binning in r and phi:
///   \<acts_surface_binning nr=\"2\" rmin=\"25\" rmax=\"100\" nphi=\"22\"
///   phimin=\"-3.1415\" phimax=\"3.1415\" \/ \>
/// - Variable binning in z:
///   \<acts_surface_binning zboundaries=\"-100,-90,90,100\" \/ \>
///
/// And 2D combinations of this are allowed.
///
/// @param variantParams [in,out] the variant parameters that will be overwritten
/// @param xmlBinning the surface binning
/// @param bname the binning base name, e.g. surface_binning, material_binning
/// @param bvals the boundary values, i.e. x,y,z,phi,r
///
void decodeBinning(dd4hep::rec::VariantParameters &variantParams,
                   const xml_comp_t &xmlBinning, const std::string &bname,
                   const std::vector<std::string> &bvals);

/// @brief This method converts the DD4hep binning into the Acts ProtoBinning
///
/// @param dd4hepElement the element which has a binning description attached
/// @param bname the binning base name, e.g. surface_binning, material_binning
///
/// @return a vector of proto binning descriptions
std::vector<Acts::Experimental::ProtoBinning> convertBinning(
    const dd4hep::DetElement &dd4hepElement, const std::string &bname);

}  // namespace DD4hepBinningHelpers
}  // namespace Acts
