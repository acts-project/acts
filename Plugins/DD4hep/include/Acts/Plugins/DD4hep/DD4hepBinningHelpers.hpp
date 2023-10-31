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
#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

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
///   \< surface_binning nr=\"2\" rmin=\"25\" rmax=\"100\" nphi=\"22\"
///   phimin=\"-3.1415\" phimax=\"3.1415\" \/ \>
/// - Variable binning in z:
///   \< surface_binning zboundaries=\"-100,-90,90,100\" \/ \>
///
/// And 2D combinations of this are allowed.
///
/// @param variantParams [in,out] the variant parameters that will be overwritten
/// @param xmlBinning the surface binning
/// @param bname the binning base name, e.g. surface_binning, material_binning
/// @param bvals the boundary values, i.e. x,y,z,phi,r
///
inline void decodeBinning(dd4hep::rec::VariantParameters &variantParams,
                          const xml_comp_t &xmlBinning,
                          const std::string &bname,
                          const std::vector<std::string> &bvals) {
  // Set the surface binninng parameter to true
  variantParams.set<int>(std::string(bname + "_dim"), bvals.size());
  for (const auto &bv : bvals) {
    // Gather the number of bins, 0 indicates variable binning
    int nBins = Acts::getAttrValueOr<int>(xmlBinning, std::string("n" + bv), 0);
    // Gather the bin expansion parameter, expansion of 0 is default
    int nExpansion =
        Acts::getAttrValueOr<int>(xmlBinning, std::string(bv + "expansion"), 0);
    variantParams.set<int>(bname + "_" + bv + "_exp", nExpansion);
    // Equidistant binning detected
    if (nBins > 0) {
      // Set the type identificatio
      variantParams.set<std::string>(bname + "_" + bv + "_type", "equidistant");
      // Set the number of bins
      variantParams.set<int>(bname + "_" + bv + "_n", nBins);
      // Set min/max paraeter
      variantParams.set<double>(
          bname + "_" + bv + "_min",
          xmlBinning.attr<double>(std::string(bv + "min").c_str()));
      variantParams.set<double>(
          bname + "_" + bv + "_max",
          xmlBinning.attr<double>(std::string(bv + "max").c_str()));
    } else {
      // Variable binning detected
      variantParams.set<std::string>(bname + "_" + bv + "_type", "variable");
      // Get the number of bins explicitly
      auto boundaries =
          xmlBinning.attr<std::string>(std::string(bv + "boundaries").c_str());
      std::string del = ",";
      int end = boundaries.find(del);
      int ib = 0;
      // Unit conversion
      double unitScalar = 1.;
      if (bv != "phi") {
        unitScalar = Acts::UnitConstants::mm / dd4hep::millimeter;
      }
      // Split and convert
      while (end != -1) {
        double bR = unitScalar * dd4hep::_toFloat(boundaries.substr(0, end));
        variantParams.set<double>(
            bname + "_" + bv + "_b" + std::to_string(ib++), bR);
        boundaries.erase(boundaries.begin(), boundaries.begin() + end + 1);
        end = boundaries.find(del);
      }
      double bR = unitScalar * std::stod(boundaries.substr(0, end));
      variantParams.set<double>(bname + "_" + bv + "_b" + std::to_string(ib),
                                bR);
      // The number of bins are needed to unpack the data
      variantParams.set<int>(bname + "_" + bv + "_n", ib);
    }
  }
}

/// @brief This method converts the DD4hep binning into the Acts ProtoBinning
///
/// @param dd4hepElement the element which has a binning description attached
/// @param bname the binning base name, e.g. surface_binning, material_binning
///
/// @return a vector of proto binning descriptions
inline std::vector<Acts::Experimental::ProtoBinning> convertBinning(
    const dd4hep::DetElement &dd4hepElement, const std::string &bname) {
  std::vector<Experimental::ProtoBinning> protoBinnings;
  for (const auto &[ab, bVal] : allowedBinnings) {
    auto type =
        getParamOr<std::string>(bname + "_" + ab + "_type", dd4hepElement, "");
    if (!type.empty()) {
      // Default binning is bound
      auto bType = Acts::detail::AxisBoundaryType::Bound;
      // Equidistant or variable binning
      detail::AxisType aType = type == "equidistant"
                                   ? detail::AxisType::Equidistant
                                   : detail::AxisType::Variable;
      int nBins = getParamOr<int>(bname + "_" + ab + "_n", dd4hepElement, 0);
      int nExpansion =
          getParamOr<int>(bname + "_" + ab + "_exp", dd4hepElement, 0);
      if (aType == detail::AxisType::Equidistant) {
        // Equidistant binning
        auto min = getParamOr<ActsScalar>(bname + "_" + ab + "_min",
                                          dd4hepElement, 0.);
        auto max = getParamOr<ActsScalar>(bname + "_" + ab + "_max",
                                          dd4hepElement, 0.);
        // Check for closed phi binning
        if (bVal == binPhi && (max - min) > 1.9 * M_PI) {
          bType = Acts::detail::AxisBoundaryType::Closed;
        }
        protoBinnings.push_back(Experimental::ProtoBinning(
            bVal, bType, min, max, nBins, nExpansion));
      } else {
        // Variable binning
        std::vector<ActsScalar> edges;
        for (int ib = 0; ib <= nBins; ++ib) {
          edges.push_back(getParamOr<ActsScalar>(
              bname + "_" + ab + "_b" + std::to_string(ib), dd4hepElement, 0.));
        }
        // Check for closed phi binning
        if (bVal == binPhi && (edges.back() - edges.front()) > 1.9 * M_PI) {
          bType = Acts::detail::AxisBoundaryType::Closed;
        }
        protoBinnings.push_back(
            Experimental::ProtoBinning(bVal, bType, edges, nExpansion));
      }
    }
  }
  return protoBinnings;
}

}  // namespace Acts
