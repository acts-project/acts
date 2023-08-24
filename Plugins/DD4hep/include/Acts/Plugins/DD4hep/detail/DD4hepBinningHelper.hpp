// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <string>
#include <vector>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Objects.h>
#include <DDRec/DetectorData.h>
#include <XML/Utilities.h>

namespace Acts {
namespace detail {

/// Helper method to decode the  binning from what would appear in the
/// xml into variant parameters, such that it can be understood in the
/// downstream processing.
///
/// This parses the dediced <surface_binning/> tag
/// - allowed/understood binnings are x,y,z,phi,r
/// - allowed/unserstood types are equidistant/variable (those are
/// auto-detected)
///
/// Example for e.g. bname = "surface_binning":
///
/// - Equidistant binning in r and phi:
///   <surface_binning nr="2" rmin="25" rmax="100" nphi="22" phimin="-3.1415"
///   phimax="3.1415"/>
/// - Variable binning in z:
///   <surface_binning zboundaries="-100,-90,90,100"/>
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
                   const std::vector<std::string> &bvals) {
  // Set the surface binninng parameter to true
  variantParams.set<bool>(bname, true);
  for (const auto &bv : bvals) {
    // Gather the number of bins, 0 indicates vairable binning
    int nBins =
        Acts::detail::getAttrValueOr<int>(xmlBinning, std::string("n" + bv), 0);
    // Gather the bin expansion parameter, expansion of 0 is default
    int nExpansion = Acts::detail::getAttrValueOr<int>(
        xmlBinning, std::string(bv + "expansion"), 0);
    variantParams.set<int>(bname + "_" + bv + "_exp", nExpansion);
    // Equidistant binning detected
    if (nBins > 0) {
      // Set the type identificatio
      variantParams.set<std::string>(bname + "_" + bv + "_type", "equidistant");
      // Set the number of bins
      variantParams.set<int>(bname + "_" + bv + "_count", nBins);
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
      // Get the number of bins explicitely
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

}  // namespace detail
}  // namespace Acts