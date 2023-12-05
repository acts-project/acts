// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepBinningHelpers.hpp"

void Acts::DD4hepBinningHelpers::decodeBinning(
    dd4hep::rec::VariantParameters &variantParams, const xml_comp_t &xmlBinning,
    const std::string &bname, const std::vector<std::string> &bvals) {
  // Set the surface binninng parameter to true
  variantParams.set<int>(std::string(bname + "_dim"), bvals.size());
  for (const auto &bv : bvals) {
    // Gather the number of bins, 0 indicates variable binning
    int nBins = Acts::getAttrValueOr<int>(xmlBinning, std::string("n" + bv), 0);
    // Gather the bin expansion parameter, expansion of 0 is default
    int nExpansion =
        Acts::getAttrValueOr<int>(xmlBinning, std::string(bv + "expansion"), 0);
    // Auto-range detection
    bool autoRange = Acts::getAttrValueOr<bool>(
        xmlBinning, std::string(bv + "autorange"), false);
    variantParams.set<bool>(bname + "_" + bv + "_autorange", autoRange);
    variantParams.set<int>(bname + "_" + bv + "_exp", nExpansion);
    // Equidistant binning detected
    if (nBins > 0) {
      // Set the type identificatio
      variantParams.set<std::string>(bname + "_" + bv + "_type", "equidistant");
      // Set the number of bins
      variantParams.set<int>(bname + "_" + bv + "_n", nBins);
      // Set min/max paraeter
      if (!autoRange) {
        variantParams.set<double>(
            bname + "_" + bv + "_min",
            xmlBinning.attr<double>(std::string(bv + "min").c_str()));
        variantParams.set<double>(
            bname + "_" + bv + "_max",
            xmlBinning.attr<double>(std::string(bv + "max").c_str()));
      }
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

std::vector<Acts::Experimental::ProtoBinning>
Acts::DD4hepBinningHelpers::convertBinning(
    const dd4hep::DetElement &dd4hepElement, const std::string &bname) {
  // Return proto binning vector
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
      // Indicate auto-range checking
      bool autoRange = getParamOr<bool>(bname + "_" + ab + "_autorange",
                                        dd4hepElement, false);
      // Equidistant binning
      if (aType == detail::AxisType::Equidistant) {
        if (autoRange) {
          protoBinnings.push_back(
              Experimental::ProtoBinning(bVal, bType, nBins, nExpansion));
        } else {
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
        }
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
