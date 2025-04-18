// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepBinningHelpers.hpp"

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
          ActsScalar minDefault = bVal == binPhi ? -M_PI : 0.;
          ActsScalar maxDefault = bVal == binPhi ? M_PI : 0.;
          auto min = getParamOr<ActsScalar>(bname + "_" + ab + "_min",
                                            dd4hepElement, minDefault);
          auto max = getParamOr<ActsScalar>(bname + "_" + ab + "_max",
                                            dd4hepElement, maxDefault);
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
