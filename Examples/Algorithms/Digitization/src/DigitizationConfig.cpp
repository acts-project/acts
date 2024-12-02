// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationConfig.hpp"

namespace ActsExamples {

std::vector<double> GeometricConfig::variances(
    const std::array<std::size_t, 2u>& csizes,
    const std::array<std::size_t, 2u>& cmins) const {
  std::vector<double> rVariances;
  for (const auto& bIndex : indices) {
    double var = 0.;
    if (varianceMap.contains(bIndex)) {
      // Try to find the variance for this cluster size
      std::size_t lsize =
          std::min(csizes[bIndex], varianceMap.at(bIndex).size());
      var = varianceMap.at(bIndex).at(lsize - 1);
    } else {
      // Pitch size ofer / sqrt(12) as error instead
      std::size_t ictr = cmins[bIndex] + csizes[bIndex] / 2;
      var = std::pow(segmentation.binningData()[bIndex].width(ictr), 2) / 12.0;
    }
    rVariances.push_back(var);
  }
  return rVariances;
}

}  // namespace ActsExamples
