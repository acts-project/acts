// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/GeometryIdentifierHooks.hpp"

namespace det {
namespace GeometryIdentifierHooks {

Acts::GeometryIdentifier stripEndcapODD(Acts::GeometryIdentifier identifier,
                                        const Acts::Surface& surface) {
  // @TODO : This implementation hardcode the parameters for the ODD,
  // ideally a more generic implementation (using the python binding would be
  // preferable)

  // Define the spliting parameters for the ODD.
  // In case of geometry change please edit those accordingly

  // volume where the spliting should be performed (strip endcaps)
  std::vector<int> volumes = {28, 30};
  // list of the radial cutoff point (smallest to largest)
  std::vector<double> radialCuts = {580};
  // list of the radial cutoff point for each volumes
  std::vector<std::vector<double>> cut = {radialCuts, radialCuts};

  auto it = std::find(volumes.begin(), volumes.end(), identifier.volume());
  if (it != volumes.end()) {
    std::vector<double> layercut = cut[it - volumes.begin()];
    Acts::Vector3 center = surface.center(Acts::GeometryContext());
    double radius = sqrt(center[0] * center[0] + center[1] * center[1]);
    identifier.setExtra(1);
    for (unsigned i = 0; i < layercut.size(); i++) {
      if (radius > layercut[i]) {
        identifier.setExtra(i + 2);
      }
    }
  }
  return identifier;
}

}  // namespace GeometryIdentifierHooks
}  // namespace det
