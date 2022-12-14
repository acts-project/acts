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

Acts::GeometryIdentifier RadiusGeometryIdentifierDecorator::decorateIdentifier(
    Acts::GeometryIdentifier identifier, const Acts::Surface& surface) const {
  auto it = volumeToRadialCuts.find(identifier.volume());
  if (it != volumeToRadialCuts.end()) {
    const std::vector<double> &layercut = volumeToRadialCuts.at(identifier.volume());
    Acts::Vector3 center = surface.center(Acts::GeometryContext());
    double radius = sqrt(center[0] * center[0] + center[1] * center[1]);
    std::cout << "HOOK " << radius << " " << it->first << "\n";
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
