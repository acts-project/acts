// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/GeoIDHooks.hpp"

namespace det {
namespace geoIDHook {

Acts::GeometryIdentifier stripEndcapODD(Acts::GeometryIdentifier identifier,
                                        const Acts::Surface& surface) {
  if (identifier.volume() == 28 || identifier.volume() == 30) {
    Acts::Vector3 center = surface.center(Acts::GeometryContext());
    double radius = sqrt(center[0] * center[0] + center[1] * center[1]);
    if (radius < 850) {
      identifier.setExtra(1);
    } else {
      identifier.setExtra(2);
    }
  }
  return identifier;
}

}  // namespace geoIDHook
}  // namespace det
