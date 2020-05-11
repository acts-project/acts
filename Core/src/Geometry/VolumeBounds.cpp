// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/VolumeBounds.hpp"

Acts::OrientedSurfaces Acts::VolumeBounds::orientedSurfaces(
    const Transform3D* transform) const noexcept(false) {
  // Get the raw surfaces first
  auto rawSurfaces = decomposeToSurfaces(transform);
  auto orientations = boundaryOrientations();

  OrientedSurfaces oSurfaces;
  oSurfaces.reserve(rawSurfaces.size());
  if (rawSurfaces.size() > orientations.size()) {
    throw std::invalid_argument(
        "VolumeBounds: not enough surface orientations given.");
  }
  //  Fill the surface vector and the associated orientataions
  for (size_t is = 0; is < rawSurfaces.size(); ++is) {
    oSurfaces.push_back(OrientedSurface(rawSurfaces[is], orientations[is]));
  }
  return oSurfaces;
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream& Acts::operator<<(std::ostream& sl, const Acts::VolumeBounds& vb) {
  return vb.toStream(sl);
}
