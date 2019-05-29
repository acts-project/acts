// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/ImpactPoint3dEstimator.hpp"

double Acts::ImpactPoint3dEstimator::calculateDistance(
    const BoundParameters& params, const Vector3D& refPos) const {
  // normalized momentum
  const Vector3D normMomentum = params.momentum().normalized();

  // calculate the path length
  double pathLength = (refPos - params.position()).dot(normMomentum);

  // Point of closest approach in 3D
  Vector3D closestApp = params.position() + pathLength * normMomentum;
  // Difference vector
  Vector3D deltaR = refPos - closestApp;

  // return distance
  return deltaR.norm();
}
