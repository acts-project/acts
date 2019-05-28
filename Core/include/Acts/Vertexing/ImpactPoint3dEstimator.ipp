// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename input_track_t>
double Acts::ImpactPoint3dEstimator<input_track_t>::calculateDistance(
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

template <typename input_track_t>
Acts::BoundParameters
Acts::ImpactPoint3dEstimator<input_track_t>::getParamsAtIP3d(
    const GeometryContext& geoCtx, const TrackAtVertex<input_track_t>& trk,
    const Vector3D& refPos) const {
  // temp TODO
  std::unique_ptr<BoundSymMatrix> covMat = std::make_unique<BoundSymMatrix>();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));
  TrackParametersBase::ParVector_t paramVec;
  return BoundParameters(geoCtx, std::move(covMat), paramVec, perigeeSurface);
}