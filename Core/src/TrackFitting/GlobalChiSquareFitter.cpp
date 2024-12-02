// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

void Acts::Experimental::updateGx2fParams(
    BoundTrackParameters& params, const Eigen::VectorXd& deltaParamsExtended,
    const std::size_t nMaterialSurfaces,
    std::unordered_map<GeometryIdentifier, ScatteringProperties>& scatteringMap,
    const std::vector<GeometryIdentifier>& geoIdVector) {
  // update params
  params.parameters() +=
      deltaParamsExtended.topLeftCorner<eBoundSize, 1>().eval();

  // update the scattering angles.
  for (std::size_t matSurface = 0; matSurface < nMaterialSurfaces;
       matSurface++) {
    const std::size_t deltaPosition = eBoundSize + 2 * matSurface;
    const GeometryIdentifier geoId = geoIdVector[matSurface];
    const auto scatteringMapId = scatteringMap.find(geoId);
    assert(scatteringMapId != scatteringMap.end() &&
           "No scattering angles found for material surface.");
    scatteringMapId->second.scatteringAngles().block<2, 1>(2, 0) +=
        deltaParamsExtended.block<2, 1>(deltaPosition, 0).eval();
  }

  return;
}

void Acts::Experimental::updateGx2fCovarianceParams(
    BoundMatrix& fullCovariancePredicted, Gx2fSystem& extendedSystem) {
  // make invertible
  for (std::size_t i = 0; i < extendedSystem.nDims(); ++i) {
    if (extendedSystem.aMatrix()(i, i) == 0.) {
      extendedSystem.aMatrix()(i, i) = 1.;
    }
  }

  visit_measurement(extendedSystem.findRequiredNdf(), [&](auto N) {
    fullCovariancePredicted.topLeftCorner<N, N>() =
        extendedSystem.aMatrix().inverse().topLeftCorner<N, N>();
  });

  return;
}
