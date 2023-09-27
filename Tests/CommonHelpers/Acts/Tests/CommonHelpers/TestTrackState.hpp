// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"

#include <random>

namespace Acts::Test {

struct TestTrackState {
  std::shared_ptr<Surface> surface;
  TestSourceLink sourceLink;
  BoundTrackParameters predicted;
  BoundTrackParameters filtered;
  BoundTrackParameters smoothed;
  BoundMatrix jacobian;
  double chi2;
  double pathLength;

  // Generate a random TestTrackState.
  //
  // @param rng Random number generator
  // @param size_t nMeasurement either 1 or 2
  template <typename rng_t>
  TestTrackState(rng_t& rng, size_t measdim)
      : surface(Surface::makeShared<PlaneSurface>(Vector3::Zero(),
                                                  Vector3::UnitZ())),
        // set bogus parameters first since they are not default-constructible
        predicted(surface, BoundVector::Zero()),
        filtered(surface, BoundVector::Zero()),
        smoothed(surface, BoundVector::Zero()),
        jacobian(BoundMatrix::Identity()),
        chi2(std::chi_squared_distribution<double>(measdim)(rng)),
        pathLength(std::uniform_real_distribution<ActsScalar>(
            1 * Acts::UnitConstants::mm, 10 * Acts::UnitConstants::mm)(rng)) {
    // set a random geometry identifier to uniquely identify each surface
    auto geoId =
        std::uniform_int_distribution<GeometryIdentifier::Value>()(rng);
    surface->assignGeometryId(geoId);

    // create source link w/ inline 1d or 2d measurement data
    if (measdim == 1u) {
      auto [par, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
      sourceLink = TestSourceLink(eBoundLoc0, par[0], cov(0, 0), geoId);
    } else if (measdim == 2u) {
      auto [par, cov] = generateParametersCovariance<ActsScalar, 2u>(rng);
      sourceLink = TestSourceLink(eBoundLoc1, eBoundQOverP, par, cov, geoId);
    } else {
      throw std::runtime_error("invalid number of measurement dimensions");
    }

    // create track parameters
    auto [trkPar, trkCov] = generateBoundParametersCovariance(rng);
    // trkPar[eBoundPhi] = 45_degree;
    // trkPar[eBoundTheta] = 90_degree;
    // trkPar[eBoundQOverP] = 5.;
    // predicted
    predicted = BoundTrackParameters(surface, trkPar, trkCov);
    // filtered, modified q/p, reduced covariance
    // trkPar[eBoundQOverP] = 10.;
    filtered = BoundTrackParameters(surface, trkPar, 0.75 * trkCov);
    // smoothed, modified q/p, further reduced covariance
    // trkPar[eBoundQOverP] = 15.;
    smoothed = BoundTrackParameters(surface, trkPar, 0.5 * trkCov);

    // propagation jacobian is identity + corrections
    for (Eigen::Index c = 0; c < jacobian.cols(); ++c) {
      for (Eigen::Index r = 0; r < jacobian.rows(); ++r) {
        jacobian(c, r) +=
            std::uniform_real_distribution<ActsScalar>(-0.125, 0.125)(rng);
      }
    }
  }
};

// Fill a TrackStateProxy with values from a TestTrackState.
//
// @param[in] pc TestTrackState with the input values
// @param[in] mask Specifies which components are used/filled
// @param[out] ts TrackStateProxy which is filled
// @param [in] measdim Dimension of the measurement
template <typename track_state_t>
void fillTrackState(const TestTrackState& pc, TrackStatePropMask mask,
                    track_state_t& ts) {
  // always set the reference surface
  ts.setReferenceSurface(pc.predicted.referenceSurface().getSharedPtr());

  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
    ts.predicted() = pc.predicted.parameters();
    assert(pc.predicted.covariance().has_value());
    ts.predictedCovariance() = *(pc.predicted.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
    ts.filtered() = pc.filtered.parameters();
    assert(pc.filtered.covariance().has_value());
    ts.filteredCovariance() = *(pc.filtered.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
    ts.smoothed() = pc.smoothed.parameters();
    assert(pc.smoothed.covariance().has_value());
    ts.smoothedCovariance() = *(pc.smoothed.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian)) {
    ts.jacobian() = pc.jacobian;
  }
  ts.chi2() = pc.chi2;
  ts.pathLength() = pc.pathLength;
  // source link defines the uncalibrated measurement
  ts.setUncalibratedSourceLink(Acts::SourceLink{pc.sourceLink});
  // create calibrated measurements from source link
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
    testSourceLinkCalibrator<VectorMultiTrajectory>(Acts::GeometryContext{},
                                                    ts);
  }
}

}  // namespace Acts::Test
