// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/BFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include <random>
#include <stdexcept>

namespace {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;

using Stepper = Acts::EigenStepper<>;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using CKF = Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother>;

struct TrackFinderFunctionImpl {
  CKF trackFinder;

  TrackFinderFunctionImpl(CKF&& f) : trackFinder(std::move(f)) {}

  ActsExamples::TrackFindingAlgorithm::TrackFinderResult operator()(
      const ActsExamples::IndexSourceLinkContainer& sourcelinks,
      const ActsExamples::TrackParametersContainer& initialParameters,
      const ActsExamples::TrackFindingAlgorithm::TrackFinderOptions& options)
      const {
    return trackFinder.findTracks(sourcelinks, initialParameters, options);
  };
};

}  // namespace

ActsExamples::TrackFindingAlgorithm::TrackFinderFunction
ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::BFieldProvider> magneticField) {
  Stepper stepper(std::move(magneticField));
  Navigator navigator(trackingGeometry);
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
  Propagator propagator(std::move(stepper), std::move(navigator));
  CKF trackFinder(std::move(propagator));

  // build the track finder functions. owns the track finder object.
  return TrackFinderFunctionImpl(std::move(trackFinder));
}
