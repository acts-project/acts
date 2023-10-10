// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
class TrackingGeometry;
}  // namespace Acts

namespace {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;

using Stepper = Acts::EigenStepper<>;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using CKF =
    Acts::CombinatorialKalmanFilter<Propagator, Acts::VectorMultiTrajectory>;

using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;

struct TrackFinderFunctionImpl
    : public ActsExamples::TrackFindingAlgorithm::TrackFinderFunction {
  CKF trackFinder;

  TrackFinderFunctionImpl(CKF&& f) : trackFinder(std::move(f)) {}

  ActsExamples::TrackFindingAlgorithm::TrackFinderResult operator()(
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFindingAlgorithm::TrackFinderOptions& options,
      TrackContainer& tracks) const override {
    return trackFinder.findTracks(initialParameters, options, tracks);
  };
};

}  // namespace

std::shared_ptr<ActsExamples::TrackFindingAlgorithm::TrackFinderFunction>
ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    const Acts::Logger& logger) {
  Stepper stepper(std::move(magneticField));
  Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("Propagator"));
  CKF trackFinder(std::move(propagator), logger.cloneWithSuffix("Finder"));

  // build the track finder functions. owns the track finder object.
  return std::make_shared<TrackFinderFunctionImpl>(std::move(trackFinder));
}
