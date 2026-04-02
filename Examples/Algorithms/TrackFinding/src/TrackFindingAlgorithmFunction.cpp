// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/MultiStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include <memory>
#include <utility>

namespace ActsExamples {

namespace {

using Stepper = Acts::SympyStepper;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using CKF = Acts::CombinatorialKalmanFilter<Propagator, TrackContainer>;

using BremStepper = Acts::MultiStepperLoop<Stepper>;
using BremPropagator = Acts::Propagator<BremStepper, Navigator>;
using BremCKF = Acts::CombinatorialKalmanFilter<BremPropagator, TrackContainer>;

template <typename CKFType>
struct TrackFinderFunctionImpl
    : public TrackFindingAlgorithm::TrackFinderFunction {
  CKFType trackFinder;

  explicit TrackFinderFunctionImpl(CKFType&& f) : trackFinder(std::move(f)) {}

  TrackFindingAlgorithm::TrackFinderResult operator()(
      const TrackParameters& initialParameters,
      const TrackFindingAlgorithm::TrackFinderOptions& options,
      TrackContainer& tracks, TrackProxy rootBranch) const override {
    return trackFinder.findTracks(initialParameters, options, tracks,
                                  rootBranch);
  }
};

}  // namespace

std::shared_ptr<TrackFindingAlgorithm::TrackFinderFunction>
TrackFindingAlgorithm::makeTrackFinderFunction(
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
  return std::make_shared<TrackFinderFunctionImpl<CKF>>(std::move(trackFinder));
}

std::shared_ptr<TrackFindingAlgorithm::TrackFinderFunction>
TrackFindingAlgorithm::makeBremTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    const Acts::Logger& logger) {
  BremStepper stepper(std::move(magneticField));
  Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  BremPropagator propagator(std::move(stepper), std::move(navigator),
                            logger.cloneWithSuffix("BremPropagator"));
  BremCKF trackFinder(std::move(propagator),
                      logger.cloneWithSuffix("BremFinder"));

  // build the track finder functions. owns the track finder object.
  return std::make_shared<TrackFinderFunctionImpl<BremCKF>>(
      std::move(trackFinder));
}

}  // namespace ActsExamples
