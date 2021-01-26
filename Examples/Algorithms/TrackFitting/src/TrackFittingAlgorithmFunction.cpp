// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

namespace {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;
using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;
using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
using DirectFitter = Acts::KalmanFitter<DirectPropagator, Updater, Smoother>;

struct TrackFitterFunctionImpl {
  Fitter trackFitter;

  TrackFitterFunctionImpl(Fitter&& f) : trackFitter(std::move(f)) {}

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::TrackFitterOptions& options)
      const {
    return trackFitter.fit(sourceLinks, initialParameters, options);
  };
};

struct DirectedFitterFunctionImpl {
  DirectFitter fitter;
  DirectedFitterFunctionImpl(DirectFitter&& f) : fitter(std::move(f)) {}

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::TrackFitterOptions& options,
      const std::vector<const Acts::Surface*>& sSequence) const {
    return fitter.fit(sourceLinks, initialParameters, options, sSequence);
  };
};

}  // namespace

ActsExamples::TrackFittingAlgorithm::TrackFitterFunction
ActsExamples::TrackFittingAlgorithm::makeTrackFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::BFieldProvider> magneticField) {
  Stepper stepper(std::move(magneticField));
  Acts::Navigator navigator(trackingGeometry);
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  return TrackFitterFunctionImpl(std::move(trackFitter));
}

ActsExamples::TrackFittingAlgorithm::DirectedTrackFitterFunction
ActsExamples::TrackFittingAlgorithm::makeTrackFitterFunction(
    std::shared_ptr<const Acts::BFieldProvider> magneticField) {
  // construct all components for the fitter
  Stepper stepper(std::move(magneticField));
  Acts::DirectNavigator navigator;
  DirectPropagator propagator(std::move(stepper), std::move(navigator));
  DirectFitter fitter(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  return DirectedFitterFunctionImpl(std::move(fitter));
}
