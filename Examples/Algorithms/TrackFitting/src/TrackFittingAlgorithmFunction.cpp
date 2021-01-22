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

template <typename track_fitter_t>
struct TrackFitterFunctionImpl {
  track_fitter_t trackFitter;

  TrackFitterFunctionImpl(track_fitter_t&& f) : trackFitter(std::move(f)) {}

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::TrackFitterOptions& options)
      const {
    return trackFitter.fit(sourceLinks, initialParameters, options);
  };
};

template <typename Fitter>
struct DirectedFitterFunctionImpl {
  Fitter fitter;
  DirectedFitterFunctionImpl(Fitter&& f) : fitter(std::move(f)) {}

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
    MagneticField magneticField) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [trackingGeometry](auto&& inputField) -> TrackFitterFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using SharedMagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::EigenStepper<>;
        using Navigator = Acts::Navigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

        // construct all components for the fitter
        auto field =
            std::make_shared<SharedMagneticField>(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter trackFitter(std::move(propagator));

        // build the fitter functions. owns the fitter object.
        return TrackFitterFunctionImpl<Fitter>(std::move(trackFitter));
      },
      std::move(magneticField));
}

ActsExamples::TrackFittingAlgorithm::DirectedTrackFitterFunction
ActsExamples::TrackFittingAlgorithm::makeTrackFitterFunction(
    MagneticField magneticField) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [](auto&& inputField) -> DirectedTrackFitterFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using SharedMagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::EigenStepper<>;
        using Navigator = Acts::DirectNavigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

        // construct all components for the fitter
        auto field =
            std::make_shared<SharedMagneticField>(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter fitter(std::move(propagator));

        // build the fitter functions. owns the fitter object.
        return DirectedFitterFunctionImpl<Fitter>(std::move(fitter));
      },
      std::move(magneticField));
}
