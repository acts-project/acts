// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/GsfFitterFunction.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

using namespace ActsExamples;

namespace {

using MultiStepper = Acts::MultiEigenStepperLoop<>;
using Propagator = Acts::Propagator<MultiStepper, Acts::Navigator>;
using Fitter =
    Acts::Experimental::GaussianSumFitter<Propagator,
                                          Acts::VectorMultiTrajectory>;
using DirectPropagator = Acts::Propagator<MultiStepper, Acts::DirectNavigator>;
using DirectFitter =
    Acts::Experimental::GaussianSumFitter<DirectPropagator,
                                          Acts::VectorMultiTrajectory>;

struct GsfFitterFunctionImpl
    : public ActsExamples::TrackFittingAlgorithm::TrackFitterFunction {
  Fitter fitter;
  DirectFitter directFitter;

  Acts::GainMatrixUpdater updater;

  std::size_t maxComponents;
  bool abortOnError;
  bool disableAllMaterialHandling;

  GsfFitterFunctionImpl(Fitter&& f, DirectFitter&& df)
      : fitter(std::move(f)), directFitter(std::move(df)) {}

  auto makeGsfOptions(
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options)
      const {
    Acts::Experimental::GsfExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &updater);

    Acts::Experimental::GsfOptions<Acts::VectorMultiTrajectory> gsfOptions{
        options.geoContext,
        options.magFieldContext,
        options.calibrationContext,
        extensions,
        options.logger,
        options.propOptions,
        &(*options.referenceSurface),
        maxComponents,
        abortOnError,
        disableAllMaterialHandling};
    gsfOptions.extensions.calibrator
        .template connect<&ActsExamples::MeasurementCalibrator::calibrate>(
            &options.calibrator.get());

    return gsfOptions;
  }

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsExamples::IndexSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options,
      std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const override {
    const auto gsfOptions = makeGsfOptions(options);
    return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                      gsfOptions, trajectory);
  }

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsExamples::IndexSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options,
      const std::vector<const Acts::Surface*>& surfaceSequence,
      std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const override {
    const auto gsfOptions = makeGsfOptions(options);
    return directFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                            initialParameters, gsfOptions, surfaceSequence,
                            trajectory);
  }
};

}  // namespace

std::shared_ptr<TrackFittingAlgorithm::TrackFitterFunction>
ActsExamples::makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    std::size_t maxComponents, bool abortOnError,
    bool disableAllMaterialHandling) {
  MultiStepper stepper(std::move(magneticField));

  // Standard fitter
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator));

  // Direct fitter
  Acts::DirectNavigator directNavigator;
  DirectPropagator directPropagator(stepper, directNavigator);
  DirectFitter directTrackFitter(std::move(directPropagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction = std::make_shared<GsfFitterFunctionImpl>(
      std::move(trackFitter), std::move(directTrackFitter));
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;

  return fitterFunction;
}
