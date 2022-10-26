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
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <filesystem>

using namespace ActsExamples;

namespace {

using MultiStepper = Acts::MultiEigenStepperLoop<>;
using Propagator = Acts::Propagator<MultiStepper, Acts::Navigator>;
using DirectPropagator = Acts::Propagator<MultiStepper, Acts::DirectNavigator>;

#if USE_CUSTOM_BETHE_HEITLER
using BHApprox = Acts::BetheHeitlerSimulatedAnnealingMinimizer
#else
using BHApprox = Acts::AtlasBetheHeitlerApprox<6, 5>;
#endif

using Fitter = Acts::GaussianSumFitter<Propagator, BHApprox, Acts::VectorMultiTrajectory>;
using DirectFitter =
    Acts::GaussianSumFitter<DirectPropagator, BHApprox, Acts::VectorMultiTrajectory>;

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
    Acts::GsfExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &updater);

    Acts::GsfOptions<Acts::VectorMultiTrajectory> gsfOptions{
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
    std::string lowBetheHeitlerPath, std::string highBetheHeitlerPath,
    std::size_t maxComponents, bool abortOnError,
    bool disableAllMaterialHandling) {
  MultiStepper stepper(std::move(magneticField));

#if USE_CUSTOM_BETHE_HEITLER
  constexpr std::size_t NComponents = 12;

  auto gen = std::make_shared<std::mt19937>(23465u);

  const auto iterations = 3000;
  const auto temperatures = []() {
    std::vector<double> t;
    for (int i = 0; i < iterations; ++i) {
      t.push_back(1 * std::exp(-0.0001 * i));
    }
    return t;
  }();

  auto next = [](std::array<double, 3 * NComponents> ps,
                 std::mt19937& generator) {
    auto val_dist = std::uniform_real_distribution{-0.5, 0.5};

    auto idx = std::uniform_int_distribution(0ul, ps.size())(generator);
    ps[idx] += ps[idx] * val_dist(generator);

    return ps;
  };

   const auto startValue = std::array<Acts::detail::GaussianComponent, NComponents>{{
        {1, 0.5, 1.e-5},
        {2, 0.99, 1.e-5},
        {2, 0.99, 1.e-5},
        {2, 0.99, 1.e-5},
        {2,.99, 1.e-5},
        {2,.99, 1.e-5},
        {2, 0.99, 1.e-5},
        {2,.99, 1.e-5},
        {2,.99, 1.e-5},
        {2, 0.99, 1.e-5},
        {2,.99, 1.e-5},
        {2,.99, 1.e-5},
    }};

  auto bhapp = Acts::BetheHeitlerSimulatedAnnealingMinimizer(
      temperatures, startValue, gen, next);
#else
  auto makeBehteHeitlerApprox = [&]() {
    if (std::filesystem::exists(lowBetheHeitlerPath) &&
        std::filesystem::exists(highBetheHeitlerPath)) {
      return Acts::AtlasBetheHeitlerApprox<6, 5>::loadFromFile(
          lowBetheHeitlerPath, highBetheHeitlerPath);
    } else {
      std::cout
          << "WARNING: Could not find files, use standard configuration\n";
      return Acts::AtlasBetheHeitlerApprox<6, 5>(Acts::bh_cdf_cmps6_order5_data,
                                                 Acts::bh_cdf_cmps6_order5_data,
                                                 true, true);
    }
  };

  auto bhapp = makeBehteHeitlerApprox();
#endif

  // Standard fitter
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator), std::move(bhapp));

  // Direct fitter
  Acts::DirectNavigator directNavigator;
  DirectPropagator directPropagator(stepper, directNavigator);
  DirectFitter directTrackFitter(std::move(directPropagator), std::move(bhapp));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction = std::make_shared<GsfFitterFunctionImpl>(
      std::move(trackFitter), std::move(directTrackFitter));
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;

  return fitterFunction;
}
