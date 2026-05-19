// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/TrackFitting/ReferenceTrajectoryBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <memory>
#include <utility>
#include <vector>

namespace {

using Stepper = Acts::SympyStepper;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using ReferenceTrajectoryBuilder =
    Acts::Experimental::ReferenceTrajectoryBuilder<Propagator,
                                                   Acts::VectorMultiTrajectory>;
using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
using DirectReferenceTrajectoryBuilder =
    Acts::Experimental::ReferenceTrajectoryBuilder<DirectPropagator,
                                                   Acts::VectorMultiTrajectory>;

using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;

using namespace ActsExamples;

struct KalmanReferenceTrajectoryFitterFunctionImpl final
    : public TrackFitterFunction {
  Propagator extrapolator;

  ReferenceTrajectoryBuilder referenceTrajectoryBuilder;
  DirectReferenceTrajectoryBuilder directReferenceTrajectoryBuilder;

  Acts::GainMatrixUpdater kfUpdater;
  Acts::MbfSmoother kfSmoother;

  bool multipleScattering = false;
  bool energyLoss = false;
  Acts::FreeToBoundCorrection freeToBoundCorrection;

  IndexSourceLink::SurfaceAccessor slSurfaceAccessor;

  std::unique_ptr<const Acts::Logger> m_smootherLogger;
  std::unique_ptr<const Acts::Logger> m_extrapolationLogger;

  KalmanReferenceTrajectoryFitterFunctionImpl(
      Propagator&& ext, ReferenceTrajectoryBuilder&& rtBuilder,
      DirectReferenceTrajectoryBuilder&& dRtBuilder,
      const Acts::TrackingGeometry& trkGeo, const Acts::Logger& logger)
      : extrapolator(std::move(ext)),
        referenceTrajectoryBuilder(std::move(rtBuilder)),
        directReferenceTrajectoryBuilder(std::move(dRtBuilder)),
        slSurfaceAccessor{trkGeo} {
    m_smootherLogger = logger.cloneWithSuffix("Smoother");
    m_extrapolationLogger = logger.cloneWithSuffix("Extrapolation");
  }

  auto makeReferenceTrajectoryBuilderOptions(
      const GeneralFitterOptions& options) const {
    Acts::Experimental::ReferenceTrajectoryBuilderOptions<
        Acts::VectorMultiTrajectory>
        rtOptions(options.geoContext, options.magFieldContext,
                  options.propOptions, nullptr);

    rtOptions.multipleScattering = multipleScattering;
    rtOptions.energyLoss = energyLoss;
    rtOptions.freeToBoundCorrection = freeToBoundCorrection;

    return rtOptions;
  }

  TrackFitterResult operator()(const std::vector<Acts::SourceLink>& sourceLinks,
                               const TrackParameters& initialParameters,
                               const GeneralFitterOptions& options,
                               const MeasurementCalibratorAdapter& calibrator,
                               TrackContainer& tracks) const override {
    const auto referenceTrajectoryBuilderOptions =
        makeReferenceTrajectoryBuilderOptions(options);

    Acts::SourceLinkSurfaceAccessor sourceLinkAccessorDelegate;
    sourceLinkAccessorDelegate
        .connect<&IndexSourceLink::SurfaceAccessor::operator()>(
            &slSurfaceAccessor);

    ReferenceTrajectoryBuilder::Calibrator calibratorDelegate;
    calibratorDelegate.connect<&MeasurementCalibratorAdapter::calibrate>(
        &calibrator);

    ReferenceTrajectoryBuilder::Updater updaterDelegate;
    updaterDelegate.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);

    auto buildResult = referenceTrajectoryBuilder.build(
        initialParameters, referenceTrajectoryBuilderOptions, tracks);
    if (!buildResult.ok()) {
      return buildResult.error();
    }
    auto track = *buildResult;

    referenceTrajectoryBuilder.attachSourceLinks(track, sourceLinks,
                                                 sourceLinkAccessorDelegate);

    referenceTrajectoryBuilder.calibrateMeasurements(options.geoContext,
                                                     options.calibrationContext,
                                                     track, calibratorDelegate);

    auto fitterResult = referenceTrajectoryBuilder.filter(
        options.geoContext, track, updaterDelegate);
    if (!fitterResult.ok()) {
      return fitterResult.error();
    }

    auto smoothRes =
        kfSmoother(options.geoContext, tracks.trackStateContainer(),
                   track.tipIndex(), *m_smootherLogger);
    if (!smoothRes.ok()) {
      return smoothRes.error();
    }

    if (!track.hasReferenceSurface() && options.referenceSurface != nullptr) {
      typename Propagator::template Options<> extrapolationOptions(
          options.geoContext, options.magFieldContext);
      auto extrapolationResult = extrapolateTrackToReferenceSurface(
          track, *options.referenceSurface, extrapolator, extrapolationOptions,
          Acts::TrackExtrapolationStrategy::first, *m_extrapolationLogger);

      if (!extrapolationResult.ok()) {
        return extrapolationResult.error();
      }
    }

    Acts::calculateTrackQuantities(track);

    return TrackFitterResult::success(track);
  }

  TrackFitterResult operator()(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const TrackParameters& initialParameters,
      const GeneralFitterOptions& options,
      const RefittingCalibrator& calibrator,
      const std::vector<const Acts::Surface*>& surfaceSequence,
      TrackContainer& tracks) const override {
    const auto referenceTrajectoryBuilderOptions =
        makeReferenceTrajectoryBuilderOptions(options);

    Acts::SourceLinkSurfaceAccessor sourceLinkAccessorDelegate;
    sourceLinkAccessorDelegate.connect<&RefittingCalibrator::accessSurface>();

    DirectReferenceTrajectoryBuilder::Calibrator calibratorDelegate;
    calibratorDelegate.connect<&RefittingCalibrator::calibrate>(&calibrator);

    DirectReferenceTrajectoryBuilder::Updater updaterDelegate;
    updaterDelegate.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);

    auto buildResult = directReferenceTrajectoryBuilder.build(
        initialParameters, referenceTrajectoryBuilderOptions, surfaceSequence,
        tracks);
    if (!buildResult.ok()) {
      return buildResult.error();
    }
    auto track = *buildResult;

    directReferenceTrajectoryBuilder.attachSourceLinks(
        track, sourceLinks, sourceLinkAccessorDelegate);

    directReferenceTrajectoryBuilder.calibrateMeasurements(
        options.geoContext, options.calibrationContext, track,
        calibratorDelegate);

    directReferenceTrajectoryBuilder.filter(options.geoContext, track,
                                            updaterDelegate);

    auto smoothRes = kfSmoother(options.geoContext,
                                tracks.trackStateContainer(), track.tipIndex());
    if (!smoothRes.ok()) {
      return smoothRes.error();
    }

    if (!track.hasReferenceSurface() && options.referenceSurface != nullptr) {
      typename Propagator::template Options<> extrapolationOptions(
          options.geoContext, options.magFieldContext);
      auto extrapolationResult = extrapolateTrackToReferenceSurface(
          track, *options.referenceSurface, extrapolator, extrapolationOptions,
          Acts::TrackExtrapolationStrategy::first);

      if (!extrapolationResult.ok()) {
        return extrapolationResult.error();
      }
    }

    Acts::calculateTrackQuantities(track);

    return TrackFitterResult::success(track);
  }
};

}  // namespace

std::shared_ptr<TrackFitterFunction>
ActsExamples::makeKalmanReferenceTrajectoryFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    const Acts::FreeToBoundCorrection& freeToBoundCorrection,
    bool useJosephFormulation, const Acts::Logger& logger) {
  const Stepper stepper(std::move(magneticField));

  const auto& geo = *trackingGeometry;
  Acts::Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(stepper, navigator,
                        logger.cloneWithSuffix("Propagator"));
  ReferenceTrajectoryBuilder referenceTrajectoryBuilder(
      std::move(propagator),
      logger.cloneWithSuffix("ReferenceTrajectoryBuilder"));

  Acts::DirectNavigator directNavigator{
      logger.cloneWithSuffix("DirectNavigator")};
  DirectPropagator directPropagator(stepper, std::move(directNavigator),
                                    logger.cloneWithSuffix("DirectPropagator"));
  DirectReferenceTrajectoryBuilder directReferenceTrajectoryBuilder(
      std::move(directPropagator),
      logger.cloneWithSuffix("DirectReferenceTrajectoryBuilder"));

  Propagator extrapolator(stepper, navigator,
                          logger.cloneWithSuffix("Extrapolator"));

  auto fitterFunction =
      std::make_shared<KalmanReferenceTrajectoryFitterFunctionImpl>(
          std::move(extrapolator), std::move(referenceTrajectoryBuilder),
          std::move(directReferenceTrajectoryBuilder), geo, logger);
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;
  fitterFunction->kfUpdater = Acts::GainMatrixUpdater(useJosephFormulation);

  return fitterFunction;
}
