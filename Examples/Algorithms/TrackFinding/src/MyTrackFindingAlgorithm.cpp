// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MyTrackFindingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <utility>

ActsExamples::MyTrackFindingAlgorithm::MyTrackFindingAlgorithm(
    Config config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("MyTrackFindingAlgorithm", level),
      m_cfg(std::move(config)),
      m_trackSelector(
          m_cfg.trackSelectorCfg.has_value()
              ? std::visit(
                    [](const auto& cfg) -> std::optional<Acts::TrackSelector> {
                      return {cfg};
                    },
                    m_cfg.trackSelectorCfg.value())
              : std::nullopt) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing initial track parameters input collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing tracks output collection");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode ActsExamples::MyTrackFindingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& sourceLinks = m_inputSourceLinks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  Acts::MeasurementSelector measSel{m_cfg.measurementSelectorCfg};

  Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
      extensions;
  extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);
  extensions.measurementSelector
      .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
          &measSel);

  IndexSourceLinkAccessor slAccessor;
  slAccessor.container = &sourceLinks;
  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

  using Propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
  using CKF =
      Acts::CombinatorialKalmanFilter<Propagator, Acts::VectorMultiTrajectory>;
  using Options =
      Acts::CombinatorialKalmanFilterOptions<IndexSourceLinkAccessor::Iterator,
                                             Acts::VectorMultiTrajectory>;

  Acts::EigenStepper<> stepper(m_cfg.magneticField);
  Acts::Navigator navigator({m_cfg.trackingGeometry},
                            logger().cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger().cloneWithSuffix("Propagator"));
  CKF ckf(std::move(propagator), logger().cloneWithSuffix("Finder"));

  Acts::PropagatorPlainOptions fwdPropOptions;
  fwdPropOptions.maxSteps = m_cfg.maxSteps;
  fwdPropOptions.direction = Acts::Direction::Forward;

  Options fwdOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
                     slAccessorDelegate, extensions, fwdPropOptions,
                     pSurface.get());
  fwdOptions.filterTargetSurface = nullptr;
  fwdOptions.smoothing = true;
  fwdOptions.smoothingTargetSurface = pSurface.get();
  fwdOptions.smoothingTargetSurfaceStrategy =
      Acts::CombinatorialKalmanFilterTargetSurfaceStrategy::first;

  Acts::PropagatorPlainOptions bwdPropOptions;
  bwdPropOptions.maxSteps = m_cfg.maxSteps;
  bwdPropOptions.direction = Acts::Direction::Backward;

  Options bwdOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
                     slAccessorDelegate, extensions, bwdPropOptions,
                     pSurface.get());
  bwdOptions.filterTargetSurface = pSurface.get();
  bwdOptions.smoothing = true;  // TODO this should not be necessary
  bwdOptions.smoothingTargetSurface = pSurface.get();
  bwdOptions.smoothingTargetSurfaceStrategy =
      Acts::CombinatorialKalmanFilterTargetSurfaceStrategy::last;

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

  auto trackContainerTemp = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainerTemp =
      std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer tracks(trackContainer, trackStateContainer);
  TrackContainer tracksTemp(trackContainerTemp, trackStateContainerTemp);

  tracks.addColumn<unsigned int>("trackGroup");
  tracksTemp.addColumn<unsigned int>("trackGroup");
  Acts::ProxyAccessor<unsigned int> seedNumber("trackGroup");

  unsigned int nSeed = 0;

  for (std::size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    // Clear trackContainerTemp and trackStateContainerTemp
    tracksTemp.clear();

    auto fwdResult =
        ckf.findTracks(initialParameters.at(iseed), fwdOptions, tracksTemp);
    nSeed++;

    if (!fwdResult.ok()) {
      ACTS_WARNING("Forward track finding failed for seed "
                   << iseed << " with error" << fwdResult.error());
      continue;
    }

    auto& fwdTracksForSeed = fwdResult.value();
    for (auto& fwdTrack : fwdTracksForSeed) {
      std::optional<Acts::VectorMultiTrajectory::TrackStateProxy> firstState;
      for (auto st : fwdTrack.trackStatesReversed()) {
        bool isMeasurement =
            st.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag);
        bool isOutlier = st.typeFlags().test(Acts::TrackStateFlag::OutlierFlag);
        // We are excluding non measurement states and outlier here. Those can
        // decrease resolution because only the smoothing corrected the very
        // first prediction as filtering is not possible.
        if (isMeasurement && !isOutlier) {
          firstState = st;
        }
      }

      Acts::BoundTrackParameters bwdInitialParameters(
          firstState->referenceSurface().getSharedPtr(),
          firstState->parameters(), firstState->covariance(),
          initialParameters.at(iseed).particleHypothesis());

      auto bwdResult =
          ckf.findTracks(bwdInitialParameters, bwdOptions, tracksTemp);

      if (!bwdResult.ok()) {
        ACTS_WARNING("Backward track finding failed for seed "
                     << iseed << " with error" << bwdResult.error());
        continue;
      }

      auto& bwdTracks = bwdResult.value();
      for (auto& bwdTrack : bwdTracks) {
        if (bwdTrack.nTrackStates() < 2) {
          continue;
        }

        bwdTrack.reverseTrackStates(true);
        // Set the seed number, this number decrease by 1 since the seed number
        // has already been updated
        seedNumber(bwdTrack) = nSeed - 1;

        auto innermostFwdState = fwdTrack.trackStatesReversed().begin();
        for (auto i = innermostFwdState;
             ++i != fwdTrack.trackStatesReversed().end();
             i = ++innermostFwdState) {
        }

        (*innermostFwdState).previous() =
            (*std::next(bwdTrack.trackStatesReversed().begin())).index();
        bwdTrack.tipIndex() = fwdTrack.tipIndex();

        Acts::calculateTrackQuantities(bwdTrack);

        if (!m_trackSelector.has_value() ||
            m_trackSelector->isValidTrack(bwdTrack)) {
          auto destProxy = tracks.getTrack(tracks.addTrack());
          destProxy.copyFrom(bwdTrack, true);
        }
      }
    }
  }

  ACTS_DEBUG("Track finding with " << tracks.size() << " track candidates.");

  auto constTrackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer));

  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
      std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTracks(ctx, std::move(constTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
