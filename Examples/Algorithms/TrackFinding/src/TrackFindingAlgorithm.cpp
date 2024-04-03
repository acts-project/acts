// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
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

#include <boost/histogram.hpp>

ActsExamples::TrackFindingAlgorithm::TrackFindingAlgorithm(
    Config config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFindingAlgorithm", level),
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

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithm::execute(
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
  Acts::MeasurementSelector measSel{m_cfg.measurementSelectorCfg};

  Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
      extensions;
  extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.measurementSelector
      .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
          &measSel);

  IndexSourceLinkAccessor slAccessor;
  slAccessor.container = &sourceLinks;
  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

  Acts::PropagatorPlainOptions firstPropOptions;
  firstPropOptions.maxSteps = m_cfg.maxSteps;
  firstPropOptions.direction = Acts::Direction::Forward;

  Acts::PropagatorPlainOptions secondPropOptions;
  secondPropOptions.maxSteps = m_cfg.maxSteps;
  secondPropOptions.direction = firstPropOptions.direction.invert();

  // Set the CombinatorialKalmanFilter options
  ActsExamples::TrackFindingAlgorithm::TrackFinderOptions firstOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, slAccessorDelegate,
      extensions, firstPropOptions);

  ActsExamples::TrackFindingAlgorithm::TrackFinderOptions secondOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, slAccessorDelegate,
      extensions, secondPropOptions);
  secondOptions.targetSurface = pSurface.get();

  Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  Acts::PropagatorOptions<Acts::ActionList<Acts::MaterialInteractor>,
                          Acts::AbortList<Acts::EndOfWorldReached>>
      extrapolationOptions(ctx.geoContext, ctx.magFieldContext);

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

  for (std::size_t iSeed = 0; iSeed < initialParameters.size(); ++iSeed) {
    // Clear trackContainerTemp and trackStateContainerTemp
    tracksTemp.clear();

    const Acts::BoundTrackParameters& firstInitialParameters =
        initialParameters.at(iSeed);

    auto firstResult =
        (*m_cfg.findTracks)(firstInitialParameters, firstOptions, tracksTemp);
    m_nTotalSeeds++;
    nSeed++;

    if (!firstResult.ok()) {
      m_nFailedSeeds++;
      ACTS_WARNING("Track finding failed for seed " << iSeed << " with error"
                                                    << firstResult.error());
      continue;
    }

    auto& firstTracksForSeed = firstResult.value();
    for (auto& firstTrack : firstTracksForSeed) {
      // TODO select track

      auto firstSmoothingResult =
          Acts::smoothTrack(ctx.geoContext, firstTrack, logger());
      if (!firstSmoothingResult.ok()) {
        m_nFailedSmoothing++;
        ACTS_ERROR("First smoothing for seed "
                   << iSeed << " and track " << firstTrack.index()
                   << " failed with error " << firstSmoothingResult.error());
        continue;
      }

      std::size_t nSecond = 0;

      // Set the seed number, this number decrease by 1 since the seed number
      // has already been updated
      seedNumber(firstTrack) = nSeed - 1;

      if (m_cfg.twoWay) {
        std::optional<Acts::VectorMultiTrajectory::TrackStateProxy> firstState;
        for (auto st : firstTrack.trackStatesReversed()) {
          bool isMeasurement =
              st.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag);
          bool isOutlier =
              st.typeFlags().test(Acts::TrackStateFlag::OutlierFlag);
          // We are excluding non measurement states and outlier here. Those can
          // decrease resolution because only the smoothing corrected the very
          // first prediction as filtering is not possible.
          if (isMeasurement && !isOutlier) {
            firstState = st;
          }
        }

        if (firstState.has_value()) {
          Acts::BoundTrackParameters secondInitialParameters(
              firstState->referenceSurface().getSharedPtr(),
              firstState->parameters(), firstState->covariance(),
              firstInitialParameters.particleHypothesis());

          auto secondResult = (*m_cfg.findTracks)(secondInitialParameters,
                                                  secondOptions, tracksTemp);

          if (!secondResult.ok()) {
            ACTS_WARNING("Second track finding failed for seed "
                         << iSeed << " with error" << secondResult.error());
          } else {
            auto firstFirstState =
                std::next(firstTrack.trackStatesReversed().begin(),
                          firstTrack.nTrackStates() - 1);

            auto& secondTracksForSeed = secondResult.value();
            for (auto& secondTrack : secondTracksForSeed) {
              // TODO select track
              if (secondTrack.nTrackStates() < 2) {
                continue;
              }

              secondTrack.reverseTrackStates(true);
              seedNumber(secondTrack) = nSeed - 1;

              (*firstFirstState).previous() =
                  (*std::next(secondTrack.trackStatesReversed().begin()))
                      .index();
              secondTrack.tipIndex() = firstTrack.tipIndex();

              Acts::calculateTrackQuantities(secondTrack);

              // TODO this extrapolation should not be necessary as the CKF will
              // target this surface anyway and could communicate the parameters
              auto secondExtrapolationResult =
                  Acts::extrapolateTrackToReferenceSurface(
                      secondTrack, *pSurface, extrapolator,
                      extrapolationOptions, m_cfg.extrapolationStrategy,
                      logger());
              if (!secondExtrapolationResult.ok()) {
                m_nFailedExtrapolation++;
                ACTS_ERROR("Second extrapolation for seed "
                           << iSeed << " and track " << secondTrack.index()
                           << " failed with error "
                           << secondExtrapolationResult.error());
                continue;
              }

              if (!m_trackSelector.has_value() ||
                  m_trackSelector->isValidTrack(secondTrack)) {
                auto destProxy = tracks.makeTrack();
                destProxy.copyFrom(secondTrack, true);
              }

              ++nSecond;
            }
          }
        }
      }

      if (nSecond == 0) {
        auto firstExtrapolationResult =
            Acts::extrapolateTrackToReferenceSurface(
                firstTrack, *pSurface, extrapolator, extrapolationOptions,
                m_cfg.extrapolationStrategy, logger());
        if (!firstExtrapolationResult.ok()) {
          m_nFailedExtrapolation++;
          ACTS_ERROR("Extrapolation for seed "
                     << iSeed << " and track " << firstTrack.index()
                     << " failed with error "
                     << firstExtrapolationResult.error());
          continue;
        }

        if (!m_trackSelector.has_value() ||
            m_trackSelector->isValidTrack(firstTrack)) {
          auto destProxy = tracks.makeTrack();
          destProxy.copyFrom(firstTrack, true);
        }
      }
    }
  }

  // Compute shared hits from all the reconstructed tracks
  if (m_cfg.computeSharedHits) {
    computeSharedHits(sourceLinks, tracks);
  }

  ACTS_DEBUG("Finalized track finding with " << tracks.size()
                                             << " track candidates.");

  m_memoryStatistics.local().hist +=
      tracks.trackStateContainer().statistics().hist;

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

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithm::finalize() {
  ACTS_INFO("TrackFindingAlgorithm statistics:");
  ACTS_INFO("- total seeds: " << m_nTotalSeeds);
  ACTS_INFO("- failed seeds: " << m_nFailedSeeds);
  ACTS_INFO("- failed smoothing: " << m_nFailedSmoothing);
  ACTS_INFO("- failed extrapolation: " << m_nFailedExtrapolation);
  ACTS_INFO("- failure ratio seeds: " << static_cast<double>(m_nFailedSeeds) /
                                             m_nTotalSeeds);

  auto memoryStatistics =
      m_memoryStatistics.combine([](const auto& a, const auto& b) {
        Acts::VectorMultiTrajectory::Statistics c;
        c.hist = a.hist + b.hist;
        return c;
      });
  std::stringstream ss;
  memoryStatistics.toStream(ss);
  ACTS_DEBUG("Track State memory statistics (averaged):\n" << ss.str());
  return ProcessCode::SUCCESS;
}
