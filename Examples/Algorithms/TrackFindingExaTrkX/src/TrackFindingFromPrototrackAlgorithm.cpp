// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingFromPrototrackAlgorithm.hpp"

#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"

#include <algorithm>
#include <ranges>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace {

using namespace ActsExamples;

struct ProtoTrackSourceLinkAccessor
    : GeometryIdMultisetAccessor<IndexSourceLink> {
  using BaseIterator = GeometryIdMultisetAccessor<IndexSourceLink>::Iterator;
  using Iterator = Acts::SourceLinkAdapterIterator<BaseIterator>;

  std::unique_ptr<const Acts::Logger> loggerPtr;
  Container protoTrackSourceLinks;
  bool onlyPrototrackMeasurements = false;

  // get the range of elements with requested geoId
  std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
    const auto& logger = *loggerPtr;

    ACTS_VERBOSE("Prototrack only mode? " << std::boolalpha
                                          << onlyPrototrackMeasurements);
    if (protoTrackSourceLinks.contains(surface.geometryId())) {
      auto [begin, end] =
          protoTrackSourceLinks.equal_range(surface.geometryId());
      ACTS_VERBOSE("Select " << std::distance(begin, end)
                             << " source-links from prototrack on "
                             << surface.geometryId());
      return {Iterator{begin}, Iterator{end}};
    } else {
      ACTS_VERBOSE("No source-links found on " << surface.geometryId());
    }

    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());

    if (onlyPrototrackMeasurements) {
      ACTS_VERBOSE("Return empty range");
      return {Iterator{begin}, Iterator{begin}};
    } else {
      ACTS_VERBOSE("Select " << std::distance(begin, end)
                             << " source-links from collection on "
                             << surface.geometryId());
      return {Iterator{begin}, Iterator{end}};
    }
  }
};

}  // namespace

namespace ActsExamples {

TrackFindingFromPrototrackAlgorithm::TrackFindingFromPrototrackAlgorithm(
    Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm(cfg.tag + "CkfFromProtoTracks", lvl), m_cfg(cfg) {
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode TrackFindingFromPrototrackAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  using Extrapolator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
  using ExtrapolatorOptions = Extrapolator::template Options<
      Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;

  Extrapolator extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  ExtrapolatorOptions extrapolationOptions(ctx.geoContext, ctx.magFieldContext);

  if (initialParameters.size() != protoTracks.size()) {
    ACTS_FATAL("Inconsistent number of parameters and prototracks");
    return ProcessCode::ABORT;
  }

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::PropagatorPlainOptions pOptions(ctx.geoContext, ctx.magFieldContext);
  pOptions.maxSteps = 500;

  PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);
  Acts::GainMatrixUpdater kfUpdater;
  Acts::MeasurementSelector measSel{m_cfg.measurementSelectorCfg};

  // The source link accessor
  ProtoTrackSourceLinkAccessor sourceLinkAccessor;
  sourceLinkAccessor.loggerPtr = logger().clone("SourceLinkAccessor");
  sourceLinkAccessor.onlyPrototrackMeasurements =
      m_cfg.onlyPrototrackMeasurements;
  sourceLinkAccessor.container = &measurements.orderedIndices();

  using TrackStateCreatorType =
      Acts::TrackStateCreator<IndexSourceLinkAccessor::Iterator,
                              TrackContainer>;
  TrackStateCreatorType trackStateCreator;
  trackStateCreator.sourceLinkAccessor
      .template connect<&ProtoTrackSourceLinkAccessor::range>(
          &sourceLinkAccessor);
  trackStateCreator.calibrator
      .connect<&MeasurementCalibratorAdapter::calibrate>(&calibrator);
  trackStateCreator.measurementSelector
      .connect<&Acts::MeasurementSelector::select<
          typename TrackContainer::TrackStateContainerBackend>>(&measSel);

  Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<
      typename TrackContainer::TrackStateContainerBackend>>(&kfUpdater);
  extensions.createTrackStates
      .template connect<&TrackStateCreatorType ::createTrackStates>(
          &trackStateCreator);

  // Set the CombinatorialKalmanFilter options
  TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, extensions,
      pOptions, pSurface.get());

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer tracks(trackContainer, trackStateContainer);

  tracks.addColumn<unsigned int>("trackGroup");
  Acts::ProxyAccessor<unsigned int> seedNumber("trackGroup");

  std::size_t nSeed = 0;
  std::size_t nFailedFit = 0, nFailedSmoothing = 0, nFailedExtrapolation = 0;

  std::size_t nMeasurementIncrease = 0;
  std::size_t nMeasurementConstant = 0;
  std::size_t nMeasurementDecrease = 0;

  std::vector<std::size_t> nTracksPerSeeds;
  nTracksPerSeeds.reserve(initialParameters.size());

  for (auto i = 0ul; i < initialParameters.size(); ++i) {
    if (m_cfg.pickSeed != -1 && static_cast<int>(i) != m_cfg.pickSeed) {
      continue;
    }
    sourceLinkAccessor.protoTrackSourceLinks.clear();

    // Fill the source links via their indices from the container
    ACTS_VERBOSE("Collect source links from proto track");
    for (const auto measIndex : protoTracks.at(i)) {
      ConstVariableBoundMeasurementProxy measurement =
          measurements.getMeasurement(measIndex);
      IndexSourceLink sourceLink(measurement.geometryId(), measIndex);
      sourceLinkAccessor.protoTrackSourceLinks.insert(sourceLink);
      ACTS_VERBOSE("- " << measIndex << " | " << measurement.geometryId());
    }

    auto rootBranch = tracks.makeTrack();
    auto result = (*m_cfg.findTracks)(initialParameters.at(i), options, tracks,
                                      rootBranch);
    nSeed++;

    if (!result.ok()) {
      nFailedFit++;
      ACTS_ERROR("Track finding failed for proto track " << i << " with error"
                                                         << result.error());
      continue;
    }

    auto& tracksForSeed = result.value();

    nTracksPerSeeds.push_back(tracksForSeed.size());

    bool anySmoothFailed = false, anyExtrapolationFailed = false;
    for (auto& track : tracksForSeed) {
      // Set the seed number, this number decrease by 1 since the seed number
      // has already been updated
      seedNumber(track) = nSeed - 1;

      auto smoothingResult = Acts::smoothTrack(ctx.geoContext, track, logger());
      if (!smoothingResult.ok()) {
        ACTS_ERROR("Smoothing for seed " << i << " failed with error "
                                         << smoothingResult.error().message());
        anySmoothFailed = true;
        continue;
      }

      auto extrapolationResult = Acts::extrapolateTrackToReferenceSurface(
          track, *pSurface, extrapolator, extrapolationOptions,
          Acts::TrackExtrapolationStrategy::firstOrLast, logger());
      if (!extrapolationResult.ok()) {
        ACTS_ERROR("Extrapolation for seed "
                   << i << " failed with error "
                   << extrapolationResult.error().message());
        ACTS_DEBUG("Assign start parameters instead");
        anyExtrapolationFailed = true;
        continue;
      }
      ACTS_VERBOSE("Prototrack extension through fit: "
                   << protoTracks.at(i).size() << " -> "
                   << track.nMeasurements());
      if (protoTracks.at(i).size() < track.nMeasurements()) {
        nMeasurementIncrease++;
      } else if (protoTracks.at(i).size() == track.nMeasurements()) {
        nMeasurementConstant++;
      } else {
        nMeasurementDecrease++;
      }
      if (track.nMeasurements() < protoTracks.at(i).size()) {
        for (auto mid : protoTracks.at(i)) {
          auto geoId = measurements.getMeasurement(mid).geometryId();
          auto found = std::find_if(
              track.trackStatesReversed().cbegin(),
              track.trackStatesReversed().cend(), [&](auto ts) {
                return ts.hasReferenceSurface()
                           ? ts.referenceSurface().geometryId() == geoId
                           : false;
              });
          if (found == track.trackStatesReversed().cend()) {
            ACTS_VERBOSE(" - " << geoId);
          } else {
            ACTS_VERBOSE(" + " << geoId);
          }
        }
      }
    }
    nFailedSmoothing += static_cast<std::size_t>(anySmoothFailed);
    nFailedExtrapolation += static_cast<std::size_t>(anyExtrapolationFailed);
  }

  ACTS_DEBUG("Tracks with measurement increase: "
             << nMeasurementIncrease << ", decrease: " << nMeasurementDecrease
             << ", constant: " << nMeasurementConstant);

  {
    std::lock_guard<std::mutex> guard(m_mutex);

    std::copy(nTracksPerSeeds.begin(), nTracksPerSeeds.end(),
              std::back_inserter(m_nTracksPerSeeds));
  }

  // TODO The computeSharedHits function is still a member function of
  // TrackFindingAlgorithm, but could also be a free function. Uncomment this
  // once this is done.
  // Compute shared hits from all the reconstructed tracks if
  // (m_cfg.computeSharedHits) {
  //   computeSharedHits(measurements, tracks);
  // }

  auto nFailed = nFailedFit + nFailedExtrapolation + nFailedSmoothing;
  ACTS_DEBUG("Event " << ctx.eventNumber << ": " << nFailed << " / " << nSeed
                      << " failed (" << ((100.f * nFailed) / nSeed) << "%)");
  ACTS_DEBUG("- fit failed: " << nFailedFit);
  ACTS_DEBUG("- smoothing failed: " << nFailedSmoothing);
  ACTS_DEBUG("- extrapolation failed: " << nFailedExtrapolation);
  ACTS_DEBUG("Finalized track finding with " << tracks.size()
                                             << " track candidates.");
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

ActsExamples::ProcessCode TrackFindingFromPrototrackAlgorithm::finalize() {
  assert(std::distance(m_nTracksPerSeeds.begin(), m_nTracksPerSeeds.end()) > 0);

  ACTS_INFO("TrackFindingFromPrototracksAlgorithm statistics:");
  namespace ba = boost::accumulators;
  using Accumulator = ba::accumulator_set<
      float, ba::features<ba::tag::sum, ba::tag::mean, ba::tag::variance>>;

  Accumulator totalAcc;
  std::ranges::for_each(m_nTracksPerSeeds,
                        [&](auto v) { totalAcc(static_cast<float>(v)); });
  ACTS_INFO("- total number tracks: " << ba::sum(totalAcc));
  ACTS_INFO("- avg tracks per seed: " << ba::mean(totalAcc) << " +- "
                                      << std::sqrt(ba::variance(totalAcc)));

  return {};
}

}  // namespace ActsExamples
