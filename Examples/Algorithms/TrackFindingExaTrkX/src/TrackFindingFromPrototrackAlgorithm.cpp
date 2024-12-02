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
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"

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

    if (protoTrackSourceLinks.contains(surface.geometryId())) {
      auto [begin, end] =
          protoTrackSourceLinks.equal_range(surface.geometryId());
      ACTS_VERBOSE("Select " << std::distance(begin, end)
                             << " source-links from prototrack on "
                             << surface.geometryId());
      return {Iterator{begin}, Iterator{end}};
    }

    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());
    ACTS_VERBOSE("Select " << std::distance(begin, end)
                           << " source-links from collection on "
                           << surface.geometryId());
    if (onlyPrototrackMeasurements) {
      return {Iterator{begin}, Iterator{begin}};
    } else {
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
  Acts::GainMatrixSmoother kfSmoother;
  Acts::MeasurementSelector measSel{m_cfg.measurementSelectorCfg};

  Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
  extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<
      typename TrackContainer::TrackStateContainerBackend>>(&kfUpdater);
  extensions.measurementSelector.connect<&Acts::MeasurementSelector::select<
      typename TrackContainer::TrackStateContainerBackend>>(&measSel);

  // The source link accessor
  ProtoTrackSourceLinkAccessor sourceLinkAccessor;
  sourceLinkAccessor.loggerPtr = logger().clone("SourceLinkAccessor");
  sourceLinkAccessor.onlyPrototrackMeasurements =
      m_cfg.onlyPrototrackMeasurements;
  sourceLinkAccessor.container = &measurements.orderedIndices();

  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&ProtoTrackSourceLinkAccessor::range>(
      &sourceLinkAccessor);

  // Set the CombinatorialKalmanFilter options
  TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, slAccessorDelegate,
      extensions, pOptions, &(*pSurface));

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer tracks(trackContainer, trackStateContainer);

  tracks.addColumn<unsigned int>("trackGroup");
  Acts::ProxyAccessor<unsigned int> seedNumber("trackGroup");

  std::size_t nSeed = 0;
  std::size_t nFailed = 0;

  std::size_t nMeasurementIncrease = 0;
  std::size_t nMeasurementConstant = 0;
  std::size_t nMeasurementDecrease = 0;

  std::vector<std::size_t> nTracksPerSeeds;
  nTracksPerSeeds.reserve(initialParameters.size());

  for (auto i = 0ul; i < initialParameters.size(); ++i) {
    sourceLinkAccessor.protoTrackSourceLinks.clear();

    // Fill the source links via their indices from the container
    for (const auto hitIndex : protoTracks.at(i)) {
      if (auto it = measurements.orderedIndices().nth(hitIndex);
          it != measurements.orderedIndices().end()) {
        sourceLinkAccessor.protoTrackSourceLinks.insert(*it);
      } else {
        ACTS_FATAL("Proto track " << i << " contains invalid hit index"
                                  << hitIndex);
        return ProcessCode::ABORT;
      }
    }

    auto rootBranch = tracks.makeTrack();
    auto result = (*m_cfg.findTracks)(initialParameters.at(i), options, tracks,
                                      rootBranch);
    nSeed++;

    if (!result.ok()) {
      nFailed++;
      ACTS_WARNING("Track finding failed for proto track " << i << " with error"
                                                           << result.error());
      continue;
    }

    auto& tracksForSeed = result.value();

    nTracksPerSeeds.push_back(tracksForSeed.size());

    for (auto& track : tracksForSeed) {
      // Set the seed number, this number decrease by 1 since the seed number
      // has already been updated
      seedNumber(track) = nSeed - 1;

#if 1
      track.parameters() = initialParameters.at(i).parameters();
      if (initialParameters.at(i).covariance()) {
        track.covariance() = *initialParameters.at(i).covariance();
      }
      track.setReferenceSurface(
          initialParameters.at(i).referenceSurface().getSharedPtr());
#else
      auto smoothingResult = Acts::smoothTrack(ctx.geoContext, track, logger());
      if (!smoothingResult.ok()) {
        ACTS_ERROR("Smoothing for seed " << i << " failed with error "
                                         << smoothingResult.error().message());
        continue;
      }

      auto extrapolationResult = Acts::extrapolateTrackToReferenceSurface(
          track, *pSurface, extrapolator, extrapolationOptions,
          Acts::TrackExtrapolationStrategy::firstOrLast, logger());
      if (!extrapolationResult.ok()) {
        ACTS_WARNING("Extrapolation for seed "
                     << i << " failed with error "
                     << extrapolationResult.error().message());
        ACTS_WARNING("Assign start parameters instead");
        continue;
      }
#endif
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
    }
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

  ACTS_INFO("Event " << ctx.eventNumber << ": " << nFailed << " / " << nSeed
                     << " failed (" << ((100.f * nFailed) / nSeed) << "%)");
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
  std::for_each(m_nTracksPerSeeds.begin(), m_nTracksPerSeeds.end(),
                [&](auto v) { totalAcc(static_cast<float>(v)); });
  ACTS_INFO("- total number tracks: " << ba::sum(totalAcc));
  ACTS_INFO("- avg tracks per seed: " << ba::mean(totalAcc) << " +- "
                                      << std::sqrt(ba::variance(totalAcc)));

  return {};
}

}  // namespace ActsExamples
