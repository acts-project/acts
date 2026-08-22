// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingGnn/TrackFindingFromProtoTracksAlgorithm.hpp"

#include "Acts/EventData/ProxyAccessor.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/TrackFinding/TrackStateCreator.hpp"

#include <algorithm>
#include <ranges>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace Acts;

namespace {

using namespace ActsExamples;

/// Creates track states from the measurements of the current proto track,
/// falling back to all measurements on surfaces the proto track does not
/// touch.
class ProtoTrackStateCreator final
    : public ActsExamples::TrackStateCreatorBase<ProtoTrackStateCreator> {
 public:
  using Container = GeometryIdMultiset<IndexSourceLink>;

  std::unique_ptr<const Logger> loggerPtr;
  const MeasurementContainer* measurements = nullptr;
  Container protoTrackSourceLinks;

  auto measurementRange(const State& state) const {
    const auto& logger = *loggerPtr;

    if (protoTrackSourceLinks.contains(state.surface->geometryId())) {
      auto [begin, end] =
          protoTrackSourceLinks.equal_range(state.surface->geometryId());
      ACTS_VERBOSE("Select " << std::distance(begin, end)
                             << " source-links from proto track on "
                             << state.surface->geometryId());
      return std::ranges::subrange(begin, end);
    }

    assert(measurements != nullptr);
    auto [begin, end] =
        measurements->orderedIndices().equal_range(state.surface->geometryId());
    ACTS_VERBOSE("Select " << std::distance(begin, end)
                           << " source-links from collection on "
                           << state.surface->geometryId());
    return std::ranges::subrange(begin, end);
  }

  ConstVariableBoundMeasurementProxy calibrate(
      const State& /*state*/, const IndexSourceLink& measurement) const {
    return measurements->getMeasurement(measurement.index());
  }
};

}  // namespace

namespace ActsExamples {

TrackFindingFromProtoTracksAlgorithm::TrackFindingFromProtoTracksAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm(cfg.tag + "CkfFromProtoTracks", std::move(logger)),
      m_cfg(cfg) {
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode TrackFindingFromProtoTracksAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  if (initialParameters.size() != protoTracks.size()) {
    ACTS_FATAL("Inconsistent number of parameters and proto tracks");
    return ProcessCode::ABORT;
  }

  // Construct a perigee surface as the target surface
  auto pSurface = Surface::makeShared<PerigeeSurface>(Vector3{0., 0., 0.});

  PropagatorPlainOptions pOptions(ctx.geoContext, ctx.magFieldContext);
  pOptions.maxSteps = 10000;

  GainMatrixUpdater kfUpdater;
  GainMatrixSmoother kfSmoother;

  ProtoTrackStateCreator trackStateCreator;
  trackStateCreator.loggerPtr = logger().clone("TrackStateCreator");
  trackStateCreator.measurements = &measurements;
  trackStateCreator.cuts = m_cfg.trackStateSelection;

  CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
  extensions.updater.connect<&GainMatrixUpdater::operator()<
      typename TrackContainer::TrackStateContainerBackend>>(&kfUpdater);
  extensions.trackStateCreator
      .template connect<&ProtoTrackStateCreator::createTrackStates>(
          &trackStateCreator);

  // Set the CombinatorialKalmanFilter options
  TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, extensions,
      pOptions, &(*pSurface));

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto trackContainer = std::make_shared<VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<VectorMultiTrajectory>();

  TrackContainer tracks(trackContainer, trackStateContainer);

  tracks.addColumn<unsigned int>("trackGroup");
  ProxyAccessor<unsigned int> seedNumber("trackGroup");

  std::size_t nSeed = 0;
  std::size_t nFailed = 0;

  std::vector<std::size_t> nTracksPerSeeds;
  nTracksPerSeeds.reserve(initialParameters.size());

  for (auto i = 0ul; i < initialParameters.size(); ++i) {
    trackStateCreator.protoTrackSourceLinks.clear();

    // Fill the source links via their indices from the container
    for (const auto hitIndex : protoTracks.at(i)) {
      if (auto it = measurements.orderedIndices().nth(hitIndex);
          it != measurements.orderedIndices().end()) {
        trackStateCreator.protoTrackSourceLinks.insert(*it);
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
    }
  }

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
  auto constTrackStateContainer = std::make_shared<ConstVectorMultiTrajectory>(
      std::move(*trackStateContainer));

  auto constTrackContainer =
      std::make_shared<ConstVectorTrackContainer>(std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTracks(ctx, std::move(constTracks));
  return ProcessCode::SUCCESS;
}

ProcessCode TrackFindingFromProtoTracksAlgorithm::finalize() {
  assert(std::distance(m_nTracksPerSeeds.begin(), m_nTracksPerSeeds.end()) > 0);

  ACTS_INFO("TrackFindingFromProtoTracksAlgorithm statistics:");
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
