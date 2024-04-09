// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/SourceLink.hpp"
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
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
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
#include <unordered_set>
#include <utility>

#include <boost/functional/hash.hpp>

namespace ActsExamples {
namespace {

/// Source link indices of the bottom, middle, top measurements.
/// In case of strip seeds only the first source link of the pair is used.
using SeedIdentifier = std::array<Index, 3>;

SeedIdentifier makeSeedIdentifier(const SimSeed& seed) {
  SeedIdentifier result;

  for (const auto& [i, sp] : Acts::enumerate(seed.sp())) {
    const Acts::SourceLink& firstSourceLink = sp->sourceLinks().front();
    result.at(i) = firstSourceLink.get<IndexSourceLink>().index();
  }

  return result;
}

template <typename Visitor>
void visitSeedIdentifiers(const TrackProxy& track, Visitor visitor) {
  std::vector<Index> sourceLinkIndices;
  sourceLinkIndices.reserve(track.nMeasurements());
  for (const auto& trackState : track.trackStatesReversed()) {
    if (!trackState.hasUncalibratedSourceLink()) {
      continue;
    }
    const Acts::SourceLink& sourceLink = trackState.getUncalibratedSourceLink();
    sourceLinkIndices.push_back(sourceLink.get<IndexSourceLink>().index());
  }

  for (std::size_t i = 0; i < sourceLinkIndices.size(); ++i) {
    for (std::size_t j = i + 1; j < sourceLinkIndices.size(); ++j) {
      for (std::size_t k = j + 1; k < sourceLinkIndices.size(); ++k) {
        // Putting them into reverse order (k, j, i) to compensate for the
        // `trackStatesReversed` above.
        visitor({sourceLinkIndices.at(k), sourceLinkIndices.at(j),
                 sourceLinkIndices.at(i)});
      }
    }
  }
}

class MeasurementSelector {
 public:
  using Traj = Acts::VectorMultiTrajectory;

  explicit MeasurementSelector(Acts::MeasurementSelector selector)
      : m_selector(std::move(selector)) {}

  void setSeed(const std::optional<SimSeed>& seed) { m_seed = seed; }

  Acts::Result<std::pair<std::vector<Traj::TrackStateProxy>::iterator,
                         std::vector<Traj::TrackStateProxy>::iterator>>
  select(std::vector<Traj::TrackStateProxy>& candidates, bool& isOutlier,
         const Acts::Logger& logger) const {
    if (m_seed.has_value()) {
      std::vector<Traj::TrackStateProxy> newCandidates;

      for (const auto& candidate : candidates) {
        if (isSeedCandidate(candidate)) {
          newCandidates.push_back(candidate);
        }
      }

      if (!newCandidates.empty()) {
        candidates = std::move(newCandidates);
      }
    }

    return m_selector.select<Acts::VectorMultiTrajectory>(candidates, isOutlier,
                                                          logger);
  }

 private:
  Acts::MeasurementSelector m_selector;
  std::optional<SimSeed> m_seed;

  bool isSeedCandidate(const Traj::TrackStateProxy& candidate) const {
    assert(candidate.hasUncalibratedSourceLink());

    const Acts::SourceLink& sourceLink = candidate.getUncalibratedSourceLink();

    for (const auto& sp : m_seed->sp()) {
      for (const auto& sl : sp->sourceLinks()) {
        if (sourceLink.get<IndexSourceLink>() == sl.get<IndexSourceLink>()) {
          return true;
        }
      }
    }

    return false;
  }
};

}  // namespace
}  // namespace ActsExamples

template <class T, std::size_t N>
struct std::hash<std::array<T, N>> {
  std::size_t operator()(const std::array<T, N>& array) const {
    std::hash<T> hasher;
    std::size_t result = 0;
    for (auto&& element : array) {
      boost::hash_combine(result, hasher(element));
    }
    return result;
  }
};

namespace ActsExamples {

TrackFindingAlgorithm::TrackFindingAlgorithm(Config config,
                                             Acts::Logging::Level level)
    : IAlgorithm("TrackFindingAlgorithm", level),
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

  if (m_cfg.seedDeduplication && m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument(
        "Missing seeds input collection. This is "
        "required for seed deduplication.");
  }
  if (m_cfg.stayOnSeed && m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument(
        "Missing seeds input collection. This is "
        "required for staying on seed.");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_inputSeeds.maybeInitialize(m_cfg.inputSeeds);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode TrackFindingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& sourceLinks = m_inputSourceLinks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);
  const SimSeedContainer* seeds = nullptr;

  if (m_inputSeeds.isInitialized()) {
    seeds = &m_inputSeeds(ctx);

    if (initialParameters.size() != seeds->size()) {
      ACTS_ERROR("Number of initial parameters and seeds do not match. "
                 << initialParameters.size() << " != " << seeds->size());
    }
  }

  std::unordered_set<SeedIdentifier> encounteredSeeds;

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);
  Acts::GainMatrixUpdater kfUpdater;
  MeasurementSelector measSel{
      Acts::MeasurementSelector(m_cfg.measurementSelectorCfg)};

  Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
      extensions;
  extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.measurementSelector.connect<&MeasurementSelector::select>(
      &measSel);

  IndexSourceLinkAccessor slAccessor;
  slAccessor.container = &sourceLinks;
  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = m_cfg.maxSteps;
  pOptions.direction = Acts::Direction::Forward;

  // Set the CombinatorialKalmanFilter options
  TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, slAccessorDelegate,
      extensions, pOptions);

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

  auto addTrack = [&](const TrackProxy& track) {
    if (m_trackSelector.has_value() && !m_trackSelector->isValidTrack(track)) {
      return;
    }

    auto destProxy = tracks.makeTrack();
    // make sure we copy track states!
    destProxy.copyFrom(track, true);

    visitSeedIdentifiers(track, [&](const SeedIdentifier& seedIdentifier) {
      encounteredSeeds.insert(seedIdentifier);
    });
  };

  for (std::size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    m_nTotalSeeds++;

    if (seeds != nullptr) {
      const SimSeed& seed = seeds->at(iseed);

      if (m_cfg.seedDeduplication) {
        SeedIdentifier seedIdentifier = makeSeedIdentifier(seed);
        if (encounteredSeeds.find(seedIdentifier) != encounteredSeeds.end()) {
          m_nDeduplicatedSeeds++;
          ACTS_VERBOSE("Skipping seed " << iseed << " due to deduplication.");
          continue;
        }
      }

      if (m_cfg.stayOnSeed) {
        measSel.setSeed(seed);
      }
    }

    // Clear trackContainerTemp and trackStateContainerTemp
    tracksTemp.clear();

    auto result =
        (*m_cfg.findTracks)(initialParameters.at(iseed), options, tracksTemp);
    nSeed++;

    if (!result.ok()) {
      m_nFailedSeeds++;
      ACTS_WARNING("Track finding failed for seed " << iseed << " with error"
                                                    << result.error());
      continue;
    }

    auto& tracksForSeed = result.value();
    for (auto& track : tracksForSeed) {
      auto smoothingResult = Acts::smoothTrack(ctx.geoContext, track, logger());
      if (!smoothingResult.ok()) {
        m_nFailedSmoothing++;
        ACTS_ERROR("Smoothing for seed "
                   << iseed << " and track " << track.index()
                   << " failed with error " << smoothingResult.error());
        continue;
      }

      auto extrapolationResult = Acts::extrapolateTrackToReferenceSurface(
          track, *pSurface, extrapolator, extrapolationOptions,
          m_cfg.extrapolationStrategy, logger());
      if (!extrapolationResult.ok()) {
        m_nFailedExtrapolation++;
        ACTS_ERROR("Extrapolation for seed "
                   << iseed << " and track " << track.index()
                   << " failed with error " << extrapolationResult.error());
        continue;
      }

      // Set the seed number, this number decrease by 1 since the seed number
      seedNumber(track) = nSeed - 1;

      addTrack(track);
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
  return ProcessCode::SUCCESS;
}

ProcessCode TrackFindingAlgorithm::finalize() {
  ACTS_INFO("TrackFindingAlgorithm statistics:");
  ACTS_INFO("- total seeds: " << m_nTotalSeeds);
  ACTS_INFO("- deduplicated seeds: " << m_nDeduplicatedSeeds);
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

}  // namespace ActsExamples
