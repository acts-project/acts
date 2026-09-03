// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFinding/TrackStateCreator.hpp"

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace ActsExamples {

namespace {

/// Creates track states for the measurements on a surface, restricted to the
/// current seed when `stayOnSeed` is set.
class TrackStateCreator final
    : public TrackStateCreatorBase<TrackStateCreator> {
 public:
  const MeasurementSubset* measurements = nullptr;

  /// Resolve the seed to measurement indices once, so the per measurement
  /// check does not have to walk space points and unpack source links.
  void setSeed(const std::optional<ConstSeedProxy>& seed) {
    m_seedIndices.clear();
    if (!seed.has_value()) {
      return;
    }
    for (const ConstSpacePointProxy sp : seed->spacePoints()) {
      for (const Acts::SourceLink& sl : sp.sourceLinks()) {
        m_seedIndices.push_back(sl.get<IndexSourceLink>().index());
      }
    }
  }

  auto measurementRange(const State& state) const {
    assert(measurements != nullptr);
    auto [begin, end] =
        measurements->orderedIndices().equal_range(state.surface->geometryId());
    return std::ranges::subrange(begin, end);
  }

  ConstVariableBoundMeasurementProxy calibrate(
      const State& /*state*/, const IndexSourceLink& measurement) const {
    // note this has to be `getMeasurement`, which takes an index into the
    // underlying container, and not `at`, which takes a position in the subset
    return measurements->getMeasurement(measurement.index());
  }

  bool hasMeasurementPreselection(const State& /*state*/) const {
    return !m_seedIndices.empty();
  }

  bool preselectMeasurement(const State& /*state*/,
                            const IndexSourceLink& measurement) const {
    return std::ranges::find(m_seedIndices, measurement.index()) !=
           m_seedIndices.end();
  }

 private:
  /// measurement indices of the current seed, at most a handful
  std::vector<Index> m_seedIndices;
};

/// Flags the seeds whose measurements are a subset of an already found track.
///
/// In case of strip seeds only the first source link of the pair is used.
class SeedCoverage {
 public:
  /// Index the seeds by their measurements.
  ///
  /// @param seeds The seeds to index.
  void index(const SeedContainer& seeds) {
    m_seedSize.assign(seeds.size(), 0);
    m_covered.assign(seeds.size(), false);
    m_counts.assign(seeds.size(), 0);

    for (std::size_t iSeed = 0; iSeed < seeds.size(); ++iSeed) {
      const ConstSeedProxy seed = seeds.at(iSeed);
      for (const SpacePointIndex spIndex : seed.spacePointIndices()) {
        const ConstSpacePointProxy sp = seeds.spacePointContainer().at(spIndex);
        if (sp.sourceLinks().empty()) {
          continue;
        }
        const Index measurement =
            sp.sourceLinks().front().get<IndexSourceLink>().index();
        m_measurementToSeeds[measurement].push_back(iSeed);
        ++m_seedSize[iSeed];
      }
    }
  }

  /// Flag every seed whose measurements the track covers.
  ///
  /// @param track The track that was found.
  void addCoverageFrom(const TrackProxy& track) {
    for (const auto& trackState : track.trackStatesReversed()) {
      if (!trackState.hasUncalibratedSourceLink()) {
        continue;
      }
      const Index measurement =
          trackState.getUncalibratedSourceLink().get<IndexSourceLink>().index();
      const auto it = m_measurementToSeeds.find(measurement);
      if (it == m_measurementToSeeds.end()) {
        continue;
      }
      for (const std::size_t iSeed : it->second) {
        if (m_counts[iSeed] == 0) {
          m_touched.push_back(iSeed);
        }
        ++m_counts[iSeed];
        if (m_counts[iSeed] >= m_seedSize[iSeed]) {
          m_covered[iSeed] = true;
        }
      }
    }

    for (const std::size_t iSeed : m_touched) {
      m_counts[iSeed] = 0;
    }
    m_touched.clear();
  }

  /// Whether a previously found track covers the seed.
  ///
  /// @param iSeed The index of the seed.
  /// @return True if the seed is covered.
  bool isCovered(std::size_t iSeed) const { return m_covered.at(iSeed); }

 private:
  std::unordered_map<Index, std::vector<std::size_t>> m_measurementToSeeds;
  std::vector<std::uint32_t> m_seedSize;
  std::vector<bool> m_covered;

  /// per track scratch, reset through `m_touched`
  std::vector<std::uint32_t> m_counts;
  std::vector<std::size_t> m_touched;
};

class BranchStopper {
 public:
  using BranchStopperResult =
      Acts::CombinatorialKalmanFilterBranchStopperResult;

  struct BranchState {
    std::size_t nPixelHoles = 0;
    std::size_t nStripHoles = 0;
  };

  static constexpr Acts::ProxyAccessor<BranchState> branchStateAccessor =
      Acts::ProxyAccessor<BranchState>(Acts::hashString("MyBranchState"));

  mutable std::atomic<std::size_t> m_nStoppedBranches{0};

  explicit BranchStopper(const TrackFindingAlgorithm::Config& config)
      : m_cfg(config) {}

  BranchStopperResult operator()(
      const TrackContainer::TrackProxy& track,
      const TrackContainer::TrackStateProxy& trackState) const {
    if (!m_cfg.trackSelectorCfg.has_value()) {
      return BranchStopperResult::Continue;
    }

    const Acts::TrackSelector::Config* singleConfig = std::visit(
        [&](const auto& config) -> const Acts::TrackSelector::Config* {
          using T = std::decay_t<decltype(config)>;
          if constexpr (std::is_same_v<T, Acts::TrackSelector::Config>) {
            return &config;
          } else if constexpr (std::is_same_v<
                                   T, Acts::TrackSelector::EtaBinnedConfig>) {
            double theta = trackState.parameters()[Acts::eBoundTheta];
            double eta = Acts::AngleHelpers::etaFromTheta(theta);
            return config.hasCuts(eta) ? &config.getCuts(eta) : nullptr;
          }
        },
        *m_cfg.trackSelectorCfg);

    if (singleConfig == nullptr) {
      ++m_nStoppedBranches;
      return BranchStopperResult::StopAndDrop;
    }

    bool tooManyHolesPS = false;
    if (!(m_cfg.pixelVolumeIds.empty() && m_cfg.stripVolumeIds.empty())) {
      auto& branchState = branchStateAccessor(track);
      // count both holes and outliers as holes for pixel/strip counts
      if (trackState.typeFlags().isHole() ||
          trackState.typeFlags().isOutlier()) {
        auto volumeId = trackState.referenceSurface().geometryId().volume();
        if (Acts::rangeContainsValue(m_cfg.pixelVolumeIds, volumeId)) {
          ++branchState.nPixelHoles;
        } else if (Acts::rangeContainsValue(m_cfg.stripVolumeIds, volumeId)) {
          ++branchState.nStripHoles;
        }
      }
      tooManyHolesPS = branchState.nPixelHoles > m_cfg.maxPixelHoles ||
                       branchState.nStripHoles > m_cfg.maxStripHoles;
    }

    bool enoughMeasurements =
        track.nMeasurements() >= singleConfig->minMeasurements;
    bool tooManyHoles =
        track.nHoles() > singleConfig->maxHoles || tooManyHolesPS;
    bool tooManyOutliers = track.nOutliers() > singleConfig->maxOutliers;
    bool tooManyHolesAndOutliers = (track.nHoles() + track.nOutliers()) >
                                   singleConfig->maxHolesAndOutliers;

    if (tooManyHoles || tooManyOutliers || tooManyHolesAndOutliers) {
      ++m_nStoppedBranches;
      return enoughMeasurements ? BranchStopperResult::StopAndKeep
                                : BranchStopperResult::StopAndDrop;
    }

    return BranchStopperResult::Continue;
  }

 private:
  const TrackFindingAlgorithm::Config& m_cfg;
};

}  // namespace

TrackFindingAlgorithm::TrackFindingAlgorithm(
    Config config, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TrackFindingAlgorithm", std::move(logger)),
      m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
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

  if (m_cfg.trackSelectorCfg.has_value()) {
    m_trackSelector = std::visit(
        [](const auto& cfg) -> std::optional<Acts::TrackSelector> {
          return Acts::TrackSelector(cfg);
        },
        m_cfg.trackSelectorCfg.value());
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_inputSeeds.maybeInitialize(m_cfg.inputSeeds);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode TrackFindingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);
  const SeedContainer* seeds = nullptr;

  if (m_inputSeeds.isInitialized()) {
    seeds = &m_inputSeeds(ctx);

    if (initialParameters.size() != seeds->size()) {
      ACTS_ERROR("Number of initial parameters and seeds do not match. "
                 << initialParameters.size() << " != " << seeds->size());
    }
  }

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::GainMatrixUpdater kfUpdater(m_cfg.useJosephFormulation);

  using Extensions = Acts::CombinatorialKalmanFilterExtensions<TrackContainer>;

  BranchStopper branchStopper(m_cfg);

  TrackStateCreator trackStateCreator;
  trackStateCreator.measurements = &measurements;
  trackStateCreator.cuts = m_cfg.trackStateSelection;

  Extensions extensions;
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<
      typename TrackContainer::TrackStateContainerBackend>>(&kfUpdater);
  extensions.branchStopper.connect<&BranchStopper::operator()>(&branchStopper);
  extensions.trackStateCreator
      .template connect<&TrackStateCreator::createTrackStates>(
          &trackStateCreator);
  extensions.mixtureReducer.connect<&Acts::reduceMixtureWithKLDistance>();

  Acts::PropagatorPlainOptions firstPropOptions(ctx.recoGeoContext,
                                                ctx.magFieldContext);
  firstPropOptions.maxSteps = m_cfg.maxSteps;
  firstPropOptions.direction = m_cfg.reverseSearch ? Acts::Direction::Backward()
                                                   : Acts::Direction::Forward();
  firstPropOptions.constrainToVolumeIds = m_cfg.constrainToVolumeIds;
  firstPropOptions.endOfWorldVolumeIds = m_cfg.endOfWorldVolumeIds;

  Acts::PropagatorPlainOptions secondPropOptions(ctx.recoGeoContext,
                                                 ctx.magFieldContext);
  secondPropOptions.maxSteps = m_cfg.maxSteps;
  secondPropOptions.direction = firstPropOptions.direction.invert();
  secondPropOptions.constrainToVolumeIds = m_cfg.constrainToVolumeIds;
  secondPropOptions.endOfWorldVolumeIds = m_cfg.endOfWorldVolumeIds;

  // Set the CombinatorialKalmanFilter options
  TrackFinderOptions firstOptions(ctx.recoGeoContext, ctx.magFieldContext,
                                  ctx.calibContext, extensions,
                                  firstPropOptions);

  firstOptions.targetSurface = m_cfg.reverseSearch ? pSurface.get() : nullptr;
  firstOptions.recordMaterialStates = m_cfg.recordMaterialStates;
  firstOptions.betheHeitlerApprox =
      std::make_shared<Acts::AtlasBetheHeitlerApprox>(
          Acts::makeDefaultBetheHeitlerApprox());

  TrackFinderOptions secondOptions(ctx.recoGeoContext, ctx.magFieldContext,
                                   ctx.calibContext, extensions,
                                   secondPropOptions);
  secondOptions.targetSurface = m_cfg.reverseSearch ? nullptr : pSurface.get();
  secondOptions.skipPrePropagationUpdate = true;
  secondOptions.recordMaterialStates = m_cfg.recordMaterialStates;
  secondOptions.betheHeitlerApprox = firstOptions.betheHeitlerApprox;

  // the brem finder cannot skip material states, so it gets its own options
  TrackFinderOptions firstBremOptions = firstOptions;
  firstBremOptions.recordMaterialStates = true;
  TrackFinderOptions secondBremOptions = secondOptions;
  secondBremOptions.recordMaterialStates = true;

  using Extrapolator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
  using ExtrapolatorOptions = Extrapolator::template Options<
      Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;

  Extrapolator extrapolator(
      Acts::SympyStepper(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  ExtrapolatorOptions extrapolationOptions(ctx.recoGeoContext,
                                           ctx.magFieldContext);
  extrapolationOptions.constrainToVolumeIds = m_cfg.constrainToVolumeIds;
  extrapolationOptions.endOfWorldVolumeIds = m_cfg.endOfWorldVolumeIds;

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

  // Note that not all backends support PODs as column types
  tracks.addColumn<BranchStopper::BranchState>("MyBranchState");
  tracksTemp.addColumn<BranchStopper::BranchState>("MyBranchState");

  tracks.addColumn<unsigned int>("trackGroup");
  tracksTemp.addColumn<unsigned int>("trackGroup");
  Acts::ProxyAccessor<unsigned int> seedNumber("trackGroup");

  unsigned int nSeed = 0;

  // A map indicating whether a seed has been discovered already
  SeedCoverage seedCoverage;

  auto addTrack = [&](const TrackProxy& track) {
    ++m_nFoundTracks;

    // trim the track if requested
    if (m_cfg.trimTracks) {
      Acts::trimTrack(track, true, true, true, true);
    }
    Acts::calculateTrackQuantities(track);

    if (m_trackSelector.has_value() && !m_trackSelector->isValidTrack(track)) {
      return;
    }

    // flag seeds which are covered by the track
    seedCoverage.addCoverageFrom(track);

    ++m_nSelectedTracks;

    auto destProxy = tracks.makeTrack();
    // make sure we copy track states!
    destProxy.copyFrom(track);
  };

  if (seeds != nullptr && m_cfg.seedDeduplication) {
    // Index the seeds for deduplication
    seedCoverage.index(*seeds);
  }

  for (std::size_t iSeed = 0; iSeed < initialParameters.size(); ++iSeed) {
    m_nTotalSeeds++;

    if (seeds != nullptr) {
      const ConstSeedProxy seed = seeds->at(iSeed);

      if (m_cfg.seedDeduplication) {
        // check if an already found track covers the seed
        if (seedCoverage.isCovered(iSeed)) {
          m_nDeduplicatedSeeds++;
          ACTS_VERBOSE("Skipping seed " << iSeed << " due to deduplication.");
          continue;
        }
      }

      if (m_cfg.stayOnSeed) {
        trackStateCreator.setSeed(seed);
      }
    }

    // Clear trackContainerTemp and trackStateContainerTemp
    tracksTemp.clear();

    const Acts::BoundTrackParameters& firstInitialParameters =
        initialParameters.at(iSeed);
    ACTS_VERBOSE("Processing seed " << iSeed << " with initial parameters "
                                    << firstInitialParameters);

    const bool useBrem = m_cfg.findTracksBrem != nullptr &&
                         firstInitialParameters.particleHypothesis() ==
                             Acts::ParticleHypothesis::electron();
    const TrackFinderFunction& findTracks =
        useBrem ? *m_cfg.findTracksBrem : *m_cfg.findTracks;
    const TrackFinderOptions& firstFindOptions =
        useBrem ? firstBremOptions : firstOptions;
    const TrackFinderOptions& secondFindOptions =
        useBrem ? secondBremOptions : secondOptions;

    auto firstRootBranch = tracksTemp.makeTrack();
    auto firstResult = findTracks(firstInitialParameters, firstFindOptions,
                                  tracksTemp, firstRootBranch);
    nSeed++;

    if (!firstResult.ok()) {
      m_nFailedSeeds++;
      ACTS_WARNING("Track finding failed for seed " << iSeed << " with error"
                                                    << firstResult.error());
      continue;
    }

    auto& firstTracksForSeed = firstResult.value();
    for (auto& firstTrack : firstTracksForSeed) {
      // TODO a copy of the track should not be necessary but is the safest way
      //      with the current EDM
      // TODO a lightweight copy without copying all the track state components
      //      might be a solution
      auto trackCandidate = tracksTemp.makeTrack();
      trackCandidate.copyFrom(firstTrack);

      Acts::Result<void> firstSmoothingResult{
          Acts::smoothTrack(ctx.recoGeoContext, trackCandidate, logger())};
      if (!firstSmoothingResult.ok()) {
        m_nFailedSmoothing++;
        ACTS_ERROR("First smoothing for seed "
                   << iSeed << " and track " << firstTrack.index()
                   << " failed with error " << firstSmoothingResult.error());
        continue;
      }

      // number of second tracks found
      std::size_t nSecond = 0;

      // Set the seed number, this number decrease by 1 since the seed number
      // has already been updated
      seedNumber(trackCandidate) = nSeed - 1;

      if (m_cfg.twoWay) {
        std::optional<Acts::VectorMultiTrajectory::TrackStateProxy>
            firstMeasurementOpt;
        for (auto trackState : trackCandidate.trackStatesReversed()) {
          // We are excluding non measurement states and outlier here. Those can
          // decrease resolution because only the smoothing corrected the very
          // first prediction as filtering is not possible.
          if (trackState.typeFlags().isMeasurement()) {
            firstMeasurementOpt = trackState;
          }
        }

        if (firstMeasurementOpt.has_value()) {
          TrackContainer::TrackStateProxy firstMeasurement{
              firstMeasurementOpt.value()};
          TrackContainer::ConstTrackStateProxy firstMeasurementConst{
              firstMeasurement};

          Acts::BoundTrackParameters secondInitialParameters =
              trackCandidate.createParametersFromState(firstMeasurementConst);

          if (!secondInitialParameters.referenceSurface().insideBounds(
                  secondInitialParameters.localPosition())) {
            m_nSkippedSecondPass++;
            ACTS_DEBUG(
                "Smoothing of first pass fit produced out-of-bounds parameters "
                "relative to the surface. Skipping second pass.");
            continue;
          }

          auto secondRootBranch = tracksTemp.makeTrack();
          secondRootBranch.copyFromWithoutStates(trackCandidate);
          auto secondResult =
              findTracks(secondInitialParameters, secondFindOptions, tracksTemp,
                         secondRootBranch);

          if (!secondResult.ok()) {
            ACTS_WARNING("Second track finding failed for seed "
                         << iSeed << " with error" << secondResult.error());
          } else {
            // store the original previous state to restore it later
            auto originalFirstMeasurementPrevious = firstMeasurement.previous();

            auto& secondTracksForSeed = secondResult.value();
            for (auto& secondTrack : secondTracksForSeed) {
              // TODO a copy of the track should not be necessary but is the
              //      safest way with the current EDM
              // TODO a lightweight copy without copying all the track state
              //      components might be a solution
              auto secondTrackCopy = tracksTemp.makeTrack();
              secondTrackCopy.copyFrom(secondTrack);

              // Note that this is only valid if there are no branches
              // We disallow this by breaking this look after a second track was
              // processed
              secondTrackCopy.reverseTrackStates(true);

              firstMeasurement.previous() =
                  secondTrackCopy.outermostTrackState().index();

              // Retain tip and stem index of the first track
              auto tipIndex = trackCandidate.tipIndex();
              auto stemIndex = trackCandidate.stemIndex();
              trackCandidate.copyFromWithoutStates(secondTrackCopy);
              trackCandidate.tipIndex() = tipIndex;
              trackCandidate.stemIndex() = stemIndex;

              // finalize the track candidate

              bool doExtrapolate = true;

              if (!m_cfg.reverseSearch) {
                // these parameters are already extrapolated by the CKF and have
                // the optimal resolution. note that we did not smooth all the
                // states.

                // only extrapolate if we did not do it already
                doExtrapolate = !trackCandidate.hasReferenceSurface();
              } else {
                // smooth the full track and extrapolate to the reference

                auto secondSmoothingResult = Acts::smoothTrack(
                    ctx.recoGeoContext, trackCandidate, logger());
                if (!secondSmoothingResult.ok()) {
                  m_nFailedSmoothing++;
                  ACTS_ERROR("Second smoothing for seed "
                             << iSeed << " and track " << secondTrack.index()
                             << " failed with error "
                             << secondSmoothingResult.error());
                  continue;
                }

                trackCandidate.reverseTrackStates(true);
              }

              if (doExtrapolate) {
                auto secondExtrapolationResult =
                    Acts::extrapolateTrackToReferenceSurface(
                        trackCandidate, *pSurface, extrapolator,
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
              }

              addTrack(trackCandidate);

              ++nSecond;
            }

            // restore the original previous state
            firstMeasurement.previous() = originalFirstMeasurementPrevious;
          }
        }
      }

      // if no second track was found, we will use only the first track
      if (nSecond == 0) {
        // restore the track to the original state
        auto tipIndex = trackCandidate.tipIndex();
        auto stemIndex = trackCandidate.stemIndex();
        trackCandidate.copyFromWithoutStates(firstTrack);
        trackCandidate.tipIndex() = tipIndex;
        trackCandidate.stemIndex() = stemIndex;

        auto firstExtrapolationResult =
            Acts::extrapolateTrackToReferenceSurface(
                trackCandidate, *pSurface, extrapolator, extrapolationOptions,
                m_cfg.extrapolationStrategy, logger());
        if (!firstExtrapolationResult.ok()) {
          m_nFailedExtrapolation++;
          ACTS_ERROR("Extrapolation for seed "
                     << iSeed << " and track " << firstTrack.index()
                     << " failed with error "
                     << firstExtrapolationResult.error());
          continue;
        }

        addTrack(trackCandidate);
      }
    }
  }

  // Compute shared hits from all the reconstructed tracks
  if (m_cfg.computeSharedHits) {
    computeSharedHits(tracks, measurements);
  }

  ACTS_DEBUG("Finalized track finding with " << tracks.size()
                                             << " track candidates.");

  m_nStoppedBranches += branchStopper.m_nStoppedBranches;

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
  ACTS_INFO("- found tracks: " << m_nFoundTracks);
  ACTS_INFO("- selected tracks: " << m_nSelectedTracks);
  ACTS_INFO("- stopped branches: " << m_nStoppedBranches);
  ACTS_INFO("- skipped second pass: " << m_nSkippedSecondPass);

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

// TODO this is somewhat duplicated in AmbiguityResolutionAlgorithm.cpp
// TODO we should make a common implementation in the core at some point
void TrackFindingAlgorithm::computeSharedHits(
    TrackContainer& tracks, const MeasurementSubset& measurements) const {
  // Compute shared hits from all the reconstructed tracks
  // Compute nSharedhits and Update ckf results
  // hit index -> list of multi traj indexes [traj, meas]

  std::vector<std::size_t> firstTrackOnTheHit(
      measurements.container().size(), std::numeric_limits<std::size_t>::max());
  std::vector<std::size_t> firstStateOnTheHit(
      measurements.container().size(), std::numeric_limits<std::size_t>::max());

  for (auto track : tracks) {
    for (auto state : track.trackStatesReversed()) {
      if (!state.typeFlags().isMeasurement()) {
        continue;
      }

      std::size_t hitIndex = state.getUncalibratedSourceLink()
                                 .template get<IndexSourceLink>()
                                 .index();

      // Check if hit not already used
      if (firstTrackOnTheHit.at(hitIndex) ==
          std::numeric_limits<std::size_t>::max()) {
        firstTrackOnTheHit.at(hitIndex) = track.index();
        firstStateOnTheHit.at(hitIndex) = state.index();
        continue;
      }

      // if already used, control if first track state has been marked
      // as shared
      std::size_t indexFirstTrack = firstTrackOnTheHit.at(hitIndex);
      std::size_t indexFirstState = firstStateOnTheHit.at(hitIndex);

      auto firstState = tracks.getTrack(indexFirstTrack)
                            .container()
                            .trackStateContainer()
                            .getTrackState(indexFirstState);
      firstState.typeFlags().setIsSharedHit();

      // Decorate this track state
      state.typeFlags().setIsSharedHit();
    }
  }
}

}  // namespace ActsExamples
