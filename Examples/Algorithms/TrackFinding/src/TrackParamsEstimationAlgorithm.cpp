// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TrackParamsEstimationAlgorithm", std::move(logger)),
      m_cfg(cfg),
      m_estimator(m_cfg, this->logger().clone()) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }

  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_inputTracks.maybeInitialize(m_cfg.inputProtoTracks);
  m_inputParticleHypotheses.maybeInitialize(m_cfg.inputParticleHypotheses);
  m_inputFreeParameters.maybeInitialize(m_cfg.inputFreeParameters);

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputSeeds.maybeInitialize(m_cfg.outputSeeds);
  m_outputTracks.maybeInitialize(m_cfg.outputProtoTracks);
}

ProcessCode TrackParamsEstimationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  auto const& seeds = m_inputSeeds(ctx);
  ACTS_VERBOSE("Read " << seeds.size() << " seeds");

  TrackParametersContainer trackParameters;
  trackParameters.reserve(seeds.size());

  SeedContainer outputSeeds;
  if (m_outputSeeds.isInitialized()) {
    outputSeeds.assignSpacePointContainer(seeds.spacePointContainer());
    outputSeeds.reserve(seeds.size());
  }

  const ProtoTrackContainer* inputTracks = nullptr;
  ProtoTrackContainer outputTracks;
  if (m_inputTracks.isInitialized() && m_outputTracks.isInitialized()) {
    const auto& inputTracksRef = m_inputTracks(ctx);
    if (seeds.size() != inputTracksRef.size()) {
      ACTS_FATAL("Inconsistent number of seeds and proto tracks");
      return ProcessCode::ABORT;
    }
    inputTracks = &inputTracksRef;
    outputTracks.reserve(seeds.size());
  }

  const std::vector<Acts::ParticleHypothesis>* inputParticleHypotheses =
      nullptr;
  if (m_inputParticleHypotheses.isInitialized()) {
    const auto& inputParticleHypothesesRef = m_inputParticleHypotheses(ctx);
    if (seeds.size() != inputParticleHypothesesRef.size()) {
      ACTS_FATAL("Inconsistent number of seeds and particle hypotheses");
      return ProcessCode::ABORT;
    }
    inputParticleHypotheses = &inputParticleHypothesesRef;
  }

  const std::vector<Acts::FreeVector>* inputFreeParameters = nullptr;
  if (m_inputFreeParameters.isInitialized()) {
    const auto& inputFreeParametersRef = m_inputFreeParameters(ctx);
    if (seeds.size() != inputFreeParametersRef.size()) {
      ACTS_FATAL("Inconsistent number of seeds and free parameters");
      return ProcessCode::ABORT;
    }
    inputFreeParameters = &inputFreeParametersRef;
  }

  SeedParamsEstimator::Cache cache = m_estimator.makeCache(ctx.magFieldContext);

  struct {
    std::size_t selection = 0;
    std::size_t fit = 0;
    std::size_t degenerate = 0;
    std::size_t transport = 0;

    std::size_t total() const {
      return selection + fit + degenerate + transport;
    }
  } skipped;

  // Loop over all found seeds to estimate track parameters
  for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto& seed = seeds[iseed];
    if (seed.spacePoints().size() < 3) {
      ACTS_WARNING("Seed " << iseed << " has less than 3 space points, skip");
      continue;
    }

    const Acts::ParticleHypothesis hypothesis =
        inputParticleHypotheses != nullptr ? inputParticleHypotheses->at(iseed)
                                           : m_cfg.particleHypothesis;

    const Acts::Result<TrackParameters, SeedParamsEstimator::Failure> estimate =
        m_estimator.estimate(
            ctx.recoGeoContext, ctx.magFieldContext, cache, seed, hypothesis,
            inputFreeParameters != nullptr ? &inputFreeParameters->at(iseed)
                                           : nullptr,
            m_cfg.parametersOnLastSpacePoint
                ? SeedParamsEstimator::Anchor::Last
                : SeedParamsEstimator::Anchor::First);
    if (!estimate.ok()) {
      using enum SeedParamsEstimator::Failure;
      switch (estimate.error()) {
        case FieldLookup:
          return ProcessCode::ABORT;
        case SpacePointSelection:
          ++skipped.selection;
          break;
        case Fit:
          ++skipped.fit;
          break;
        case Degenerate:
          ++skipped.degenerate;
          break;
        case Transport:
          ++skipped.transport;
          break;
        case MissingSurface:
        case WeakField:
          break;
      }
      ACTS_DEBUG("Seed " << iseed << " has no estimate, skip");
      continue;
    }

    const TrackParameters& trackParams =
        trackParameters.emplace_back(*estimate);
    ACTS_VERBOSE("Estimated track parameters: " << trackParams);
    if (m_outputSeeds.isInitialized()) {
      auto newSp = outputSeeds.createSeed();
      // TODO copy shorthand
      newSp.assignSpacePointIndices(seed.spacePointIndices());
      newSp.quality() = seed.quality();
      newSp.vertexZ() = seed.vertexZ();
    }
    if (m_outputTracks.isInitialized() && inputTracks != nullptr) {
      outputTracks.push_back(inputTracks->at(iseed));
    }
  }

  ACTS_DEBUG("Estimated " << trackParameters.size() << " track parameters from "
                          << seeds.size() << " seeds");
  if (skipped.total() > 0) {
    ACTS_DEBUG(
        "Skipped " << skipped.selection
                   << " seeds without a space point selection, " << skipped.fit
                   << " without a fit, " << skipped.degenerate
                   << " with a degenerate estimate and " << skipped.transport
                   << " without a transport to the surface");
  }

  m_outputTrackParameters(ctx, std::move(trackParameters));
  if (m_outputSeeds.isInitialized()) {
    m_outputSeeds(ctx, std::move(outputSeeds));
  }

  if (m_outputTracks.isInitialized()) {
    m_outputTracks(ctx, std::move(outputTracks));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
