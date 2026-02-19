// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <cstddef>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("TrackParamsEstimationAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (!m_cfg.magneticField) {
    throw std::invalid_argument("Missing magnetic field");
  }

  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_inputTracks.maybeInitialize(m_cfg.inputProtoTracks);
  m_inputParticleHypotheses.maybeInitialize(m_cfg.inputParticleHypotheses);

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

  auto bCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);

  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};

  // Loop over all found seeds to estimate track parameters
  for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto& seed = seeds[iseed];
    if (seed.spacePoints().size() < 3) {
      ACTS_WARNING("Seed " << iseed << " has less than 3 space points");
      continue;
    } else if (seed.spacePoints().size() > 3) {
      ACTS_WARNING(
          "Seed "
          << iseed
          << " has more than 3 space points, only the first 3 will be used");
    }

    // Get the bottom space point and its reference surface
    const ConstSpacePointProxy bottomSp = seed.spacePoints()[0];
    const ConstSpacePointProxy middleSp = seed.spacePoints()[1];
    const ConstSpacePointProxy topSp = seed.spacePoints()[2];
    if (bottomSp.sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in the space point");
      continue;
    }

    const Acts::Vector3 bottomSpVec{bottomSp.x(), bottomSp.y(), bottomSp.z()};
    const Acts::Vector3 middleSpVec{middleSp.x(), middleSp.y(), middleSp.z()};
    const Acts::Vector3 topSpVec{topSp.x(), topSp.y(), topSp.z()};

    const Acts::SourceLink& bottomSourceLink = bottomSp.sourceLinks()[0];
    const Acts::Surface* bottomSurface = surfaceAccessor(bottomSourceLink);
    if (bottomSurface == nullptr) {
      ACTS_WARNING(
          "Surface from source link is not found in the tracking geometry");
      continue;
    }

    // Get the magnetic field at the bottom space point
    auto fieldRes = m_cfg.magneticField->getField(bottomSpVec, bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return ProcessCode::ABORT;
    }
    const Acts::Vector3 field = *fieldRes;

    if (field.norm() < m_cfg.bFieldMin) {
      ACTS_WARNING("Magnetic field at seed " << iseed << " is too small "
                                             << field.norm());
      continue;
    }

    // Estimate the track parameters from seed
    Acts::Result<Acts::BoundVector> boundParams =
        Acts::estimateTrackParamsFromSeed(
            ctx.geoContext, *bottomSurface, bottomSpVec,
            std::isnan(bottomSp.time()) ? 0.0 : bottomSp.time(), middleSpVec,
            topSpVec, field);
    if (!boundParams.ok()) {
      ACTS_WARNING("Failed to estimate track parameters from seed: "
                   << boundParams.error().message());
      continue;
    }

    Acts::EstimateTrackParamCovarianceConfig config{
        .initialSigmas =
            Eigen::Map<const Acts::BoundVector>{m_cfg.initialSigmas.data()},
        .initialSigmaQoverPt = m_cfg.initialSigmaQoverPt,
        .initialSigmaPtRel = m_cfg.initialSigmaPtRel,
        .initialVarInflation = Eigen::Map<const Acts::BoundVector>{
            m_cfg.initialVarInflation.data()}};

    const Acts::BoundMatrix cov = Acts::estimateTrackParamCovariance(
        config, *boundParams, !std::isnan(bottomSp.time()));

    const Acts::ParticleHypothesis hypothesis =
        inputParticleHypotheses != nullptr ? inputParticleHypotheses->at(iseed)
                                           : m_cfg.particleHypothesis;

    const TrackParameters& trackParams = trackParameters.emplace_back(
        bottomSurface->getSharedPtr(), *boundParams, cov, hypothesis);
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

  ACTS_DEBUG("Estimated " << trackParameters.size() << " track parameters");

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
