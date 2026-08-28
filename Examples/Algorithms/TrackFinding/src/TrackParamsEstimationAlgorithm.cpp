// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <ostream>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

namespace {

/// Moves the estimate onto the surface of its first space point. That is at
/// most the distance of a space point from its own module, so there is nothing
/// in between to navigate.
using SeedPropagator =
    Acts::Propagator<Acts::SympyStepper, Acts::VoidNavigator>;
using SeedPropagatorOptions = SeedPropagator::Options<>;

/// Fit a helix through the given positions and express it at the first of them.
///
/// @param positions the global positions in track order, innermost first
/// @param weights relative weights of the positions, empty for uniform. Only
///        the fit of more than three positions uses them.
/// @param bField the magnetic field vector at the first position
/// @param t0 the time at the first position
/// @param refineIterations geometric refinement iterations of the circle fit
///
/// @return the free parameters at the first position
Acts::Result<Acts::FreeVector> estimateFreeParams(
    std::span<const Acts::Vector3> positions, std::span<const double> weights,
    const Acts::Vector3& bField, double t0, std::size_t refineIterations) {
  // three positions determine a helix exactly, more are fitted
  if (positions.size() == 3) {
    return Acts::Result<Acts::FreeVector>::success(
        Acts::estimateTrackParamsFromSeed(positions[0], t0, positions[1],
                                          positions[2], bField));
  }
  return Acts::estimateTrackParamsFromSpacePoints(positions, bField, t0,
                                                  refineIterations, weights);
}

/// Express free parameters on a surface.
///
/// The parameters are propagated to the surface, which is a no-op if they are
/// already on it. A space point does not have to be on the surface of its
/// source link: a strip space point can sit between the two modules it was
/// built from, and either of them can be behind it.
///
/// @param propagator the propagator to transport with
/// @param options the propagation options, direction is overridden
/// @param freeParams the free parameters to express
/// @param surface the surface to express them on
/// @param hypothesis the particle hypothesis of the parameters
///
/// @return the bound parameters on the surface
Acts::Result<Acts::BoundVector> transportToSurface(
    const SeedPropagator& propagator, const SeedPropagatorOptions& options,
    const Acts::FreeVector& freeParams, const Acts::Surface& surface,
    const Acts::ParticleHypothesis& hypothesis) {
  const Acts::Vector3 direction = freeParams.segment<3>(Acts::eFreeDir0);

  // only for the direction to propagate in, so the bounds do not matter here
  const Acts::Intersection3D intersection =
      surface
          .intersect(options.geoContext, freeParams.segment<3>(Acts::eFreePos0),
                     direction)
          .closest();
  if (!intersection.isValid()) {
    return Acts::Result<Acts::BoundVector>::failure(
        Acts::PropagatorError::Failure);
  }

  SeedPropagatorOptions surfaceOptions = options;
  surfaceOptions.direction =
      Acts::Direction::fromScalarZeroAsPositive(intersection.pathLength());

  const Acts::BoundTrackParameters start =
      Acts::BoundTrackParameters::createCurvilinear(
          freeParams.segment<4>(Acts::eFreePos0), direction,
          freeParams[Acts::eFreeQOverP], std::nullopt, hypothesis);

  auto result = propagator.propagate(start, surface, surfaceOptions);
  if (!result.ok()) {
    return Acts::Result<Acts::BoundVector>::failure(result.error());
  }
  return Acts::Result<Acts::BoundVector>::success(
      result->endParameters.value().parameters());
}

}  // namespace

TrackParamsEstimationAlgorithm::SpacePointWeight
TrackParamsEstimationAlgorithm::inverseRadiusPowerWeight(double exponent) {
  return [exponent](const Acts::Vector3& position) {
    const double r = Acts::VectorHelpers::perp(position);
    if (!(r > 0)) {
      return 1.;
    }
    return std::pow(r, -exponent);
  };
}

TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TrackParamsEstimationAlgorithm", std::move(logger)),
      m_cfg(cfg) {
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

  const SpacePointContainer& spacePoints = seeds.spacePointContainer();

  const SeedPropagator propagator(Acts::SympyStepper(m_cfg.magneticField),
                                  Acts::VoidNavigator(),
                                  logger().cloneWithSuffix("Propagator"));
  const SeedPropagatorOptions propagatorOptions(ctx.recoGeoContext,
                                                ctx.magFieldContext);

  // reused across seeds by the helix fit
  std::vector<Acts::Vector3> positions;
  std::vector<double> weights;

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

    const std::optional<std::vector<SpacePointIndex>> selected =
        selectSeedSpacePoints(spacePoints, seed.spacePointIndices(),
                              m_cfg.spacePointSelection,
                              m_cfg.minTransverseDistance);
    if (!selected.has_value()) {
      ACTS_DEBUG("Seed " << iseed << " has no space point selection, skip");
      ++skipped.selection;
      continue;
    }

    // Get the bottom space point and its reference surface
    const ConstSpacePointProxy bottomSp = spacePoints.at((*selected)[0]);
    if (bottomSp.sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in the space point");
      continue;
    }

    const Acts::Vector3 bottomSpVec{bottomSp.x(), bottomSp.y(), bottomSp.z()};

    const Acts::SourceLink& bottomSourceLink = bottomSp.sourceLinks()[0];
    const Acts::Surface* bottomSurface = surfaceAccessor(bottomSourceLink);
    if (bottomSurface == nullptr) {
      ACTS_WARNING(
          "Surface from source link is not found in the tracking geometry");
      continue;
    }

    // Get the magnetic field at the bottom space point
    const auto fieldRes = m_cfg.magneticField->getField(bottomSpVec, bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return ProcessCode::ABORT;
    }
    const Acts::Vector3& field = *fieldRes;

    if (field.norm() < m_cfg.bFieldMin) {
      ACTS_WARNING("Magnetic field at seed " << iseed << " is too small "
                                             << field.norm());
      continue;
    }

    positions.clear();
    weights.clear();
    for (const SpacePointIndex index : *selected) {
      const ConstSpacePointProxy sp = spacePoints.at(index);
      positions.emplace_back(sp.x(), sp.y(), sp.z());
      if (m_cfg.spacePointWeight) {
        weights.push_back(m_cfg.spacePointWeight(positions.back()));
      }
    }

    const double t0 = std::isnan(bottomSp.time()) ? 0.0 : bottomSp.time();

    const Acts::Result<Acts::FreeVector> freeParams = estimateFreeParams(
        positions, weights, field, t0, m_cfg.geometricRefineIterations);
    if (!freeParams.ok()) {
      ACTS_DEBUG("Seed " << iseed << " could not be fitted: "
                         << freeParams.error().message());
      ++skipped.fit;
      continue;
    }
    // Degenerate space points, e.g. a bottom and a middle space point at the
    // same transverse position, make the estimate not a number rather than
    // merely wrong. A straight line fit instead leaves q/p at zero, an
    // infinite momentum, which makes the time part of the covariance transport
    // not a number further down the line.
    if (!freeParams->allFinite() || (*freeParams)[Acts::eFreeQOverP] == 0) {
      ACTS_DEBUG("Seed " << iseed << " has a degenerate estimate, skip");
      ++skipped.degenerate;
      continue;
    }

    const Acts::ParticleHypothesis hypothesis =
        inputParticleHypotheses != nullptr ? inputParticleHypotheses->at(iseed)
                                           : m_cfg.particleHypothesis;

    const Acts::Result<Acts::BoundVector> boundParams = transportToSurface(
        propagator, propagatorOptions, *freeParams, *bottomSurface, hypothesis);
    if (!boundParams.ok()) {
      ACTS_DEBUG("Seed " << iseed
                         << " could not be transported to the surface of its "
                            "first space point: "
                         << boundParams.error().message());
      ++skipped.transport;
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
