// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedParamsEstimator.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <cmath>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

namespace {

// the step goes from a space point to its own module, or from given free
// parameters across the seed, so the material on the way is ignored
using SeedPropagator =
    Acts::Propagator<Acts::SympyStepper, Acts::VoidNavigator>;
using SeedPropagatorOptions = SeedPropagator::Options<>;

Acts::Result<Acts::FreeVector> estimateFreeParams(
    std::span<const Acts::Vector3> positions, const Acts::Vector3& bField,
    double t0, std::span<const double> weights, std::size_t refineIterations) {
  // three positions determine a helix exactly, more are fitted
  if (positions.size() == 3) {
    return Acts::Result<Acts::FreeVector>::success(
        Acts::estimateTrackParamsFromSeed(positions[0], t0, positions[1],
                                          positions[2], bField));
  }

  return Acts::estimateTrackParamsFromSpacePoints(positions, bField, t0,
                                                  refineIterations, weights);
}

// a space point does not have to lie on the surface of its own measurement
Acts::Result<Acts::BoundVector> transportToSurface(
    const SeedPropagator& propagator, const SeedPropagatorOptions& options,
    const Acts::FreeVector& freeParams, const Acts::Surface& surface,
    const Acts::ParticleHypothesis& hypothesis) {
  const Acts::Vector3 direction = freeParams.segment<3>(Acts::eFreeDir0);

  // only for the direction to propagate in, so the bounds do not matter
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

SeedParamsEstimator::SpacePointWeight
SeedParamsEstimator::inverseRadiusPowerWeight(double exponent) {
  return [exponent](const Acts::Vector3& position) {
    const double r = Acts::VectorHelpers::perp(position);
    if (r <= 0) {
      return 1.;
    }
    return std::pow(r, -exponent);
  };
}

SeedParamsEstimator::SeedParamsEstimator(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(std::move(cfg)),
      m_logger(logger != nullptr
                   ? std::move(logger)
                   : Acts::getDefaultLogger("SeedParamsEstimator",
                                            Acts::Logging::INFO)) {
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (!m_cfg.magneticField) {
    throw std::invalid_argument("Missing magnetic field");
  }
  m_propagatorLogger = m_logger->cloneWithSuffix("Propagator");
}

SeedParamsEstimator::Cache SeedParamsEstimator::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  return Cache{.field = m_cfg.magneticField->makeCache(mctx)};
}

Acts::Result<TrackParameters, SeedParamsEstimator::Failure>
SeedParamsEstimator::estimate(const Acts::GeometryContext& gctx,
                              const Acts::MagneticFieldContext& mctx,
                              Cache& cache, const ConstSeedProxy& seed,
                              const Acts::ParticleHypothesis& hypothesis,
                              const Acts::FreeVector* freeParams,
                              Anchor anchor) const {
  using Result = Acts::Result<TrackParameters, Failure>;

  const SpacePointContainer& spacePoints =
      seed.container().spacePointContainer();

  // the fit runs from the anchor, so for the last space point it runs over
  // the seed in reverse and the selection picks from the outer end
  std::span<const SpacePointIndex> candidates = seed.spacePointIndices();
  if (anchor == Anchor::Last) {
    cache.candidates.assign(candidates.rbegin(), candidates.rend());
    candidates = cache.candidates;
  }
  std::span<const SpacePointIndex> selected = candidates;

  // given free parameters replace the fit, unless they have no momentum
  const bool hasFreeParams =
      freeParams != nullptr && (*freeParams)[Acts::eFreeQOverP] != 0;

  // the triplet selections pick three space points
  if (m_cfg.spacePointSelection != SeedSpacePointSelection::All) {
    const std::optional<std::array<SpacePointIndex, 3>> selectedTriplet =
        selectSeedSpacePoints(spacePoints, candidates,
                              m_cfg.spacePointSelection,
                              m_cfg.minTransverseDistance);
    // a seed with given free parameters keeps all of its space points
    if (!selectedTriplet.has_value() && !hasFreeParams) {
      return Result::failure(Failure::SpacePointSelection);
    }
    if (selectedTriplet.has_value()) {
      cache.triplet = *selectedTriplet;
      selected = cache.triplet;
    }
  }

  // the fit starts at the anchor space point
  const ConstSpacePointProxy anchorSp = spacePoints.at(selected.front());
  if (anchorSp.sourceLinks().empty()) {
    ACTS_WARNING("Missing source link in the space point");
    return Result::failure(Failure::MissingSurface);
  }
  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};
  const Acts::Surface* anchorSurface =
      surfaceAccessor(anchorSp.sourceLinks()[0]);
  if (anchorSurface == nullptr) {
    ACTS_WARNING(
        "Surface from source link is not found in the tracking geometry");
    return Result::failure(Failure::MissingSurface);
  }

  const Acts::Vector3 anchorSpVec{anchorSp.x(), anchorSp.y(), anchorSp.z()};

  // Get the magnetic field at the anchor space point
  const auto fieldRes = m_cfg.magneticField->getField(anchorSpVec, cache.field);
  if (!fieldRes.ok()) {
    ACTS_ERROR("Field lookup error: " << fieldRes.error());
    return Result::failure(Failure::FieldLookup);
  }
  const Acts::Vector3& field = *fieldRes;

  if (field.norm() < m_cfg.bFieldMin) {
    ACTS_WARNING("Magnetic field at seed is too small " << field.norm());
    return Result::failure(Failure::WeakField);
  }

  cache.positions.clear();
  for (const SpacePointIndex index : selected) {
    const ConstSpacePointProxy sp = spacePoints.at(index);
    cache.positions.emplace_back(sp.x(), sp.y(), sp.z());
  }

  cache.weights.clear();
  if (m_cfg.spacePointWeight) {
    for (const Acts::Vector3& position : cache.positions) {
      cache.weights.emplace_back(m_cfg.spacePointWeight(position));
    }
  }

  auto timeOf = [](const ConstSpacePointProxy& sp) {
    return std::isnan(sp.time()) ? 0.0 : sp.time();
  };
  const double t0 = timeOf(anchorSp);

  // a fit against the flight direction sees the opposite direction and charge
  auto fit = [&]() {
    Acts::Result<Acts::FreeVector> res =
        estimateFreeParams(cache.positions, field, t0, cache.weights,
                           m_cfg.geometricRefineIterations);
    if (res.ok() && anchor == Anchor::Last) {
      res->segment<3>(Acts::eFreeDir0) *= -1;
      (*res)[Acts::eFreeQOverP] *= -1;
    }
    return res;
  };

  Acts::Result<Acts::FreeVector> free =
      hasFreeParams ? Acts::Result<Acts::FreeVector>::success(*freeParams)
                    : fit();
  // given free parameters lie at the inner end of the seed
  if (hasFreeParams) {
    (*free)[Acts::eFreeTime] = timeOf(
        anchor == Anchor::First ? anchorSp : spacePoints.at(candidates.back()));
  }
  if (!free.ok()) {
    ACTS_DEBUG("Seed could not be fitted: " << free.error().message());
    return Result::failure(Failure::Fit);
  }
  // degenerate space points make the estimate not a number, a straight line
  // fit leaves q/p at zero, which the covariance transport cannot take
  if (!free->allFinite() || (*free)[Acts::eFreeQOverP] == 0) {
    return Result::failure(Failure::Degenerate);
  }

  const SeedPropagator propagator(Acts::SympyStepper(m_cfg.magneticField),
                                  Acts::VoidNavigator(), m_propagatorLogger);
  const SeedPropagatorOptions propagatorOptions(gctx, mctx);

  Acts::Result<Acts::BoundVector> boundParams = transportToSurface(
      propagator, propagatorOptions, *free, *anchorSurface, hypothesis);
  // given free parameters which do not reach the surface give way to the fit
  if (!boundParams.ok() && hasFreeParams) {
    free = fit();
    if (free.ok() && free->allFinite() && (*free)[Acts::eFreeQOverP] != 0) {
      boundParams = transportToSurface(propagator, propagatorOptions, *free,
                                       *anchorSurface, hypothesis);
    }
  }
  if (!boundParams.ok()) {
    ACTS_DEBUG(
        "Seed could not be transported to the surface of its anchor "
        "space point: "
        << boundParams.error().message());
    return Result::failure(Failure::Transport);
  }

  const Acts::EstimateTrackParamCovarianceConfig covConfig{
      .initialSigmas =
          Eigen::Map<const Acts::BoundVector>{m_cfg.initialSigmas.data()},
      .initialSigmaQoverPt = m_cfg.initialSigmaQoverPt,
      .initialSigmaPtRel = m_cfg.initialSigmaPtRel,
      .initialVarInflation =
          Eigen::Map<const Acts::BoundVector>{m_cfg.initialVarInflation.data()},
      .noTimeVarInflation = m_cfg.noTimeVarInflation};

  const Acts::BoundMatrix cov = Acts::estimateTrackParamCovariance(
      covConfig, *boundParams, !std::isnan(anchorSp.time()));

  return Result::success(TrackParameters(anchorSurface->getSharedPtr(),
                                         *boundParams, cov, hypothesis));
}

}  // namespace ActsExamples
