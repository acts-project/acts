// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/ParametersOnSurface.hpp"

#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <utility>

std::optional<Acts::BoundTrackParameters>
ActsExamples::truthParametersOnSurface(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    Index measurementIndex, const SimParticle& particle,
    const SimHitContainer& simHits,
    const MeasurementSimHitsMap& measurementSimHitsMap,
    const Acts::Logger& logger) {
  using Acts::VectorHelpers::phi;
  using Acts::VectorHelpers::theta;

  using enum Acts::BoundIndices;

  const auto indices =
      makeRange(measurementSimHitsMap.equal_range(measurementIndex));
  if (indices.empty()) {
    ACTS_WARNING("No truth hits associated to measurement " << measurementIndex
                                                            << " found");
    return std::nullopt;
  }

  const auto [truthLocal, truthPos4, truthUnitDir] =
      averageSimHits(gctx, surface, simHits, indices, logger);

  // position, direction, and time are averaged over the hits above. the
  // momentum is not: an average over hits of different particles is not the
  // momentum of any of them. take the first hit instead, which exists because
  // the range was checked to be non-empty.
  const auto simHitIdx0 = indices.begin()->second;
  const auto& simHit0 = *simHits.nth(simHitIdx0);
  const auto momentum = simHit0.momentum4Before().segment<3>(Acts::eMom0);

  Acts::BoundVector params = Acts::BoundVector::Zero();
  params[eBoundLoc0] = truthLocal[Acts::ePos0];
  params[eBoundLoc1] = truthLocal[Acts::ePos1];
  params[eBoundPhi] = phi(truthUnitDir);
  params[eBoundTheta] = theta(truthUnitDir);
  params[eBoundQOverP] =
      particle.hypothesis().qOverP(momentum.norm(), particle.charge());
  params[eBoundTime] = truthPos4[Acts::eTime];

  return Acts::BoundTrackParameters(surface.getSharedPtr(), params,
                                    std::nullopt, particle.hypothesis());
}

std::optional<Acts::BoundTrackParameters> ActsExamples::recoParametersOnSurface(
    const ConstTrackStateProxy& state,
    std::optional<TrackParameterType> parameterType,
    const Acts::ParticleHypothesis& hypothesis) {
  using enum TrackParameterType;

  if (!state.hasReferenceSurface()) {
    return std::nullopt;
  }

  const auto stateParameters =
      [&]() -> std::optional<std::pair<Acts::BoundVector, Acts::BoundMatrix>> {
    if (!parameterType.has_value()) {
      if (!state.hasSmoothed() && !state.hasFiltered() &&
          !state.hasPredicted()) {
        return std::nullopt;
      }
      // the choice is the proxy's: `parameters()` returns smoothed, else
      // filtered, else predicted. `Unbiased` is not among them and has to be
      // requested explicitly.
      return std::pair(state.parameters(), state.covariance());
    }
    if (parameterType == Predicted && state.hasPredicted()) {
      return std::pair(state.predicted(), state.predictedCovariance());
    }
    if (parameterType == Filtered && state.hasFiltered()) {
      return std::pair(state.filtered(), state.filteredCovariance());
    }
    if (parameterType == Smoothed && state.hasSmoothed()) {
      return std::pair(state.smoothed(), state.smoothedCovariance());
    }
    if (parameterType == Unbiased && state.hasSmoothed() &&
        state.hasProjector() && state.hasCalibrated()) {
      return Acts::calculateUnbiasedParametersCovariance(
          Acts::AnyConstTrackStateProxy{state});
    }
    return std::nullopt;
  }();

  if (!stateParameters.has_value()) {
    return std::nullopt;
  }

  return Acts::BoundTrackParameters(state.referenceSurface().getSharedPtr(),
                                    stateParameters->first,
                                    stateParameters->second, hypothesis);
}

bool ActsExamples::parametersUseOwnMeasurement(
    TrackParameterType parameterType) {
  using enum TrackParameterType;

  return parameterType == Filtered || parameterType == Smoothed;
}

std::optional<ActsExamples::MeasurementResidual>
ActsExamples::measurementResidual(const ConstTrackStateProxy& state,
                                  const Acts::BoundTrackParameters& parameters,
                                  TrackParameterType parameterType) {
  if (!state.hasCalibrated() || !state.hasProjector()) {
    return std::nullopt;
  }

  const Acts::VariableBoundSubspaceHelper subspace =
      state.projectorSubspaceHelper();
  if (subspace.empty()) {
    return std::nullopt;
  }

  // the projector only selects bound indices, so `H x` is `x(subspace)` and
  // `H P H^T` is `P(subspace, subspace)`. no matrix products needed.
  const Acts::DynamicVector measurement = state.effectiveCalibrated();
  const Acts::DynamicMatrix measurementCovariance =
      state.effectiveCalibratedCovariance();

  const Acts::BoundVector& parameterVector = parameters.parameters();
  const Acts::BoundMatrix parameterCovariance =
      parameters.covariance().value_or(Acts::BoundMatrix::Zero());

  // parameters that used the measurement are pulled towards it, so their
  // covariance is correlated with the measurement and has to be subtracted
  // instead of added
  const double sign = parametersUseOwnMeasurement(parameterType) ? -1. : 1.;

  const std::size_t size = subspace.size();
  Acts::DynamicVector residual(size);
  Acts::DynamicMatrix covariance = measurementCovariance;
  for (std::size_t i = 0; i < size; ++i) {
    // the reference is subtracted from the reconstructed value, matching the
    // convention of the truth residuals
    residual[i] = parameterVector[subspace[i]] - measurement[i];

    for (std::size_t j = 0; j < size; ++j) {
      covariance(i, j) += sign * parameterCovariance(subspace[i], subspace[j]);
    }
  }

  return MeasurementResidual{subspace, subspace.expandVector(residual),
                             subspace.expandMatrix(covariance)};
}
