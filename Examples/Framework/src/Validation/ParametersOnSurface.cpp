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

  // momentum averaging makes even less sense than averaging position and
  // direction. use the first momentum
  // we assume that the indices are within valid ranges so we do not need to
  // check their validity again.
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
      // best available parameters, i.e. smoothed, filtered, or predicted
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
      return Acts::calculateUnbiasedParametersCovariance(state);
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
