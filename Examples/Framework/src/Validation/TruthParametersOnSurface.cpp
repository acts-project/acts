// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TruthParametersOnSurface.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/Utilities/Range.hpp"

namespace ActsExamples {

std::optional<TruthParametersOnSurface> truthParametersOnSurface(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    Index measurementIndex, double truthCharge, const SimHitContainer& simHits,
    const MeasurementSimHitsMap& measurementSimHitsMap,
    const Acts::Logger& logger) {
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
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

  TruthParametersOnSurface truth;
  truth.params[eBoundLoc0] = truthLocal[Acts::ePos0];
  truth.params[eBoundLoc1] = truthLocal[Acts::ePos1];
  truth.params[eBoundPhi] = phi(truthUnitDir);
  truth.params[eBoundTheta] = theta(truthUnitDir);
  truth.params[eBoundQOverP] = truthCharge / momentum.norm();
  truth.params[eBoundTime] = truthPos4[Acts::eTime];
  truth.eta = eta(truthUnitDir);
  truth.phi = truth.params[eBoundPhi];
  truth.pt = perp(momentum);
  truth.charge = static_cast<int>(truthCharge);
  return truth;
}

}  // namespace ActsExamples
