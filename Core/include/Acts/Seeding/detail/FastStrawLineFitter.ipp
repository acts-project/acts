// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"

#include "Acts/Utilities/Enumerate.hpp"
namespace Acts::Experimental::detail {

template <CompositeSpacePointContainer StrawCont_t>
std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const StrawCont_t& measurements,
    const std::vector<std::int32_t>& signs) const {
  if (measurements.size() != signs.size()) {
    ACTS_WARNING("Not all measurements are associated with a drift sign");
    return std::nullopt;
  }

  return fit(fillAuxiliaries(measurements, signs));
}

template <CompositeSpacePointContainer StrawCont_t>
FastStrawLineFitter::FitAuxiliaries FastStrawLineFitter::fillAuxiliaries(
    const StrawCont_t& measurements,
    const std::vector<std::int32_t>& signs) const {
  FitAuxiliaries auxVars{};
  std::vector<double> invCovs(signs.size(), -1.);

  /// Calculate first the center of gravity
  std::uint32_t nValid{0};
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    if (!strawMeas->isStraw()) {
      ACTS_WARNING("The measurement is not a straw");
      continue;
    }
    const double cov =
        strawMeas->covariance()[toUnderlying(ResidualIdx::bending)];
    if (cov < std::numeric_limits<double>::epsilon()) {
      ACTS_WARNING("The covariance (" << cov
                                      << ") of the measurement is invalid.");
      continue;
    }
    auto& invCov = (invCovs[sIdx] = 1. / cov);
    auxVars.covNorm += invCov;
    auxVars.centerOfGrav += invCov * strawMeas->localPosition();
    ++nValid;
  }
  if (nValid < 3) {
    ACTS_WARNING(
        "At least 3 measurements are required to perform the straw line ift");
    return auxVars;
  }
  auxVars.covNorm = 1. / auxVars.covNorm;
  auxVars.centerOfGrav *= auxVars.covNorm;

  // Now calculate the fit constants
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    const Vector pos = strawMeas->localPosition() - auxVars.centerOfGrav;
    const double y = pos.dot(strawMeas->toNextSensor());
    const double z = pos.dot(strawMeas->planeNormal());
    const double r = strawMeas->driftRadius();
    const auto& invCov = invCovs[sIdx];
    const double sInvCov = invCov * signs[sIdx];

    auxVars.T_zzyy += invCov * (Acts::pow(z, 2) - Acts::pow(y, 2));
    auxVars.T_yz += invCov * z * y;

    auxVars.T_rz += sInvCov * z * r;
    auxVars.T_ry += sInvCov * y * r;
    auxVars.fitY0 += sInvCov * r;
  }
  auxVars.fitY0 *= auxVars.covNorm;

  return auxVars;
}

}  // namespace Acts::Experimental::detail
