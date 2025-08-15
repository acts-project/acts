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

  auto result = fit(fillAuxiliaries(measurements, signs));
  if (!result) {
    return std::nullopt;
  }
  /// Calculate the chi2
  const double cosTheta{std::cos(result->theta)},
      sinTheta{std::sin(result->theta)};
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    if (!strawMeas->isStraw()) {
      continue;
    }
    const double cov =
        strawMeas->covariance()[toUnderlying(ResidualIdx::bending)];
    if (cov < std::numeric_limits<double>::epsilon()) {
      continue;
    }
    const Vector& pos = strawMeas->localPosition();
    const double y = pos.dot(strawMeas->toNextSensor());
    const double z = pos.dot(strawMeas->planeNormal());
    const double dist = Acts::abs((y - result->y0) * cosTheta - z * sinTheta);
    ACTS_VERBOSE("chi2 calculation -  Distance straw ("
                 << y << ", " << z << "), r: " << strawMeas->driftRadius()
                 << " - track: " << dist);
    result->chi2 += Acts::pow(dist - strawMeas->driftRadius(), 2) / cov;
  }
  ACTS_VERBOSE("Overall chi2: " << result->chi2 << ", nDoF: " << result->nDoF
                                << ", redChi2: "
                                << (result->chi2 / result->nDoF));
  return result;
}

template <CompositeSpacePointContainer StrawCont_t>
FastStrawLineFitter::FitAuxiliaries FastStrawLineFitter::fillAuxiliaries(
    const StrawCont_t& measurements,
    const std::vector<std::int32_t>& signs) const {
  FitAuxiliaries auxVars{};
  std::vector<double> invCovs(signs.size(), -1.);
  Vector centerOfGravity{Vector::Zero()};

  /// Calculate first the center of gravity

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
    centerOfGravity += invCov * strawMeas->localPosition();
    ++auxVars.nDoF;
  }
  if (auxVars.nDoF < 3) {
    ACTS_WARNING(
        "At least 3 measurements are required to perform the straw line ift");
    return auxVars;
  }
  /// Reduce the number of degrees of freedom by 2 to account for the two free
  /// parameters
  auxVars.nDoF -= 2u;
  auxVars.covNorm = 1. / auxVars.covNorm;
  centerOfGravity *= auxVars.covNorm;

  // Now calculate the fit constants
  bool centerSet{false};
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    const auto& invCov = invCovs[sIdx];
    /// The invalid measurements were marked
    if (invCov < 0.) {
      continue;
    }
    if (!centerSet) {
      auxVars.centerY = centerOfGravity.dot(strawMeas->toNextSensor());
      auxVars.centerZ = centerOfGravity.dot(strawMeas->planeNormal());
      centerSet = true;
    }
    const Vector pos = strawMeas->localPosition() - centerOfGravity;
    const double y = pos.dot(strawMeas->toNextSensor());
    const double z = pos.dot(strawMeas->planeNormal());
    const double r = strawMeas->driftRadius();

    auxVars.T_zzyy += invCov * (Acts::pow(z, 2) - Acts::pow(y, 2));
    auxVars.T_yz += invCov * z * y;
    const double sInvCov = -invCov * signs[sIdx];
    auxVars.T_rz += sInvCov * z * r;
    auxVars.T_ry += sInvCov * y * r;
    auxVars.fitY0 += sInvCov * r;
  }
  auxVars.fitY0 *= auxVars.covNorm;

  return auxVars;
}

}  // namespace Acts::Experimental::detail
