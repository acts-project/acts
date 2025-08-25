// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <format>
namespace Acts::Experimental::detail {

template <CompositeSpacePointContainer StrawCont_t>
std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const StrawCont_t& measurements, const std::vector<int>& signs) const {
  if (measurements.size() != signs.size()) {
    ACTS_WARNING(
        __func__ << "() - " << __LINE__
                 << ": Not all measurements are associated with a drift sign");
    return std::nullopt;
  }

  auto result = fit(fillAuxiliaries(measurements, signs));
  if (!result) {
    return std::nullopt;
  }
  /// Calculate the chi2
  calcPostFitChi2(measurements, *result);
  return result;
}

template <CompositeSpacePointContainer StrawCont_t>
void FastStrawLineFitter::calcPostFitChi2(const StrawCont_t& measurements,
                                          FitResult& result) const {
  const TrigonomHelper angles{result.theta};
  result.chi2 = 0.;
  for (const auto& strawMeas : measurements) {
    result.chi2 += chi2Term(angles, result.y0, *strawMeas);
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Overall chi2: "
                      << result.chi2 << ", nDoF: " << result.nDoF
                      << ", redChi2: " << (result.chi2 / result.nDoF));
}

template <CompositeSpacePoint Point_t>
double FastStrawLineFitter::chi2Term(const TrigonomHelper& angle,
                                     const double y0, const Point_t& strawMeas,
                                     std::optional<double> r) const {
  if (!strawMeas.isStraw()) {
    return 0.;
  }
  const double cov = strawMeas.covariance()[s_covIdx];
  if (cov < std::numeric_limits<double>::epsilon()) {
    return 0.;
  }
  const Vector& pos = strawMeas.localPosition();
  const double y = pos.dot(strawMeas.toNextSensor());
  const double z = pos.dot(strawMeas.planeNormal());
  const double dist = Acts::abs((y - y0) * angle.cosTheta - z * angle.sinTheta);
  ACTS_VERBOSE(__func__ << "() - " << __LINE__ << ": Distance straw (" << y
                        << ", " << z
                        << "), r: " << r.value_or(strawMeas.driftRadius())
                        << " - track: " << dist);
  return Acts::pow(dist - r.value_or(strawMeas.driftRadius()), 2) / cov;
}

template <CompositeSpacePointContainer StrawCont_t>
FastStrawLineFitter::FitAuxiliaries FastStrawLineFitter::fillAuxiliaries(
    const StrawCont_t& measurements, const std::vector<int>& signs) const {
  FitAuxiliaries auxVars{};
  auxVars.invCovs.resize(signs.size(), -1.);
  Vector centerOfGravity{Vector::Zero()};

  /// Calculate first the center of gravity
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    if (!strawMeas->isStraw()) {
      ACTS_WARNING(__func__ << "() - " << __LINE__
                            << ": The measurement is not a straw");
      continue;
    }
    const double cov = strawMeas->covariance()[s_covIdx];
    if (cov < std::numeric_limits<double>::epsilon()) {
      ACTS_WARNING(__func__ << "() - " << __LINE__ << ": The covariance ("
                            << cov << ") of the measurement is invalid.");
      continue;
    }
    auto& invCov = (auxVars.invCovs[sIdx] = 1. / cov);
    auxVars.covNorm += invCov;
    centerOfGravity += invCov * strawMeas->localPosition();
    ++auxVars.nDoF;
  }
  if (auxVars.nDoF < 3) {
    std::stringstream sstr{};
    for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
      sstr << " --- " << (sIdx + 1) << ") "
           << toString(strawMeas->localPosition())
           << ", r: " << strawMeas->driftRadius()
           << ", weight: " << auxVars.invCovs[sIdx] << std::endl;
    }
    ACTS_WARNING(__func__ << "() - " << __LINE__
                          << ": At least 3 measurements are required to "
                             "perform the straw line fit\n"
                          << sstr.str());
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
    const auto& invCov = auxVars.invCovs[sIdx];
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

    auxVars.T_zzyy += invCov * (Acts::square(z) - Acts::square(y));
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
