// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <utility>

namespace Acts {

/// @addtogroup track_finding
/// @{

/// The non-Gaussian tail of multiple scattering in one material slab, as a
/// second Gaussian component of the projected scattering angle.
struct ScatteringTailComponent {
  /// Probability of the tail component
  double weight{};
  /// Width of the tail component in units of the Gaussian core width
  double widthScale{};
};

/// Model of the scattering tail of one material slab.
using ScatteringTailModel = Delegate<ScatteringTailComponent(
    const MaterialSlab& slab, const ParticleHypothesis& particleHypothesis,
    double qOverP)>;

/// Scattering tail in silicon, fitted to the general mixture model of Fatras
/// for Z = 14 and beta = 1.
///
/// @param slab the material slab, only its thickness in X0 is used
/// @return the tail component
ScatteringTailComponent siliconScatteringTail(
    const MaterialSlab& slab, const ParticleHypothesis& /*particleHypothesis*/,
    double /*qOverP*/);

/// The scattering tail since the last measurement as the second component of
/// the predicted parameters: `(1 - weight) N(C) + weight N(C + covariance)`.
struct ScatteringTailPrior {
  /// Prior probability of the tail component
  double weight{};
  /// Covariance of the tail component in excess of the core prediction C
  BoundMatrix covariance = BoundMatrix::Zero();
};

/// The chi2 of a measurement under the two-component prediction,
/// `chi2_core - 2 softplus(a)` with `a` the log likelihood ratio of the tail
/// over the core. It is below the core chi2 and close to it for a small
/// residual, so the core chi2 cuts keep their meaning.
///
/// @param residual the residual to the core prediction
/// @param residualCovariance the residual covariance of the core prediction
/// @param tailCovariance the projected covariance excess of the tail
/// @param logPriorRatio `ln(weight / (1 - weight))` of the prior
/// @param rejectAbove returns `{rejectAbove, -inf}` without transcendental
///        functions if a lower bound of the chi2 is above this value
/// @return the chi2 and `a`, the posterior tail weight is `1 / (1 + exp(-a))`
template <int N>
std::pair<double, double> scatteringTailChi2(
    const Vector<N>& residual, const SquareMatrix<N>& residualCovariance,
    const SquareMatrix<N>& tailCovariance, double logPriorRatio,
    double rejectAbove = std::numeric_limits<double>::infinity()) {
  const SquareMatrix<N> tailCovarianceSum = residualCovariance + tailCovariance;
  const double chi2Core = residual.dot(residualCovariance.inverse() * residual);
  const double chi2Tail = residual.dot(tailCovarianceSum.inverse() * residual);
  // a <= aMax because det core <= det tail, and softplus(a) <= max(a, 0) + ln 2
  const double aMax = logPriorRatio + 0.5 * (chi2Core - chi2Tail);
  if (chi2Core - 2 * (std::max(aMax, 0.) + std::numbers::ln2) >= rejectAbove) {
    return {rejectAbove, -std::numeric_limits<double>::infinity()};
  }
  const double a = aMax + 0.5 * std::log(residualCovariance.determinant() /
                                         tailCovarianceSum.determinant());
  const double softplus =
      a > 0 ? a + std::log1p(std::exp(-a)) : std::log1p(std::exp(a));
  return {chi2Core - 2 * softplus, a};
}

/// The chi2 of a calibrated track state under the two-component prediction.
///
/// @param trackState the track state with the core prediction
/// @param prior the scattering tail
/// @param rejectAbove see the overload above
/// @return the chi2 and `a`, see the overload above
template <typename track_state_t>
std::pair<double, double> scatteringTailChi2(
    const track_state_t& trackState, const ScatteringTailPrior& prior,
    double rejectAbove = std::numeric_limits<double>::infinity()) {
  return visit_measurement(
      trackState.calibratedSize(), [&](auto N) -> std::pair<double, double> {
        constexpr std::size_t kSize = decltype(N)::value;
        const auto projector = trackState.projectorSubspaceIndices();
        const FixedBoundSubspaceHelper<kSize> subspace(
            std::span<const std::uint8_t, kSize>(projector.begin(),
                                                 projector.begin() + kSize));
        const typename TrackStateTraits<kSize, true>::Calibrated calibrated{
            trackState.effectiveCalibrated().data()};
        const typename TrackStateTraits<kSize, true>::CalibratedCovariance
            calibratedCovariance{
                trackState.effectiveCalibratedCovariance().data()};
        return scatteringTailChi2<static_cast<int>(kSize)>(
            calibrated - subspace.projectVector(trackState.predicted()),
            calibratedCovariance +
                subspace.projectMatrix(trackState.predictedCovariance()),
            subspace.projectMatrix(prior.covariance),
            std::log(prior.weight / (1 - prior.weight)), rejectAbove);
      });
}

/// @}

}  // namespace Acts
