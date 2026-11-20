// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/TrackFinding/ScatteringTail.hpp"

#include <cmath>
#include <optional>

namespace Acts::detail {

/// The scattering tail of the material since the last measurement of the CKF.
///
/// Each scatterer i adds `eps_i D_i`, with `D_i = (k_i^2 - 1) Q_i` the excess
/// of its tail over its core noise `Q_i`. At most one scatterer is assumed in
/// its tail, so the tail covariance is `sum eps_i D_i / sum eps_i`.
class ScatteringTailState {
 public:
  /// Forget all scatterers.
  void reset() {
    m_sum.setZero();
    m_weightSum = 0;
    m_logNoTail = 0;
  }

  /// Transport with the jacobian of one segment.
  void transport(const BoundMatrix& jacobian) {
    if (m_weightSum > 0) {
      m_sum = jacobian * m_sum * jacobian.transpose();
    }
  }

  /// Add one scatterer with its core noise in phi and theta.
  void add(const ScatteringTailComponent& tail, double variancePhi,
           double varianceTheta) {
    if (tail.weight <= 0 || tail.widthScale <= 1) {
      return;
    }
    const double excess = tail.weight * (tail.widthScale * tail.widthScale - 1);
    m_sum(eBoundPhi, eBoundPhi) += excess * variancePhi;
    m_sum(eBoundTheta, eBoundTheta) += excess * varianceTheta;
    m_weightSum += tail.weight;
    m_logNoTail += std::log1p(-tail.weight);
  }

  /// The tail prior, unless it adds less than @p gate times the predicted
  /// variance to both local coordinates.
  std::optional<ScatteringTailPrior> prior(
      const BoundMatrix& predictedCovariance, double gate) const {
    const auto below = [&](BoundIndices i) {
      return m_sum(i, i) < gate * m_weightSum * predictedCovariance(i, i);
    };
    if (m_weightSum == 0 || (below(eBoundLoc0) && below(eBoundLoc1))) {
      return std::nullopt;
    }
    return ScatteringTailPrior{.weight = -std::expm1(m_logNoTail),
                               .covariance = m_sum / m_weightSum};
  }

 private:
  BoundMatrix m_sum = BoundMatrix::Zero();
  double m_weightSum = 0;
  double m_logNoTail = 0;
};

}  // namespace Acts::detail
