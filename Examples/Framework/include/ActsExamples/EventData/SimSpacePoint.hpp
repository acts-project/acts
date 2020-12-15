// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "ActsExamples/EventData/Index.hpp"

#include <cmath>

namespace ActsExamples {

/// Space point representation of a measurement suitable for track seeding.
class SimSpacePoint {
 public:
  template <typename position_t, typename covariance_t>
  SimSpacePoint(const Eigen::MatrixBase<position_t>& pos,
                const Eigen::MatrixBase<covariance_t>& cov,
                Index measurementIndex)
      : m_x(pos[Acts::ePos0]),
        m_y(pos[Acts::ePos1]),
        m_z(pos[Acts::ePos2]),
        m_r(std::hypot(m_x, m_y)),
        //      r = sqrt(x² + y²)
        //  dr/dx = (1 / sqrt(x² + y²)) * 2 * x = 2 * x / r
        // var(r) = (dr/dx)² * var(x) + (dr/dy)² * var(y)
        //        = (4 / r²) * (x² * var(x) + y*2 * var(y))
        m_varianceR(
            (4 / (m_r * m_r)) *
            (m_x * m_x * static_cast<float>(cov(Acts::ePos0, Acts::ePos0)) +
             m_y * m_y * static_cast<float>(cov(Acts::ePos1, Acts::ePos1)))),
        m_varianceZ(cov(Acts::ePos2, Acts::ePos2)),
        m_measurementIndex(measurementIndex) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(position_t, 3);
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(covariance_t, 3, 3);
  }

  constexpr float x() const { return m_x; }
  constexpr float y() const { return m_y; }
  constexpr float z() const { return m_z; }
  constexpr float r() const { return m_r; }
  constexpr float varianceR() const { return m_varianceR; }
  constexpr float varianceZ() const { return m_varianceZ; }

  constexpr Index measurementIndex() const { return m_measurementIndex; }

 private:
  // Global position
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  // VarianceR/Z of the global position
  float m_varianceR;
  float m_varianceZ;
  // Index of the corresponding measurement
  Index m_measurementIndex;
};

constexpr bool operator==(const SimSpacePoint& lhs, const SimSpacePoint& rhs) {
  // TODO would it be sufficient to check just the index under the assumption
  //   that the same measurement index always produces the same space point?
  // no need to check r since it is fully defined by x/y
  return (lhs.measurementIndex() == rhs.measurementIndex()) and
         (lhs.x() == rhs.x()) and (lhs.y() == rhs.y()) and
         (lhs.z() == rhs.z()) and (lhs.varianceR() == rhs.varianceR()) and
         (lhs.varianceZ() == rhs.varianceZ());
}

}  // namespace ActsExamples
