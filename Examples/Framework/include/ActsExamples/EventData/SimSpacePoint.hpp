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
#include <vector>

namespace ActsExamples {

/// Space point representation of a measurement suitable for track seeding.
class SimSpacePoint {
 public:
  /// Construct the space point from global position and selected variances.
  ///
  /// @tparam position_t Input position type
  /// @param pos Global position
  /// @param varRho Measurement variance of the global transverse distance
  /// @param varZ Measurement variance of the global longitudinal position
  /// @param measurementIndex Index of the underlying measurement
  /// @param topHalfStripLength half of the length of the top strip
  /// @param bottomHalfStripLength half of the length of the bottom strip
  /// @param topStripDirection direction of the top strip
  /// @param bottomStripDirection direction of the bottom strip
  /// @param stripCenterDistance distance between the center of the two strips
  /// @param topStripCenterPosition position of the center of the top strip
  /// @param validDoubleMeasurementDetails boolean to check if double measurements are valid
  template <typename position_t>
  SimSpacePoint(const Eigen::MatrixBase<position_t>& pos, float varRho,
                float varZ, Index measurementIndex,
                const float topHalfStripLength,
                const float bottomHalfStripLength,
                const Acts::Vector3 topStripDirection,
                const Acts::Vector3 bottomStripDirection,
                const Acts::Vector3 stripCenterDistance,
                const Acts::Vector3 topStripCenterPosition)
      : m_x(pos[Acts::ePos0]),
        m_y(pos[Acts::ePos1]),
        m_z(pos[Acts::ePos2]),
        m_rho(std::hypot(m_x, m_y)),
        m_varianceRho(varRho),
        m_varianceZ(varZ),
        m_measurementIndex(measurementIndex),
        m_topHalfStripLength(topHalfStripLength),
        m_bottomHalfStripLength(bottomHalfStripLength),
        m_topStripDirection(topStripDirection),
        m_bottomStripDirection(bottomStripDirection),
        m_stripCenterDistance(stripCenterDistance),
        m_topStripCenterPosition(topStripCenterPosition),
        m_validDoubleMeasurementDetails(true) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(position_t, 3);
  }

  /// Construct the space point from global position and selected variances.
  ///
  /// @tparam position_t Input position type
  /// @param pos Global position
  /// @param varRho Measurement variance of the global transverse distance
  /// @param varZ Measurement variance of the global longitudinal position
  /// @param measurementIndex Index of the underlying measurement
  template <typename position_t>
  SimSpacePoint(const Eigen::MatrixBase<position_t>& pos, float varRho,
                float varZ, Index measurementIndex)
      : m_x(pos[Acts::ePos0]),
        m_y(pos[Acts::ePos1]),
        m_z(pos[Acts::ePos2]),
        m_rho(std::hypot(m_x, m_y)),
        m_varianceRho(varRho),
        m_varianceZ(varZ),
        m_measurementIndex(measurementIndex) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(position_t, 3);
  }

  constexpr float x() const { return m_x; }
  constexpr float y() const { return m_y; }
  constexpr float z() const { return m_z; }
  constexpr float r() const { return m_rho; }
  constexpr float varianceR() const { return m_varianceRho; }
  constexpr float varianceZ() const { return m_varianceZ; }

  constexpr Index measurementIndex() const { return m_measurementIndex; }

  constexpr float topHalfStripLength() const { return m_topHalfStripLength; }
  constexpr float bottomHalfStripLength() const {
    return m_bottomHalfStripLength;
  }
  Acts::Vector3 topStripDirection() const { return m_topStripDirection; }
  Acts::Vector3 bottomStripDirection() const { return m_bottomStripDirection; }
  Acts::Vector3 stripCenterDistance() const { return m_stripCenterDistance; }
  Acts::Vector3 topStripCenterPosition() const {
    return m_topStripCenterPosition;
  }
  constexpr bool validDoubleMeasurementDetails() const {
    return m_validDoubleMeasurementDetails;
  }

 private:
  // Global position
  float m_x;
  float m_y;
  float m_z;
  float m_rho;
  // Variance in rho/z of the global coordinates
  float m_varianceRho;
  float m_varianceZ;
  // Index of the corresponding measurement
  Index m_measurementIndex;

  // half of the length of the top strip
  float m_topHalfStripLength = 0;
  // half of the length of the bottom strip
  float m_bottomHalfStripLength = 0;
  // direction of the top strip
  Acts::Vector3 m_topStripDirection = {0, 0, 0};
  // direction of the bottom strip
  Acts::Vector3 m_bottomStripDirection = {0, 0, 0};
  // distance between the center of the two strips
  Acts::Vector3 m_stripCenterDistance = {0, 0, 0};
  // position of the center of the bottom strip
  Acts::Vector3 m_topStripCenterPosition = {0, 0, 0};
  bool m_validDoubleMeasurementDetails = false;
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

/// Container of space points.
using SimSpacePointContainer = std::vector<SimSpacePoint>;

}  // namespace ActsExamples
