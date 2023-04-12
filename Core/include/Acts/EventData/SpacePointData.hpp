// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <limits>
#include <vector>

namespace Acts {

/// @class SpacePointData
/// This class contains auxiliary and mutable data associated to the
/// external space points provided by the customers
/// These variables are used internally by the seeding algorithm, that
/// reads and updates them
/// The variables collected here are also dynamic variables only present
/// for strip space points
class SpacePointData {
 public:
  /// @brief Default constructor
  SpacePointData() = default;

  /// No copies
  SpacePointData(const SpacePointData& other) = delete;
  SpacePointData& operator=(const SpacePointData& other) = delete;

  /// @brief Move operations
  SpacePointData(SpacePointData&& other) noexcept = default;
  SpacePointData& operator=(SpacePointData&& other) noexcept = default;

  /// @brief Destructor
  ~SpacePointData() = default;

  /// @brief Getters
  const float& quality(std::size_t idx) const;
  const float& deltaR(std::size_t idx) const;

  /// @brief Setters
  void setQuality(std::size_t idx, const float& value);
  void setDeltaR(std::size_t idx, const float& value);

  /// @brief Resize vectors
  void resize(std::size_t n, bool resizeDynamic = false);

  /// @brief clear vectors
  void clear();

  ///
  bool hasDynamicVariable() const { return not m_topHalfStripLength.empty(); }

  const float& getTopHalfStripLength(std::size_t idx) const {
    return m_topHalfStripLength[idx];
  }

  const float& getBottomHalfStripLength(std::size_t idx) const {
    return m_bottomHalfStripLength[idx];
  }

  const Acts::Vector3& getTopStripDirection(std::size_t idx) const {
    return m_topStripDirection[idx];
  }

  const Acts::Vector3& getBottomStripDirection(std::size_t idx) const {
    return m_bottomStripDirection[idx];
  }

  const Acts::Vector3& getStripCenterDistance(std::size_t idx) const {
    return m_stripCenterDistance[idx];
  }

  const Acts::Vector3& getTopStripCenterPosition(std::size_t idx) const {
    return m_topStripCenterPosition[idx];
  }

  void setTopHalfStripLength(std::size_t idx, const float& value) {
    m_topHalfStripLength[idx] = value;
  }

  void setBottomHalfStripLength(std::size_t idx, const float& value) {
    m_bottomHalfStripLength[idx] = value;
  }

  void setTopStripDirection(std::size_t idx, const Acts::Vector3& value) {
    m_topStripDirection[idx] = value;
  }

  void setBottomStripDirection(std::size_t idx, const Acts::Vector3& value) {
    m_bottomStripDirection[idx] = value;
  }

  void setStripCenterDistance(std::size_t idx, const Acts::Vector3& value) {
    m_stripCenterDistance[idx] = value;
  }

  void setTopStripCenterPosition(std::size_t idx, const Acts::Vector3& value) {
    m_topStripCenterPosition[idx] = value;
  }

 private:
  /// Mutable variables
  std::vector<float> m_quality{};
  std::vector<float> m_deltaR{};

  /// dynamic variables
  std::vector<float> m_topHalfStripLength{};
  std::vector<float> m_bottomHalfStripLength{};
  std::vector<Acts::Vector3> m_topStripDirection{};
  std::vector<Acts::Vector3> m_bottomStripDirection{};
  std::vector<Acts::Vector3> m_stripCenterDistance{};
  std::vector<Acts::Vector3> m_topStripCenterPosition{};
};

inline const float& SpacePointData::quality(std::size_t idx) const {
  return m_quality[idx];
}

inline const float& SpacePointData::deltaR(std::size_t idx) const {
  return m_deltaR[idx];
}

inline void SpacePointData::setQuality(std::size_t idx, const float& value) {
  if (value > m_quality[idx]) {
    m_quality[idx] = value;
  }
}

inline void SpacePointData::setDeltaR(std::size_t idx, const float& value) {
  m_deltaR[idx] = value;
}

inline void SpacePointData::resize(std::size_t n, bool resizeDynamic) {
  clear();

  m_quality.resize(n, -std::numeric_limits<float>::infinity());
  m_deltaR.resize(n, 0.);

  if (resizeDynamic) {
    m_topHalfStripLength.resize(n, 0.);
    m_bottomHalfStripLength.resize(n, 0.);
    m_topStripDirection.resize(n, {0, 0, 0});
    m_bottomStripDirection.resize(n, {0, 0, 0});
    m_stripCenterDistance.resize(n, {0, 0, 0});
    m_topStripCenterPosition.resize(n, {0, 0, 0});
  }
}

inline void SpacePointData::clear() {
  // mutable variables
  m_quality.clear();
  m_deltaR.clear();
  // dynamicvariables
  m_topHalfStripLength.clear();
  m_bottomHalfStripLength.clear();
  m_topStripDirection.clear();
  m_bottomStripDirection.clear();
  m_stripCenterDistance.clear();
  m_topStripCenterPosition.clear();
}

}  // namespace Acts
