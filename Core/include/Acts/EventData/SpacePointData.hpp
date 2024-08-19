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
  float quality(std::size_t idx) const;
  float deltaR(std::size_t idx) const;

  /// @brief Setters
  void setQuality(std::size_t idx, const float value);
  void setDeltaR(std::size_t idx, const float value);

  /// @brief Resize vectors
  void resize(std::size_t n, bool resizeDynamic = false);

  /// @brief clear vectors
  void clear();

  ///
  bool hasDynamicVariable() const { return !m_topStripVector.empty(); }

  const Acts::Vector3& getTopStripVector(std::size_t idx) const {
    return m_topStripVector[idx];
  }

  const Acts::Vector3& getBottomStripVector(std::size_t idx) const {
    return m_bottomStripVector[idx];
  }

  const Acts::Vector3& getStripCenterDistance(std::size_t idx) const {
    return m_stripCenterDistance[idx];
  }

  const Acts::Vector3& getTopStripCenterPosition(std::size_t idx) const {
    return m_topStripCenterPosition[idx];
  }

  void setTopStripVector(std::size_t idx, const Acts::Vector3& value) {
    m_topStripVector[idx] = value;
  }

  void setBottomStripVector(std::size_t idx, const Acts::Vector3& value) {
    m_bottomStripVector[idx] = value;
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
  std::vector<Acts::Vector3> m_topStripVector{};
  std::vector<Acts::Vector3> m_bottomStripVector{};
  std::vector<Acts::Vector3> m_stripCenterDistance{};
  std::vector<Acts::Vector3> m_topStripCenterPosition{};
};

inline float SpacePointData::quality(std::size_t idx) const {
  return m_quality[idx];
}

inline float SpacePointData::deltaR(std::size_t idx) const {
  return m_deltaR[idx];
}

inline void SpacePointData::setQuality(std::size_t idx, const float value) {
  if (value > m_quality[idx]) {
    m_quality[idx] = value;
  }
}

inline void SpacePointData::setDeltaR(std::size_t idx, const float value) {
  m_deltaR[idx] = value;
}

inline void SpacePointData::resize(std::size_t n, bool resizeDynamic) {
  clear();

  m_quality.resize(n, -std::numeric_limits<float>::infinity());
  m_deltaR.resize(n, 0.);

  if (resizeDynamic) {
    m_topStripVector.resize(n, {0, 0, 0});
    m_bottomStripVector.resize(n, {0, 0, 0});
    m_stripCenterDistance.resize(n, {0, 0, 0});
    m_topStripCenterPosition.resize(n, {0, 0, 0});
  }
}

inline void SpacePointData::clear() {
  // mutable variables
  m_quality.clear();
  m_deltaR.clear();
  // dynamicvariables
  m_topStripVector.clear();
  m_bottomStripVector.clear();
  m_stripCenterDistance.clear();
  m_topStripCenterPosition.clear();
}

}  // namespace Acts
