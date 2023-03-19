// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <vector>
#include <unordered_map>
#include <iostream>

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
  ~SpacePointData() {
    std::cout << "Num of mutable variables: " << m_quality.size() << "\n";
    std::cout << "Num of dynamic variables: " << m_topHalfStripLength.size() << "\n";
    std::cout << "     * " << static_cast<float>(m_topHalfStripLength.size()) / static_cast<float>(m_quality.size()) * 100. << "\n"; 
  }

  /// @brief Getters
  const float& quality(std::size_t idx) const;
  const float& deltaR(std::size_t idx) const;

  /// @brief Setters
  void setQuality(std::size_t idx, const float& value);
  void setDeltaR(std::size_t idx, const float& value);

  /// @brief Reserve memory
  void reserve(std::size_t n);

  /// @brief Resize vectors
  void resize(std::size_t n);

  /// @brief clear vectors
  void clear();

  ///
  bool hasDynamicVariable(std::size_t idx) const
  { return m_topHalfStripLength.find(idx) != m_topHalfStripLength.end(); }

  
  const float& getTopHalfStripLength(std::size_t idx) const
  { return m_topHalfStripLength.at(idx); }

  const float& getBottomHalfStripLength(std::size_t idx) const
  { return m_bottomHalfStripLength.at(idx); }

  const Acts::Vector3& getTopStripDirection(std::size_t idx) const
  { return m_topStripDirection.at(idx); }

  const Acts::Vector3& getBottomStripDirection(std::size_t idx) const
  { return m_bottomStripDirection.at(idx); }

  const Acts::Vector3& getStripCenterDistance(std::size_t idx) const
  { return m_stripCenterDistance.at(idx); }

  const Acts::Vector3& getTopStripCenterPosition(std::size_t idx) const
  { return m_topStripCenterPosition.at(idx); }


  void setTopHalfStripLength(std::size_t idx, const float& value)
  { m_topHalfStripLength[idx] = value; }
  
  void setBottomHalfStripLength(std::size_t idx, const float& value)
  { m_bottomHalfStripLength[idx] = value; }
  
  void setTopStripDirection(std::size_t idx, const Acts::Vector3& value)
  { m_topStripDirection[idx] = value; }
  
  void setBottomStripDirection(std::size_t idx, const Acts::Vector3& value)
  { m_bottomStripDirection[idx] = value; }
  
  void setStripCenterDistance(std::size_t idx, const Acts::Vector3& value)
  { m_stripCenterDistance[idx] = value; }
  
  void setTopStripCenterPosition(std::size_t idx, const Acts::Vector3& value)
  { m_topStripCenterPosition[idx] = value; }
  
  
 private:
  /// Mutable variables
  std::vector<float> m_quality{};
  std::vector<float> m_deltaR{};

  /// dynamic variables
  std::unordered_map<std::size_t, float> m_topHalfStripLength;
  std::unordered_map<std::size_t, float> m_bottomHalfStripLength;
  std::unordered_map<std::size_t, Acts::Vector3> m_topStripDirection;
  std::unordered_map<std::size_t, Acts::Vector3> m_bottomStripDirection;
  std::unordered_map<std::size_t, Acts::Vector3> m_stripCenterDistance;
  std::unordered_map<std::size_t, Acts::Vector3> m_topStripCenterPosition;
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

inline void SpacePointData::reserve(std::size_t n) {
  m_quality.reserve(n);
  m_deltaR.reserve(n);
}

inline void SpacePointData::resize(std::size_t n) {
  m_quality.resize(n, -std::numeric_limits<float>::infinity());
  m_deltaR.resize(n, 0.);
}

inline void SpacePointData::clear() {
  m_quality.clear();
  m_deltaR.clear();
}

}  // namespace Acts
