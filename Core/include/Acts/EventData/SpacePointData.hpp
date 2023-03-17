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

  /// @brief Reserve memory
  void reserve(std::size_t n);

  /// @brief Resize vectors
  void resize(std::size_t n);

  /// @brief clear vectors
  void clear();

 private:
  std::vector<float> m_quality{};
  std::vector<float> m_deltaR{};
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
