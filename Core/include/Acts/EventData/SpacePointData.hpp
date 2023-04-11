// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
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
  bool hasDynamicVariable() const;

  std::any component(Acts::HashedString key, std::size_t n) const;
  
  const float& getTopHalfStripLength(std::size_t idx) const;
  const float& getBottomHalfStripLength(std::size_t idx) const;
  const Acts::Vector3& getTopStripDirection(std::size_t idx) const;
  const Acts::Vector3& getBottomStripDirection(std::size_t idx) const;
  const Acts::Vector3& getStripCenterDistance(std::size_t idx) const;
  const Acts::Vector3& getTopStripCenterPosition(std::size_t idx) const;

  void setTopHalfStripLength(std::size_t idx, const float& value);
  void setBottomHalfStripLength(std::size_t idx, const float& value);
  void setTopStripDirection(std::size_t idx, const Acts::Vector3& value);
  void setBottomStripDirection(std::size_t idx, const Acts::Vector3& value);
  void setStripCenterDistance(std::size_t idx, const Acts::Vector3& value);
  void setTopStripCenterPosition(std::size_t idx, const Acts::Vector3& value);

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

}  // namespace Acts

#include "Acts/EventData/SpacePointData.ipp"
