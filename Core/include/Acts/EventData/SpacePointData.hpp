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
  const float& x(const std::size_t& idx) const;
  const float& y(const std::size_t& idx) const;
  const float& z(const std::size_t& idx) const;
  const float& radius(const std::size_t& idx) const;
  const float& phi(const std::size_t& idx) const;
  const float& varianceZ(const std::size_t& idx) const;
  const float& varianceR(const std::size_t& idx) const;
  
  const float& quality(const std::size_t& idx) const;
  const float& deltaR(const std::size_t& idx) const;

  /// @brief Setters
  void setX(const std::size_t& idx, const float& value);
  void setY(const std::size_t& idx, const float& value);
  void setZ(const std::size_t& idx, const float& value);
  void setRadius(const std::size_t& idx, const float& value);
  void setPhi(const std::size_t& idx, const float& value);
  void setVarianceZ(const std::size_t& idx, const float& value);
  void setVarianceR(const std::size_t& idx, const float& value);
  
  void setQuality(const std::size_t& idx, const float& value);
  void setDeltaR(const std::size_t& idx, const float& value);

  /// @brief Resize vectors
  void resize(const std::size_t& n, bool resizeDynamic = false);

  /// @brief clear vectors
  void clear();

  ///
  bool hasDynamicVariable() const;

  std::any component(Acts::HashedString key, const std::size_t& n) const;
  
  const float& getTopHalfStripLength(const std::size_t& idx) const;
  const float& getBottomHalfStripLength(const std::size_t& idx) const;
  const Acts::Vector3& getTopStripDirection(const std::size_t& idx) const;
  const Acts::Vector3& getBottomStripDirection(const std::size_t& idx) const;
  const Acts::Vector3& getStripCenterDistance(const std::size_t& idx) const;
  const Acts::Vector3& getTopStripCenterPosition(const std::size_t& idx) const;

  void setTopHalfStripLength(const std::size_t& idx, const float& value);
  void setBottomHalfStripLength(const std::size_t& idx, const float& value);
  void setTopStripDirection(const std::size_t& idx, const Acts::Vector3& value);
  void setBottomStripDirection(const std::size_t& idx, const Acts::Vector3& value);
  void setStripCenterDistance(const std::size_t& idx, const Acts::Vector3& value);
  void setTopStripCenterPosition(const std::size_t& idx, const Acts::Vector3& value);

 private:
  /// base variables
  std::vector<float> m_x{};
  std::vector<float> m_y{};
  std::vector<float> m_z{};
  std::vector<float> m_radius{};
  std::vector<float> m_phi{};
  std::vector<float> m_varianceR{};
  std::vector<float> m_varianceZ{};
  
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
