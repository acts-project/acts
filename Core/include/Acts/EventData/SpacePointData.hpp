// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <vector>

namespace Acts {

/// @class SpacePointData
/// This class contains auxiliary data associated to the
/// external spacepoints provided by the customers
/// These variables are used internally by the seeding algorithm, that
/// reads them
/// The variables collected here are also dynamic variables only present
/// for strip spacepoints
class SpacePointData {
 public:
  /// @brief Default constructor
  SpacePointData() = default;

  /// No copies
  SpacePointData(const SpacePointData& other) = delete;
  SpacePointData& operator=(const SpacePointData& other) = delete;

  /// @brief Move operations
  /// @param other SpacePointData object to move from
  SpacePointData(SpacePointData&& other) noexcept = default;
  /// Move assignment operator
  /// @param other SpacePointData object to move from
  /// @return Reference to this object
  SpacePointData& operator=(SpacePointData&& other) noexcept = default;

  /// @brief Destructor
  ~SpacePointData() = default;

  /// @param idx Index of the spacepoint
  /// @return X coordinate value
  float x(const std::size_t idx) const;
  /// Get y coordinate of spacepoint
  /// @param idx Index of the spacepoint
  /// @return Y coordinate value
  float y(const std::size_t idx) const;
  /// Get z coordinate of spacepoint
  /// @param idx Index of the spacepoint
  /// @return Z coordinate value
  float z(const std::size_t idx) const;
  /// Get radial distance of spacepoint
  /// @param idx Index of the spacepoint
  /// @return Radial distance value
  float radius(const std::size_t idx) const;
  /// Get azimuthal angle of spacepoint
  /// @param idx Index of the spacepoint
  /// @return Azimuthal angle value
  float phi(const std::size_t idx) const;
  /// Get z coordinate variance of spacepoint
  /// @param idx Index of the spacepoint
  /// @return Z variance value
  float varianceZ(const std::size_t idx) const;
  /// Get radial variance of spacepoint
  /// @param idx Index of the spacepoint
  /// @return Radial variance value
  float varianceR(const std::size_t idx) const;

  /// @param idx Index of the spacepoint to modify
  /// @param value New x coordinate value to set
  void setX(const std::size_t idx, const float value);
  /// Set y coordinate of spacepoint
  /// @param idx Index of the spacepoint to modify
  /// @param value New y coordinate value to set
  void setY(const std::size_t idx, const float value);
  /// Set z coordinate of spacepoint
  /// @param idx Index of the spacepoint to modify
  /// @param value New z coordinate value to set
  void setZ(const std::size_t idx, const float value);
  /// Set radial distance of spacepoint
  /// @param idx Index of the spacepoint to modify
  /// @param value New radial distance value to set
  void setRadius(const std::size_t idx, const float value);
  /// Set azimuthal angle of spacepoint
  /// @param idx Index of the spacepoint to modify
  /// @param value New azimuthal angle value to set
  void setPhi(const std::size_t idx, const float value);
  /// Set z coordinate variance of spacepoint
  /// @param idx Index of the spacepoint to modify
  /// @param value New z variance value to set
  void setVarianceZ(const std::size_t idx, const float value);
  /// Set radial variance of spacepoint
  /// @param idx Index of the spacepoint to modify
  /// @param value New radial variance value to set
  void setVarianceR(const std::size_t idx, const float value);

  /// @brief Resize vectors
  /// @param n New size for the data vectors
  /// @param resizeDynamic Whether to resize dynamic data containers
  void resize(const std::size_t n, bool resizeDynamic = false);

  /// @brief clear vectors
  void clear();

  /// Check if spacepoint data has dynamic variables
  /// @return True if dynamic variables (strip data) are present
  bool hasDynamicVariable() const;

  /// Get top strip vector for strip spacepoints
  /// @param idx Index of the spacepoint
  /// @return Reference to top strip vector
  const Acts::Vector3& topStripVector(const std::size_t idx) const;
  /// Get bottom strip vector for strip spacepoints
  /// @param idx Index of the spacepoint
  /// @return Reference to bottom strip vector
  const Acts::Vector3& bottomStripVector(const std::size_t idx) const;
  /// Get strip center distance vector for strip spacepoints
  /// @param idx Index of the spacepoint
  /// @return Reference to strip center distance vector
  const Acts::Vector3& stripCenterDistance(const std::size_t idx) const;
  /// Get top strip center position for strip spacepoints
  /// @param idx Index of the spacepoint
  /// @return Reference to top strip center position vector
  const Acts::Vector3& topStripCenterPosition(const std::size_t idx) const;

  /// Set top strip vector for strip spacepoints
  /// @param idx Index of the spacepoint to modify
  /// @param value New top strip vector to set
  void setTopStripVector(const std::size_t idx, const Acts::Vector3& value);
  /// Set bottom strip vector for strip spacepoints
  /// @param idx Index of the spacepoint to modify
  /// @param value New bottom strip vector to set
  void setBottomStripVector(const std::size_t idx, const Acts::Vector3& value);
  /// Set strip center distance vector for strip spacepoints
  /// @param idx Index of the spacepoint to modify
  /// @param value New strip center distance vector to set
  void setStripCenterDistance(const std::size_t idx,
                              const Acts::Vector3& value);
  /// Set top strip center position for strip spacepoints
  /// @param idx Index of the spacepoint to modify
  /// @param value New top strip center position vector to set
  void setTopStripCenterPosition(const std::size_t idx,
                                 const Acts::Vector3& value);

 private:
  /// base variables
  std::vector<float> m_x{};
  std::vector<float> m_y{};
  std::vector<float> m_z{};
  std::vector<float> m_radius{};
  std::vector<float> m_phi{};
  std::vector<float> m_varianceR{};
  std::vector<float> m_varianceZ{};

  /// dynamic variables
  std::vector<Acts::Vector3> m_topStripVector{};
  std::vector<Acts::Vector3> m_bottomStripVector{};
  std::vector<Acts::Vector3> m_stripCenterDistance{};
  std::vector<Acts::Vector3> m_topStripCenterPosition{};
};

}  // namespace Acts

#include "Acts/EventData/SpacePointData.ipp"
