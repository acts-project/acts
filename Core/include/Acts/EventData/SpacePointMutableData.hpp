// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Acts {

/// @class SpacePointMutableData
/// This class contains mutable data associated to the
/// external space points provided by the customers
/// These variables are used mainly internally by the seeding algorithm, that
/// reads and updates them for seed selection purposes.
/// The quality is also accessed after the seeding for an additional selection
/// round on the candidates
class SpacePointMutableData {
 public:
  /// @brief Default constructor
  SpacePointMutableData() = default;

  /// No copies
  SpacePointMutableData(const SpacePointMutableData& other) = delete;
  SpacePointMutableData& operator=(const SpacePointMutableData& other) = delete;

  /// @brief Move operations
  /// @param other Source object to move from
  SpacePointMutableData(SpacePointMutableData&& other) noexcept = default;
  /// Move assignment operator
  /// @param other Source object to move from
  /// @return Reference to this object after move assignment
  SpacePointMutableData& operator=(SpacePointMutableData&& other) noexcept =
      default;

  /// @brief Destructor
  ~SpacePointMutableData() = default;

  /// @brief Getters
  /// @param idx Index of the space point
  /// @return Quality value at the given index
  float quality(const std::size_t idx) const;
  /// Get deltaR value for space point
  /// @param idx Index of the space point
  /// @return DeltaR value at the given index
  float deltaR(const std::size_t idx) const;

  /// @brief Setters
  /// @param idx Index of the space point
  /// @param value Quality value to set
  void setQuality(const std::size_t idx, const float value);
  /// Set the delta R value for a specific space point.
  /// @param idx Index of the space point
  /// @param value Delta R value to set
  void setDeltaR(const std::size_t idx, const float value);

  /// @brief Resize vectors
  /// @param n New size for the vectors
  void resize(const std::size_t n);

  /// @brief clear vectors
  void clear();

 private:
  /// Variables
  std::vector<float> m_quality{};
  std::vector<float> m_deltaR{};
};

}  // namespace Acts

#include "Acts/EventData/SpacePointMutableData.ipp"
