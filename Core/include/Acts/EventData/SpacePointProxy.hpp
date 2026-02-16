// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace Acts {

/// Proxy giving read access to space point data in a container.
template <typename container_t>
class SpacePointProxy {
 public:
  /// Container type
  using ContainerType = container_t;
  /// Value type of the container
  using ValueType = typename ContainerType::ValueType;

  // Never take the ownership of the container
  SpacePointProxy(container_t&& container, std::size_t index) = delete;
  /// Construct from container reference and index
  /// @param container Reference to the space point container
  /// @param index Index of the space point in the container
  SpacePointProxy(const container_t& container, std::size_t index);

  /// Access the external space point
  /// @return Const reference to the external space point
  const ValueType& externalSpacePoint() const;
  /// Get the index of this space point
  /// @return Index of this space point
  std::size_t index() const;

  /// Get the x coordinate
  /// @return x coordinate
  float x() const;
  /// Get the y coordinate
  /// @return y coordinate
  float y() const;
  /// Get the z coordinate
  /// @return z coordinate
  float z() const;
  /// Get the azimuthal angle
  /// @return Azimuthal angle phi
  float phi() const;
  /// Get the radius
  /// @return Radius in the transverse plane
  float radius() const;
  /// Get the variance in r
  /// @return Variance in r
  float varianceR() const;
  /// Get the variance in z
  /// @return Variance in z
  float varianceZ() const;

  /// Get the top strip vector
  /// @return Top strip vector
  const Acts::Vector3& topStripVector() const;
  /// Get the bottom strip vector
  /// @return Bottom strip vector
  const Acts::Vector3& bottomStripVector() const;
  /// Get the strip center distance
  /// @return Strip center distance vector
  const Acts::Vector3& stripCenterDistance() const;
  /// Get the top strip center position
  /// @return Top strip center position vector
  const Acts::Vector3& topStripCenterPosition() const;

 private:
  const container_t& container() const;

  const container_t* m_container{nullptr};
  std::size_t m_index{0ul};
};

}  // namespace Acts

#include "Acts/EventData/SpacePointProxy.ipp"
