// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>
#include <optional>
#include <span>

#include <Eigen/Core>

namespace Acts::Experimental {

class SpacePointContainer2;

/// A proxy class for accessing individual space points.
template <bool read_only>
class SpacePointProxy2 {
 public:
  /// Indicates whether this space point proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SpacePointIndex2;

  using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

  /// Constructs a space point proxy for the given container and index.
  /// @param container The container holding the space point.
  /// @param index The index of the space point in the container.
  SpacePointProxy2(ContainerType &container, IndexType index) noexcept
      : m_container(&container), m_index(index) {}

  /// Copy construct a space point proxy.
  /// @param other The space point proxy to copy.
  SpacePointProxy2(const SpacePointProxy2 &other) noexcept = default;

  /// Copy construct a mutable space point proxy.
  /// @param other The mutable space point proxy to copy.
  explicit SpacePointProxy2(const SpacePointProxy2<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Gets the container holding the space point.
  /// @return A reference to the container holding the space point.
  SpacePointContainer2 &container() noexcept { return *m_container; }
  /// Gets the container holding the space point.
  /// @return A const reference to the container holding the space point.
  const SpacePointContainer2 &container() const noexcept {
    return *m_container;
  }
  /// Gets the index of the space point in the container.
  /// @return The index of the space point in the container.
  IndexType index() const noexcept { return m_index; }

  /// Mutable access to the source links of the space point.
  /// @return A mutable span of source links associated with the space point.
  std::span<SourceLink> sourceLinks() noexcept
    requires(!ReadOnly)
  {
    return m_container->sourceLinks(m_index);
  }
  /// Mutable access to the x coordinate of the space point.
  /// @return A mutable reference to the x coordinate of the space point.
  float &x() noexcept
    requires(!ReadOnly)
  {
    return m_container->x(m_index);
  }
  /// Mutable access to the y coordinate of the space point.
  /// @return A mutable reference to the y coordinate of the space point.
  float &y() noexcept
    requires(!ReadOnly)
  {
    return m_container->y(m_index);
  }
  /// Mutable access to the z coordinate of the space point.
  /// @return A mutable reference to the z coordinate of the space point.
  float &z() noexcept
    requires(!ReadOnly)
  {
    return m_container->z(m_index);
  }

  /// Mutable access to the extra r coordinate of the space point.
  /// @return A mutable reference to the r coordinate of the space point.
  float &r() noexcept
    requires(!ReadOnly)
  {
    return m_container->r(m_index);
  }
  /// Mutable access to the extra phi coordinate of the space point.
  /// @return A mutable reference to the phi coordinate of the space point.
  float &phi() noexcept
    requires(!ReadOnly)
  {
    return m_container->phi(m_index);
  }
  /// Mutable access to the extra time information of the space point.
  /// @return A mutable reference to the time information of the space point.
  std::optional<float> &time() noexcept
    requires(!ReadOnly)
  {
    return m_container->time(m_index);
  }
  /// Mutable access to the variance in Z direction of the space point.
  /// @return A mutable reference to the variance in Z direction of the space point.
  float &varianceZ() noexcept
    requires(!ReadOnly)
  {
    return m_container->varianceZ(m_index);
  }
  /// Mutable access to the variance in R direction of the space point.
  /// @return A mutable reference to the variance in R direction of the space point.
  float &varianceR() noexcept
    requires(!ReadOnly)
  {
    return m_container->varianceR(m_index);
  }
  /// Mutable access to the extra top strip vector of the space point.
  /// @return A mutable reference to the top strip vector of the space point.
  Eigen::Vector3f &topStripVector() noexcept
    requires(!ReadOnly)
  {
    return m_container->topStripVector(m_index);
  }
  /// Mutable access to the extra bottom strip vector of the space point.
  /// @return A mutable reference to the bottom strip vector of the space point.
  Eigen::Vector3f &bottomStripVector() noexcept
    requires(!ReadOnly)
  {
    return m_container->bottomStripVector(m_index);
  }
  /// Mutable access to the extra strip center distance of the space point.
  /// @return A mutable reference to the strip center distance of the space point.
  Eigen::Vector3f &stripCenterDistance() noexcept
    requires(!ReadOnly)
  {
    return m_container->stripCenterDistance(m_index);
  }
  /// Mutable access to the extra top strip center of the space point.
  /// @return A mutable reference to the top strip center of the space point.
  Eigen::Vector3f &topStripCenter() noexcept
    requires(!ReadOnly)
  {
    return m_container->topStripCenter(m_index);
  }
  /// Mutable access to the copy from index of the space point.
  /// @return A mutable reference to the copy from index of the space point.
  std::size_t &copyFromIndex() noexcept
    requires(!ReadOnly)
  {
    return m_container->copyFromIndex(m_index);
  }

  /// Mutable access to the extra column of data for the space point.
  /// @param column The extra column to access.
  /// @return A mutable reference to the value in the extra column for the space
  ///         point.
  template <typename column_proxy>
  typename column_proxy::ValueType &extra(column_proxy column) noexcept
    requires(!ReadOnly)
  {
    return m_container->extra(column, m_index);
  }

  /// Const access to the x coordinate of the space point.
  /// @return The x coordinate of the space point.
  float x() const noexcept { return m_container->x(m_index); }
  /// Const access to the y coordinate of the space point.
  /// @return The y coordinate of the space point.
  float y() const noexcept { return m_container->y(m_index); }
  /// Const access to the z coordinate of the space point.
  /// @return The z coordinate of the space point.
  float z() const noexcept { return m_container->z(m_index); }

  /// Const access to the r coordinate of the space point.
  /// @return The r coordinate of the space point.
  float r() const noexcept { return m_container->r(m_index); }
  /// Const access to the phi coordinate of the space point.
  /// @return The phi coordinate of the space point.
  float phi() const noexcept { return m_container->phi(m_index); }
  /// Const access to the time information of the space point.
  /// @return An optional containing the time information of the space point, or
  std::optional<float> time() const noexcept {
    return m_container->time(m_index);
  }
  /// Const access to the variance in Z direction of the space point.
  /// @return The variance in Z direction of the space point.
  float varianceZ() const noexcept { return m_container->varianceZ(m_index); }
  /// Const access to the variance in R direction of the space point.
  /// @return The variance in R direction of the space point.
  float varianceR() const noexcept { return m_container->varianceR(m_index); }
  /// Const access to the extra top strip vector of the space point.
  /// @return A const reference to the top strip vector of the space point.
  const Eigen::Vector3f &topStripVector() const noexcept {
    return m_container->topStripVector(m_index);
  }
  /// Const access to the extra bottom strip vector of the space point.
  /// @return A const reference to the bottom strip vector of the space point.
  const Eigen::Vector3f &bottomStripVector() const noexcept {
    return m_container->bottomStripVector(m_index);
  }
  /// Const access to the extra strip center distance of the space point.
  /// @return A const reference to the strip center distance of the space point.
  const Eigen::Vector3f &stripCenterDistance() const noexcept {
    return m_container->stripCenterDistance(m_index);
  }
  /// Const access to the extra top strip center of the space point.
  /// @return A const reference to the top strip center of the space point.
  const Eigen::Vector3f &topStripCenter() const noexcept {
    return m_container->topStripCenter(m_index);
  }
  /// Const access to the copy from index of the space point.
  /// @return A const reference to the copy from index of the space point.
  std::size_t copyFromIndex() const noexcept {
    return m_container->copyFromIndex(m_index);
  }

  /// Const access to the extra column of data for the space point.
  /// @param column The extra column to access.
  /// @return A const reference to the value in the extra column for the space
  ///         point.
  template <typename column_proxy>
  const typename column_proxy::ValueType &extra(
      const column_proxy &column) const noexcept {
    return m_container->extra(column, m_index);
  }

 private:
  ContainerType *m_container{};
  IndexType m_index{};
};

}  // namespace Acts::Experimental
