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
#include <span>

namespace Acts {

class SpacePointContainer2;
template <typename T, bool read_only>
class SpacePointColumnProxy;
template <typename T>
using MutableSpacePointColumnProxy = SpacePointColumnProxy<T, false>;
template <typename T>
using ConstSpacePointColumnProxy = SpacePointColumnProxy<T, true>;

/// A proxy class for accessing individual spacepoints.
template <bool read_only>
class SpacePointProxy2 {
 public:
  /// Indicates whether this spacepoint proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Type alias for spacepoint index type
  using Index = SpacePointIndex2;

  /// Type alias for container type (const if read-only)
  using Container = const_if_t<ReadOnly, SpacePointContainer2>;

  /// Constructs a spacepoint proxy for the given container and index.
  /// @param container The container holding the spacepoint.
  /// @param index The index of the spacepoint in the container.
  SpacePointProxy2(Container &container, Index index) noexcept
      : m_container(&container), m_index(index) {}

  /// Copy construct a spacepoint proxy.
  /// @param other The spacepoint proxy to copy.
  SpacePointProxy2(const SpacePointProxy2 &other) noexcept = default;

  /// Copy construct a mutable spacepoint proxy.
  /// @param other The mutable spacepoint proxy to copy.
  explicit SpacePointProxy2(const SpacePointProxy2<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Returns a const proxy of the spacepoint.
  /// @return A const proxy of the spacepoint.
  SpacePointProxy2<true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_index};
  }

  /// Gets the container holding the spacepoint.
  /// @return A reference to the container holding the spacepoint.
  SpacePointContainer2 &container() const noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the spacepoint.
  /// @return A const reference to the container holding the spacepoint.
  const SpacePointContainer2 &container() const noexcept {
    return *m_container;
  }
  /// Gets the index of the spacepoint in the container.
  /// @return The index of the spacepoint in the container.
  Index index() const noexcept { return m_index; }

  /// Assigns source links to the spacepoint.
  /// @param sourceLinks A span of source links to assign to the spacepoint.
  /// @throws std::out_of_range if the index is out of range.
  /// @throws std::logic_error if no source links column is available.
  /// @throws std::logic_error if source links are already assigned to the spacepoint.
  void assignSourceLinks(std::span<const SourceLink> sourceLinks) const
    requires(!ReadOnly)
  {
    m_container->assignSourceLinks(m_index, sourceLinks);
  }

  /// Mutable access to the source links of the spacepoint.
  /// @return A mutable span of source links associated with the spacepoint.
  std::span<SourceLink> sourceLinks() const noexcept
    requires(!ReadOnly)
  {
    return m_container->sourceLinks(m_index);
  }
  /// Mutable access to the x coordinate of the spacepoint.
  /// @return A mutable reference to the x coordinate of the spacepoint.
  float &x() const noexcept
    requires(!ReadOnly)
  {
    return m_container->x(m_index);
  }
  /// Mutable access to the y coordinate of the spacepoint.
  /// @return A mutable reference to the y coordinate of the spacepoint.
  float &y() const noexcept
    requires(!ReadOnly)
  {
    return m_container->y(m_index);
  }
  /// Mutable access to the z coordinate of the spacepoint.
  /// @return A mutable reference to the z coordinate of the spacepoint.
  float &z() const noexcept
    requires(!ReadOnly)
  {
    return m_container->z(m_index);
  }
  /// Mutable access to the r coordinate of the spacepoint.
  /// @return A mutable reference to the r coordinate of the spacepoint.
  float &r() const noexcept
    requires(!ReadOnly)
  {
    return m_container->r(m_index);
  }
  /// Mutable access to the phi coordinate of the spacepoint.
  /// @return A mutable reference to the phi coordinate of the spacepoint.
  float &phi() const noexcept
    requires(!ReadOnly)
  {
    return m_container->phi(m_index);
  }
  /// Mutable access to the time information of the spacepoint.
  /// @return A mutable reference to the time information of the spacepoint.
  float &time() const noexcept
    requires(!ReadOnly)
  {
    return m_container->time(m_index);
  }
  /// Mutable access to the variance in Z direction of the spacepoint.
  /// @return A mutable reference to the variance in Z direction of the spacepoint.
  float &varianceZ() const noexcept
    requires(!ReadOnly)
  {
    return m_container->varianceZ(m_index);
  }
  /// Mutable access to the variance in R direction of the spacepoint.
  /// @return A mutable reference to the variance in R direction of the spacepoint.
  float &varianceR() const noexcept
    requires(!ReadOnly)
  {
    return m_container->varianceR(m_index);
  }
  /// Mutable access to the top strip vector of the spacepoint.
  /// @return A mutable reference to the top strip vector of the spacepoint.
  std::array<float, 3> &topStripVector() const noexcept
    requires(!ReadOnly)
  {
    return m_container->topStripVector(m_index);
  }
  /// Mutable access to the bottom strip vector of the spacepoint.
  /// @return A mutable reference to the bottom strip vector of the spacepoint.
  std::array<float, 3> &bottomStripVector() const noexcept
    requires(!ReadOnly)
  {
    return m_container->bottomStripVector(m_index);
  }
  /// Mutable access to the strip center distance of the spacepoint.
  /// @return A mutable reference to the strip center distance of the spacepoint.
  std::array<float, 3> &stripCenterDistance() const noexcept
    requires(!ReadOnly)
  {
    return m_container->stripCenterDistance(m_index);
  }
  /// Mutable access to the top strip center of the spacepoint.
  /// @return A mutable reference to the top strip center of the spacepoint.
  std::array<float, 3> &topStripCenter() const noexcept
    requires(!ReadOnly)
  {
    return m_container->topStripCenter(m_index);
  }
  /// Mutable access to the copy from index of the spacepoint.
  /// @return A mutable reference to the copy from index of the spacepoint.
  SpacePointIndex2 &copyFromIndex() const noexcept
    requires(!ReadOnly)
  {
    return m_container->copyFromIndex(m_index);
  }
  /// @brief Get mutable reference to XY coordinates of the spacepoint
  /// @return Mutable reference to array containing [x, y] coordinates
  std::array<float, 2> &xy() const noexcept
    requires(!ReadOnly)
  {
    return m_container->xy(m_index);
  }
  /// @brief Get mutable reference to ZR coordinates of the spacepoint
  /// @return Mutable reference to array containing [z, r] coordinates
  std::array<float, 2> &zr() const noexcept
    requires(!ReadOnly)
  {
    return m_container->zr(m_index);
  }
  /// @brief Get mutable reference to XYZ coordinates of the spacepoint
  /// @return Mutable reference to array containing [x, y, z] coordinates
  std::array<float, 3> &xyz() const noexcept
    requires(!ReadOnly)
  {
    return m_container->xyz(m_index);
  }
  /// @brief Get mutable reference to XYZR coordinates of the spacepoint
  /// @return Mutable reference to array containing [x, y, z, r] coordinates
  std::array<float, 4> &xyzr() const noexcept
    requires(!ReadOnly)
  {
    return m_container->xyzr(m_index);
  }
  /// @brief Get mutable reference to ZR coordinate variances
  /// @return Mutable reference to array containing [var_z, var_r] variances
  std::array<float, 2> &varianceZR() const noexcept
    requires(!ReadOnly)
  {
    return m_container->varianceZR(m_index);
  }

  /// Mutable access to an extra column of data for the spacepoint.
  /// @param column The extra column proxy to access.
  /// @return A mutable reference to the value in the extra column for the spacepoint.
  template <typename T>
  T &extra(MutableSpacePointColumnProxy<T> &column) const noexcept {
    return column[m_index];
  }

  /// Const access to the source links of the spacepoint.
  /// @return A const span of source links associated with the spacepoint.
  std::span<const SourceLink> sourceLinks() const noexcept {
    return m_container->sourceLinks(m_index);
  }
  /// Const access to the x coordinate of the spacepoint.
  /// @return The x coordinate of the spacepoint.
  float x() const noexcept { return m_container->x(m_index); }
  /// Const access to the y coordinate of the spacepoint.
  /// @return The y coordinate of the spacepoint.
  float y() const noexcept { return m_container->y(m_index); }
  /// Const access to the z coordinate of the spacepoint.
  /// @return The z coordinate of the spacepoint.
  float z() const noexcept { return m_container->z(m_index); }
  /// Const access to the r coordinate of the spacepoint.
  /// @return The r coordinate of the spacepoint.
  float r() const noexcept { return m_container->r(m_index); }
  /// Const access to the phi coordinate of the spacepoint.
  /// @return The phi coordinate of the spacepoint.
  float phi() const noexcept { return m_container->phi(m_index); }
  /// Const access to the time information of the spacepoint.
  /// @return An optional containing the time information of the spacepoint, or
  float time() const noexcept { return m_container->time(m_index); }
  /// Const access to the variance in Z direction of the spacepoint.
  /// @return The variance in Z direction of the spacepoint.
  float varianceZ() const noexcept { return m_container->varianceZ(m_index); }
  /// Const access to the variance in R direction of the spacepoint.
  /// @return The variance in R direction of the spacepoint.
  float varianceR() const noexcept { return m_container->varianceR(m_index); }
  /// Const access to the top strip vector of the spacepoint.
  /// @return A const reference to the top strip vector of the spacepoint.
  const std::array<float, 3> &topStripVector() const noexcept {
    return m_container->topStripVector(m_index);
  }
  /// Const access to the bottom strip vector of the spacepoint.
  /// @return A const reference to the bottom strip vector of the spacepoint.
  const std::array<float, 3> &bottomStripVector() const noexcept {
    return m_container->bottomStripVector(m_index);
  }
  /// Const access to the strip center distance of the spacepoint.
  /// @return A const reference to the strip center distance of the spacepoint.
  const std::array<float, 3> &stripCenterDistance() const noexcept {
    return m_container->stripCenterDistance(m_index);
  }
  /// Const access to the top strip center of the spacepoint.
  /// @return A const reference to the top strip center of the spacepoint.
  const std::array<float, 3> &topStripCenter() const noexcept {
    return m_container->topStripCenter(m_index);
  }
  /// Const access to the copy from index of the spacepoint.
  /// @return A const reference to the copy from index of the spacepoint.
  SpacePointIndex2 copyFromIndex() const noexcept {
    return m_container->copyFromIndex(m_index);
  }
  /// @brief Get const reference to XY coordinates of the spacepoint
  /// @return Const reference to array containing [x, y] coordinates
  const std::array<float, 2> &xy() const noexcept {
    return m_container->xy(m_index);
  }
  /// @brief Get const reference to ZR coordinates of the spacepoint
  /// @return Const reference to array containing [z, r] coordinates
  const std::array<float, 2> &zr() const noexcept {
    return m_container->zr(m_index);
  }
  /// @brief Get const reference to XYZ coordinates of the spacepoint
  /// @return Const reference to array containing [x, y, z] coordinates
  const std::array<float, 3> &xyz() const noexcept {
    return m_container->xyz(m_index);
  }
  /// @brief Get const reference to XYZR coordinates of the spacepoint
  /// @return Const reference to array containing [x, y, z, r] coordinates
  const std::array<float, 4> &xyzr() const noexcept {
    return m_container->xyzr(m_index);
  }
  /// @brief Get const reference to ZR coordinate variances
  /// @return Const reference to array containing [var_z, var_r] variances
  const std::array<float, 2> &varianceZR() const noexcept {
    return m_container->varianceZR(m_index);
  }

  /// Const access to an extra column of data for the spacepoint.
  /// @param column The extra column proxy to access.
  /// @return A const reference to the value in the extra column for the spacepoint.
  template <typename T>
  const T &extra(const ConstSpacePointColumnProxy<T> &column) const noexcept {
    return column[m_index];
  }

  /// Returns the resolved index of the spacepoint.
  /// This resolves the index if the spacepoint was copied from another index.
  /// @return The resolved index of the spacepoint.
  SpacePointIndex2 resolvedIndex() const noexcept {
    return m_container->resolvedIndex(m_index);
  }

 private:
  Container *m_container{};
  Index m_index{};
};

}  // namespace Acts
