// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointColumns.hpp"
#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>
#include <span>

namespace Acts {

class SpacePointContainer;
template <typename T, bool read_only>
class SpacePointColumnProxy;
template <typename T>
using MutableSpacePointColumnProxy = SpacePointColumnProxy<T, false>;
template <typename T>
using ConstSpacePointColumnProxy = SpacePointColumnProxy<T, true>;

/// A proxy class for accessing individual space points.
template <bool read_only>
class SpacePointProxy {
 public:
  /// Indicates whether this space point proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Type alias for space point index type
  using Index = SpacePointIndex;

  /// Type alias for container type (const if read-only)
  using Container = const_if_t<ReadOnly, SpacePointContainer>;

  /// Constructs a space point proxy for the given container and index.
  /// @param container The container holding the space point.
  /// @param index The index of the space point in the container.
  SpacePointProxy(Container &container, Index index) noexcept
      : m_container(&container), m_index(index) {}

  /// Copy construct a space point proxy.
  /// @param other The space point proxy to copy.
  SpacePointProxy(const SpacePointProxy &other) noexcept = default;

  /// Copy construct a mutable space point proxy.
  /// @param other The mutable space point proxy to copy.
  explicit SpacePointProxy(const SpacePointProxy<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Copy assign a space point proxy.
  /// @param other The space point proxy to copy.
  /// @return Reference to this space point proxy after assignment.
  SpacePointProxy &operator=(const SpacePointProxy &other) noexcept = default;

  /// Copy assign a mutable space point proxy.
  /// @param other The mutable space point proxy to copy.
  /// @return Reference to this space point proxy after assignment.
  SpacePointProxy &operator=(const SpacePointProxy<false> &other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Move assign a space point proxy.
  /// @param other The space point proxy to move.
  /// @return Reference to this space point proxy after assignment.
  SpacePointProxy &operator=(SpacePointProxy &&other) noexcept = default;

  /// Move assign a mutable space point proxy.
  /// @param other The mutable space point proxy to move.
  /// @return Reference to this space point proxy after assignment.
  SpacePointProxy &operator=(SpacePointProxy<false> &&other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Returns a const proxy of the space point.
  /// @return A const proxy of the space point.
  SpacePointProxy<true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_index};
  }

  /// Gets the container holding the space point.
  /// @return A reference to the container holding the space point.
  SpacePointContainer &container() const noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the space point.
  /// @return A const reference to the container holding the space point.
  const SpacePointContainer &container() const noexcept { return *m_container; }
  /// Gets the index of the space point in the container.
  /// @return The index of the space point in the container.
  Index index() const noexcept { return m_index; }

  /// Assigns source links to the space point.
  /// @param sourceLinks A span of source links to assign to the space point.
  /// @throws std::out_of_range if the index is out of range.
  /// @throws std::logic_error if no source links column is available.
  /// @throws std::logic_error if source links are already assigned to the space point.
  void assignSourceLinks(std::span<const SourceLink> sourceLinks) const
    requires(!ReadOnly)
  {
    if (m_index >= m_container->size()) {
      throw std::out_of_range("Index out of range in SpacePointContainer");
    }
    if (!m_container->m_sourceLinkOffsetColumn.has_value() ||
        !m_container->m_sourceLinkCountColumn.has_value()) {
      throw std::logic_error("No source links column available");
    }
    if (accessImpl(m_container->m_sourceLinkCountColumn) != 0) {
      throw std::logic_error(
          "Source links already assigned to the space point");
    }

    accessImpl(m_container->m_sourceLinkOffsetColumn) =
        static_cast<std::uint32_t>(m_container->m_sourceLinks.size());
    accessImpl(m_container->m_sourceLinkCountColumn) =
        static_cast<std::uint8_t>(sourceLinks.size());
    m_container->m_sourceLinks.insert(m_container->m_sourceLinks.end(),
                                      sourceLinks.begin(), sourceLinks.end());
  }

  /// Mutable access to the `source links` of the space point.
  /// @return A mutable span of `source links` associated with the space point.
  std::span<SourceLink> sourceLinks() const noexcept
    requires(!ReadOnly)
  {
    const std::size_t offset =
        accessImpl(m_container->m_sourceLinkOffsetColumn);
    const std::size_t count = accessImpl(m_container->m_sourceLinkCountColumn);
    return std::span<SourceLink>(m_container->m_sourceLinks.data() + offset,
                                 count);
  }
  /// Mutable access to the `copied from index` of the space point.
  /// @return A mutable reference to the `copied from index` of the space point.
  SpacePointIndex &copiedFromIndex() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_copiedFromIndexColumn);
  }
  /// Mutable access to the `x` coordinate of the space point.
  /// @return A mutable reference to the `x` coordinate of the space point.
  float &x() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_xColumn);
  }
  /// Mutable access to the `y` coordinate of the space point.
  /// @return A mutable reference to the `y` coordinate of the space point.
  float &y() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_yColumn);
  }
  /// Mutable access to the `z` coordinate of the space point.
  /// @return A mutable reference to the `z` coordinate of the space point.
  float &z() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_zColumn);
  }
  /// Mutable access to the `r` coordinate of the space point.
  /// @return A mutable reference to the `r` coordinate of the space point.
  float &r() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_rColumn);
  }
  /// Mutable access to the `phi` coordinate of the space point.
  /// @return A mutable reference to the `phi` coordinate of the space point.
  float &phi() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_phiColumn);
  }
  /// Mutable access to the `time` information of the space point.
  /// @return A mutable reference to the `time` information of the space point.
  float &time() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_timeColumn);
  }
  /// Mutable access to the `variance z` of the space point.
  /// @return A mutable reference to the `variance z` of the space point.
  float &varianceZ() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_varianceZColumn);
  }
  /// Mutable access to the `variance r` of the space point.
  /// @return A mutable reference to the `variance r` of the space point.
  float &varianceR() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_varianceRColumn);
  }
  /// Mutable access to the `variance t` of the space point.
  /// @return A mutable reference to the `variance t` of the space point.
  float &varianceT() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_varianceTColumn);
  }
  /// Mutable access to the `outer strip calibration details` of the space
  /// point.
  /// @return A mutable reference to the `outer strip calibration details` of the space point.
  OuterStripSpacePointCalibrationDetails &outerStripCalibrationDetails()
      const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_outerStripCalibrationDetailsColumn);
  }
  /// Mutable access to the `XY` coordinates of the space point
  /// @return A mutable reference to array containing `[x, y]` coordinates
  std::array<float, 2> &xy() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_xyColumn);
  }
  /// Mutable access to the `ZR` coordinates of the space point
  /// @return A mutable reference to array containing `[z, r]` coordinates
  std::array<float, 2> &zr() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_zrColumn);
  }
  /// Mutable access to the `XYZ` coordinates of the space point
  /// @return A mutable reference to array containing `[x, y, z]` coordinates
  std::array<float, 3> &xyz() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_xyzColumn);
  }
  /// Mutable access to the `XYZR` coordinates of the space point
  /// @return A mutable reference to array containing `[x, y, z, r]` coordinates
  std::array<float, 4> &xyzr() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_xyzrColumn);
  }
  /// Mutable access to the `ZR` coordinate variances of the space point
  /// @return A mutable reference to array containing `[var_z, var_r]` variances
  std::array<float, 2> &varianceZR() const noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_varianceZRColumn);
  }

  /// Mutable access to an extra column of data for the space point.
  /// @param column The extra column proxy to access.
  /// @return A mutable reference to the value in the extra column for the space point.
  template <typename T>
  T &extra(MutableSpacePointColumnProxy<T> &column) const noexcept {
    return column[m_index];
  }

  /// Const access to the `source links` of the space point.
  /// @return A const span of `source links` associated with the space point.
  std::span<const SourceLink> sourceLinks() const noexcept {
    const std::size_t offset =
        accessImpl(m_container->m_sourceLinkOffsetColumn);
    const std::size_t count = accessImpl(m_container->m_sourceLinkCountColumn);
    return std::span<const SourceLink>(
        m_container->m_sourceLinks.data() + offset, count);
  }
  /// Const access to the `copied from index` of the space point.
  /// @return A const reference to the `copied from index` of the space point.
  SpacePointIndex copiedFromIndex() const noexcept {
    return accessImpl(m_container->m_copiedFromIndexColumn);
  }
  /// Const access to the `x` coordinate of the space point.
  /// @return The `x` coordinate of the space point.
  float x() const noexcept { return accessImpl(m_container->m_xColumn); }
  /// Const access to the `y` coordinate of the space point.
  /// @return The `y` coordinate of the space point.
  float y() const noexcept { return accessImpl(m_container->m_yColumn); }
  /// Const access to the `z` coordinate of the space point.
  /// @return The `z` coordinate of the space point.
  float z() const noexcept { return accessImpl(m_container->m_zColumn); }
  /// Const access to the `r` coordinate of the space point.
  /// @return The `r` coordinate of the space point.
  float r() const noexcept { return accessImpl(m_container->m_rColumn); }
  /// Const access to the `phi` coordinate of the space point.
  /// @return The `phi` coordinate of the space point.
  float phi() const noexcept { return accessImpl(m_container->m_phiColumn); }
  /// Const access to the `time` information of the space point.
  /// @return An optional containing the `time` information of the space point, or
  float time() const noexcept { return accessImpl(m_container->m_timeColumn); }
  /// Const access to the `variance z` of the space point.
  /// @return Value of the `variance z` of the space point.
  float varianceZ() const noexcept {
    return accessImpl(m_container->m_varianceZColumn);
  }
  /// Const access to the `variance r` of the space point.
  /// @return Value of the `variance r` of the space point.
  float varianceR() const noexcept {
    return accessImpl(m_container->m_varianceRColumn);
  }
  /// Const access to the `variance t` of the space point.
  /// @return Value of the `variance t` of the space point.
  float varianceT() const noexcept {
    return accessImpl(m_container->m_varianceTColumn);
  }
  /// Const access to the `outer strip calibration details` of the space point.
  /// @return A const reference to the `outer strip calibration details` of the space point.
  const OuterStripSpacePointCalibrationDetails &outerStripCalibrationDetails()
      const noexcept {
    return accessImpl(m_container->m_outerStripCalibrationDetailsColumn);
  }
  /// Const access to the `XY` coordinates of the space point
  /// @return A const reference to array containing `[x, y]` coordinates
  const std::array<float, 2> &xy() const noexcept {
    return accessImpl(m_container->m_xyColumn);
  }
  /// Const access to the `ZR` coordinates of the space point
  /// @return A const reference to array containing `[z, r]` coordinates
  const std::array<float, 2> &zr() const noexcept {
    return accessImpl(m_container->m_zrColumn);
  }
  /// Const access to the `XYZ` coordinates of the space point
  /// @return A const reference to array containing `[x, y, z]` coordinates
  const std::array<float, 3> &xyz() const noexcept {
    return accessImpl(m_container->m_xyzColumn);
  }
  /// Const access to the `XYZR` coordinates of the space point
  /// @return A const reference to array containing `[x, y, z, r]` coordinates
  const std::array<float, 4> &xyzr() const noexcept {
    return accessImpl(m_container->m_xyzrColumn);
  }
  /// Const access to the `ZR` coordinate variances
  /// @return A const reference to array containing `[var_z, var_r]` variances
  const std::array<float, 2> &varianceZR() const noexcept {
    return accessImpl(m_container->m_varianceZRColumn);
  }

  /// Const access to an extra column of data for the space point.
  /// @param column The extra column proxy to access.
  /// @return A const reference to the value in the extra column for the space point.
  template <typename T>
  const T &extra(const ConstSpacePointColumnProxy<T> &column) const noexcept {
    assert(m_index < column.size() && "Index out of bounds");
    return column[m_index];
  }

  /// Copies the specified columns from another space point to this space point.
  /// @param other The space point proxy to copy from.
  /// @param columnsToCopy The columns to copy from the other space point.
  template <bool other_read_only>
  void copyFrom(const SpacePointProxy<other_read_only> &other,
                SpacePointColumns columnsToCopy) const
    requires(!ReadOnly)
  {
    m_container->copyFrom(m_index, other.container(), other.index(),
                          columnsToCopy);
  }

 private:
  Container *m_container{nullptr};
  Index m_index{0};

  template <typename OptColumn>
  auto &accessImpl(OptColumn &&column) const {
    assert(column.has_value() && "Column does not exist");
    assert(m_index < column->size() && "Index out of bounds");
    return column->proxy(*m_container)[m_index];
  }
};

}  // namespace Acts
