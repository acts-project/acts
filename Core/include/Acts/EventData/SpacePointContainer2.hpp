// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointColumnProxy2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/detail/SpacePointContainer2Column.hpp"
#include "Acts/Utilities/EnumBitwiseOperators.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "Acts/Utilities/detail/ContainerRange.hpp"
#include "Acts/Utilities/detail/ContainerSubset.hpp"

#include <cassert>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Acts {

/// Sentinel value for spacepoints without timing information
static constexpr float NoTime = std::numeric_limits<float>::quiet_NaN();

class SpacePointContainer2;
template <bool read_only>
class SpacePointProxy2;
/// Mutable proxy to a spacepoint allowing modification
using MutableSpacePointProxy2 = SpacePointProxy2<false>;
/// Const proxy to a spacepoint for read-only access
using ConstSpacePointProxy2 = SpacePointProxy2<true>;

/// Enumeration of available columns for spacepoint data storage
enum class SpacePointColumns : std::uint32_t {
  None = 0,  ///< No columns

  SourceLinks = 1 << 0,           ///< Source link information
  X = 1 << 1,                     ///< X coordinate
  Y = 1 << 2,                     ///< Y coordinate
  Z = 1 << 3,                     ///< Z coordinate
  R = 1 << 4,                     ///< Radial coordinate
  Phi = 1 << 5,                   ///< Azimuthal angle
  Time = 1 << 6,                  ///< Time information
  VarianceZ = 1 << 7,             ///< Variance in Z direction
  VarianceR = 1 << 8,             ///< Variance in radial direction
  TopStripVector = 1 << 9,        ///< Vector for the top strip
  BottomStripVector = 1 << 10,    ///< Vector for the bottom strip
  StripCenterDistance = 1 << 11,  ///< Distance to the strip center
  TopStripCenter = 1 << 12,       ///< Center of the top strip
  CopyFromIndex = 1 << 13,        ///< Copy from index

  // packed columns for performance reasons
  XY = 1 << 14,          ///< X and Y coordinates
  ZR = 1 << 15,          ///< Z and R coordinates
  XYZ = 1 << 16,         ///< X, Y, and Z coordinates
  XYZR = 1 << 17,        ///< X, Y, Z, and R coordinates
  VarianceZR = 1 << 18,  ///< Variance in Z and R directions

  /// All strip-related columns
  Strip =
      TopStripVector | BottomStripVector | StripCenterDistance | TopStripCenter,
};

/// Enable bitwise operators for SpacePointColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SpacePointColumns);

/// A container for spacepoints, which can hold additional columns of data
/// and allows for efficient access to spacepoints and their associated source
/// links. Individual spacepoints are addressed via index. A proxy object
/// simplifies the handling.
class SpacePointContainer2 {
 public:
  /// Type alias for spacepoint index in container
  using Index = SpacePointIndex2;
  /// Type alias for range of spacepoint indices
  using IndexRange = SpacePointIndexRange2;
  /// Type alias for subset of spacepoint indices
  using IndexSubset = SpacePointIndexSubset2;
  /// Type alias for mutable spacepoint proxy
  using MutableProxy = MutableSpacePointProxy2;
  /// Type alias for const spacepoint proxy
  using ConstProxy = ConstSpacePointProxy2;

  /// Constructs and empty spacepoint container.
  /// @param columns The columns to create in the container.
  explicit SpacePointContainer2(
      SpacePointColumns columns = SpacePointColumns::None) noexcept;

  /// Constructs a copy of the given spacepoint container.
  /// @param other The spacepoint container to copy.
  SpacePointContainer2(const SpacePointContainer2 &other) noexcept;

  /// Move constructs a spacepoint container.
  /// @param other The spacepoint container to move.
  SpacePointContainer2(SpacePointContainer2 &&other) noexcept;

  /// Detructs the spacepoint container.
  ~SpacePointContainer2() noexcept = default;

  /// Assignment operator for copying a spacepoint container.
  /// @param other The spacepoint container to copy.
  /// @return A reference to this spacepoint container.
  SpacePointContainer2 &operator=(const SpacePointContainer2 &other) noexcept;

  /// Move assignment operator for a spacepoint container.
  /// @param other The spacepoint container to move.
  /// @return A reference to this spacepoint container.
  SpacePointContainer2 &operator=(SpacePointContainer2 &&other) noexcept;

  /// Returns the number of spacepoints in the container.
  /// @return The number of spacepoints in the container.
  std::uint32_t size() const noexcept { return m_size; }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of spacepoints.
  /// This will reserve space for the source links and other columns as well.
  /// @param size The number of spacepoints to reserve space for.
  /// @param averageSourceLinks The average number of source links per spacepoint.
  void reserve(std::uint32_t size, float averageSourceLinks = 1) noexcept;

  /// Clears the container, removing all spacepoints and columns.
  void clear() noexcept;

  /// Creates a new spacepoint at the end of the container.
  /// @return A mutable proxy to the newly created spacepoint.
  MutableProxy createSpacePoint() noexcept;

  /// Creates additional columns. This will create the columns if they do not
  /// already exist.
  /// @param columns The columns to create.
  void createColumns(SpacePointColumns columns) noexcept;

  /// Drops the specified columns from the container.
  /// This will only drop columns if they exist.
  /// @param columns The columns to drop.
  void dropColumns(SpacePointColumns columns) noexcept;

  /// Checks if the container has the given Columns.
  /// @param columns The Columns to check for.
  /// @return True if the container has all the specified Columns, false
  ///         otherwise.
  bool hasColumns(SpacePointColumns columns) const noexcept {
    return (m_knownColumns & columns) == columns;
  }

  /// Creates a new column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created column.
  /// @throws std::runtime_error if a column with the same name already exists.
  /// @throws std::runtime_error if the column name is reserved.
  template <typename T>
  MutableSpacePointColumnProxy<T> createColumn(const std::string &name) {
    return createColumnImpl<ColumnHolder<T>>(name);
  }

  /// Drops the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @throws std::runtime_error if the column does not exist.
  /// @throws std::runtime_error if the column name is reserved.
  void dropColumn(const std::string &name);

  /// Checks if an Column with the given name exists.
  /// @param name The name of the column.
  /// @return True if the column exists, false otherwise.
  bool hasColumn(const std::string &name) const noexcept {
    return m_namedColumns.contains(name);
  }

  /// Returns a mutable reference to the Column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the Column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  ConstSpacePointColumnProxy<T> column(const std::string &name) const {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a mutable proxy to the x coordinate column.
  /// @return A mutable proxy to the x coordinate column.
  MutableSpacePointColumnProxy<float> xColumn() noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    return m_xColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the y coordinate column.
  /// @return A mutable proxy to the y coordinate column.
  MutableSpacePointColumnProxy<float> yColumn() noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    return m_yColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the z coordinate column.
  /// @return A mutable proxy to the z coordinate column.
  MutableSpacePointColumnProxy<float> zColumn() noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    return m_zColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the r coordinate column.
  /// @return A mutable proxy to the r coordinate column.
  MutableSpacePointColumnProxy<float> rColumn() noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    return m_rColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the phi coordinate column.
  /// @return A mutable proxy to the phi coordinate column.
  MutableSpacePointColumnProxy<float> phiColumn() noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    return m_phiColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the time column.
  /// @return A mutable proxy to the time column.
  MutableSpacePointColumnProxy<float> timeColumn() noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    return m_timeColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the variance in Z direction column.
  /// @return A mutable proxy to the variance in Z direction column.
  MutableSpacePointColumnProxy<float> varianceZColumn() noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    return m_varianceZColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the variance in R direction column.
  /// @return A mutable proxy to the variance in R direction column.
  MutableSpacePointColumnProxy<float> varianceRColumn() noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    return m_varianceRColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the top strip vector column.
  /// @return A mutable proxy to the top strip vector column.
  MutableSpacePointColumnProxy<std::array<float, 3>>
  topStripVectorColumn() noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    return m_topStripVectorColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the bottom strip vector column.
  /// @return A mutable proxy to the bottom strip vector column.
  MutableSpacePointColumnProxy<std::array<float, 3>>
  bottomStripVectorColumn() noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    return m_bottomStripVectorColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the strip center distance column.
  /// @return A mutable proxy to the strip center distance column.
  MutableSpacePointColumnProxy<std::array<float, 3>>
  stripCenterDistanceColumn() noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    return m_stripCenterDistanceColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the top strip center column.
  /// @return A mutable proxy to the top strip center column.
  MutableSpacePointColumnProxy<std::array<float, 3>>
  topStripCenterColumn() noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    return m_topStripCenterColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the copy from index column.
  /// @return A mutable proxy to the copy from index column.
  MutableSpacePointColumnProxy<SpacePointIndex2>
  copyFromIndexColumn() noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    return m_copyFromIndexColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the xy coordinates column.
  /// @return A mutable proxy to the xy coordinates column.
  MutableSpacePointColumnProxy<std::array<float, 2>> xyColumn() noexcept {
    assert(m_xyColumn.has_value() && "Column 'xy' does not exist");
    return m_xyColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the zr coordinates column.
  /// @return A mutable proxy to the zr coordinates column.
  MutableSpacePointColumnProxy<std::array<float, 2>> zrColumn() noexcept {
    assert(m_zrColumn.has_value() && "Column 'zr' does not exist");
    return m_zrColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the xyz coordinates column.
  /// @return A mutable proxy to the xyz coordinates column.
  MutableSpacePointColumnProxy<std::array<float, 3>> xyzColumn() noexcept {
    assert(m_xyzColumn.has_value() && "Column 'xyz' does not exist");
    return m_xyzColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the xyzr coordinates column.
  /// @return A mutable proxy to the xyzr coordinates column.
  MutableSpacePointColumnProxy<std::array<float, 4>> xyzrColumn() noexcept {
    assert(m_xyzrColumn.has_value() && "Column 'xyzr' does not exist");
    return m_xyzrColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the variance zr column.
  /// @return A mutable proxy to the variance zr column.
  MutableSpacePointColumnProxy<std::array<float, 2>>
  varianceZRColumn() noexcept {
    assert(m_varianceZRColumn.has_value() &&
           "Column 'varianceZR' does not exist");
    return m_varianceZRColumn->proxy(*this);
  }

  /// Returns a const proxy to the x coordinate column.
  /// @return A const proxy to the x coordinate column.
  ConstSpacePointColumnProxy<float> xColumn() const noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    return m_xColumn->proxy(*this);
  }
  /// Returns a const proxy to the y coordinate column.
  /// @return A const proxy to the y coordinate column.
  ConstSpacePointColumnProxy<float> yColumn() const noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    return m_yColumn->proxy(*this);
  }
  /// Returns a const proxy to the z coordinate column.
  /// @return A const proxy to the z coordinate column.
  ConstSpacePointColumnProxy<float> zColumn() const noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    return m_zColumn->proxy(*this);
  }
  /// Returns a const proxy to the r coordinate column.
  /// @return A const proxy to the r coordinate column.
  ConstSpacePointColumnProxy<float> rColumn() const noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    return m_rColumn->proxy(*this);
  }
  /// Returns a const proxy to the phi coordinate column.
  /// @return A const proxy to the phi coordinate column.
  ConstSpacePointColumnProxy<float> phiColumn() const noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    return m_phiColumn->proxy(*this);
  }
  /// Returns a const proxy to the time column.
  /// @return A const proxy to the time column.
  ConstSpacePointColumnProxy<float> timeColumn() const noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    return m_timeColumn->proxy(*this);
  }
  /// Returns a const proxy to the variance in Z direction column.
  /// @return A const proxy to the variance in Z direction column.
  ConstSpacePointColumnProxy<float> varianceZColumn() const noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    return m_varianceZColumn->proxy(*this);
  }
  /// Returns a const proxy to the variance in R direction column.
  /// @return A const proxy to the variance in R direction column.
  ConstSpacePointColumnProxy<float> varianceRColumn() const noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    return m_varianceRColumn->proxy(*this);
  }
  /// Returns a const proxy to the top strip vector column.
  /// @return A const proxy to the top strip vector column.
  ConstSpacePointColumnProxy<std::array<float, 3>> topStripVectorColumn()
      const noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    return m_topStripVectorColumn->proxy(*this);
  }
  /// Returns a const proxy to the bottom strip vector column.
  /// @return A const proxy to the bottom strip vector column.
  ConstSpacePointColumnProxy<std::array<float, 3>> bottomStripVectorColumn()
      const noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    return m_bottomStripVectorColumn->proxy(*this);
  }
  /// Returns a const proxy to the strip center distance column.
  /// @return A const proxy to the strip center distance column.
  ConstSpacePointColumnProxy<std::array<float, 3>> stripCenterDistanceColumn()
      const noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    return m_stripCenterDistanceColumn->proxy(*this);
  }
  /// Returns a const proxy to the top strip center column.
  /// @return A const proxy to the top strip center column.
  ConstSpacePointColumnProxy<std::array<float, 3>> topStripCenterColumn()
      const noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    return m_topStripCenterColumn->proxy(*this);
  }
  /// Returns a const proxy to the copy from index column.
  /// @return A const proxy to the copy from index column.
  ConstSpacePointColumnProxy<SpacePointIndex2> copyFromIndexColumn()
      const noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    return m_copyFromIndexColumn->proxy(*this);
  }
  /// Returns a const proxy to the xy coordinates column.
  /// @return A const proxy to the xy coordinates column.
  ConstSpacePointColumnProxy<std::array<float, 2>> xyColumn() const noexcept {
    assert(m_xyColumn.has_value() && "Column 'xy' does not exist");
    return m_xyColumn->proxy(*this);
  }
  /// Returns a const proxy to the zr coordinates column.
  /// @return A const proxy to the zr coordinates column.
  ConstSpacePointColumnProxy<std::array<float, 2>> zrColumn() const noexcept {
    assert(m_zrColumn.has_value() && "Column 'zr' does not exist");
    return m_zrColumn->proxy(*this);
  }
  /// Returns a const proxy to the xyz coordinates column.
  /// @return A const proxy to the xyz coordinates column.
  ConstSpacePointColumnProxy<std::array<float, 3>> xyzColumn() const noexcept {
    assert(m_xyzColumn.has_value() && "Column 'xyz' does not exist");
    return m_xyzColumn->proxy(*this);
  }
  /// Returns a const proxy to the xyzr coordinates column.
  /// @return A const proxy to the xyzr coordinates column.
  ConstSpacePointColumnProxy<std::array<float, 4>> xyzrColumn() const noexcept {
    assert(m_xyzrColumn.has_value() && "Column 'xyzr' does not exist");
    return m_xyzrColumn->proxy(*this);
  }
  /// Returns a const proxy to the variance zr column.
  /// @return A const proxy to the variance zr column.
  ConstSpacePointColumnProxy<std::array<float, 2>> varianceZRColumn()
      const noexcept {
    assert(m_varianceZRColumn.has_value() &&
           "Column 'varianceZR' does not exist");
    return m_varianceZRColumn->proxy(*this);
  }

  /// Returns a mutable proxy to the spacepoint at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the spacepoint to access.
  /// @return A mutable proxy to the spacepoint at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the spacepoint at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the spacepoint to access.
  /// @return A const proxy to the spacepoint at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the spacepoint at the given index.
  /// @param index The index of the spacepoint to access.
  /// @return A mutable proxy to the spacepoint at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the spacepoint at the given index.
  /// @param index The index of the spacepoint to access.
  /// @return A const proxy to the spacepoint at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Assigns source links to the spacepoint at the given index.
  /// @param index The index of the spacepoint to assign source links to.
  /// @param sourceLinks A span of source links to assign to the spacepoint.
  /// @throws std::out_of_range if the index is out of range.
  /// @throws std::logic_error if no source links column is available.
  /// @throws std::logic_error if source links are already assigned to the spacepoint.
  void assignSourceLinks(Index index, std::span<const SourceLink> sourceLinks);

  /// Mutable access to the source links at the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the source link at the given index.
  std::span<SourceLink> sourceLinks(Index index) {
    assert(m_sourceLinkOffsetColumn.has_value() &&
           m_sourceLinkCountColumn.has_value() &&
           "Column 'sourceLinks' does not exist");
    assert(index < m_sourceLinkOffsetColumn->size() &&
           index < m_sourceLinkCountColumn->size() && "Index out of bounds");
    return std::span<SourceLink>(
        m_sourceLinks.data() + m_sourceLinkOffsetColumn->proxy(*this)[index],
        m_sourceLinkCountColumn->proxy(*this)[index]);
  }
  /// Mutable access to the x coordinate of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the x coordinate of the spacepoint.
  float &x(Index index) noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    assert(index < m_xColumn->size() && "Index out of bounds");
    return m_xColumn->proxy(*this)[index];
  }
  /// Mutable access to the y coordinate of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the y coordinate of the spacepoint.
  float &y(Index index) noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    assert(index < m_yColumn->size() && "Index out of bounds");
    return m_yColumn->proxy(*this)[index];
  }
  /// Mutable access to the z coordinate of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the z coordinate of the spacepoint.
  float &z(Index index) noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    assert(index < m_zColumn->size() && "Index out of bounds");
    return m_zColumn->proxy(*this)[index];
  }
  /// Mutable access to the r coordinate of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the r coordinate of the spacepoint.
  float &r(Index index) noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    assert(index < m_rColumn->size() && "Index out of bounds");
    return m_rColumn->proxy(*this)[index];
  }
  /// Mutable access to the phi coordinate of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the phi coordinate of the spacepoint.
  float &phi(Index index) noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    assert(index < m_phiColumn->size() && "Index out of bounds");
    return m_phiColumn->proxy(*this)[index];
  }
  /// Mutable access to the time information of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the time information of the spacepoint.
  float &time(Index index) noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    assert(index < m_timeColumn->size() && "Index out of bounds");
    return m_timeColumn->proxy(*this)[index];
  }
  /// Mutable access to the variance in Z direction of the spacepoint at
  /// the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the variance in Z direction of the spacepoint.
  float &varianceZ(Index index) noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    assert(index < m_varianceZColumn->size() && "Index out of bounds");
    return m_varianceZColumn->proxy(*this)[index];
  }
  /// Mutable access to the variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the variance in R direction of the spacepoint.
  float &varianceR(Index index) noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    assert(index < m_varianceRColumn->size() && "Index out of bounds");
    return m_varianceRColumn->proxy(*this)[index];
  }
  /// Mutable access to the top strip vector of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the top strip vector of the spacepoint.
  std::array<float, 3> &topStripVector(Index index) noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    assert(index < m_topStripVectorColumn->size() && "Index out of bounds");
    return m_topStripVectorColumn->proxy(*this)[index];
  }
  /// Mutable access to the bottom strip vector of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the bottom strip vector of the spacepoint.
  std::array<float, 3> &bottomStripVector(Index index) noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    assert(index < m_bottomStripVectorColumn->size() && "Index out of bounds");
    return m_bottomStripVectorColumn->proxy(*this)[index];
  }
  /// Mutable access to the strip center distance of the spacepoint at
  /// the given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the strip center distance of the spacepoint.
  std::array<float, 3> &stripCenterDistance(Index index) noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    assert(index < m_stripCenterDistanceColumn->size() &&
           "Index out of bounds");
    return m_stripCenterDistanceColumn->proxy(*this)[index];
  }
  /// Mutable access to the top strip center of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the top strip center of the spacepoint.
  std::array<float, 3> &topStripCenter(Index index) noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    assert(index < m_topStripCenterColumn->size() && "Index out of bounds");
    return m_topStripCenterColumn->proxy(*this)[index];
  }
  /// Mutable access to the copy from index of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the copy from index of the spacepoint.
  SpacePointIndex2 &copyFromIndex(Index index) noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    assert(index < m_copyFromIndexColumn->size() && "Index out of bounds");
    return m_copyFromIndexColumn->proxy(*this)[index];
  }
  /// Mutable access to the xy coordinates of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the xy coordinates of the spacepoint.
  std::array<float, 2> &xy(Index index) noexcept {
    assert(m_xyColumn.has_value() && "Column 'xy' does not exist");
    assert(index < m_xyColumn->size() && "Index out of bounds");
    return m_xyColumn->proxy(*this)[index];
  }
  /// Mutable access to the zr coordinates of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the zr coordinates of the spacepoint.
  std::array<float, 2> &zr(Index index) noexcept {
    assert(m_zrColumn.has_value() && "Column 'zr' does not exist");
    assert(index < m_zrColumn->size() && "Index out of bounds");
    return m_zrColumn->proxy(*this)[index];
  }
  /// Mutable access to the xyz coordinates of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the xyz coordinates of the spacepoint.
  std::array<float, 3> &xyz(Index index) noexcept {
    assert(m_xyzColumn.has_value() && "Column 'xyz' does not exist");
    assert(index < m_xyzColumn->size() && "Index out of bounds");
    return m_xyzColumn->proxy(*this)[index];
  }
  /// Mutable access to the xyzr coordinates of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the xyzr coordinates of the spacepoint.
  std::array<float, 4> &xyzr(Index index) noexcept {
    assert(m_xyzrColumn.has_value() && "Column 'xyzr' does not exist");
    assert(index < m_xyzrColumn->size() && "Index out of bounds");
    return m_xyzrColumn->proxy(*this)[index];
  }
  /// Mutable access to the variance zr of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A mutable reference to the variance zr of the spacepoint.
  std::array<float, 2> &varianceZR(Index index) noexcept {
    assert(m_varianceZRColumn.has_value() &&
           "Column 'varianceZR' does not exist");
    assert(index < m_varianceZRColumn->size() && "Index out of bounds");
    return m_varianceZRColumn->proxy(*this)[index];
  }

  /// Const access to the source links at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const span to the source links at the given index.
  std::span<const SourceLink> sourceLinks(Index index) const noexcept {
    assert(m_sourceLinkOffsetColumn.has_value() &&
           m_sourceLinkCountColumn.has_value() &&
           "Column 'sourceLinks' does not exist");
    assert(index < m_sourceLinkOffsetColumn->size() &&
           index < m_sourceLinkCountColumn->size() && "Index out of bounds");
    return std::span<const SourceLink>(
        m_sourceLinks.data() + m_sourceLinkOffsetColumn->proxy(*this)[index],
        m_sourceLinkCountColumn->proxy(*this)[index]);
  }
  /// Const access to the x coordinate of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the x coordinate of the spacepoint.
  float x(Index index) const noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    assert(index < m_xColumn->size() && "Index out of bounds");
    return m_xColumn->proxy(*this)[index];
  }
  /// Const access to the y coordinate of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the y coordinate of the spacepoint.
  float y(Index index) const noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    assert(index < m_yColumn->size() && "Index out of bounds");
    return m_yColumn->proxy(*this)[index];
  }
  /// Const access to the z coordinate of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the z coordinate of the spacepoint.
  float z(Index index) const noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    assert(index < m_zColumn->size() && "Index out of bounds");
    return m_zColumn->proxy(*this)[index];
  }
  /// Const access to the r coordinate of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the r coordinate of the spacepoint.
  float r(Index index) const noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    assert(index < m_rColumn->size() && "Index out of bounds");
    return m_rColumn->proxy(*this)[index];
  }
  /// Const access to the phi coordinate of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the phi coordinate of the spacepoint.
  float phi(Index index) const noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    assert(index < m_phiColumn->size() && "Index out of bounds");
    return m_phiColumn->proxy(*this)[index];
  }
  /// Const access to the time information of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the time information of the spacepoint.
  float time(Index index) const noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    assert(index < m_timeColumn->size() && "Index out of bounds");
    return m_timeColumn->proxy(*this)[index];
  }
  /// Const access to the variance in Z direction of the spacepoint at
  /// the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the variance in Z direction of the spacepoint.
  float varianceZ(Index index) const noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    assert(index < m_varianceZColumn->size() && "Index out of bounds");
    return m_varianceZColumn->proxy(*this)[index];
  }
  /// Const access to the variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the variance in R direction of the spacepoint.
  float varianceR(Index index) const noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    assert(index < m_varianceRColumn->size() && "Index out of bounds");
    return m_varianceRColumn->proxy(*this)[index];
  }
  /// Const access to the top strip vector of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the top strip vector of the spacepoint.
  const std::array<float, 3> &topStripVector(Index index) const noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    assert(index < m_topStripVectorColumn->size() && "Index out of bounds");
    return m_topStripVectorColumn->proxy(*this)[index];
  }
  /// Const access to the bottom strip vector of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the bottom strip vector of the spacepoint.
  const std::array<float, 3> &bottomStripVector(Index index) const noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    assert(index < m_bottomStripVectorColumn->size() && "Index out of bounds");
    return m_bottomStripVectorColumn->proxy(*this)[index];
  }
  /// Const access to the strip center distance of the spacepoint at
  /// the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the strip center distance of the spacepoint.
  const std::array<float, 3> &stripCenterDistance(Index index) const noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    assert(index < m_stripCenterDistanceColumn->size() &&
           "Index out of bounds");
    return m_stripCenterDistanceColumn->proxy(*this)[index];
  }
  /// Const access to the top strip center of the spacepoint at the
  /// given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the top strip center of the spacepoint.
  const std::array<float, 3> &topStripCenter(Index index) const noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    assert(index < m_topStripCenterColumn->size() && "Index out of bounds");
    return m_topStripCenterColumn->proxy(*this)[index];
  }
  /// Const access to the copy from index of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the copy from index of the spacepoint.
  SpacePointIndex2 copyFromIndex(Index index) const noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    assert(index < m_copyFromIndexColumn->size() && "Index out of bounds");
    return m_copyFromIndexColumn->proxy(*this)[index];
  }
  /// Const access to the xy coordinates of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the xy coordinates of the spacepoint.
  const std::array<float, 2> &xy(Index index) const noexcept {
    assert(m_xyColumn.has_value() && "Column 'xy' does not exist");
    assert(index < m_xyColumn->size() && "Index out of bounds");
    return m_xyColumn->proxy(*this)[index];
  }
  /// Const access to the zr coordinates of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the zr coordinates of the spacepoint.
  const std::array<float, 2> &zr(Index index) const noexcept {
    assert(m_zrColumn.has_value() && "Column 'zr' does not exist");
    assert(index < m_zrColumn->size() && "Index out of bounds");
    return m_zrColumn->proxy(*this)[index];
  }
  /// Const access to the xyz coordinates of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the xyz coordinates of the spacepoint.
  const std::array<float, 3> &xyz(Index index) const noexcept {
    assert(m_xyzColumn.has_value() && "Column 'xyz' does not exist");
    assert(index < m_xyzColumn->size() && "Index out of bounds");
    return m_xyzColumn->proxy(*this)[index];
  }
  /// Const access to the xyzr coordinates of the spacepoint at the given
  /// index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the xyzr coordinates of the spacepoint.
  const std::array<float, 4> &xyzr(Index index) const noexcept {
    assert(m_xyzrColumn.has_value() && "Column 'xyzr' does not exist");
    assert(index < m_xyzrColumn->size() && "Index out of bounds");
    return m_xyzrColumn->proxy(*this)[index];
  }
  /// Const access to the variance zr of the spacepoint at the given index.
  /// @param index The index of the spacepoint.
  /// @return A const reference to the variance zr of the spacepoint.
  const std::array<float, 2> &varianceZR(Index index) const noexcept {
    assert(m_varianceZRColumn.has_value() &&
           "Column 'varianceZR' does not exist");
    assert(index < m_varianceZRColumn->size() && "Index out of bounds");
    return m_varianceZRColumn->proxy(*this)[index];
  }

  /// Resolves the index to the actual index in the container.
  /// If the copyFromIndex column is set, it will return the index from that
  /// column. Otherwise, it will return the index itself.
  /// @param index The index to resolve.
  /// @return The resolved index.
  SpacePointIndex2 resolvedIndex(Index index) const noexcept {
    if (m_copyFromIndexColumn.has_value()) {
      return this->copyFromIndex(index);
    }
    return index;
  }

  /// Type alias for template iterator over spacepoints in container
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      SpacePointContainer2,
      std::conditional_t<read_only, ConstSpacePointProxy2,
                         MutableSpacePointProxy2>,
      Index, read_only>;

  /// Type alias for mutable iterator over spacepoints
  using iterator = Iterator<false>;
  /// Type alias for const iterator over spacepoints
  using const_iterator = Iterator<true>;

  /// @brief Returns mutable iterator to the beginning of the container
  /// @return Mutable iterator pointing to the first spacepoint
  iterator begin() noexcept { return iterator(*this, 0); }
  /// @brief Returns mutable iterator to the end of the container
  /// @return Mutable iterator pointing past the last spacepoint
  iterator end() noexcept { return iterator(*this, size()); }

  /// @brief Returns const iterator to the beginning of the container
  /// @return Const iterator pointing to the first spacepoint
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /// @brief Returns const iterator to the end of the container
  /// @return Const iterator pointing past the last spacepoint
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

  /// Range facade over contiguous index spans.
  template <bool read_only>
  class Range
      : public Acts::detail::ContainerRange<Range<read_only>, Range<true>,
                                            SpacePointContainer2, Index,
                                            read_only> {
   public:
    /// Base class type
    using Base =
        Acts::detail::ContainerRange<Range<read_only>, Range<true>,
                                     SpacePointContainer2, Index, read_only>;

    using Base::Base;

    /// Zip this range with additional columns
    /// @param columns Additional columns to zip with the range
    /// @return Zipped range with additional columns
    template <typename... Ts>
    auto zip(const ConstSpacePointColumnProxy<Ts> &...columns) const noexcept {
      return Base::container().zip(Base::range(), columns...);
    }
  };
  /// Type alias for mutable range of spacepoints
  using MutableRange = Range<false>;
  /// Type alias for const range of spacepoints
  using ConstRange = Range<true>;

  /// Creates a range of spacepoints from the given index range.
  /// @param range The index range to create the range from.
  /// @return A mutable range of spacepoints.
  MutableRange range(const IndexRange &range) noexcept {
    return MutableRange(*this, range);
  }
  /// Creates a range of spacepoints from the given index range.
  /// @param range The index range to create the range from.
  /// @return A const range of spacepoints.
  ConstRange range(const IndexRange &range) const noexcept {
    return ConstRange(*this, range);
  }

  /// Subset facade over arbitrary index sets.
  template <bool read_only>
  class Subset : public Acts::detail::ContainerSubset<
                     Subset<read_only>, Subset<true>, SpacePointContainer2,
                     std::conditional_t<read_only, ConstSpacePointProxy2,
                                        MutableSpacePointProxy2>,
                     SpacePointIndex2, read_only> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<
        Subset<read_only>, Subset<true>, SpacePointContainer2,
        std::conditional_t<read_only, ConstSpacePointProxy2,
                           MutableSpacePointProxy2>,
        SpacePointIndex2, read_only>;

    using Base::Base;

    /// Zip this subset with additional columns
    /// @param columns Additional columns to zip with the subset
    /// @return Zipped subset with additional columns
    template <typename... Ts>
    auto zip(const ConstSpacePointColumnProxy<Ts> &...columns) const noexcept {
      return Base::container().zip(Base::subset(), columns...);
    }
  };
  /// Type alias for mutable subset of spacepoints
  using MutableSubset = Subset<false>;
  /// Type alias for const subset of spacepoints
  using ConstSubset = Subset<true>;

  /// Creates a mutable subset of spacepoints from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A mutable subset of spacepoints.
  MutableSubset subset(const IndexSubset &subset) noexcept {
    return MutableSubset(*this, subset);
  }
  /// Creates a const subset of spacepoints from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A const subset of spacepoints.
  ConstSubset subset(const IndexSubset &subset) const noexcept {
    return ConstSubset(*this, subset);
  }

  /// Creates a zipped mutable range of spacepoint data from the given columns.
  /// @param columns The columns to zip.
  /// @return A zipped mutable range of spacepoint data.
  template <typename... Ts>
  auto zip(const MutableSpacePointColumnProxy<Ts> &...columns) noexcept {
    return Acts::zip(std::ranges::iota_view<Index, Index>(0, size()),
                     columns.data()...);
  }
  /// Creates a zipped const range of spacepoint data from the given columns.
  /// @param columns The columns to zip.
  /// @return A zipped const range of spacepoint data.
  template <typename... Ts>
  auto zip(const ConstSpacePointColumnProxy<Ts> &...columns) const noexcept {
    return Acts::zip(std::ranges::iota_view<Index, Index>(0, size()),
                     columns.data()...);
  }

  /// Creates a zipped mutable range of spacepoint data from the given columns.
  /// @param range The index range to zip.
  /// @param columns The columns to zip.
  /// @return A zipped mutable range of spacepoint data.
  template <typename... Ts>
  auto zip(const IndexRange &range,
           const MutableSpacePointColumnProxy<Ts> &...columns) noexcept {
    return Acts::zip(
        std::ranges::iota_view<Index, Index>(range.first, range.second),
        columns.data().subspan(range.first, range.second - range.first)...);
  }
  /// Creates a zipped const range of spacepoint data from the given columns.
  /// @param range The index range to create the zipped range from.
  /// @param columns The columns to zip.
  /// @return A zipped const range of spacepoint data.
  template <typename... Ts>
  auto zip(const IndexRange &range,
           const ConstSpacePointColumnProxy<Ts> &...columns) const noexcept {
    return Acts::zip(
        std::ranges::iota_view<Index, Index>(range.first, range.second),
        columns.data().subspan(range.first, range.second - range.first)...);
  }

  /// @brief Create a zipped range over subset indices and mutable column data
  /// @tparam Ts Column data types to zip with indices
  /// @param subset Index subset to iterate over
  /// @param columns Mutable column proxies to zip with indices
  /// @return Zipped range for iteration over indices and column data
  template <typename... Ts>
  auto zip(const IndexSubset &subset,
           const MutableSpacePointColumnProxy<Ts> &...columns) noexcept {
    return Acts::zip(subset, columns.subset(subset)...);
  }

  /// @brief Create a zipped range over subset indices and const column data
  /// @tparam Ts Column data types to zip with indices
  /// @param subset Index subset to iterate over
  /// @param columns Const column proxies to zip with indices
  /// @return Const zipped range for iteration over indices and column data
  template <typename... Ts>
  auto zip(const IndexSubset &subset,
           const ConstSpacePointColumnProxy<Ts> &...columns) const noexcept {
    return Acts::zip(subset, columns.subset(subset)...);
  }

 private:
  using ColumnHolderBase = detail::sp::ColumnHolderBase;
  template <typename T>
  using ColumnHolder = detail::sp::ColumnHolder<T>;

  std::uint32_t m_size{0};

  std::unordered_map<
      std::string,
      std::pair<ColumnHolderBase *, std::unique_ptr<ColumnHolderBase>>,
      std::hash<std::string_view>, std::equal_to<>>
      m_namedColumns;
  SpacePointColumns m_knownColumns{SpacePointColumns::None};

  std::vector<SourceLink> m_sourceLinks;

  std::optional<ColumnHolder<SpacePointIndex2>> m_sourceLinkOffsetColumn;
  std::optional<ColumnHolder<std::uint8_t>> m_sourceLinkCountColumn;

  std::optional<ColumnHolder<float>> m_xColumn;
  std::optional<ColumnHolder<float>> m_yColumn;
  std::optional<ColumnHolder<float>> m_zColumn;

  // cylindrical coordinates
  std::optional<ColumnHolder<float>> m_rColumn;
  std::optional<ColumnHolder<float>> m_phiColumn;
  // time information
  std::optional<ColumnHolder<float>> m_timeColumn;
  // covariance information
  std::optional<ColumnHolder<float>> m_varianceZColumn;
  std::optional<ColumnHolder<float>> m_varianceRColumn;
  // strip information
  std::optional<ColumnHolder<std::array<float, 3>>> m_topStripVectorColumn;
  std::optional<ColumnHolder<std::array<float, 3>>> m_bottomStripVectorColumn;
  std::optional<ColumnHolder<std::array<float, 3>>> m_stripCenterDistanceColumn;
  std::optional<ColumnHolder<std::array<float, 3>>> m_topStripCenterColumn;
  // copy information
  std::optional<ColumnHolder<SpacePointIndex2>> m_copyFromIndexColumn;

  std::optional<ColumnHolder<std::array<float, 2>>> m_xyColumn;
  std::optional<ColumnHolder<std::array<float, 2>>> m_zrColumn;
  std::optional<ColumnHolder<std::array<float, 3>>> m_xyzColumn;
  std::optional<ColumnHolder<std::array<float, 4>>> m_xyzrColumn;
  std::optional<ColumnHolder<std::array<float, 2>>> m_varianceZRColumn;

  static auto knownColumnMasks() noexcept {
    using enum SpacePointColumns;
    return std::tuple(SourceLinks, SourceLinks, X, Y, Z, R, Phi, Time,
                      VarianceZ, VarianceR, TopStripVector, BottomStripVector,
                      StripCenterDistance, TopStripCenter, CopyFromIndex, XY,
                      ZR, XYZ, XYZR, VarianceZR);
  }

  static auto knownColumnNames() noexcept {
    return std::tuple("sourceLinkOffset", "sourceLinkCount", "x", "y", "z", "r",
                      "phi", "time", "varianceZ", "varianceR", "topStripVector",
                      "bottomStripVector", "stripCenterDistance",
                      "topStripCenter", "copyFromIndex", "xy", "zr", "xyz",
                      "xyzr", "varianceZR");
  }

  static auto knownColumnDefaults() noexcept {
    return std::tuple(
        SpacePointIndex2{0}, std::uint8_t{0}, float{0}, float{0}, float{0},
        float{0}, float{0}, float{NoTime}, float{0}, float{0},
        std::array<float, 3>{0, 0, 0}, std::array<float, 3>{0, 0, 0},
        std::array<float, 3>{0, 0, 0}, std::array<float, 3>{0, 0, 0},
        SpacePointIndex2{0}, std::array<float, 2>{0, 0},
        std::array<float, 2>{0, 0}, std::array<float, 3>{0, 0, 0},
        std::array<float, 4>{0, 0, 0, 0}, std::array<float, 2>{0, 0});
  }

  auto knownColumns() & noexcept {
    return std::tie(m_sourceLinkOffsetColumn, m_sourceLinkCountColumn,
                    m_xColumn, m_yColumn, m_zColumn, m_rColumn, m_phiColumn,
                    m_timeColumn, m_varianceZColumn, m_varianceRColumn,
                    m_topStripVectorColumn, m_bottomStripVectorColumn,
                    m_stripCenterDistanceColumn, m_topStripCenterColumn,
                    m_copyFromIndexColumn, m_xyColumn, m_zrColumn, m_xyzColumn,
                    m_xyzrColumn, m_varianceZRColumn);
  }
  auto knownColumns() const & noexcept {
    return std::tie(m_sourceLinkOffsetColumn, m_sourceLinkCountColumn,
                    m_xColumn, m_yColumn, m_zColumn, m_rColumn, m_phiColumn,
                    m_timeColumn, m_varianceZColumn, m_varianceRColumn,
                    m_topStripVectorColumn, m_bottomStripVectorColumn,
                    m_stripCenterDistanceColumn, m_topStripCenterColumn,
                    m_copyFromIndexColumn, m_xyColumn, m_zrColumn, m_xyzColumn,
                    m_xyzrColumn, m_varianceZRColumn);
  }
  auto knownColumns() && noexcept {
    return std::tuple(
        std::move(m_sourceLinkOffsetColumn), std::move(m_sourceLinkCountColumn),
        std::move(m_xColumn), std::move(m_yColumn), std::move(m_zColumn),
        std::move(m_rColumn), std::move(m_phiColumn), std::move(m_timeColumn),
        std::move(m_varianceZColumn), std::move(m_varianceRColumn),
        std::move(m_topStripVectorColumn), std::move(m_bottomStripVectorColumn),
        std::move(m_stripCenterDistanceColumn),
        std::move(m_topStripCenterColumn), std::move(m_copyFromIndexColumn),
        std::move(m_xyColumn), std::move(m_zrColumn), std::move(m_xyzColumn),
        std::move(m_xyzrColumn), std::move(m_varianceZRColumn));
  }

  void copyColumns(const SpacePointContainer2 &other);
  void moveColumns(SpacePointContainer2 &other) noexcept;

  static bool reservedColumn(const std::string &name) noexcept;

  template <typename Holder>
  auto createColumnImpl(const std::string &name) {
    if (reservedColumn(name)) {
      throw std::runtime_error("Column name is reserved: " + name);
    }
    if (hasColumn(name)) {
      throw std::runtime_error("Column already exists: " + name);
    }
    auto holder = std::make_unique<Holder>();
    holder->resize(size());
    auto proxy = holder->proxy(*this);
    m_namedColumns.try_emplace(name,
                               std::pair{holder.get(), std::move(holder)});
    return proxy;
  }

  template <typename Holder>
  auto columnImpl(const std::string &name) const {
    auto it = m_namedColumns.find(name);
    if (it == m_namedColumns.end()) {
      throw std::runtime_error("Column not found: " + name);
    }
    auto &holder = dynamic_cast<Holder &>(*it->second.first);
    return holder.proxy();
  }
};

}  // namespace Acts

#include "Acts/EventData/SpacePointContainer2.ipp"
