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
#include "Acts/EventData/SpacePointColumns.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/detail/SpacePointContainer2Column.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "Acts/Utilities/detail/ContainerRange.hpp"
#include "Acts/Utilities/detail/ContainerSubset.hpp"

#include <cassert>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Acts {

/// Sentinel value for space points without timing information
static constexpr float NoTime = std::numeric_limits<float>::quiet_NaN();

class SpacePointContainer2;
template <bool read_only>
class SpacePointProxy2;
/// Mutable proxy to a space point allowing modification
using MutableSpacePointProxy2 = SpacePointProxy2<false>;
/// Const proxy to a space point for read-only access
using ConstSpacePointProxy2 = SpacePointProxy2<true>;

/// A container for space points, which can hold additional columns of data
/// and allows for efficient access to space points and their associated source
/// links. Individual space points are addressed via index. A proxy object
/// simplifies the handling.
class SpacePointContainer2 {
 public:
  /// Type alias for space point index in container
  using Index = SpacePointIndex2;
  /// Type alias for range of space point indices
  using IndexRange = SpacePointIndexRange2;
  /// Type alias for subset of space point indices
  using IndexSubset = SpacePointIndexSubset2;
  /// Type alias for mutable space point proxy
  using MutableProxy = MutableSpacePointProxy2;
  /// Type alias for const space point proxy
  using ConstProxy = ConstSpacePointProxy2;
  /// Type alias for mutable column proxy
  template <typename T>
  using MutableColumnProxy = MutableSpacePointColumnProxy<T>;
  /// Type alias for const column proxy
  template <typename T>
  using ConstColumnProxy = ConstSpacePointColumnProxy<T>;

  /// Constructs and empty space point container.
  /// @param columns The columns to create in the container.
  explicit SpacePointContainer2(
      SpacePointColumns columns = SpacePointColumns::None) noexcept;

  /// Constructs a copy of the given space point container.
  /// @param other The space point container to copy.
  SpacePointContainer2(const SpacePointContainer2 &other) noexcept;

  /// Move constructs a space point container.
  /// @param other The space point container to move.
  SpacePointContainer2(SpacePointContainer2 &&other) noexcept;

  /// Detructs the space point container.
  ~SpacePointContainer2() noexcept = default;

  /// Assignment operator for copying a space point container.
  /// @param other The space point container to copy.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(const SpacePointContainer2 &other) noexcept;

  /// Move assignment operator for a space point container.
  /// @param other The space point container to move.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(SpacePointContainer2 &&other) noexcept;

  /// Returns the number of space points in the container.
  /// @return The number of space points in the container.
  std::uint32_t size() const noexcept { return m_size; }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of space points.
  /// This will reserve space for the source links and other columns as well.
  /// @param size The number of space points to reserve space for.
  /// @param averageSourceLinks The average number of source links per space point.
  void reserve(std::uint32_t size, float averageSourceLinks = 1) noexcept;

  /// Clears the container, removing all space points and columns.
  void clear() noexcept;

  /// Creates a new space point at the end of the container.
  /// @return A mutable proxy to the newly created space point.
  MutableProxy createSpacePoint() noexcept;

  /// Copies the specified columns from another spacepoint to this spacepoint
  /// @param index The index of the spacepoint to copy to in this container.
  /// @param otherContainer The space point container to copy from.
  /// @param otherIndex The index of the spacepoint to copy from in the other container.
  /// @param columnsToCopy The columns to copy from the other spacepoint.
  void copyFrom(Index index, const SpacePointContainer2 &otherContainer,
                Index otherIndex, SpacePointColumns columnsToCopy);

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
  MutableColumnProxy<T> createColumn(const std::string &name) {
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
    return m_allColumns.contains(name);
  }

  /// Returns a mutable reference to the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  MutableColumnProxy<T> column(const std::string &name) {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a const reference to the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A const reference to the column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  ConstColumnProxy<T> column(const std::string &name) const {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a mutable proxy to the x coordinate column.
  /// @return A mutable proxy to the x coordinate column.
  MutableColumnProxy<float> xColumn() noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    return m_xColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the y coordinate column.
  /// @return A mutable proxy to the y coordinate column.
  MutableColumnProxy<float> yColumn() noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    return m_yColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the z coordinate column.
  /// @return A mutable proxy to the z coordinate column.
  MutableColumnProxy<float> zColumn() noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    return m_zColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the r coordinate column.
  /// @return A mutable proxy to the r coordinate column.
  MutableColumnProxy<float> rColumn() noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    return m_rColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the phi coordinate column.
  /// @return A mutable proxy to the phi coordinate column.
  MutableColumnProxy<float> phiColumn() noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    return m_phiColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the time column.
  /// @return A mutable proxy to the time column.
  MutableColumnProxy<float> timeColumn() noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    return m_timeColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the variance in Z direction column.
  /// @return A mutable proxy to the variance in Z direction column.
  MutableColumnProxy<float> varianceZColumn() noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    return m_varianceZColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the variance in R direction column.
  /// @return A mutable proxy to the variance in R direction column.
  MutableColumnProxy<float> varianceRColumn() noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    return m_varianceRColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `top strip vector` column.
  /// @return A mutable proxy to the `top strip vector` column.
  MutableColumnProxy<std::array<float, 3>> topStripVectorColumn() noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    return m_topStripVectorColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `bottom strip vector` column.
  /// @return A mutable proxy to the `bottom strip vector` column.
  MutableColumnProxy<std::array<float, 3>> bottomStripVectorColumn() noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    return m_bottomStripVectorColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `strip center distance` column.
  /// @return A mutable proxy to the `strip center distance` column.
  MutableColumnProxy<std::array<float, 3>>
  stripCenterDistanceColumn() noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    return m_stripCenterDistanceColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `top strip center` column.
  /// @return A mutable proxy to the `top strip center` column.
  MutableColumnProxy<std::array<float, 3>> topStripCenterColumn() noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    return m_topStripCenterColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `copy from index` column.
  /// @return A mutable proxy to the `copy from index` column.
  MutableColumnProxy<SpacePointIndex2> copyFromIndexColumn() noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    return m_copyFromIndexColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `xy` coordinates column.
  /// @return A mutable proxy to the `xy` coordinates column.
  MutableColumnProxy<std::array<float, 2>> xyColumn() noexcept {
    assert(m_xyColumn.has_value() && "Column 'xy' does not exist");
    return m_xyColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `zr` coordinates column.
  /// @return A mutable proxy to the `zr` coordinates column.
  MutableColumnProxy<std::array<float, 2>> zrColumn() noexcept {
    assert(m_zrColumn.has_value() && "Column 'zr' does not exist");
    return m_zrColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `xyz` coordinates column.
  /// @return A mutable proxy to the `xyz` coordinates column.
  MutableColumnProxy<std::array<float, 3>> xyzColumn() noexcept {
    assert(m_xyzColumn.has_value() && "Column 'xyz' does not exist");
    return m_xyzColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `xyzr` coordinates column.
  /// @return A mutable proxy to the `xyzr` coordinates column.
  MutableColumnProxy<std::array<float, 4>> xyzrColumn() noexcept {
    assert(m_xyzrColumn.has_value() && "Column 'xyzr' does not exist");
    return m_xyzrColumn->proxy(*this);
  }
  /// Returns a mutable proxy to the `variance zr` column.
  /// @return A mutable proxy to the `variance zr` column.
  MutableColumnProxy<std::array<float, 2>> varianceZRColumn() noexcept {
    assert(m_varianceZRColumn.has_value() &&
           "Column 'varianceZR' does not exist");
    return m_varianceZRColumn->proxy(*this);
  }

  /// Returns a const proxy to the x coordinate column.
  /// @return A const proxy to the x coordinate column.
  ConstColumnProxy<float> xColumn() const noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    return m_xColumn->proxy(*this);
  }
  /// Returns a const proxy to the y coordinate column.
  /// @return A const proxy to the y coordinate column.
  ConstColumnProxy<float> yColumn() const noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    return m_yColumn->proxy(*this);
  }
  /// Returns a const proxy to the z coordinate column.
  /// @return A const proxy to the z coordinate column.
  ConstColumnProxy<float> zColumn() const noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    return m_zColumn->proxy(*this);
  }
  /// Returns a const proxy to the r coordinate column.
  /// @return A const proxy to the r coordinate column.
  ConstColumnProxy<float> rColumn() const noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    return m_rColumn->proxy(*this);
  }
  /// Returns a const proxy to the phi coordinate column.
  /// @return A const proxy to the phi coordinate column.
  ConstColumnProxy<float> phiColumn() const noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    return m_phiColumn->proxy(*this);
  }
  /// Returns a const proxy to the time column.
  /// @return A const proxy to the time column.
  ConstColumnProxy<float> timeColumn() const noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    return m_timeColumn->proxy(*this);
  }
  /// Returns a const proxy to the variance in Z direction column.
  /// @return A const proxy to the variance in Z direction column.
  ConstColumnProxy<float> varianceZColumn() const noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    return m_varianceZColumn->proxy(*this);
  }
  /// Returns a const proxy to the variance in R direction column.
  /// @return A const proxy to the variance in R direction column.
  ConstColumnProxy<float> varianceRColumn() const noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    return m_varianceRColumn->proxy(*this);
  }
  /// Returns a const proxy to the `top strip vector` column.
  /// @return A const proxy to the `top strip vector` column.
  ConstColumnProxy<std::array<float, 3>> topStripVectorColumn() const noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    return m_topStripVectorColumn->proxy(*this);
  }
  /// Returns a const proxy to the `bottom strip vector` column.
  /// @return A const proxy to the `bottom strip vector` column.
  ConstColumnProxy<std::array<float, 3>> bottomStripVectorColumn()
      const noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    return m_bottomStripVectorColumn->proxy(*this);
  }
  /// Returns a const proxy to the `strip center distance` column.
  /// @return A const proxy to the `strip center distance` column.
  ConstColumnProxy<std::array<float, 3>> stripCenterDistanceColumn()
      const noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    return m_stripCenterDistanceColumn->proxy(*this);
  }
  /// Returns a const proxy to the `top strip center` column.
  /// @return A const proxy to the `top strip center` column.
  ConstColumnProxy<std::array<float, 3>> topStripCenterColumn() const noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    return m_topStripCenterColumn->proxy(*this);
  }
  /// Returns a const proxy to the `copy from index` column.
  /// @return A const proxy to the `copy from index` column.
  ConstColumnProxy<SpacePointIndex2> copyFromIndexColumn() const noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    return m_copyFromIndexColumn->proxy(*this);
  }
  /// Returns a const proxy to the `xy` coordinates column.
  /// @return A const proxy to the `xy` coordinates column.
  ConstColumnProxy<std::array<float, 2>> xyColumn() const noexcept {
    assert(m_xyColumn.has_value() && "Column 'xy' does not exist");
    return m_xyColumn->proxy(*this);
  }
  /// Returns a const proxy to the `zr` coordinates column.
  /// @return A const proxy to the `zr` coordinates column.
  ConstColumnProxy<std::array<float, 2>> zrColumn() const noexcept {
    assert(m_zrColumn.has_value() && "Column 'zr' does not exist");
    return m_zrColumn->proxy(*this);
  }
  /// Returns a const proxy to the `xyz` coordinates column.
  /// @return A const proxy to the `xyz` coordinates column.
  ConstColumnProxy<std::array<float, 3>> xyzColumn() const noexcept {
    assert(m_xyzColumn.has_value() && "Column 'xyz' does not exist");
    return m_xyzColumn->proxy(*this);
  }
  /// Returns a const proxy to the `xyzr` coordinates column.
  /// @return A const proxy to the `xyzr` coordinates column.
  ConstColumnProxy<std::array<float, 4>> xyzrColumn() const noexcept {
    assert(m_xyzrColumn.has_value() && "Column 'xyzr' does not exist");
    return m_xyzrColumn->proxy(*this);
  }
  /// Returns a const proxy to the `variance zr` column.
  /// @return A const proxy to the `variance zr` column.
  ConstColumnProxy<std::array<float, 2>> varianceZRColumn() const noexcept {
    assert(m_varianceZRColumn.has_value() &&
           "Column 'varianceZR' does not exist");
    return m_varianceZRColumn->proxy(*this);
  }

  /// Returns a mutable proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Type alias for template iterator over space points in container
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      SpacePointContainer2,
      std::conditional_t<read_only, ConstProxy, MutableProxy>, Index,
      read_only>;

  /// Type alias for mutable iterator over space points
  using iterator = Iterator<false>;
  /// Type alias for const iterator over space points
  using const_iterator = Iterator<true>;

  /// @brief Returns mutable iterator to the beginning of the container
  /// @return Mutable iterator pointing to the first space point
  iterator begin() noexcept { return iterator(*this, 0); }
  /// @brief Returns mutable iterator to the end of the container
  /// @return Mutable iterator pointing past the last space point
  iterator end() noexcept { return iterator(*this, size()); }

  /// @brief Returns const iterator to the beginning of the container
  /// @return Const iterator pointing to the first space point
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /// @brief Returns const iterator to the end of the container
  /// @return Const iterator pointing past the last space point
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
    auto zip(const ConstColumnProxy<Ts> &...columns) const noexcept {
      return Base::container().zip(Base::range(), columns...);
    }
  };
  /// Type alias for mutable range of space points
  using MutableRange = Range<false>;
  /// Type alias for const range of space points
  using ConstRange = Range<true>;

  /// Creates a range of space points from the given index range.
  /// @param range The index range to create the range from.
  /// @return A mutable range of space points.
  MutableRange range(const IndexRange &range) noexcept {
    return MutableRange(*this, range);
  }
  /// Creates a range of space points from the given index range.
  /// @param range The index range to create the range from.
  /// @return A const range of space points.
  ConstRange range(const IndexRange &range) const noexcept {
    return ConstRange(*this, range);
  }

  /// Subset facade over arbitrary index sets.
  template <bool read_only>
  class Subset : public Acts::detail::ContainerSubset<
                     Subset<read_only>, Subset<true>, SpacePointContainer2,
                     std::conditional_t<read_only, ConstProxy, MutableProxy>,
                     SpacePointIndex2, read_only> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<
        Subset<read_only>, Subset<true>, SpacePointContainer2,
        std::conditional_t<read_only, ConstProxy, MutableProxy>,
        SpacePointIndex2, read_only>;

    using Base::Base;

    /// Zip this subset with additional columns
    /// @param columns Additional columns to zip with the subset
    /// @return Zipped subset with additional columns
    template <typename... Ts>
    auto zip(const ConstColumnProxy<Ts> &...columns) const noexcept {
      return Base::container().zip(Base::subset(), columns...);
    }
  };
  /// Type alias for mutable subset of space points
  using MutableSubset = Subset<false>;
  /// Type alias for const subset of space points
  using ConstSubset = Subset<true>;

  /// Creates a mutable subset of space points from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A mutable subset of space points.
  MutableSubset subset(const IndexSubset &subset) noexcept {
    return MutableSubset(*this, subset);
  }
  /// Creates a const subset of space points from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A const subset of space points.
  ConstSubset subset(const IndexSubset &subset) const noexcept {
    return ConstSubset(*this, subset);
  }

  /// Creates a zipped mutable range of space point data from the given columns.
  /// @param columns The columns to zip.
  /// @return A zipped mutable range of space point data.
  template <typename... Ts>
  auto zip(const MutableColumnProxy<Ts> &...columns) noexcept {
    return Acts::zip(std::ranges::iota_view<Index, Index>(0, size()),
                     columns.data()...);
  }
  /// Creates a zipped const range of space point data from the given columns.
  /// @param columns The columns to zip.
  /// @return A zipped const range of space point data.
  template <typename... Ts>
  auto zip(const ConstColumnProxy<Ts> &...columns) const noexcept {
    return Acts::zip(std::ranges::iota_view<Index, Index>(0, size()),
                     columns.data()...);
  }

  /// Creates a zipped mutable range of space point data from the given columns.
  /// @param range The index range to zip.
  /// @param columns The columns to zip.
  /// @return A zipped mutable range of space point data.
  template <typename... Ts>
  auto zip(const IndexRange &range,
           const MutableColumnProxy<Ts> &...columns) noexcept {
    return Acts::zip(
        std::ranges::iota_view<Index, Index>(range.first, range.second),
        columns.data().subspan(range.first, range.second - range.first)...);
  }
  /// Creates a zipped const range of space point data from the given columns.
  /// @param range The index range to create the zipped range from.
  /// @param columns The columns to zip.
  /// @return A zipped const range of space point data.
  template <typename... Ts>
  auto zip(const IndexRange &range,
           const ConstColumnProxy<Ts> &...columns) const noexcept {
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
           const MutableColumnProxy<Ts> &...columns) noexcept {
    return Acts::zip(subset, columns.subset(subset)...);
  }

  /// @brief Create a zipped range over subset indices and const column data
  /// @tparam Ts Column data types to zip with indices
  /// @param subset Index subset to iterate over
  /// @param columns Const column proxies to zip with indices
  /// @return Const zipped range for iteration over indices and column data
  template <typename... Ts>
  auto zip(const IndexSubset &subset,
           const ConstColumnProxy<Ts> &...columns) const noexcept {
    return Acts::zip(subset, columns.subset(subset)...);
  }

 private:
  template <bool>
  friend class SpacePointProxy2;

  using ColumnHolderBase = detail::sp::ColumnHolderBase;
  template <typename T>
  using ColumnHolder = detail::sp::ColumnHolder<T>;

  std::uint32_t m_size{0};

  std::unordered_map<std::string, ColumnHolderBase *> m_allColumns;
  SpacePointColumns m_knownColumns{SpacePointColumns::None};
  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_dynamicColumns;

  std::vector<SourceLink> m_sourceLinks;

  std::optional<ColumnHolder<std::uint32_t>> m_sourceLinkOffsetColumn;
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
                      StripCenterDistance, TopStripCenter, CopyFromIndex,
                      PackedXY, PackedZR, PackedXYZ, PackedXYZR,
                      PackedVarianceZR);
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
        std::uint32_t{0}, std::uint8_t{0}, float{0}, float{0}, float{0},
        float{0}, float{0}, float{NoTime}, float{0}, float{0},
        std::array<float, 3>{0, 0, 0}, std::array<float, 3>{0, 0, 0},
        std::array<float, 3>{0, 0, 0}, std::array<float, 3>{0, 0, 0},
        std::uint32_t{0}, std::array<float, 2>{0, 0},
        std::array<float, 2>{0, 0}, std::array<float, 3>{0, 0, 0},
        std::array<float, 4>{0, 0, 0, 0}, std::array<float, 2>{0, 0});
  }

  template <typename Self>
  static auto knownColumns(Self &&self) noexcept {
    return std::tie(self.m_sourceLinkOffsetColumn, self.m_sourceLinkCountColumn,
                    self.m_xColumn, self.m_yColumn, self.m_zColumn,
                    self.m_rColumn, self.m_phiColumn, self.m_timeColumn,
                    self.m_varianceZColumn, self.m_varianceRColumn,
                    self.m_topStripVectorColumn, self.m_bottomStripVectorColumn,
                    self.m_stripCenterDistanceColumn,
                    self.m_topStripCenterColumn, self.m_copyFromIndexColumn,
                    self.m_xyColumn, self.m_zrColumn, self.m_xyzColumn,
                    self.m_xyzrColumn, self.m_varianceZRColumn);
  }
  auto knownColumns() & noexcept { return knownColumns(*this); }
  auto knownColumns() const & noexcept { return knownColumns(*this); }
  auto knownColumns() && noexcept { return knownColumns(*this); }

  void copyColumns(const SpacePointContainer2 &other);
  void moveColumns(SpacePointContainer2 &other) noexcept;

  static bool reservedColumn(const std::string &name) noexcept;

  template <typename Holder>
  MutableColumnProxy<typename Holder::Value> createColumnImpl(
      const std::string &name) {
    if (reservedColumn(name)) {
      throw std::runtime_error("Column name is reserved: " + name);
    }
    if (hasColumn(name)) {
      throw std::runtime_error("Column already exists: " + name);
    }
    auto holder = std::make_unique<Holder>();
    holder->resize(size());
    auto proxy = holder->proxy(*this);
    m_allColumns.try_emplace(name, holder.get());
    m_dynamicColumns.try_emplace(name, std::move(holder));
    return proxy;
  }

  template <typename Holder, typename Self>
  static auto columnImpl(Self &&self, const std::string &name) {
    const auto it = self.m_allColumns.find(name);
    if (it == self.m_allColumns.end()) {
      throw std::runtime_error("Column not found: " + name);
    }
    auto &holder = dynamic_cast<Holder &>(*it->second);
    return holder.proxy(self);
  }

  template <typename Holder>
  MutableColumnProxy<typename Holder::Value> columnImpl(
      const std::string &name) {
    return columnImpl<Holder>(*this, name);
  }

  template <typename Holder>
  ConstColumnProxy<typename Holder::Value> columnImpl(
      const std::string &name) const {
    return columnImpl<Holder>(*this, name);
  }
};

}  // namespace Acts

#include "Acts/EventData/SpacePointContainer2.ipp"
