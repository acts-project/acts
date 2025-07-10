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
#include "Acts/Utilities/EnumBitwiseOperators.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <cassert>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace Acts::Experimental {

static constexpr float NoTime = std::numeric_limits<float>::quiet_NaN();

template <bool read_only>
class SpacePointProxy2;
using MutableSpacePointProxy2 = SpacePointProxy2<false>;
using ConstSpacePointProxy2 = SpacePointProxy2<true>;

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

  /// Default set of columns
  Default = SourceLinks | X | Y | Z,
  /// All strip-related columns
  Strip =
      TopStripVector | BottomStripVector | StripCenterDistance | TopStripCenter,
};

ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SpacePointColumns);

/// Additional column of data that can be added to the space point container.
/// The column is indexed by the space point index.
template <typename T>
class SpacePointColumnProxy {
 public:
  using Value = T;
  using Container = std::vector<Value>;

  explicit SpacePointColumnProxy(const Container &container)
      : m_container(&container) {}
  SpacePointColumnProxy(const SpacePointColumnProxy &other) = default;
  SpacePointColumnProxy(SpacePointColumnProxy &&other) noexcept = default;
  SpacePointColumnProxy &operator=(const SpacePointColumnProxy &other) =
      default;
  SpacePointColumnProxy &operator=(SpacePointColumnProxy &&other) noexcept =
      default;

  std::span<Value> data() {
    return std::span<Value>(container().data(), container().size());
  }
  std::span<const Value> data() const {
    return std::span<const Value>(container().data(), container().size());
  }

 private:
  const Container *m_container{};

  Container &container() { return const_cast<Container &>(*m_container); }
  const Container &container() const { return *m_container; }

  friend class SpacePointContainer2;
};

class SpacePointColumnHolderBase {
 public:
  virtual ~SpacePointColumnHolderBase() = default;

  virtual std::unique_ptr<SpacePointColumnHolderBase> copy() const = 0;

  virtual std::size_t size() const = 0;
  virtual void reserve(std::size_t size) = 0;
  virtual void resize(std::size_t size) = 0;
  virtual void clear() = 0;
  virtual void emplace_back() = 0;
};

template <typename T>
class SpacePointColumnHolder final : public SpacePointColumnHolderBase {
 public:
  using Value = T;
  using Container = std::vector<Value>;
  using Proxy = SpacePointColumnProxy<Value>;

  SpacePointColumnHolder() = default;
  explicit SpacePointColumnHolder(Value defaultValue)
      : m_default(std::move(defaultValue)) {}

  Proxy proxy() const { return Proxy(m_data); }

  std::unique_ptr<SpacePointColumnHolderBase> copy() const override {
    return std::make_unique<SpacePointColumnHolder<T>>(*this);
  }

  std::size_t size() const override { return m_data.size(); }
  void reserve(std::size_t size) override { m_data.reserve(size); }
  void clear() override { m_data.clear(); }
  void resize(std::size_t size) override { m_data.resize(size, m_default); }
  void emplace_back() override { m_data.emplace_back(m_default); }

 private:
  Value m_default{};
  Container m_data;
};

/// A container for space points, which can hold additional columns of data
/// and allows for efficient access to space points and their associated source
/// links. Individual space points are addressed via index. A proxy object
/// simplifies the handling.
class SpacePointContainer2 {
 public:
  using Index = SpacePointIndex2;
  using IndexRange = SpacePointIndexRange2;
  using IndexSubset = SpacePointSubset2;
  using MutableProxy = MutableSpacePointProxy2;
  using ConstProxy = ConstSpacePointProxy2;

  /// Constructs and empty space point container.
  /// /// @param columns The columns to create in the container.
  explicit SpacePointContainer2(
      SpacePointColumns columns = SpacePointColumns::Default) noexcept;

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
  SpacePointColumnProxy<T> createColumn(const std::string &name) {
    return createColumnImpl<SpacePointColumnHolder<T>>(name);
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
  SpacePointColumnProxy<T> column(const std::string &name) const {
    return columnImpl<SpacePointColumnHolder<T>>(name);
  }

  /// Returns a proxy to the x coordinate column.
  /// @return A proxy to the r coordinate column.
  SpacePointColumnProxy<float> xColumn() const noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    return m_rColumn->proxy();
  }
  /// Returns a proxy to the y coordinate column.
  /// @return A proxy to the y coordinate column.
  SpacePointColumnProxy<float> yColumn() const noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    return m_yColumn->proxy();
  }
  /// Returns a proxy to the z coordinate column.
  /// @return A proxy to the z coordinate column.
  SpacePointColumnProxy<float> zColumn() const noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    return m_zColumn->proxy();
  }
  /// Returns a proxy to the r coordinate column.
  /// @return A proxy to the r coordinate column.
  SpacePointColumnProxy<float> rColumn() const noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    return m_rColumn->proxy();
  }
  /// Returns a proxy to the phi coordinate column.
  /// @return A proxy to the phi coordinate column.
  SpacePointColumnProxy<float> phiColumn() const noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    return m_phiColumn->proxy();
  }
  /// Returns a proxy to the time column.
  /// @return A proxy to the time column.
  SpacePointColumnProxy<float> timeColumn() const noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    return m_timeColumn->proxy();
  }
  /// Returns a proxy to the variance in Z direction column.
  /// @return A proxy to the variance in Z direction column.
  SpacePointColumnProxy<float> varianceZColumn() const noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    return m_varianceZColumn->proxy();
  }
  /// Returns a proxy to the variance in R direction column.
  /// @return A proxy to the variance in R direction column.
  SpacePointColumnProxy<float> varianceRColumn() const noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    return m_varianceRColumn->proxy();
  }
  /// Returns a proxy to the top strip vector column.
  /// @return A proxy to the top strip vector column.
  SpacePointColumnProxy<Eigen::Vector3f> topStripVectorColumn() const noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    return m_topStripVectorColumn->proxy();
  }
  /// Returns a proxy to the bottom strip vector column.
  /// @return A proxy to the bottom strip vector column.
  SpacePointColumnProxy<Eigen::Vector3f> bottomStripVectorColumn()
      const noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    return m_bottomStripVectorColumn->proxy();
  }
  /// Returns a proxy to the strip center distance column.
  /// @return A proxy to the strip center distance column.
  SpacePointColumnProxy<Eigen::Vector3f> stripCenterDistanceColumn()
      const noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    return m_stripCenterDistanceColumn->proxy();
  }
  /// Returns a proxy to the top strip center column.
  /// @return A proxy to the top strip center column.
  SpacePointColumnProxy<Eigen::Vector3f> topStripCenterColumn() const noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    return m_topStripCenterColumn->proxy();
  }
  /// Returns a proxy to the copy from index column.
  /// @return A proxy to the copy from index column.
  SpacePointColumnProxy<std::size_t> copyFromIndexColumn() const noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    return m_copyFromIndexColumn->proxy();
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

  /// Assigns source links to the space point at the given index.
  /// @param index The index of the space point to assign source links to.
  /// @param sourceLinks A span of source links to assign to the space point.
  /// @throws std::out_of_range if the index is out of range.
  /// @throws std::logic_error if no source links column is available.
  /// @throws std::logic_error if source links are already assigned to the space point.
  void assignSourceLinks(Index index, std::span<const SourceLink> sourceLinks);

  /// Mutable access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the source link at the given index.
  std::span<SourceLink> sourceLinks(Index index) {
    assert(m_sourceLinkOffsetColumn.has_value() &&
           m_sourceLinkCountColumn.has_value() &&
           "Column 'sourceLinks' does not exist");
    assert(index < m_sourceLinkOffsetColumn->size() &&
           index < m_sourceLinkCountColumn->size() && "Index out of bounds");
    return std::span<SourceLink>(
        m_sourceLinks.data() + entry(m_sourceLinkOffsetColumn->proxy(), index),
        entry(m_sourceLinkCountColumn->proxy(), index));
  }
  /// Mutable access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the x coordinate of the space point.
  float &x(Index index) noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    assert(index < m_xColumn->size() && "Index out of bounds");
    return entry(m_xColumn->proxy(), index);
  }
  /// Mutable access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the y coordinate of the space point.
  float &y(Index index) noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    assert(index < m_yColumn->size() && "Index out of bounds");
    return entry(m_yColumn->proxy(), index);
  }
  /// Mutable access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the z coordinate of the space point.
  float &z(Index index) noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    assert(index < m_zColumn->size() && "Index out of bounds");
    return entry(m_zColumn->proxy(), index);
  }
  /// Mutable access to the r coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the r coordinate of the space point.
  float &r(Index index) noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    assert(index < m_rColumn->size() && "Index out of bounds");
    return entry(m_rColumn->proxy(), index);
  }
  /// Mutable access to the phi coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the phi coordinate of the space point.
  float &phi(Index index) noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    assert(index < m_phiColumn->size() && "Index out of bounds");
    return entry(m_phiColumn->proxy(), index);
  }
  /// Mutable access to the time information of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the time information of the space point.
  float &time(Index index) noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    assert(index < m_timeColumn->size() && "Index out of bounds");
    return entry(m_timeColumn->proxy(), index);
  }
  /// Mutable access to the variance in Z direction of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the variance in Z direction of the space point.
  float &varianceZ(Index index) noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    assert(index < m_varianceZColumn->size() && "Index out of bounds");
    return entry(m_varianceZColumn->proxy(), index);
  }
  /// Mutable access to the variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the variance in R direction of the space point.
  float &varianceR(Index index) noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    assert(index < m_varianceRColumn->size() && "Index out of bounds");
    return entry(m_varianceRColumn->proxy(), index);
  }
  /// Mutable access to the top strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the top strip vector of the space point.
  Eigen::Vector3f &topStripVector(Index index) noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    assert(index < m_topStripVectorColumn->size() && "Index out of bounds");
    return entry(m_topStripVectorColumn->proxy(), index);
  }
  /// Mutable access to the bottom strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the bottom strip vector of the space point.
  Eigen::Vector3f &bottomStripVector(Index index) noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    assert(index < m_bottomStripVectorColumn->size() && "Index out of bounds");
    return entry(m_bottomStripVectorColumn->proxy(), index);
  }
  /// Mutable access to the strip center distance of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the strip center distance of the space point.
  Eigen::Vector3f &stripCenterDistance(Index index) noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    assert(index < m_stripCenterDistanceColumn->size() &&
           "Index out of bounds");
    return entry(m_stripCenterDistanceColumn->proxy(), index);
  }
  /// Mutable access to the top strip center of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the top strip center of the space point.
  Eigen::Vector3f &topStripCenter(Index index) noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    assert(index < m_topStripCenterColumn->size() && "Index out of bounds");
    return entry(m_topStripCenterColumn->proxy(), index);
  }
  /// Mutable access to the copy from index of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the copy from index of the space point.
  std::size_t &copyFromIndex(Index index) noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    assert(index < m_copyFromIndexColumn->size() && "Index out of bounds");
    return entry(m_copyFromIndexColumn->proxy(), index);
  }

  /// Const access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A const span to the source links at the given index.
  std::span<const SourceLink> sourceLinks(Index index) const noexcept {
    assert(m_sourceLinkOffsetColumn.has_value() &&
           m_sourceLinkCountColumn.has_value() &&
           "Column 'sourceLinks' does not exist");
    assert(index < m_sourceLinkOffsetColumn->size() &&
           index < m_sourceLinkCountColumn->size() && "Index out of bounds");
    return std::span<const SourceLink>(
        m_sourceLinks.data() + entry(m_sourceLinkOffsetColumn->proxy(), index),
        entry(m_sourceLinkCountColumn->proxy(), index));
  }
  /// Const access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the x coordinate of the space point.
  float x(Index index) const noexcept {
    assert(m_xColumn.has_value() && "Column 'x' does not exist");
    assert(index < m_xColumn->size() && "Index out of bounds");
    return entry(m_xColumn->proxy(), index);
  }
  /// Const access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the y coordinate of the space point.
  float y(Index index) const noexcept {
    assert(m_yColumn.has_value() && "Column 'y' does not exist");
    assert(index < m_yColumn->size() && "Index out of bounds");
    return entry(m_yColumn->proxy(), index);
  }
  /// Const access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the z coordinate of the space point.
  float z(Index index) const noexcept {
    assert(m_zColumn.has_value() && "Column 'z' does not exist");
    assert(index < m_zColumn->size() && "Index out of bounds");
    return entry(m_zColumn->proxy(), index);
  }
  /// Const access to the r coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the r coordinate of the space point.
  float r(Index index) const noexcept {
    assert(m_rColumn.has_value() && "Column 'r' does not exist");
    assert(index < m_rColumn->size() && "Index out of bounds");
    return entry(m_rColumn->proxy(), index);
  }
  /// Const access to the phi coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the phi coordinate of the space point.
  float phi(Index index) const noexcept {
    assert(m_phiColumn.has_value() && "Column 'phi' does not exist");
    assert(index < m_phiColumn->size() && "Index out of bounds");
    return entry(m_phiColumn->proxy(), index);
  }
  /// Const access to the time information of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the time information of the space point.
  float time(Index index) const noexcept {
    assert(m_timeColumn.has_value() && "Column 'time' does not exist");
    assert(index < m_timeColumn->size() && "Index out of bounds");
    return entry(m_timeColumn->proxy(), index);
  }
  /// Const access to the variance in Z direction of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the variance in Z direction of the space point.
  float varianceZ(Index index) const noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Column 'varianceZ' does not exist");
    assert(index < m_varianceZColumn->size() && "Index out of bounds");
    return entry(m_varianceZColumn->proxy(), index);
  }
  /// Const access to the variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the variance in R direction of the space point.
  float varianceR(Index index) const noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Column 'varianceR' does not exist");
    assert(index < m_varianceRColumn->size() && "Index out of bounds");
    return entry(m_varianceRColumn->proxy(), index);
  }
  /// Const access to the top strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the top strip vector of the space point.
  const Eigen::Vector3f &topStripVector(Index index) const noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Column 'topStripVector' does not exist");
    assert(index < m_topStripVectorColumn->size() && "Index out of bounds");
    return entry(m_topStripVectorColumn->proxy(), index);
  }
  /// Const access to the bottom strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the bottom strip vector of the space point.
  const Eigen::Vector3f &bottomStripVector(Index index) const noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Column 'bottomStripVector' does not exist");
    assert(index < m_bottomStripVectorColumn->size() && "Index out of bounds");
    return entry(m_bottomStripVectorColumn->proxy(), index);
  }
  /// Const access to the strip center distance of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the strip center distance of the space point.
  const Eigen::Vector3f &stripCenterDistance(Index index) const noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Column 'stripCenterDistance' does not exist");
    assert(index < m_stripCenterDistanceColumn->size() &&
           "Index out of bounds");
    return entry(m_stripCenterDistanceColumn->proxy(), index);
  }
  /// Const access to the top strip center of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the top strip center of the space point.
  const Eigen::Vector3f &topStripCenter(Index index) const noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Column 'topStripCenter' does not exist");
    assert(index < m_topStripCenterColumn->size() && "Index out of bounds");
    return entry(m_topStripCenterColumn->proxy(), index);
  }
  /// Const access to the copy from index of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the copy from index of the space point.
  std::size_t copyFromIndex(Index index) const noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Column 'copyFromIndex' does not exist");
    assert(index < m_copyFromIndexColumn->size() && "Index out of bounds");
    return entry(m_copyFromIndexColumn->proxy(), index);
  }

  /// Returns the entry for the given column and index.
  /// @param column The column to access.
  /// @param index The index of the space point.
  /// @return A mutable reference to the entry in the column at the given index.
  template <typename T>
  T &entry(SpacePointColumnProxy<T> column, Index index) noexcept {
    return column.data()[index];
  }

  /// Returns the entry for the given column and index.
  /// @param column The column to access.
  /// @param index The index of the space point.
  /// @return A const reference to the entry in the column at the given index.
  template <typename T>
  const T &entry(const SpacePointColumnProxy<T> &column,
                 Index index) const noexcept {
    return column.data()[index];
  }

  template <bool read_only>
  class Iterator {
   public:
    static constexpr bool ReadOnly = read_only;

    using Container = const_if_t<ReadOnly, SpacePointContainer2>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = SpacePointProxy2<ReadOnly>;
    using difference_type = std::ptrdiff_t;

    Iterator() noexcept = default;
    Iterator(Container &container, Index index) noexcept
        : m_container(&container), m_index(index) {}

    Iterator &operator++() noexcept {
      ++m_index;
      return *this;
    }
    Iterator operator++(int) noexcept {
      Iterator tmp(*this);
      ++(*this);
      return tmp;
    }

    value_type operator*() const noexcept {
      return value_type(*m_container, m_index);
    }

   private:
    Container *m_container{};
    Index m_index{};

    friend bool operator==(const Iterator &a, const Iterator &b) noexcept {
      return a.m_index == b.m_index;
    }
  };
  using iterator = Iterator<false>;
  using const_iterator = Iterator<true>;

  iterator begin() noexcept { return iterator(*this, 0); }
  iterator end() noexcept { return iterator(*this, size()); }

  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

  template <bool read_only>
  class Range {
   public:
    static constexpr bool ReadOnly = read_only;
    using Container = const_if_t<ReadOnly, SpacePointContainer2>;
    using Iterator = Iterator<read_only>;

    Range(Container &container, const IndexRange &range) noexcept
        : m_container(&container), m_range(range) {}

    std::size_t size() const noexcept { return m_range.second - m_range.first; }
    bool empty() const noexcept { return size() == 0; }

    Iterator begin() const noexcept {
      return Iterator(*m_container, m_range.first);
    }
    Iterator end() const noexcept {
      return Iterator(*m_container, m_range.second);
    }

   private:
    Container *m_container{};
    IndexRange m_range{};
  };
  using MutableRange = Range<false>;
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

  template <bool read_only>
  class Subset {
   public:
    static constexpr bool ReadOnly = read_only;
    using Container = const_if_t<ReadOnly, SpacePointContainer2>;

    class Iterator {
     public:
      using Container = const_if_t<ReadOnly, SpacePointContainer2>;
      using SubsetIterator = IndexSubset::iterator;

      using iterator_category = std::forward_iterator_tag;
      using value_type = SpacePointProxy2<ReadOnly>;
      using difference_type = std::ptrdiff_t;

      Iterator() noexcept = default;
      Iterator(Container &container, SubsetIterator iterator) noexcept
          : m_container(&container), m_iterator(iterator) {}

      Iterator &operator++() noexcept {
        ++m_iterator;
        return *this;
      }
      Iterator operator++(int) noexcept {
        Iterator tmp(*this);
        ++(*this);
        return tmp;
      }

      value_type operator*() const noexcept {
        return value_type(*m_container, *m_iterator);
      }

     private:
      Container *m_container{};
      SubsetIterator m_iterator{};

      friend bool operator==(const Iterator &a, const Iterator &b) noexcept {
        return a.m_iterator == b.m_iterator;
      }
    };
    using iterator = Iterator;

    Subset(Container &container, const IndexSubset &subset) noexcept
        : m_container(&container), m_subset(subset) {}

    std::size_t size() const noexcept { return m_subset.size(); }
    bool empty() const noexcept { return size() == 0; }

    iterator begin() const noexcept {
      return iterator(*m_container, m_subset.begin());
    }
    iterator end() const noexcept {
      return iterator(*m_container, m_subset.end());
    }

   private:
    Container *m_container{};
    IndexSubset m_subset{};
  };
  using MutableSubset = Subset<false>;
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

  class IndexIteratorRange {
   public:
    class Iterator {
     public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = SpacePointIndex2;
      using difference_type = std::ptrdiff_t;

      Iterator() noexcept = default;
      explicit Iterator(SpacePointIndex2 index) noexcept : m_index{index} {}

      Iterator &operator++() noexcept {
        ++m_index;
        return *this;
      }
      Iterator operator++(int) noexcept {
        Iterator tmp(*this);
        ++(*this);
        return tmp;
      }

      value_type operator*() const noexcept { return m_index; }

     private:
      SpacePointIndex2 m_index{0};

      friend bool operator==(const Iterator &a, const Iterator &b) noexcept {
        return a.m_index == b.m_index;
      }
    };
    using iterator = Iterator;

    explicit IndexIteratorRange(SpacePointIndexRange2 range) noexcept
        : m_range(range) {}

    std::size_t size() const noexcept { return m_range.second - m_range.first; }
    bool empty() const noexcept { return size() == 0; }

    iterator begin() const noexcept { return iterator(m_range.first); }
    iterator end() const noexcept { return iterator(m_range.second); }

   private:
    SpacePointIndexRange2 m_range{};
  };

  /// Creates a zipped mutable range of space point data from the given columns.
  /// @param columns The columns to zip.
  /// @return A zipped mutable range of space point data.
  template <typename... Columns>
  auto zip(const Columns &...columns) noexcept {
    return Acts::zip(IndexIteratorRange({0, size()}), columns.data()...);
  }
  /// Creates a zipped const range of space point data from the given columns.
  /// @param columns The columns to zip.
  /// @return A zipped const range of space point data.
  template <typename... Columns>
  auto zip(const Columns &...columns) const noexcept {
    return Acts::zip(IndexIteratorRange({0, size()}), columns.data()...);
  }

 private:
  std::uint32_t m_size{0};

  std::unordered_map<std::string,
                     std::pair<SpacePointColumnHolderBase *,
                               std::unique_ptr<SpacePointColumnHolderBase>>>
      m_namedColumns;
  SpacePointColumns m_knownColumns{SpacePointColumns::None};

  std::vector<SourceLink> m_sourceLinks;

  std::optional<SpacePointColumnHolder<SpacePointIndex2>>
      m_sourceLinkOffsetColumn;
  std::optional<SpacePointColumnHolder<std::uint8_t>> m_sourceLinkCountColumn;

  std::optional<SpacePointColumnHolder<float>> m_xColumn;
  std::optional<SpacePointColumnHolder<float>> m_yColumn;
  std::optional<SpacePointColumnHolder<float>> m_zColumn;

  // cylindrical coordinates
  std::optional<SpacePointColumnHolder<float>> m_rColumn;
  std::optional<SpacePointColumnHolder<float>> m_phiColumn;
  // time information
  std::optional<SpacePointColumnHolder<float>> m_timeColumn;
  // covariance information
  std::optional<SpacePointColumnHolder<float>> m_varianceZColumn;
  std::optional<SpacePointColumnHolder<float>> m_varianceRColumn;
  // strip information
  std::optional<SpacePointColumnHolder<Eigen::Vector3f>> m_topStripVectorColumn;
  std::optional<SpacePointColumnHolder<Eigen::Vector3f>>
      m_bottomStripVectorColumn;
  std::optional<SpacePointColumnHolder<Eigen::Vector3f>>
      m_stripCenterDistanceColumn;
  std::optional<SpacePointColumnHolder<Eigen::Vector3f>> m_topStripCenterColumn;
  // copy information
  std::optional<SpacePointColumnHolder<std::size_t>> m_copyFromIndexColumn;

  static auto knownColumnMaks() noexcept {
    return std::tuple(
        SpacePointColumns::SourceLinks, SpacePointColumns::SourceLinks,
        SpacePointColumns::X, SpacePointColumns::Y, SpacePointColumns::Z,
        SpacePointColumns::R, SpacePointColumns::Phi, SpacePointColumns::Time,
        SpacePointColumns::VarianceZ, SpacePointColumns::VarianceR,
        SpacePointColumns::TopStripVector, SpacePointColumns::BottomStripVector,
        SpacePointColumns::StripCenterDistance,
        SpacePointColumns::TopStripCenter, SpacePointColumns::CopyFromIndex);
  }

  static auto knownColumnNames() noexcept {
    return std::tuple("sourceLinkOffset", "sourceLinkCount", "x", "y", "z", "r",
                      "phi", "time", "varianceZ", "varianceR", "topStripVector",
                      "bottomStripVector", "stripCenterDistance",
                      "topStripCenter", "copyFromIndex");
  }

  static auto knownColumnDefaults() noexcept {
    return std::tuple(SpacePointIndex2{0}, std::uint8_t{0}, float{0}, float{0},
                      float{0}, float{0}, float{0}, float{NoTime}, float{0},
                      float{0}, Eigen::Vector3f{0, 0, 0},
                      Eigen::Vector3f{0, 0, 0}, Eigen::Vector3f{0, 0, 0},
                      Eigen::Vector3f{0, 0, 0}, std::size_t{0});
  }

  auto knownColumns() & noexcept {
    return std::tie(m_sourceLinkOffsetColumn, m_sourceLinkCountColumn,
                    m_xColumn, m_yColumn, m_zColumn, m_rColumn, m_phiColumn,
                    m_timeColumn, m_varianceZColumn, m_varianceRColumn,
                    m_topStripVectorColumn, m_bottomStripVectorColumn,
                    m_stripCenterDistanceColumn, m_topStripCenterColumn,
                    m_copyFromIndexColumn);
  }
  auto knownColumns() const & noexcept {
    return std::tie(m_sourceLinkOffsetColumn, m_sourceLinkCountColumn,
                    m_xColumn, m_yColumn, m_zColumn, m_rColumn, m_phiColumn,
                    m_timeColumn, m_varianceZColumn, m_varianceRColumn,
                    m_topStripVectorColumn, m_bottomStripVectorColumn,
                    m_stripCenterDistanceColumn, m_topStripCenterColumn,
                    m_copyFromIndexColumn);
  }
  auto knownColumns() && noexcept {
    return std::tuple(
        std::move(m_sourceLinkOffsetColumn), std::move(m_sourceLinkCountColumn),
        std::move(m_xColumn), std::move(m_yColumn), std::move(m_zColumn),
        std::move(m_rColumn), std::move(m_phiColumn), std::move(m_timeColumn),
        std::move(m_varianceZColumn), std::move(m_varianceRColumn),
        std::move(m_topStripVectorColumn), std::move(m_bottomStripVectorColumn),
        std::move(m_stripCenterDistanceColumn),
        std::move(m_topStripCenterColumn), std::move(m_copyFromIndexColumn));
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
    auto proxy = holder->proxy();
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

}  // namespace Acts::Experimental

#include "Acts/EventData/SpacePointContainer2.ipp"
