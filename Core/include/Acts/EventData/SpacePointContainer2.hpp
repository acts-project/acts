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

#include <cassert>
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

enum class SpacePointKnownExtraColumn : std::uint32_t {
  None = 0,  ///< No extra columns

  R = 1 << 0,                    ///< Radial coordinate
  Phi = 1 << 1,                  ///< Azimuthal angle
  Time = 1 << 2,                 ///< Time information
  VarianceZ = 1 << 3,            ///< Variance in Z direction
  VarianceR = 1 << 4,            ///< Variance in radial direction
  TopStripVector = 1 << 5,       ///< Vector for the top strip
  BottomStripVector = 1 << 6,    ///< Vector for the bottom strip
  StripCenterDistance = 1 << 7,  ///< Distance to the strip center
  TopStripCenter = 1 << 8,       ///< Center of the top strip
  CopyFromIndex = 1 << 9,        ///< Copy from index

  /// All strip-related columns
  Strip =
      TopStripVector | BottomStripVector | StripCenterDistance | TopStripCenter
};

ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SpacePointKnownExtraColumn);

/// Additional column of data that can be added to the space point container.
/// The column is indexed by the space point index.
template <typename T>
class SpacePointExtraColumnProxy {
 public:
  using ValueType = T;
  using ContainerType = std::vector<ValueType>;

  explicit SpacePointExtraColumnProxy(const ContainerType &data)
      : m_data(&data) {}
  SpacePointExtraColumnProxy(const SpacePointExtraColumnProxy &other) = default;
  SpacePointExtraColumnProxy(SpacePointExtraColumnProxy &&other) noexcept =
      default;
  SpacePointExtraColumnProxy &operator=(
      const SpacePointExtraColumnProxy &other) = default;
  SpacePointExtraColumnProxy &operator=(
      SpacePointExtraColumnProxy &&other) noexcept = default;

 private:
  const ContainerType *m_data{};

  ContainerType &data() { return const_cast<ContainerType &>(*m_data); }
  const ContainerType &data() const { return *m_data; }

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
class SpacePointExtraColumnHolder final : public SpacePointColumnHolderBase {
 public:
  using ValueType = T;
  using ContainerType = std::vector<ValueType>;
  using ProxyType = SpacePointExtraColumnProxy<ValueType>;

  SpacePointExtraColumnHolder() = default;
  explicit SpacePointExtraColumnHolder(ValueType defaultValue)
      : m_default(std::move(defaultValue)) {}

  ProxyType proxy() const { return ProxyType(m_data); }

  std::unique_ptr<SpacePointColumnHolderBase> copy() const override {
    return std::make_unique<SpacePointExtraColumnHolder<T>>(*this);
  }

  std::size_t size() const override { return m_data.size(); }
  void reserve(std::size_t size) override { m_data.reserve(size); }
  void clear() override { m_data.clear(); }
  void resize(std::size_t size) override { m_data.resize(size, m_default); }
  void emplace_back() override { m_data.emplace_back(m_default); }

 private:
  ValueType m_default{};
  ContainerType m_data;
};

/// A container for space points, which can hold additional columns of data
/// and allows for efficient access to space points and their associated source
/// links. Individual space points are addressed via index. A proxy object
/// simplifies the handling.
class SpacePointContainer2 {
 public:
  using IndexType = SpacePointIndex2;
  using IndexRangeType = SpacePointIndexRange2;
  using MutableProxyType = MutableSpacePointProxy2;
  using ConstProxyType = ConstSpacePointProxy2;

  /// Constructs and empty space point container.
  SpacePointContainer2() noexcept = default;

  /// Constructs a copy of the given space point container.
  /// The extra columns are copied as well.
  /// @param other The space point container to copy.
  SpacePointContainer2(const SpacePointContainer2 &other) noexcept;

  /// Move constructs a space point container.
  /// The extra columns are moved as well.
  /// @param other The space point container to move.
  SpacePointContainer2(SpacePointContainer2 &&other) noexcept;

  /// Detructs the space point container.
  ~SpacePointContainer2() noexcept = default;

  /// Assignment operator for copying a space point container.
  /// The extra columns are copied as well.
  /// @param other The space point container to copy.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(const SpacePointContainer2 &other) noexcept;

  /// Move assignment operator for a space point container.
  /// The extra columns are moved as well.
  /// @param other The space point container to move.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(SpacePointContainer2 &&other) noexcept;

  /// Returns the number of space points in the container.
  /// @return The number of space points in the container.
  std::size_t size() const noexcept { return m_sourceLinkOffsets.size(); }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of space points.
  /// This will reserve space for the source links and the extra columns as
  /// well.
  /// @param size The number of space points to reserve space for.
  /// @param averageSourceLinks The average number of source links per space point.
  void reserve(std::size_t size, float averageSourceLinks = 1) noexcept;

  /// Clears the container, removing all space points and extra columns.
  void clear() noexcept;

  /// Emplaces a new space point with the given source links and coordinates.
  /// This will create a new space point at the end of the container.
  /// @param sourceLinks The source links associated with the space point.
  /// @param x The x coordinate of the space point.
  /// @param y The y coordinate of the space point.
  /// @param z The z coordinate of the space point.
  /// @return A mutable proxy to the newly created space point.
  MutableProxyType createSpacePoint(std::span<const SourceLink> sourceLinks,
                                    float x, float y, float z) noexcept;

  /// Returns a mutable proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxyType at(IndexType index);
  /// Returns a const proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxyType at(IndexType index) const;

  /// Returns a mutable proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  MutableProxyType operator[](IndexType index) noexcept;
  /// Returns a const proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  ConstProxyType operator[](IndexType index) const noexcept;

  /// Mutable access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the source link at the given index.
  std::span<SourceLink> sourceLinks(IndexType index) {
    assert(index < m_x.size() && "Index out of bounds");
    return std::span<SourceLink>(
        m_sourceLinks.data() + m_sourceLinkOffsets[index],
        m_sourceLinkCounts[index]);
  }
  /// Mutable access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the x coordinate of the space point.
  float &x(IndexType index) noexcept {
    assert(index < m_y.size() && "Index out of bounds");
    return m_x[index];
  }
  /// Mutable access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the y coordinate of the space point.
  float &y(IndexType index) noexcept {
    assert(index < m_y.size() && "Index out of bounds");
    return m_y[index];
  }
  /// Mutable access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the z coordinate of the space point.
  float &z(IndexType index) noexcept {
    assert(index < m_z.size() && "Index out of bounds");
    return m_z[index];
  }

  /// Mutable access to the extra r coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the r coordinate of the space point.
  float &r(IndexType index) noexcept {
    assert(m_rColumn.has_value() && "Extra column 'r' does not exist");
    assert(index < m_rColumn->size() && "Index out of bounds");
    return extra(m_rColumn->proxy(), index);
  }
  /// Mutable access to the extra phi coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the phi coordinate of the space point.
  float &phi(IndexType index) noexcept {
    assert(m_phiColumn.has_value() && "Extra column 'phi' does not exist");
    assert(index < m_phiColumn->size() && "Index out of bounds");
    return extra(m_phiColumn->proxy(), index);
  }
  /// Mutable access to the extra time information of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the time information of the space point.
  float &time(IndexType index) noexcept {
    assert(m_timeColumn.has_value() && "Extra column 'time' does not exist");
    assert(index < m_timeColumn->size() && "Index out of bounds");
    return extra(m_timeColumn->proxy(), index);
  }
  /// Mutable access to the extra variance in Z direction of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the variance in Z direction of the space point.
  float &varianceZ(IndexType index) noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Extra column 'varianceZ' does not exist");
    assert(index < m_varianceZColumn->size() && "Index out of bounds");
    return extra(m_varianceZColumn->proxy(), index);
  }
  /// Mutable access to the extra variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the variance in R direction of the space point.
  float &varianceR(IndexType index) noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Extra column 'varianceR' does not exist");
    assert(index < m_varianceRColumn->size() && "Index out of bounds");
    return extra(m_varianceRColumn->proxy(), index);
  }
  /// Mutable access to the extra top strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the top strip vector of the space point.
  Eigen::Vector3f &topStripVector(IndexType index) noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Extra column 'topStripVector' does not exist");
    assert(index < m_topStripVectorColumn->size() && "Index out of bounds");
    return extra(m_topStripVectorColumn->proxy(), index);
  }
  /// Mutable access to the extra bottom strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the bottom strip vector of the space point.
  Eigen::Vector3f &bottomStripVector(IndexType index) noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Extra column 'bottomStripVector' does not exist");
    assert(index < m_bottomStripVectorColumn->size() && "Index out of bounds");
    return extra(m_bottomStripVectorColumn->proxy(), index);
  }
  /// Mutable access to the extra strip center distance of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the strip center distance of the space point.
  Eigen::Vector3f &stripCenterDistance(IndexType index) noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Extra column 'stripCenterDistance' does not exist");
    assert(index < m_stripCenterDistanceColumn->size() &&
           "Index out of bounds");
    return extra(m_stripCenterDistanceColumn->proxy(), index);
  }
  /// Mutable access to the extra top strip center of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the top strip center of the space point.
  Eigen::Vector3f &topStripCenter(IndexType index) noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Extra column 'topStripCenter' does not exist");
    assert(index < m_topStripCenterColumn->size() && "Index out of bounds");
    return extra(m_topStripCenterColumn->proxy(), index);
  }
  /// Mutable access to the copy from index of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the copy from index of the space point.
  std::size_t &copyFromIndex(IndexType index) noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Extra column 'copyFromIndex' does not exist");
    assert(index < m_copyFromIndexColumn->size() && "Index out of bounds");
    return extra(m_copyFromIndexColumn->proxy(), index);
  }

  /// Const access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A const span to the source links at the given index.
  std::span<const SourceLink> sourceLinks(IndexType index) const noexcept {
    assert(index < m_sourceLinkCounts.size() && "Index out of bounds");
    assert(index < m_sourceLinkOffsets.size() && "Index out of bounds");
    return std::span<const SourceLink>(
        m_sourceLinks.data() + m_sourceLinkOffsets[index],
        m_sourceLinkCounts[index]);
  }
  /// Const access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the x coordinate of the space point.
  float x(IndexType index) const noexcept {
    assert(index < m_x.size() && "Index out of bounds");
    return m_x[index];
  }
  /// Const access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the y coordinate of the space point.
  float y(IndexType index) const noexcept {
    assert(index < m_y.size() && "Index out of bounds");
    return m_y[index];
  }
  /// Const access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the z coordinate of the space point.
  float z(IndexType index) const noexcept {
    assert(index < m_z.size() && "Index out of bounds");
    return m_z[index];
  }

  /// Const access to the extra r coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the r coordinate of the space point.
  float r(IndexType index) const noexcept {
    assert(m_rColumn.has_value() && "Extra column 'r' does not exist");
    assert(index < m_rColumn->size() && "Index out of bounds");
    return extra(m_rColumn->proxy(), index);
  }
  /// Const access to the extra phi coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the phi coordinate of the space point.
  float phi(IndexType index) const noexcept {
    assert(m_phiColumn.has_value() && "Extra column 'phi' does not exist");
    assert(index < m_phiColumn->size() && "Index out of bounds");
    return extra(m_phiColumn->proxy(), index);
  }
  /// Const access to the extra time information of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the time information of the space point.
  float time(IndexType index) const noexcept {
    assert(m_timeColumn.has_value() && "Extra column 'time' does not exist");
    assert(index < m_timeColumn->size() && "Index out of bounds");
    return extra(m_timeColumn->proxy(), index);
  }
  /// Const access to the extra variance in Z direction of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the variance in Z direction of the space point.
  float varianceZ(IndexType index) const noexcept {
    assert(m_varianceZColumn.has_value() &&
           "Extra column 'varianceZ' does not exist");
    assert(index < m_varianceZColumn->size() && "Index out of bounds");
    return extra(m_varianceZColumn->proxy(), index);
  }
  /// Const access to the extra variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the variance in R direction of the space point.
  float varianceR(IndexType index) const noexcept {
    assert(m_varianceRColumn.has_value() &&
           "Extra column 'varianceR' does not exist");
    assert(index < m_varianceRColumn->size() && "Index out of bounds");
    return extra(m_varianceRColumn->proxy(), index);
  }
  /// Const access to the extra top strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the top strip vector of the space point.
  const Eigen::Vector3f &topStripVector(IndexType index) const noexcept {
    assert(m_topStripVectorColumn.has_value() &&
           "Extra column 'topStripVector' does not exist");
    assert(index < m_topStripVectorColumn->size() && "Index out of bounds");
    return extra(m_topStripVectorColumn->proxy(), index);
  }
  /// Const access to the extra bottom strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the bottom strip vector of the space point.
  const Eigen::Vector3f &bottomStripVector(IndexType index) const noexcept {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Extra column 'bottomStripVector' does not exist");
    assert(index < m_bottomStripVectorColumn->size() && "Index out of bounds");
    return extra(m_bottomStripVectorColumn->proxy(), index);
  }
  /// Const access to the extra strip center distance of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the strip center distance of the space point.
  const Eigen::Vector3f &stripCenterDistance(IndexType index) const noexcept {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Extra column 'stripCenterDistance' does not exist");
    assert(index < m_stripCenterDistanceColumn->size() &&
           "Index out of bounds");
    return extra(m_stripCenterDistanceColumn->proxy(), index);
  }
  /// Const access to the extra top strip center of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the top strip center of the space point.
  const Eigen::Vector3f &topStripCenter(IndexType index) const noexcept {
    assert(m_topStripCenterColumn.has_value() &&
           "Extra column 'topStripCenter' does not exist");
    assert(index < m_topStripCenterColumn->size() && "Index out of bounds");
    return extra(m_topStripCenterColumn->proxy(), index);
  }
  /// Const access to the copy from index of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the copy from index of the space point.
  std::size_t copyFromIndex(IndexType index) const noexcept {
    assert(m_copyFromIndexColumn.has_value() &&
           "Extra column 'copyFromIndex' does not exist");
    assert(index < m_copyFromIndexColumn->size() && "Index out of bounds");
    return extra(m_copyFromIndexColumn->proxy(), index);
  }

  template <typename T>
  T &extra(SpacePointExtraColumnProxy<T> column, IndexType index) noexcept {
    return column.data()[index];
  }

  template <typename T>
  const T &extra(const SpacePointExtraColumnProxy<T> &column,
                 IndexType index) const noexcept {
    return column.data()[index];
  }

  /// Creates extra columns based on the specified known extra columns.
  /// This will only create the columns if they do not already exist and fill
  /// them according to the size of the container with default values.
  /// @param columns The known extra columns to create.
  void createExtraColumns(SpacePointKnownExtraColumn columns) noexcept;

  /// Drops the specified extra columns from the container.
  /// This will only drop columns if they exist.
  /// @param columns The extra columns to drop.
  void dropExtraColumns(SpacePointKnownExtraColumn columns) noexcept;

  /// Checks if the container has the given extra columns.
  /// @param columns The extra columns to check for.
  /// @return True if the container has all the specified extra columns, false
  ///         otherwise.
  bool hasExtraColumns(SpacePointKnownExtraColumn columns) const noexcept {
    return (m_knownExtraColumns & columns) == columns;
  }

  /// Returns a proxy to the extra r coordinate column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra r coordinate column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<float> rColumn() const {
    if (!m_rColumn.has_value()) {
      throw std::runtime_error("Extra column 'r' does not exist");
    }
    return m_rColumn->proxy();
  }
  /// Returns a proxy to the extra phi coordinate column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra phi coordinate column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<float> phiColumn() const {
    if (!m_phiColumn.has_value()) {
      throw std::runtime_error("Extra column 'phi' does not exist");
    }
    return m_phiColumn->proxy();
  }
  /// Returns a proxy to the extra time column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra time column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<float> timeColumn() const {
    if (!m_timeColumn.has_value()) {
      throw std::runtime_error("Extra column 'time' does not exist");
    }
    return m_timeColumn->proxy();
  }
  /// Returns a proxy to the extra variance in Z direction column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra variance in Z direction column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<float> varianceZColumn() const {
    if (!m_varianceZColumn.has_value()) {
      throw std::runtime_error("Extra column 'varianceZ' does not exist");
    }
    return m_varianceZColumn->proxy();
  }
  /// Returns a proxy to the extra variance in R direction column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra variance in R direction column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<float> varianceRColumn() const {
    if (!m_varianceRColumn.has_value()) {
      throw std::runtime_error("Extra column 'varianceR' does not exist");
    }
    return m_varianceRColumn->proxy();
  }
  /// Returns a proxy to the extra top strip vector column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra top strip vector column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<Eigen::Vector3f> topStripVectorColumn() const {
    if (!m_topStripVectorColumn.has_value()) {
      throw std::runtime_error("Extra column 'topStripVector' does not exist");
    }
    return m_topStripVectorColumn->proxy();
  }
  /// Returns a proxy to the extra bottom strip vector column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra bottom strip vector column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<Eigen::Vector3f> bottomStripVectorColumn() const {
    if (!m_bottomStripVectorColumn.has_value()) {
      throw std::runtime_error(
          "Extra column 'bottomStripVector' does not exist");
    }
    return m_bottomStripVectorColumn->proxy();
  }
  /// Returns a proxy to the extra strip center distance column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra strip center distance column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<Eigen::Vector3f> stripCenterDistanceColumn()
      const {
    if (!m_stripCenterDistanceColumn.has_value()) {
      throw std::runtime_error(
          "Extra column 'stripCenterDistance' does not exist");
    }
    return m_stripCenterDistanceColumn->proxy();
  }
  /// Returns a proxy to the extra top strip center column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra top strip center column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<Eigen::Vector3f> topStripCenterColumn() const {
    if (!m_topStripCenterColumn.has_value()) {
      throw std::runtime_error("Extra column 'topStripCenter' does not exist");
    }
    return m_topStripCenterColumn->proxy();
  }
  /// Returns a proxy to the extra copy from index column.
  /// If the column does not exist, an exception is thrown.
  /// @return A proxy to the extra copy from index column.
  /// @throws std::runtime_error if the column does not exist.
  SpacePointExtraColumnProxy<std::size_t> copyFromIndexColumn() const {
    if (!m_copyFromIndexColumn.has_value()) {
      throw std::runtime_error("Extra column 'copyFromIndex' does not exist");
    }
    return m_copyFromIndexColumn->proxy();
  }

  /// Creates a new column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created column.
  /// @throws std::runtime_error if a column with the same name already exists.
  template <typename T>
  SpacePointExtraColumnProxy<T> createExtraColumn(const std::string &name) {
    return createExtraColumnImpl<SpacePointExtraColumnHolder<T>>(name);
  }

  /// Drops the extra column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @throws std::runtime_error if the column does not exist.
  void dropExtraColumn(const std::string &name);

  /// Checks if an extra column with the given name exists.
  /// @param name The name of the column.
  /// @return True if the column exists, false otherwise.
  bool hasExtraColumn(const std::string &name) const noexcept;

  /// Returns a mutable reference to the extra column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the extra column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  SpacePointExtraColumnProxy<T> extraColumn(const std::string &name) const {
    return extraColumnImpl<SpacePointExtraColumnHolder<T>>(name);
  }

  template <bool read_only>
  class Iterator {
   public:
    static constexpr bool ReadOnly = read_only;

    using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = SpacePointProxy2<ReadOnly>;
    using difference_type = std::ptrdiff_t;

    Iterator() noexcept = default;
    Iterator(ContainerType &container, IndexType index) noexcept
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
    ContainerType *m_container{};
    IndexType m_index{};

    friend bool operator==(const Iterator &a, const Iterator &b) noexcept {
      return a.m_index == b.m_index && a.m_container == b.m_container;
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
    using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

    using iterator = Iterator<read_only>;
    using const_iterator = Iterator<true>;

    Range(ContainerType &container, const IndexRangeType &range) noexcept
        : m_container(&container), m_range(range) {}

    std::size_t size() const noexcept { return m_range.second - m_range.first; }
    bool empty() const noexcept { return size() == 0; }

    iterator begin() const noexcept {
      return iterator(*m_container, m_range.first);
    }
    iterator end() const noexcept {
      return iterator(*m_container, m_range.second);
    }

    const_iterator cbegin() const noexcept {
      return const_iterator(*m_container, m_range.first);
    }
    const_iterator cend() const noexcept {
      return const_iterator(*m_container, m_range.second);
    }

   private:
    ContainerType *m_container{};
    IndexRangeType m_range{};
  };
  using MutableRange = Range<false>;
  using ConstRange = Range<true>;

  MutableRange range(const IndexRangeType &range) noexcept {
    return MutableRange(*this, range);
  }
  ConstRange range(const IndexRangeType &range) const noexcept {
    return ConstRange(*this, range);
  }

 private:
  std::vector<float> m_x;
  std::vector<float> m_y;
  std::vector<float> m_z;
  std::vector<std::size_t> m_sourceLinkOffsets;
  std::vector<std::uint8_t> m_sourceLinkCounts;
  std::vector<SourceLink> m_sourceLinks;

  // known extra columns
  SpacePointKnownExtraColumn m_knownExtraColumns{
      SpacePointKnownExtraColumn::None};
  // cylindrical coordinates
  std::optional<SpacePointExtraColumnHolder<float>> m_rColumn;
  std::optional<SpacePointExtraColumnHolder<float>> m_phiColumn;
  // time information
  std::optional<SpacePointExtraColumnHolder<float>> m_timeColumn;
  // covariance information
  std::optional<SpacePointExtraColumnHolder<float>> m_varianceZColumn;
  std::optional<SpacePointExtraColumnHolder<float>> m_varianceRColumn;
  // strip information
  std::optional<SpacePointExtraColumnHolder<Eigen::Vector3f>>
      m_topStripVectorColumn;
  std::optional<SpacePointExtraColumnHolder<Eigen::Vector3f>>
      m_bottomStripVectorColumn;
  std::optional<SpacePointExtraColumnHolder<Eigen::Vector3f>>
      m_stripCenterDistanceColumn;
  std::optional<SpacePointExtraColumnHolder<Eigen::Vector3f>>
      m_topStripCenterColumn;
  // copy information
  std::optional<SpacePointExtraColumnHolder<std::size_t>> m_copyFromIndexColumn;

  std::unordered_map<std::string, std::unique_ptr<SpacePointColumnHolderBase>>
      m_namedExtraColumns;

  std::vector<SpacePointColumnHolderBase *> m_extraColumns;

  auto knownExtraColumns() & noexcept {
    return std::tie(m_rColumn, m_phiColumn, m_timeColumn, m_varianceZColumn,
                    m_varianceRColumn, m_topStripVectorColumn,
                    m_bottomStripVectorColumn, m_stripCenterDistanceColumn,
                    m_topStripCenterColumn, m_copyFromIndexColumn);
  }
  auto knownExtraColumns() const & noexcept {
    return std::tie(m_rColumn, m_phiColumn, m_timeColumn, m_varianceZColumn,
                    m_varianceRColumn, m_topStripVectorColumn,
                    m_bottomStripVectorColumn, m_stripCenterDistanceColumn,
                    m_topStripCenterColumn, m_copyFromIndexColumn);
  }
  auto knownExtraColumns() && noexcept {
    return std::tuple(
        std::move(m_rColumn), std::move(m_phiColumn), std::move(m_timeColumn),
        std::move(m_varianceZColumn), std::move(m_varianceRColumn),
        std::move(m_topStripVectorColumn), std::move(m_bottomStripVectorColumn),
        std::move(m_stripCenterDistanceColumn),
        std::move(m_topStripCenterColumn), std::move(m_copyFromIndexColumn));
  }

  void copyExtraColumns(const SpacePointContainer2 &other);
  void moveExtraColumns(SpacePointContainer2 &other) noexcept;

  void initializeExtraColumns() noexcept;

  template <typename Holder>
  auto createExtraColumnImpl(const std::string &name) {
    if (hasExtraColumn(name)) {
      throw std::runtime_error("Extra column already exists: " + name);
    }
    auto holder = std::make_unique<Holder>();
    holder->resize(size());
    auto proxy = holder->proxy();
    m_extraColumns.push_back(holder.get());
    m_namedExtraColumns[name] = std::move(holder);
    return proxy;
  }
  template <typename Holder>
  auto extraColumnImpl(const std::string &name) const {
    auto it = m_namedExtraColumns.find(name);
    if (it == m_namedExtraColumns.end()) {
      throw std::runtime_error("Extra column not found: " + name);
    }
    auto &holder = dynamic_cast<Holder &>(*it->second);
    return holder.proxy();
  }
};

}  // namespace Acts::Experimental

#include "Acts/EventData/SpacePointContainer2.ipp"
