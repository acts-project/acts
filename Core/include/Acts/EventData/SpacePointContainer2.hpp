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
  SpacePointContainer2() = default;

  /// Constructs a copy of the given space point container.
  /// The extra columns are copied as well.
  /// @param other The space point container to copy.
  SpacePointContainer2(const SpacePointContainer2 &other);

  /// Move constructs a space point container.
  /// The extra columns are moved as well.
  /// @param other The space point container to move.
  SpacePointContainer2(SpacePointContainer2 &&other) noexcept = default;

  /// Detructs the space point container.
  ~SpacePointContainer2() = default;

  /// Assignment operator for copying a space point container.
  /// The extra columns are copied as well.
  /// @param other The space point container to copy.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(const SpacePointContainer2 &other);

  /// Move assignment operator for a space point container.
  /// The extra columns are moved as well.
  /// @param other The space point container to move.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(SpacePointContainer2 &&other) noexcept =
      default;

  /// Returns the number of space points in the container.
  /// @return The number of space points in the container.
  std::size_t size() const { return m_sourceLinkOffsets.size(); }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const { return size() == 0; }

  /// Reserves space for the given number of space points.
  /// This will reserve space for the source links and the extra columns as
  /// well.
  /// @param size The number of space points to reserve space for.
  /// @param averageSourceLinks The average number of source links per space point.
  void reserve(std::size_t size, float averageSourceLinks = 1);

  /// Clears the container, removing all space points and extra columns.
  void clear();

  /// Emplaces a new space point with the given source links and coordinates.
  /// This will create a new space point at the end of the container.
  /// @param sourceLinks The source links associated with the space point.
  /// @param x The x coordinate of the space point.
  /// @param y The y coordinate of the space point.
  /// @param z The z coordinate of the space point.
  /// @return A mutable proxy to the newly created space point.
  MutableProxyType createSpacePoint(std::span<const SourceLink> sourceLinks,
                                    float x, float y, float z);

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

  /// Mutable access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the source link at the given index.
  std::span<SourceLink> sourceLinks(IndexType index) {
    assert(index < m_xyz.size() && "Index out of bounds");
    return std::span<SourceLink>(
        m_sourceLinks.data() + m_sourceLinkOffsets[index],
        m_sourceLinkCounts[index]);
  }
  /// Mutable access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the x coordinate of the space point.
  float &x(IndexType index) {
    assert(index < m_xyz.size() && "Index out of bounds");
    return m_xyz[index * 3];
  }
  /// Mutable access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the y coordinate of the space point.
  float &y(IndexType index) {
    assert(index < m_xyz.size() && "Index out of bounds");
    return m_xyz[index * 3 + 1];
  }
  /// Mutable access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the z coordinate of the space point.
  float &z(IndexType index) {
    assert(index < m_xyz.size() && "Index out of bounds");
    return m_xyz[index * 3 + 2];
  }

  /// Mutable access to the extra r coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the r coordinate of the space point.
  float &r(IndexType index) {
    assert(m_rColumn.has_value() && "Extra column 'r' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_rColumn->proxy(), index);
  }
  /// Mutable access to the extra phi coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the phi coordinate of the space point.
  float &phi(IndexType index) {
    assert(m_phiColumn.has_value() && "Extra column 'phi' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_phiColumn->proxy(), index);
  }
  /// Mutable access to the extra time information of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the time information of the space point.
  float &time(IndexType index) {
    assert(m_timeColumn.has_value() && "Extra column 'time' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_timeColumn->proxy(), index);
  }
  /// Mutable access to the extra variance in Z direction of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the variance in Z direction of the space point.
  float &varianceZ(IndexType index) {
    assert(m_varianceZColumn.has_value() &&
           "Extra column 'varianceZ' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceZColumn->proxy(), index);
  }
  /// Mutable access to the extra variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the variance in R direction of the space point.
  float &varianceR(IndexType index) {
    assert(m_varianceRColumn.has_value() &&
           "Extra column 'varianceR' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceRColumn->proxy(), index);
  }
  /// Mutable access to the extra top strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the top strip vector of the space point.
  Eigen::Vector3f &topStripVector(IndexType index) {
    assert(m_topStripVectorColumn.has_value() &&
           "Extra column 'topStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripVectorColumn->proxy(), index);
  }
  /// Mutable access to the extra bottom strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the bottom strip vector of the space point.
  Eigen::Vector3f &bottomStripVector(IndexType index) {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Extra column 'bottomStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_bottomStripVectorColumn->proxy(), index);
  }
  /// Mutable access to the extra strip center distance of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the strip center distance of the space point.
  Eigen::Vector3f &stripCenterDistance(IndexType index) {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Extra column 'stripCenterDistance' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_stripCenterDistanceColumn->proxy(), index);
  }
  /// Mutable access to the extra top strip center of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the top strip center of the space point.
  Eigen::Vector3f &topStripCenter(IndexType index) {
    assert(m_topStripCenterColumn.has_value() &&
           "Extra column 'topStripCenter' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripCenterColumn->proxy(), index);
  }
  /// Mutable access to the copy from index of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the copy from index of the space point.
  std::size_t &copyFromIndex(IndexType index) {
    assert(m_copyFromIndexColumn.has_value() &&
           "Extra column 'copyFromIndex' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_copyFromIndexColumn->proxy(), index);
  }

  /// Const access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A const span to the source links at the given index.
  std::span<const SourceLink> sourceLinks(IndexType index) const {
    assert(index < m_xyz.size() && "Index out of bounds");
    return std::span<const SourceLink>(
        m_sourceLinks.data() + m_sourceLinkOffsets[index],
        m_sourceLinkCounts[index]);
  }
  /// Const access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the x coordinate of the space point.
  float x(IndexType index) const {
    assert(index < m_xyz.size() && "Index out of bounds");
    return m_xyz[index * 3];
  }
  /// Const access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the y coordinate of the space point.
  float y(IndexType index) const {
    assert(index < m_xyz.size() && "Index out of bounds");
    return m_xyz[index * 3 + 1];
  }
  /// Const access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the z coordinate of the space point.
  float z(IndexType index) const {
    assert(index < m_xyz.size() && "Index out of bounds");
    return m_xyz[index * 3 + 2];
  }

  /// Const access to the extra r coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the r coordinate of the space point.
  float r(IndexType index) const {
    assert(m_rColumn.has_value() && "Extra column 'r' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_rColumn->proxy(), index);
  }
  /// Const access to the extra phi coordinate of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the phi coordinate of the space point.
  float phi(IndexType index) const {
    assert(m_phiColumn.has_value() && "Extra column 'phi' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_phiColumn->proxy(), index);
  }
  /// Const access to the extra time information of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the time information of the space point.
  float time(IndexType index) const {
    assert(m_timeColumn.has_value() && "Extra column 'time' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_timeColumn->proxy(), index);
  }
  /// Const access to the extra variance in Z direction of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the variance in Z direction of the space point.
  float varianceZ(IndexType index) const {
    assert(m_varianceZColumn.has_value() &&
           "Extra column 'varianceZ' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceZColumn->proxy(), index);
  }
  /// Const access to the extra variance in R direction of the space
  /// point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the variance in R direction of the space point.
  float varianceR(IndexType index) const {
    assert(m_varianceRColumn.has_value() &&
           "Extra column 'varianceR' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceRColumn->proxy(), index);
  }
  /// Const access to the extra top strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the top strip vector of the space point.
  const Eigen::Vector3f &topStripVector(IndexType index) const {
    assert(m_topStripVectorColumn.has_value() &&
           "Extra column 'topStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripVectorColumn->proxy(), index);
  }
  /// Const access to the extra bottom strip vector of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the bottom strip vector of the space point.
  const Eigen::Vector3f &bottomStripVector(IndexType index) const {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Extra column 'bottomStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_bottomStripVectorColumn->proxy(), index);
  }
  /// Const access to the extra strip center distance of the space point at
  /// the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the strip center distance of the space point.
  const Eigen::Vector3f &stripCenterDistance(IndexType index) const {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Extra column 'stripCenterDistance' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_stripCenterDistanceColumn->proxy(), index);
  }
  /// Const access to the extra top strip center of the space point at the
  /// given index.
  /// @param index The index of the space point.
  /// @return A const reference to the top strip center of the space point.
  const Eigen::Vector3f &topStripCenter(IndexType index) const {
    assert(m_topStripCenterColumn.has_value() &&
           "Extra column 'topStripCenter' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripCenterColumn->proxy(), index);
  }
  /// Const access to the copy from index of the space point at the given
  /// index.
  /// @param index The index of the space point.
  /// @return A const reference to the copy from index of the space point.
  std::size_t copyFromIndex(IndexType index) const {
    assert(m_copyFromIndexColumn.has_value() &&
           "Extra column 'copyFromIndex' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_copyFromIndexColumn->proxy(), index);
  }

  template <typename T>
  T &extra(SpacePointExtraColumnProxy<T> column, IndexType index) {
    return column.data()[index];
  }

  template <typename T>
  const T &extra(const SpacePointExtraColumnProxy<T> &column,
                 IndexType index) const {
    return column.data()[index];
  }

  /// Creates extra columns based on the specified known extra columns.
  /// This will create the columns if they do not already exist.
  /// @param columns The known extra columns to create.
  void createExtraColumns(SpacePointKnownExtraColumn columns);

  /// Checks if the container has the given extra columns.
  /// @param columns The extra columns to check for.
  /// @return True if the container has all the specified extra columns, false
  ///         otherwise.
  bool hasExtraColumns(SpacePointKnownExtraColumn columns) const;

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

  /// Checks if an extra column with the given name exists.
  /// @param name The name of the column.
  /// @return True if the column exists, false otherwise.
  bool hasExtraColumn(const std::string &name) const {
    return m_namedExtraColumns.find(name) != m_namedExtraColumns.end();
  }

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

    Iterator() = default;
    Iterator(ContainerType &container, IndexType index)
        : m_container(&container), m_index(index) {}

    Iterator &operator++() {
      ++m_index;
      return *this;
    }
    Iterator operator++(int) {
      Iterator tmp(*this);
      ++(*this);
      return tmp;
    }

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};

    friend bool operator==(const Iterator &a, const Iterator &b) {
      return a.m_index == b.m_index && a.m_container == b.m_container;
    }
    friend bool operator!=(const Iterator &a, const Iterator &b) {
      return !(a == b);
    }
  };
  using iterator = Iterator<false>;
  using const_iterator = Iterator<true>;

  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, size()); }

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

  template <bool read_only>
  class Range {
   public:
    static constexpr bool ReadOnly = read_only;
    using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

    using iterator = Iterator<read_only>;
    using const_iterator = Iterator<true>;

    Range(ContainerType &container, const IndexRangeType &range)
        : m_container(&container), m_range(range) {}

    std::size_t size() const { return m_range.second - m_range.first; }
    bool empty() const { return size() == 0; }

    iterator begin() const { return iterator(*m_container, m_range.first); }
    iterator end() const { return iterator(*m_container, m_range.second); }

    const_iterator cbegin() const {
      return const_iterator(*m_container, m_range.first);
    }
    const_iterator cend() const {
      return const_iterator(*m_container, m_range.second);
    }

   private:
    ContainerType *m_container{};
    IndexRangeType m_range{};
  };
  using MutableRange = Range<false>;
  using ConstRange = Range<true>;

  MutableRange range(const IndexRangeType &range) {
    return MutableRange(*this, range);
  }
  ConstRange range(const IndexRangeType &range) const {
    return ConstRange(*this, range);
  }

 private:
  std::vector<float> m_xyz;
  std::vector<std::size_t> m_sourceLinkOffsets;
  std::vector<std::uint8_t> m_sourceLinkCounts;
  std::vector<SourceLink> m_sourceLinks;

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

  template <typename Holder>
  auto createExtraColumnImpl(const std::string &name) {
    auto it = m_namedExtraColumns.find(name);
    if (it != m_namedExtraColumns.end()) {
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
