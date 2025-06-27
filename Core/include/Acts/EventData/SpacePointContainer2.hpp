// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
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

using SpacePointIndex2 = std::size_t;
using SpacePointIndexRange2 = std::pair<SpacePointIndex2, SpacePointIndex2>;

template <bool read_only>
class SpacePointProxy2;
using MutableSpacePointProxy2 = SpacePointProxy2<false>;
using ConstSpacePointProxy2 = SpacePointProxy2<true>;

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

  enum KnownExtraColumn : std::uint32_t {
    R = 1,                      ///< Radial coordinate
    Phi = 2,                    ///< Azimuthal angle
    Time = 4,                   ///< Time information
    VarianceZ = 8,              ///< Variance in Z direction
    VarianceR = 16,             ///< Variance in radial direction
    TopStripVector = 32,        ///< Vector for the top strip
    BottomStripVector = 64,     ///< Vector for the bottom strip
    StripCenterDistance = 128,  ///< Distance to the strip center
    TopStripCenter = 256,       ///< Center of the top strip

    /// All strip-related columns
    Strip = TopStripVector | BottomStripVector | StripCenterDistance |
            TopStripCenter
  };

  /// Constructs and empty space point container.
  SpacePointContainer2() = default;

  /// Constructs a copy of the given space point container.
  /// The extra columns are copied as well.
  /// @param other The space point container to copy.
  SpacePointContainer2(const SpacePointContainer2 &other)
      : m_xyz(other.m_xyz),
        m_sourceLinkOffsets(other.m_sourceLinkOffsets),
        m_sourceLinkCounts(other.m_sourceLinkCounts),
        m_sourceLinks(other.m_sourceLinks),
        m_rColumn(other.m_rColumn),
        m_phiColumn(other.m_phiColumn),
        m_timeColumn(other.m_timeColumn),
        m_varianceZColumn(other.m_varianceZColumn),
        m_varianceRColumn(other.m_varianceRColumn),
        m_topStripVectorColumn(other.m_topStripVectorColumn),
        m_bottomStripVectorColumn(other.m_bottomStripVectorColumn),
        m_stripCenterDistanceColumn(other.m_stripCenterDistanceColumn),
        m_topStripCenterColumn(other.m_topStripCenterColumn) {
    for (const auto &[name, column] : other.m_namedExtraColumns) {
      m_namedExtraColumns.emplace(name, column->copy());
    }
  }

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
  SpacePointContainer2 &operator=(const SpacePointContainer2 &other) {
    if (this == &other) {
      return *this;
    }

    m_xyz = other.m_xyz;
    m_sourceLinkOffsets = other.m_sourceLinkOffsets;
    m_sourceLinkCounts = other.m_sourceLinkCounts;
    m_sourceLinks = other.m_sourceLinks;

    m_extraColumns.clear();
    for (const auto &[name, column] : other.m_namedExtraColumns) {
      m_namedExtraColumns.emplace(name, column->copy());
    }
    return *this;
  }

  /// Move assignment operator for a space point container.
  /// The extra columns are moved as well.
  /// @param other The space point container to move.
  /// @return A reference to this space point container.
  SpacePointContainer2 &operator=(SpacePointContainer2 &&other) noexcept =
      default;

  /// Returns the number of space points in the container.
  /// @return The number of space points in the container.
  std::size_t size() const { return m_xyz.size(); }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const { return size() == 0; }

  /// Reserves space for the given number of space points.
  /// This will reserve space for the source links and the extra columns as
  /// well.
  /// @param size The number of space points to reserve space for.
  /// @param averageSourceLinks The average number of source links per space point.
  void reserve(std::size_t size, float averageSourceLinks = 1) {
    m_xyz.reserve(size * 3);
    m_sourceLinkOffsets.reserve(size);
    m_sourceLinkCounts.reserve(size);
    m_sourceLinks.reserve(static_cast<std::size_t>(size * averageSourceLinks));

    for (auto &column : m_extraColumns) {
      column->reserve(size);
    }
  }
  /// Clears the container, removing all space points and extra columns.
  void clear() {
    m_xyz.clear();
    m_sourceLinkOffsets.clear();
    m_sourceLinkCounts.clear();
    m_sourceLinks.clear();

    for (auto &column : m_extraColumns) {
      column->clear();
    }
  }

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

  float &r(IndexType index) {
    assert(m_rColumn.has_value() && "Extra column 'r' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_rColumn->proxy(), index);
  }
  float &phi(IndexType index) {
    assert(m_phiColumn.has_value() && "Extra column 'phi' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_phiColumn->proxy(), index);
  }
  std::optional<float> &time(IndexType index) {
    assert(m_timeColumn.has_value() && "Extra column 'time' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_timeColumn->proxy(), index);
  }
  float &varianceZ(IndexType index) {
    assert(m_varianceZColumn.has_value() &&
           "Extra column 'varianceZ' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceZColumn->proxy(), index);
  }
  float &varianceR(IndexType index) {
    assert(m_varianceRColumn.has_value() &&
           "Extra column 'varianceR' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceRColumn->proxy(), index);
  }
  Eigen::Vector3f &topStripVector(IndexType index) {
    assert(m_topStripVectorColumn.has_value() &&
           "Extra column 'topStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripVectorColumn->proxy(), index);
  }
  Eigen::Vector3f &bottomStripVector(IndexType index) {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Extra column 'bottomStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_bottomStripVectorColumn->proxy(), index);
  }
  Eigen::Vector3f &stripCenterDistance(IndexType index) {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Extra column 'stripCenterDistance' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_stripCenterDistanceColumn->proxy(), index);
  }
  Eigen::Vector3f &topStripCenter(IndexType index) {
    assert(m_topStripCenterColumn.has_value() &&
           "Extra column 'topStripCenter' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripCenterColumn->proxy(), index);
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

  float r(IndexType index) const {
    assert(m_rColumn.has_value() && "Extra column 'r' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_rColumn->proxy(), index);
  }
  float phi(IndexType index) const {
    assert(m_phiColumn.has_value() && "Extra column 'phi' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_phiColumn->proxy(), index);
  }
  std::optional<float> time(IndexType index) const {
    assert(m_timeColumn.has_value() && "Extra column 'time' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_timeColumn->proxy(), index);
  }
  float varianceZ(IndexType index) const {
    assert(m_varianceZColumn.has_value() &&
           "Extra column 'varianceZ' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceZColumn->proxy(), index);
  }
  float varianceR(IndexType index) const {
    assert(m_varianceRColumn.has_value() &&
           "Extra column 'varianceR' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_varianceRColumn->proxy(), index);
  }
  const Eigen::Vector3f &topStripVector(IndexType index) const {
    assert(m_topStripVectorColumn.has_value() &&
           "Extra column 'topStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripVectorColumn->proxy(), index);
  }
  const Eigen::Vector3f &bottomStripVector(IndexType index) const {
    assert(m_bottomStripVectorColumn.has_value() &&
           "Extra column 'bottomStripVector' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_bottomStripVectorColumn->proxy(), index);
  }
  const Eigen::Vector3f &stripCenterDistance(IndexType index) const {
    assert(m_stripCenterDistanceColumn.has_value() &&
           "Extra column 'stripCenterDistance' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_stripCenterDistanceColumn->proxy(), index);
  }
  const Eigen::Vector3f &topStripCenter(IndexType index) const {
    assert(m_topStripCenterColumn.has_value() &&
           "Extra column 'topStripCenter' does not exist");
    assert(index < m_xyz.size() && "Index out of bounds");
    return extra(m_topStripCenterColumn->proxy(), index);
  }

  /// Additional column of data that can be added to the space point container.
  /// The column is indexed by the space point index.
  template <typename T>
  class ExtraColumnProxy {
   public:
    using ValueType = T;
    using ContainerType = std::vector<ValueType>;

    ExtraColumnProxy(const ExtraColumnProxy &other) = default;
    ExtraColumnProxy(ExtraColumnProxy &&other) noexcept = default;
    ExtraColumnProxy &operator=(const ExtraColumnProxy &other) = default;
    ExtraColumnProxy &operator=(ExtraColumnProxy &&other) noexcept = default;

   private:
    const ContainerType *m_data;

    explicit ExtraColumnProxy(const ContainerType &data) : m_data(&data) {}

    ContainerType &data() { return const_cast<ContainerType &>(*m_data); }
    const ContainerType &data() const { return *m_data; }

    friend class SpacePointContainer2;
  };

  template <typename T>
  T &extra(ExtraColumnProxy<T> column, IndexType index) {
    return column.data()[index];
  }

  template <typename T>
  const T &extra(const ExtraColumnProxy<T> &column, IndexType index) const {
    return column.data()[index];
  }

  /// Creates extra columns based on the specified known extra columns.
  /// This will create the columns if they do not already exist.
  /// @param columns The known extra columns to create.
  void createExtraColumns(std::uint32_t columns) {
    if ((columns & KnownExtraColumn::R) != 0 && !m_rColumn.has_value()) {
      m_rColumn.emplace();
      m_extraColumns.push_back(&*m_rColumn);
    }
    if ((columns & KnownExtraColumn::Phi) != 0 && !m_phiColumn.has_value()) {
      m_phiColumn.emplace();
      m_extraColumns.push_back(&*m_phiColumn);
    }
    if ((columns & KnownExtraColumn::Time) != 0 && !m_timeColumn.has_value()) {
      m_timeColumn.emplace();
      m_extraColumns.push_back(&*m_timeColumn);
    }
    if ((columns & KnownExtraColumn::VarianceZ) != 0 &&
        !m_varianceZColumn.has_value()) {
      m_varianceZColumn.emplace();
      m_extraColumns.push_back(&*m_varianceZColumn);
    }
    if ((columns & KnownExtraColumn::VarianceR) != 0 &&
        !m_varianceRColumn.has_value()) {
      m_varianceRColumn.emplace();
      m_extraColumns.push_back(&*m_varianceRColumn);
    }
    if ((columns & KnownExtraColumn::TopStripVector) != 0 &&
        !m_topStripVectorColumn.has_value()) {
      m_topStripVectorColumn.emplace();
      m_extraColumns.push_back(&*m_topStripVectorColumn);
    }
    if ((columns & KnownExtraColumn::BottomStripVector) != 0 &&
        !m_bottomStripVectorColumn.has_value()) {
      m_bottomStripVectorColumn.emplace();
      m_extraColumns.push_back(&*m_bottomStripVectorColumn);
    }
    if ((columns & KnownExtraColumn::StripCenterDistance) != 0 &&
        !m_stripCenterDistanceColumn.has_value()) {
      m_stripCenterDistanceColumn.emplace();
      m_extraColumns.push_back(&*m_stripCenterDistanceColumn);
    }
    if ((columns & KnownExtraColumn::TopStripCenter) != 0 &&
        !m_topStripCenterColumn.has_value()) {
      m_topStripCenterColumn.emplace();
      m_extraColumns.push_back(&*m_topStripCenterColumn);
    }
  }

  /// Checks if the container has the given extra columns.
  /// @param columns The extra columns to check for.
  /// @return True if the container has all the specified extra columns, false
  ///         otherwise.
  bool hasExtraColumns(std::uint32_t columns) const {
    if ((columns & KnownExtraColumn::R) != 0 && !m_rColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::Phi) != 0 && !m_phiColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::Time) != 0 && !m_timeColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::VarianceZ) != 0 &&
        !m_varianceZColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::VarianceR) != 0 &&
        !m_varianceRColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::TopStripVector) != 0 &&
        !m_topStripVectorColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::BottomStripVector) != 0 &&
        !m_bottomStripVectorColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::StripCenterDistance) != 0 &&
        !m_stripCenterDistanceColumn.has_value()) {
      return false;
    }
    if ((columns & KnownExtraColumn::TopStripCenter) != 0 &&
        !m_topStripCenterColumn.has_value()) {
      return false;
    }
    return true;
  }

  ExtraColumnProxy<float> rColumn() {
    if (!m_rColumn.has_value()) {
      throw std::runtime_error("Extra column 'r' does not exist");
    }
    return m_rColumn->proxy();
  }
  ExtraColumnProxy<float> phiColumn() {
    if (!m_phiColumn.has_value()) {
      throw std::runtime_error("Extra column 'phi' does not exist");
    }
    return m_phiColumn->proxy();
  }
  ExtraColumnProxy<std::optional<float>> timeColumn() {
    if (!m_timeColumn.has_value()) {
      throw std::runtime_error("Extra column 'time' does not exist");
    }
    return m_timeColumn->proxy();
  }
  ExtraColumnProxy<float> varianceZColumn() {
    if (!m_varianceZColumn.has_value()) {
      throw std::runtime_error("Extra column 'varianceZ' does not exist");
    }
    return m_varianceZColumn->proxy();
  }
  ExtraColumnProxy<float> varianceRColumn() {
    if (!m_varianceRColumn.has_value()) {
      throw std::runtime_error("Extra column 'varianceR' does not exist");
    }
    return m_varianceRColumn->proxy();
  }
  ExtraColumnProxy<Eigen::Vector3f> topStripVectorColumn() {
    if (!m_topStripVectorColumn.has_value()) {
      throw std::runtime_error("Extra column 'topStripVector' does not exist");
    }
    return m_topStripVectorColumn->proxy();
  }
  ExtraColumnProxy<Eigen::Vector3f> bottomStripVectorColumn() {
    if (!m_bottomStripVectorColumn.has_value()) {
      throw std::runtime_error(
          "Extra column 'bottomStripVector' does not exist");
    }
    return m_bottomStripVectorColumn->proxy();
  }
  ExtraColumnProxy<Eigen::Vector3f> stripCenterDistanceColumn() {
    if (!m_stripCenterDistanceColumn.has_value()) {
      throw std::runtime_error(
          "Extra column 'stripCenterDistance' does not exist");
    }
    return m_stripCenterDistanceColumn->proxy();
  }
  ExtraColumnProxy<Eigen::Vector3f> topStripCenterColumn() {
    if (!m_topStripCenterColumn.has_value()) {
      throw std::runtime_error("Extra column 'topStripCenter' does not exist");
    }
    return m_topStripCenterColumn->proxy();
  }

  /// Creates a new column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created column.
  /// @throws std::runtime_error if a column with the same name already exists.
  template <typename T>
  ExtraColumnProxy<T> createExtraColumn(const std::string &name) {
    return createExtraColumnImpl<ExtraColumnHolder<T>>(name);
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
  ExtraColumnProxy<T> extraColumn(const std::string &name) const {
    return extraColumnImpl<ExtraColumnHolder<T>>(name);
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

    bool operator==(const Iterator &other) const {
      return m_index == other.m_index && m_container == other.m_container;
    }
    bool operator!=(const Iterator &other) const { return !(*this == other); }

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
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
  class ColumnHolderBase {
   public:
    virtual ~ColumnHolderBase() = default;

    virtual std::unique_ptr<ColumnHolderBase> copy() const = 0;

    virtual void reserve(std::size_t size) = 0;
    virtual void resize(std::size_t size) = 0;
    virtual void clear() = 0;
    virtual void emplace_back() = 0;
  };
  template <typename T>
  class ExtraColumnHolder final : public ColumnHolderBase {
   public:
    using ValueType = T;
    using ContainerType = std::vector<ValueType>;
    using ProxyType = ExtraColumnProxy<ValueType>;

    std::unique_ptr<ColumnHolderBase> copy() const final {
      return std::make_unique<ExtraColumnHolder<T>>(*this);
    }
    ProxyType proxy() const { return ProxyType(m_data); }

    void reserve(std::size_t size) final { m_data.reserve(size); }
    void clear() final { m_data.clear(); }
    void resize(std::size_t size) final { m_data.resize(size); }
    void emplace_back() final { m_data.emplace_back(); }

   private:
    ContainerType m_data;
  };

  std::vector<float> m_xyz;
  std::vector<std::size_t> m_sourceLinkOffsets;
  std::vector<std::uint8_t> m_sourceLinkCounts;
  std::vector<SourceLink> m_sourceLinks;

  // cylindrical coordinates
  std::optional<ExtraColumnHolder<float>> m_rColumn;
  std::optional<ExtraColumnHolder<float>> m_phiColumn;
  // time information
  std::optional<ExtraColumnHolder<std::optional<float>>> m_timeColumn;
  // covariance information
  std::optional<ExtraColumnHolder<float>> m_varianceZColumn;
  std::optional<ExtraColumnHolder<float>> m_varianceRColumn;
  // strip information
  std::optional<ExtraColumnHolder<Eigen::Vector3f>> m_topStripVectorColumn;
  std::optional<ExtraColumnHolder<Eigen::Vector3f>> m_bottomStripVectorColumn;
  std::optional<ExtraColumnHolder<Eigen::Vector3f>> m_stripCenterDistanceColumn;
  std::optional<ExtraColumnHolder<Eigen::Vector3f>> m_topStripCenterColumn;

  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_namedExtraColumns;

  std::vector<ColumnHolderBase *> m_extraColumns;

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
  SpacePointProxy2(ContainerType &container, IndexType index)
      : m_container(&container), m_index(index) {}

  /// Copy construct a space point proxy.
  /// @param other The space point proxy to copy.
  SpacePointProxy2(const SpacePointProxy2 &other) = default;

  /// Copy construct a mutable space point proxy.
  /// @param other The mutable space point proxy to copy.
  SpacePointProxy2(const SpacePointProxy2<false> &other)
    requires(ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  /// Gets the container holding the space point.
  /// @return A reference to the container holding the space point.
  SpacePointContainer2 &container() { return *m_container; }
  /// Gets the container holding the space point.
  /// @return A const reference to the container holding the space point.
  const SpacePointContainer2 &container() const { return *m_container; }
  /// Gets the index of the space point in the container.
  /// @return The index of the space point in the container.
  IndexType index() const { return m_index; }

  /// Mutable access to the source links of the space point.
  /// @return A mutable span of source links associated with the space point.
  std::span<SourceLink> sourceLinks()
    requires(!ReadOnly)
  {
    return m_container->sourceLinks(m_index);
  }
  /// Mutable access to the x coordinate of the space point.
  /// @return A mutable reference to the x coordinate of the space point.
  float &x()
    requires(!ReadOnly)
  {
    return m_container->x(m_index);
  }
  /// Mutable access to the y coordinate of the space point.
  /// @return A mutable reference to the y coordinate of the space point.
  float &y()
    requires(!ReadOnly)
  {
    return m_container->y(m_index);
  }
  /// Mutable access to the z coordinate of the space point.
  /// @return A mutable reference to the z coordinate of the space point.
  float &z()
    requires(!ReadOnly)
  {
    return m_container->z(m_index);
  }

  float &r()
    requires(!ReadOnly)
  {
    return m_container->r(m_index);
  }
  float &phi()
    requires(!ReadOnly)
  {
    return m_container->phi(m_index);
  }
  std::optional<float> &time()
    requires(!ReadOnly)
  {
    return m_container->time(m_index);
  }
  float &varianceZ()
    requires(!ReadOnly)
  {
    return m_container->varianceZ(m_index);
  }
  float &varianceR()
    requires(!ReadOnly)
  {
    return m_container->varianceR(m_index);
  }
  Eigen::Vector3f &topStripVector()
    requires(!ReadOnly)
  {
    return m_container->topStripVector(m_index);
  }
  Eigen::Vector3f &bottomStripVector()
    requires(!ReadOnly)
  {
    return m_container->bottomStripVector(m_index);
  }
  Eigen::Vector3f &stripCenterDistance()
    requires(!ReadOnly)
  {
    return m_container->stripCenterDistance(m_index);
  }
  Eigen::Vector3f &topStripCenter()
    requires(!ReadOnly)
  {
    return m_container->topStripCenter(m_index);
  }

  /// Mutable access to the extra column of data for the space point.
  /// @param column The extra column to access.
  /// @return A mutable reference to the value in the extra column for the space
  ///         point.
  template <typename column_proxy>
  typename column_proxy::ValueType &extra(column_proxy column)
    requires(!ReadOnly)
  {
    return m_container->extra(column, m_index);
  }

  /// Const access to the x coordinate of the space point.
  /// @return The x coordinate of the space point.
  float x() const { return m_container->x(m_index); }
  /// Const access to the y coordinate of the space point.
  /// @return The y coordinate of the space point.
  float y() const { return m_container->y(m_index); }
  /// Const access to the z coordinate of the space point.
  /// @return The z coordinate of the space point.
  float z() const { return m_container->z(m_index); }

  float r() const { return m_container->r(m_index); }
  float phi() const { return m_container->phi(m_index); }
  std::optional<float> time() const { return m_container->time(m_index); }
  float varianceZ() const { return m_container->varianceZ(m_index); }
  float varianceR() const { return m_container->varianceR(m_index); }
  const Eigen::Vector3f &topStripVector() const {
    return m_container->topStripVector(m_index);
  }
  const Eigen::Vector3f &bottomStripVector() const {
    return m_container->bottomStripVector(m_index);
  }
  const Eigen::Vector3f &stripCenterDistance() const {
    return m_container->stripCenterDistance(m_index);
  }
  const Eigen::Vector3f &topStripCenter() const {
    return m_container->topStripCenter(m_index);
  }

  /// Const access to the extra column of data for the space point.
  /// @param column The extra column to access.
  /// @return A const reference to the value in the extra column for the space
  ///         point.
  template <typename column_proxy>
  const typename column_proxy::ValueType &extra(
      const column_proxy &column) const {
    return m_container->extra(column, m_index);
  }

 private:
  ContainerType *m_container{};
  IndexType m_index{};
};

inline MutableSpacePointProxy2 SpacePointContainer2::createSpacePoint(
    std::span<const SourceLink> sourceLinks, float x, float y, float z) {
  m_xyz.push_back(x);
  m_xyz.push_back(y);
  m_xyz.push_back(z);
  m_sourceLinkOffsets.push_back(m_sourceLinks.size());
  m_sourceLinkCounts.push_back(static_cast<std::uint8_t>(sourceLinks.size()));
  m_sourceLinks.insert(m_sourceLinks.end(), sourceLinks.begin(),
                       sourceLinks.end());

  for (auto &column : m_extraColumns) {
    column->emplace_back();
  }

  return MutableProxyType(*this, size() - 1);
}

inline MutableSpacePointProxy2 SpacePointContainer2::at(IndexType index) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  return MutableProxyType(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::at(IndexType index) const {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  return ConstProxyType(*this, index);
}

}  // namespace Acts::Experimental
