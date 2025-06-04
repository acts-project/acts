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

#include <iterator>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>

namespace Acts {

using SpacePointIndex2 = std::size_t;
using SpacePointIndexRange2 = std::pair<SpacePointIndex2, SpacePointIndex2>;

template <bool read_only>
class SpacePointProxy2;
using MutableSpacePointProxy2 = SpacePointProxy2<false>;
using ConstSpacePointProxy2 = SpacePointProxy2<true>;

/// A container for space points, which can hold additional columns of data
/// (both dense and sparse) and allows for efficient access to space points
/// and their associated source links. Individual space points are addressed
/// via index. A proxy object simplifies the handling.
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
  SpacePointContainer2(const SpacePointContainer2 &other)
      : m_entries(other.m_entries),
        m_xyz(other.m_xyz),
        m_sourceLinks(other.m_sourceLinks) {
    for (const auto &column : other.m_extraColumns) {
      m_extraColumns[column.first] = column.second->copy();
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

    m_entries = other.m_entries;
    m_xyz = other.m_xyz;
    m_sourceLinks = other.m_sourceLinks;

    m_extraColumns.clear();
    for (const auto &column : other.m_extraColumns) {
      m_extraColumns[column.first] = column.second->copy();
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
  std::size_t size() const { return m_entries.size(); }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const { return size() == 0; }

  /// Reserves space for the given number of space points.
  /// This will reserve space for the source links and the extra columns as
  /// well.
  /// @param size The number of space points to reserve space for.
  void reserve(std::size_t size) {
    m_entries.reserve(size);
    m_xyz.reserve(size * 3);
    m_sourceLinks.reserve(size);

    for (auto &column : m_extraColumns) {
      column.second->reserve(size);
    }
  }
  /// Clears the container, removing all space points and extra columns.
  void clear() {
    m_entries.clear();
    m_xyz.clear();
    m_sourceLinks.clear();

    for (auto &column : m_extraColumns) {
      column.second->clear();
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
    return std::span<SourceLink>(
        m_sourceLinks.data() + m_entries[index].sourceLinkOffset,
        m_entries[index].sourceLinkCount);
  }
  /// Mutable access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the x coordinate of the space point.
  float &x(IndexType index) { return m_xyz[index * 3]; }
  /// Mutable access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the y coordinate of the space point.
  float &y(IndexType index) { return m_xyz[index * 3 + 1]; }
  /// Mutable access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A mutable reference to the z coordinate of the space point.
  float &z(IndexType index) { return m_xyz[index * 3 + 2]; }

  /// Const access to the source links at the given index.
  /// @param index The index of the space point.
  /// @return A const span to the source links at the given index.
  std::span<const SourceLink> sourceLinks(IndexType index) const {
    return std::span<const SourceLink>(
        m_sourceLinks.data() + m_entries[index].sourceLinkOffset,
        m_entries[index].sourceLinkCount);
  }
  /// Const access to the x coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the x coordinate of the space point.
  float x(IndexType index) const { return m_xyz[index * 3]; }
  /// Const access to the y coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the y coordinate of the space point.
  float y(IndexType index) const { return m_xyz[index * 3 + 1]; }
  /// Const access to the z coordinate of the space point at the given index.
  /// @param index The index of the space point.
  /// @return A const reference to the z coordinate of the space point.
  float z(IndexType index) const { return m_xyz[index * 3 + 2]; }

  /// Additional dense column of data that can be added to the space point
  /// container. The column is indexed by the space point index.
  template <typename T>
  class DenseColumn {
   public:
    using ValueType = T;
    using ContainerType = std::vector<ValueType>;

    /// Mutable access to the column data.
    /// @param index The index of the space point.
    /// @return A mutable reference to the value at the given index.
    T &at(IndexType index) { return m_data[index]; }
    /// Const access to the column data.
    /// @param index The index of the space point.
    /// @return A const reference to the value at the given index.
    const T &at(IndexType index) const { return m_data[index]; }

    auto begin() { return m_data.begin(); }
    auto end() { return m_data.end(); }
    auto begin() const { return m_data.begin(); }
    auto end() const { return m_data.end(); }
    auto cbegin() const { return m_data.cbegin(); }
    auto cend() const { return m_data.cend(); }

   private:
    ContainerType m_data;

    friend class SpacePointContainer2;
  };

  /// Creates a new dense column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created dense column.
  /// @throws std::runtime_error if a column with the same name already exists.
  template <typename T>
  DenseColumn<T> &createDenseExtraColumn(const std::string &name) {
    return createExtraColumn<DenseColumnHolder<T>>(name);
  }

  /// Returns a mutable reference to the dense extra column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the dense extra column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  DenseColumn<T> &denseExtraColumn(const std::string &name) {
    return extraColumn<DenseColumnHolder<T>>(name);
  }
  /// Returns a const reference to the dense extra column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A const reference to the dense extra column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  const DenseColumn<T> &denseExtraColumn(const std::string &name) const {
    return extraColumn<DenseColumnHolder<T>>(name);
  }

  /// Additional sparse column of data that can be added to the space point
  /// container. The column is indexed by the space point index, but only
  /// contains values for space points that have a value set for that index.
  template <typename T>
  class SparseColumn {
   public:
    using ValueType = T;
    using ContainerType = std::unordered_map<IndexType, ValueType>;

    /// Returns the number of elements in the sparse column.
    /// @return The number of elements in the sparse column.
    std::size_t size() const { return m_data.size(); }
    /// Checks if the sparse column is empty.
    /// @return True if the sparse column is empty, false otherwise.
    [[nodiscard]]
    bool empty() const {
      return m_data.empty();
    }
    /// Reserves space for the given number of elements in the sparse column.
    /// @param size The number of elements to reserve space for.
    void reserve(std::size_t size) { m_data.reserve(size); }
    /// Clears the sparse column, removing all elements.
    void clear() { m_data.clear(); }

    /// Checks if the sparse column has a value for the given index.
    /// @param index The index to check.
    /// @return True if the sparse column has a value for the given index,
    ///         false otherwise.
    [[nodiscard]]
    bool has(IndexType index) const {
      return m_data.find(index) != m_data.end();
    }

    /// Mutable access to the value at the given index.
    /// @param index The index of the value to access.
    /// @return A mutable reference to the value at the given index.
    T &at(IndexType index) { return m_data[index]; }
    /// Const access to the value at the given index.
    /// @param index The index of the value to access.
    /// @return A const reference to the value at the given index.
    const T &at(IndexType index) const { return m_data[index]; }

    /// Inserts a value at the given index. If the index already exists, the
    /// value is overwritten.
    /// @param index The index to insert the value at.
    /// @param value The value to insert.
    void insert(IndexType index, const T &value) { m_data[index] = value; }
    /// Inserts a value at the given index, moving the value into the column.
    /// If the index already exists, the value is overwritten.
    /// @param index The index to insert the value at.
    /// @param value The value to insert, moved into the column.
    void insert(IndexType index, T &&value) {
      m_data[index] = std::move(value);
    }

    /// Erases the value at the given index from the sparse column.
    /// @param index The index of the value to erase.
    void erase(IndexType index) { m_data.erase(index); }

    auto begin() { return m_data.begin(); }
    auto end() { return m_data.end(); }
    auto begin() const { return m_data.begin(); }
    auto end() const { return m_data.end(); }
    auto cbegin() const { return m_data.cbegin(); }
    auto cend() const { return m_data.cend(); }

    /// Finds the value at the given index in the sparse column.
    /// @param index The index of the value to find.
    /// @return An iterator to the found value, or end() if not found.
    auto find(IndexType index) { return m_data.find(index); }
    /// Finds the value at the given index in the sparse column.
    /// @param index The index of the value to find.
    /// @return A const iterator to the found value, or end() if not found.
    auto find(IndexType index) const { return m_data.find(index); }

   private:
    ContainerType m_data;

    friend class SpacePointContainer2;
  };

  /// Creates a new sparse column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created sparse column.
  /// @throws std::runtime_error if a column with the same name already exists.
  template <typename T>
  SparseColumn<T> &createSparseExtraColumn(const std::string &name) {
    return createExtraColumn<SparseColumnHolder<T>>(name);
  }

  /// Returns a mutable reference to the sparse extra column with the given
  /// name. If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the sparse extra column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  SparseColumn<T> &sparseExtraColumn(const std::string &name) {
    return extraColumn<SparseColumnHolder<T>>(name);
  }
  /// Returns a const reference to the sparse extra column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A const reference to the sparse extra column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  const SparseColumn<T> &sparseExtraColumn(const std::string &name) const {
    return extraColumn<SparseColumnHolder<T>>(name);
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
  struct Entry {
    std::size_t sourceLinkOffset{};
    std::size_t sourceLinkCount{};

    Entry() = default;
    Entry(std::size_t sourceLinkOffset_, std::size_t sourceLinkCount_)
        : sourceLinkOffset(sourceLinkOffset_),
          sourceLinkCount(sourceLinkCount_) {}
  };

  struct ColumnHolderBase {
    virtual ~ColumnHolderBase() = default;

    virtual std::unique_ptr<ColumnHolderBase> copy() const = 0;

    virtual void reserve(std::size_t size) = 0;
    virtual void resize(std::size_t size) = 0;
    virtual void clear() = 0;
    virtual void emplace_back() = 0;
  };
  template <typename T>
  struct DenseColumnHolder final : public ColumnHolderBase {
    DenseColumn<T> column;

    std::unique_ptr<ColumnHolderBase> copy() const final {
      return std::make_unique<DenseColumnHolder<T>>(*this);
    }

    void reserve(std::size_t size) final { column.m_data.reserve(size); }
    void clear() final { column.m_data.clear(); }
    void resize(std::size_t size) final { column.m_data.resize(size); }
    void emplace_back() final { column.m_data.emplace_back(); }
  };
  template <typename T>
  struct SparseColumnHolder final : public ColumnHolderBase {
    SparseColumn<T> column;

    std::unique_ptr<ColumnHolderBase> copy() const final {
      return std::make_unique<SparseColumnHolder<T>>(*this);
    }

    void reserve(std::size_t size) final { column.reserve(size); }
    void clear() final { column.clear(); }
    void resize([[maybe_unused]] std::size_t size)
        final { /* No-op for sparse columns */ }
    void emplace_back() final { /* No-op for sparse columns */ }
  };

  std::vector<Entry> m_entries;
  std::vector<float> m_xyz;
  std::vector<SourceLink> m_sourceLinks;

  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_extraColumns;

  template <typename Holder>
  auto &createExtraColumn(const std::string &name) {
    auto it = m_extraColumns.find(name);
    if (it != m_extraColumns.end()) {
      throw std::runtime_error("Extra column already exists: " + name);
    }
    auto holder = std::make_unique<Holder>();
    auto &result = holder->column;
    m_extraColumns[name] = std::move(holder);
    return result;
  }
  template <typename Holder>
  auto &extraColumn(const std::string &name) {
    auto it = m_extraColumns.find(name);
    if (it == m_extraColumns.end()) {
      throw std::runtime_error("Extra column not found: " + name);
    }
    auto holder = dynamic_cast<Holder &>(*it->second);
    return holder.column;
  }
  template <typename Holder>
  auto &extraColumn(const std::string &name) const {
    auto it = m_extraColumns.find(name);
    if (it == m_extraColumns.end()) {
      throw std::runtime_error("Extra column not found: " + name);
    }
    auto holder = dynamic_cast<const Holder &>(*it->second);
    return holder.column;
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

  /// Mutable access to the extra column of data for the space point.
  /// @param column The extra column to access.
  /// @return A mutable reference to the value in the extra column for the space
  ///         point.
  template <typename column_type>
  column_type::ValueType &extra(column_type &column)
    requires(!ReadOnly)
  {
    return column.at(m_index);
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

  /// Const access to the extra column of data for the space point.
  /// @param column The extra column to access.
  /// @return A const reference to the value in the extra column for the space
  ///         point.
  template <typename column_type>
  const column_type::ValueType &extra(const column_type &column) const {
    return column.at(m_index);
  }

 private:
  ContainerType *m_container{};
  IndexType m_index{};
};

inline MutableSpacePointProxy2 SpacePointContainer2::createSpacePoint(
    std::span<const SourceLink> sourceLinks, float x, float y, float z) {
  m_entries.emplace_back<std::size_t, std::size_t>(m_sourceLinks.size(),
                                                   sourceLinks.size());
  m_xyz.push_back(x);
  m_xyz.push_back(y);
  m_xyz.push_back(z);
  m_sourceLinks.insert(m_sourceLinks.end(), sourceLinks.begin(),
                       sourceLinks.end());

  for (auto &column : m_extraColumns) {
    column.second->emplace_back();
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

}  // namespace Acts
