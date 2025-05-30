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
#include <string>
#include <unordered_map>
#include <utility>

namespace Acts {

using SpacePointIndex2 = std::size_t;
using SpacePointIndexRange2 = std::pair<SpacePointIndex2, SpacePointIndex2>;

class SpacePointContainer2;

template <bool read_only>
class SpacePointProxy2;
using MutableSpacePointProxy2 = SpacePointProxy2<false>;
using ConstSpacePointProxy2 = SpacePointProxy2<true>;

template <typename T>
class SpacePointColumn2 {
 public:
  using IndexType = SpacePointIndex2;
  using ValueType = T;
  using ContainerType = std::vector<ValueType>;

  T &at(IndexType index) { return m_data.at(index); }

  const T &at(IndexType index) const { return m_data.at(index); }

 private:
  ContainerType m_data;

  friend class SpacePointContainer2;
};

class SpacePointContainer2 {
 public:
  using IndexType = SpacePointIndex2;
  using RangeType = SpacePointIndexRange2;
  using MutableProxyType = MutableSpacePointProxy2;
  using ConstProxyType = ConstSpacePointProxy2;

  SpacePointContainer2() = default;

  SpacePointContainer2(const SpacePointContainer2 &other)
      : m_entries(other.m_entries),
        m_xyz(other.m_xyz),
        m_sourceLinks(other.m_sourceLinks) {
    for (const auto &column : other.m_extraColumns) {
      m_extraColumns[column.first] = column.second->copy();
    }
  }

  SpacePointContainer2(SpacePointContainer2 &&other) noexcept = default;

  ~SpacePointContainer2() = default;

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

  SpacePointContainer2 &operator=(SpacePointContainer2 &&other) noexcept =
      default;

  std::size_t size() const { return m_entries.size(); }
  bool empty() const { return size() == 0; }

  void reserve(std::size_t size) {
    m_entries.reserve(size);
    m_xyz.reserve(size * 3);
    m_sourceLinks.reserve(size);

    for (auto &column : m_extraColumns) {
      column.second->reserve(size);
    }
  }
  void clear() {
    m_entries.clear();
    m_xyz.clear();
    m_sourceLinks.clear();

    for (auto &column : m_extraColumns) {
      column.second->clear();
    }
  }

  MutableSpacePointProxy2 createSpacePoint(SourceLink sourceLink, float x,
                                           float y, float z);

  MutableProxyType at(IndexType index);

  ConstProxyType at(IndexType index) const;

  SourceLink &sourceLink(IndexType index, std::size_t sourceLinkIndex) {
    return m_sourceLinks[m_entries[index].sourceLinkOffset + sourceLinkIndex];
  }
  float &x(IndexType index) { return m_xyz[index * 3]; }
  float &y(IndexType index) { return m_xyz[index * 3 + 1]; }
  float &z(IndexType index) { return m_xyz[index * 3 + 2]; }

  std::size_t sourceLinkCount(IndexType index) const {
    return m_entries[index].sourceLinkCount;
  }
  const SourceLink &sourceLink(IndexType index,
                               std::size_t sourceLinkIndex) const {
    return m_sourceLinks[m_entries[index].sourceLinkOffset + sourceLinkIndex];
  }
  float x(IndexType index) const { return m_xyz[index * 3]; }
  float y(IndexType index) const { return m_xyz[index * 3 + 1]; }
  float z(IndexType index) const { return m_xyz[index * 3 + 2]; }

  template <typename T>
  SpacePointColumn2<T> &createExtraColumn(const std::string &name) {
    auto it = m_extraColumns.find(name);
    if (it != m_extraColumns.end()) {
      throw std::runtime_error("Extra column already exists: " + name);
    }
    auto holder = std::make_unique<ColumnHolder<T>>();
    holder->resize(size());
    auto &result = holder->column;
    m_extraColumns[name] = std::move(holder);
    return result;
  }

  template <typename T>
  SpacePointColumn2<T> &extraColumn(const std::string &name) {
    auto it = m_extraColumns.find(name);
    if (it == m_extraColumns.end()) {
      throw std::runtime_error("Extra column not found: " + name);
    }
    auto holder = dynamic_cast<ColumnHolder<T> &>(*it->second);
    return holder.column;
  }
  template <typename T>
  const SpacePointColumn2<T> &extraColumn(const std::string &name) const {
    auto it = m_extraColumns.find(name);
    if (it == m_extraColumns.end()) {
      throw std::runtime_error("Extra column not found: " + name);
    }
    auto holder = dynamic_cast<const ColumnHolder<T> &>(*it->second);
    return holder.column;
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

    Range(ContainerType &container, const RangeType &range)
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
    RangeType m_range{};
  };
  using MutableRange = Range<false>;
  using ConstRange = Range<true>;

  MutableRange range(const RangeType &range) {
    return MutableRange(*this, range);
  }
  ConstRange range(const RangeType &range) const {
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
  struct ColumnHolder final : public ColumnHolderBase {
    SpacePointColumn2<T> column;

    std::unique_ptr<ColumnHolderBase> copy() const final {
      return std::make_unique<ColumnHolder<T>>(*this);
    }

    void reserve(std::size_t size) final { column.m_data.reserve(size); }
    void clear() final { column.m_data.clear(); }
    void resize(std::size_t size) final { column.m_data.resize(size); }
    void emplace_back() final { column.m_data.emplace_back(); }
  };

  std::vector<Entry> m_entries;
  std::vector<float> m_xyz;
  std::vector<SourceLink> m_sourceLinks;

  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_extraColumns;
};

template <bool read_only>
class SpacePointProxy2 {
 public:
  /// Indicates whether this spacepoint proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SpacePointIndex2;
  using RangeType = SpacePointIndexRange2;

  using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

  SpacePointProxy2(ContainerType &container, IndexType index)
      : m_container(&container), m_index(index) {}

  SpacePointProxy2(const SpacePointProxy2 &other) = default;

  SpacePointProxy2(const SpacePointProxy2<false> &other)
    requires(ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  SpacePointContainer2 &container() { return *m_container; }

  const SpacePointContainer2 &container() const { return *m_container; }
  SpacePointIndex2 index() const { return m_index; }

  float &x()
    requires(!ReadOnly)
  {
    return m_container->x(m_index);
  }
  float &y()
    requires(!ReadOnly)
  {
    return m_container->y(m_index);
  }
  float &z()
    requires(!ReadOnly)
  {
    return m_container->z(m_index);
  }

  template <typename column_proxy_type>
  column_proxy_type::ValueType &extra(column_proxy_type &column)
    requires(!ReadOnly)
  {
    return column.at(m_index);
  }

  float x() const { return m_container->x(m_index); }
  float y() const { return m_container->y(m_index); }
  float z() const { return m_container->z(m_index); }

  template <typename column_proxy_type>
  const column_proxy_type::ValueType &extra(
      const column_proxy_type &column) const {
    return column.at(m_index);
  }

  template <bool read_only_it>
  class SourceLinkIterator {
   public:
    static constexpr bool ReadOnly = read_only_it;
    using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = const_if_t<ReadOnly, SourceLink>;
    using reference = value_type &;
    using pointer = value_type *;
    using difference_type = std::ptrdiff_t;

    SourceLinkIterator() = default;
    SourceLinkIterator(ContainerType &container, IndexType index,
                       std::size_t sourceLinkIndex)
        : m_container(&container),
          m_index(index),
          m_sourceLinkIndex(sourceLinkIndex) {}

    SourceLinkIterator &operator++() {
      ++m_sourceLinkIndex;
      return *this;
    }
    SourceLinkIterator operator++(int) {
      SourceLinkIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    bool operator==(const SourceLinkIterator &other) const {
      return m_index == other.m_index &&
             m_sourceLinkIndex == other.m_sourceLinkIndex &&
             m_container == other.m_container;
    }
    bool operator!=(const SourceLinkIterator &other) const {
      return !(*this == other);
    }

    reference operator*() const {
      return m_container->sourceLink(m_index, m_sourceLinkIndex);
    }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
    std::size_t m_sourceLinkIndex{};
  };

  template <bool read_only_range>
  class SourceLinkRange {
   public:
    static constexpr bool ReadOnly = read_only_range;
    using ContainerType = const_if_t<ReadOnly, SpacePointContainer2>;

    using iterator = SourceLinkIterator<read_only_range>;

    SourceLinkRange(ContainerType &container, IndexType index)
        : m_container(&container), m_index(index) {}

    std::size_t size() const { return m_container->sourceLinkCount(m_index); }
    bool empty() const { return size() == 0; }

    SourceLink &operator[](std::size_t index)
      requires(!ReadOnly)
    {
      return m_container->sourceLink(m_index, index);
    }
    const SourceLink &operator[](std::size_t index) const {
      return m_container->sourceLink(m_index, index);
    }

    iterator begin() const { return iterator(*m_container, m_index, 0); }
    iterator end() const {
      return iterator(*m_container, m_index,
                      m_container->sourceLinkCount(m_index));
    }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
  };

  SourceLinkRange<false> sourceLinks()
    requires(!ReadOnly)
  {
    return SourceLinkRange<false>(*m_container, m_index);
  }
  SourceLinkRange<true> sourceLinks() const {
    return SourceLinkRange<true>(*m_container, m_index);
  }

 private:
  ContainerType *m_container{};
  IndexType m_index{};
};

inline MutableSpacePointProxy2 SpacePointContainer2::createSpacePoint(
    SourceLink sourceLink, float x, float y, float z) {
  m_entries.emplace_back<std::size_t, std::size_t>(m_sourceLinks.size(), 1);
  m_xyz.push_back(x);
  m_xyz.push_back(y);
  m_xyz.push_back(z);
  m_sourceLinks.emplace_back(std::move(sourceLink));

  for (auto &column : m_extraColumns) {
    column.second->emplace_back();
  }

  return MutableSpacePointProxy2(*this, size() - 1);
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
