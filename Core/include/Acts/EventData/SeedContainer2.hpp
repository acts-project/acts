// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"

#include <limits>
#include <span>

namespace Acts {

using SeedIndex2 = std::size_t;

template <bool read_only>
class SeedProxy2;

using MutableSeedProxy2 = SeedProxy2<false>;
using ConstSeedProxy2 = SeedProxy2<true>;

class SeedContainer2 {
 public:
  using IndexType = SeedIndex2;
  using MutableProxyType = MutableSeedProxy2;
  using ConstProxyType = ConstSeedProxy2;

  std::size_t size() const { return m_entries.size(); }
  bool empty() const { return size() == 0; }

  void reserve(std::size_t size) {
    m_entries.reserve(size);
    m_spacePoints.reserve(size * 3);
  }
  void clear() {
    m_entries.clear();
    m_spacePoints.clear();
  }

  MutableProxyType createSeed(std::span<const SpacePointIndex2> spacePoints);

  MutableProxyType at(IndexType index);
  ConstProxyType at(IndexType index) const;

  float &quality(IndexType index) { return m_entries[index].quality; }
  float &vertexZ(IndexType index) { return m_entries[index].vertexZ; }

  std::span<const std::size_t> spacePointIndices(IndexType index) const {
    return std::span<const std::size_t>(
        m_spacePoints.data() + m_entries[index].spacePointOffset,
        m_spacePoints.data() + m_entries[index].spacePointOffset +
            m_entries[index].seedSize);
  }
  float quality(IndexType index) const { return m_entries[index].quality; }
  float vertexZ(IndexType index) const { return m_entries[index].vertexZ; }

  template <bool read_only>
  class SeedIterator {
   public:
    static constexpr bool ReadOnly = read_only;

    using ContainerType = const_if_t<ReadOnly, SeedContainer2>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = SeedProxy2<ReadOnly>;
    using difference_type = std::ptrdiff_t;

    SeedIterator() = default;
    SeedIterator(ContainerType &container, IndexType index)
        : m_container(&container), m_index(index) {}

    SeedIterator &operator++() {
      ++m_index;
      return *this;
    }
    SeedIterator operator++(int) {
      SeedIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    bool operator==(const SeedIterator &other) const {
      return m_index == other.m_index && m_container == other.m_container;
    }
    bool operator!=(const SeedIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
  };
  using iterator = SeedIterator<false>;
  using const_iterator = SeedIterator<true>;

  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, size()); }

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

 private:
  struct Entry {
    std::size_t seedSize{};
    std::size_t spacePointOffset{};
    float quality{-std::numeric_limits<float>::infinity()};
    float vertexZ{};
  };

  std::vector<Entry> m_entries{};
  std::vector<SpacePointIndex2> m_spacePoints{};
};

template <bool read_only>
class SeedProxy2 {
 public:
  /// Indicates whether this spacepoint proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SeedIndex2;

  using ContainerType = const_if_t<ReadOnly, SeedContainer2>;

  SeedProxy2(ContainerType &container, IndexType index)
      : m_container{&container}, m_index{index} {}

  SeedProxy2(const SeedProxy2 &other) = default;

  explicit SeedProxy2(const SeedProxy2<false> &other)
    requires(ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  SeedContainer2 &container() { return *m_container; }

  const SeedContainer2 &container() const { return *m_container; }
  IndexType index() const { return m_index; }

  std::size_t size() const {
    return m_container->spacePointIndices(m_index).size();
  }
  bool empty() const { return size() == 0; }

  float &quality()
    requires(!ReadOnly)
  {
    return m_container->quality(m_index);
  }
  float &vertexZ()
    requires(!ReadOnly)
  {
    return m_container->vertexZ(m_index);
  }

  std::span<const std::size_t> spacePointIndices() const {
    return m_container->spacePointIndices(m_index);
  }

  class SpacePointIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ConstSpacePointProxy2;
    using difference_type = std::ptrdiff_t;

    SpacePointIterator() = default;
    SpacePointIterator(const SpacePointContainer2 &spacePointContainer,
                       const std::size_t *indexPointer)
        : m_spacePointContainer{&spacePointContainer},
          m_indexPointer{indexPointer} {}

    SpacePointIterator &operator++() {
      ++m_indexPointer;
      return *this;
    }

    SpacePointIterator operator++(int) {
      SpacePointIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const SpacePointIterator &other) const {
      return m_indexPointer == other.m_indexPointer;
    }
    bool operator!=(const SpacePointIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const {
      return m_spacePointContainer->at(*m_indexPointer);
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    const std::size_t *m_indexPointer{nullptr};
  };

  class SpacePointRange {
   public:
    SpacePointRange(const SpacePointContainer2 &spacePointContainer,
                    std::span<const std::size_t> spacePointIndices)
        : m_spacePointContainer{&spacePointContainer},
          m_spacePointIndices{spacePointIndices} {}

    std::size_t size() const { return m_spacePointIndices.size(); }
    bool empty() const { return size() == 0; }

    ConstSpacePointProxy2 operator[](std::size_t index) const {
      return m_spacePointContainer->at(m_spacePointIndices[index]);
    }

    SpacePointIterator begin() const {
      return SpacePointIterator(*m_spacePointContainer,
                                m_spacePointIndices.data());
    }
    SpacePointIterator end() const {
      return SpacePointIterator(
          *m_spacePointContainer,
          m_spacePointIndices.data() + m_spacePointIndices.size());
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    std::span<const std::size_t> m_spacePointIndices;
  };

  SpacePointRange spacePoints(
      const SpacePointContainer2 &spacePointContainer) const {
    return SpacePointRange(spacePointContainer,
                           m_container->spacePointIndices(m_index));
  }

  float quality() const { return m_container->quality(m_index); }
  float vertexZ() const { return m_container->vertexZ(m_index); }

 private:
  ContainerType *m_container{nullptr};
  IndexType m_index{0};
};

inline MutableSeedProxy2 SeedContainer2::createSeed(
    std::span<const SpacePointIndex2> spacePoints) {
  m_entries.emplace_back(spacePoints.size(), m_spacePoints.size(),
                         -std::numeric_limits<float>::infinity(), 0.f);
  m_spacePoints.insert(m_spacePoints.end(), spacePoints.begin(),
                       spacePoints.end());
  return MutableProxyType(*this, size() - 1);
}

inline MutableSeedProxy2 SeedContainer2::at(IndexType index) {
  return MutableProxyType(*this, index);
}

inline ConstSeedProxy2 SeedContainer2::at(IndexType index) const {
  return ConstProxyType(*this, index);
}

}  // namespace Acts
