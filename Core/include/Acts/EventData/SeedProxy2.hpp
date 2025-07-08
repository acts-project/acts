// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>

namespace Acts::Experimental {

class SeedContainer2;

/// A proxy class for accessing individual seeds.
template <bool read_only>
class SeedProxy2 {
 public:
  /// Indicates whether this seed proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SeedIndex2;

  using ContainerType = const_if_t<ReadOnly, SeedContainer2>;

  /// Constructs a seed proxy for the given container and index.
  /// @param container The container holding the seed.
  /// @param index The index of the seed in the container.
  SeedProxy2(ContainerType &container, IndexType index) noexcept
      : m_container{&container}, m_index{index} {}

  /// Copy construct a seed proxy.
  /// @param other The seed proxy to copy.
  SeedProxy2(const SeedProxy2 &other) noexcept = default;

  /// Copy construct a mutable seed proxy.
  /// @param other The mutable seed proxy to copy.
  explicit SeedProxy2(const SeedProxy2<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Gets the container holding the seed.
  /// @return A reference to the container holding the seed.
  SeedContainer2 &container() noexcept
    requires ReadOnly
  {
    return *m_container;
  }
  /// Gets the container holding the seed.
  /// @return A const reference to the container holding the seed.
  const SeedContainer2 &container() const noexcept { return *m_container; }
  /// Gets the index of the seed in the container.
  /// @return The index of the seed in the container.
  IndexType index() const noexcept { return m_index; }

  /// Returns the size of the seed, i.e., the number of space points
  /// associated with it.
  /// @return The number of space points in the seed.
  [[nodiscard]] std::size_t size() const noexcept {
    return m_container->spacePointIndices(m_index).size();
  }
  /// Checks if the seed is empty, i.e., has no space points associated with it.
  /// @return True if the seed is empty, false otherwise.
  [[nodiscard]]
  bool empty() const noexcept {
    return size() == 0;
  }

  /// Mutable access to the space point indices of the seed.
  /// @return A mutable span of space point indices associated with the seed.
  std::span<SpacePointIndex2> spacePointIndices() noexcept
    requires(!ReadOnly)
  {
    return m_container->spacePointIndices(m_index);
  }
  /// Mutable access to the quality of the seed.
  /// @return A mutable reference to the quality of the seed.
  float &quality() noexcept
    requires(!ReadOnly)
  {
    return m_container->quality(m_index);
  }
  /// Mutable access to the vertex Z coordinate of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed.
  float &vertexZ() noexcept
    requires(!ReadOnly)
  {
    return m_container->vertexZ(m_index);
  }

  /// Const access to the space point indices of the seed.
  /// @return A span of space point indices associated with the seed.
  ///         This span is read-only and cannot be modified.
  std::span<const SpacePointIndex2> spacePointIndices() const noexcept {
    return m_container->spacePointIndices(m_index);
  }
  /// Const access to the quality of the seed.
  /// @return The quality of the seed.
  float quality() const noexcept { return m_container->quality(m_index); }
  /// Const access to the vertex Z coordinate of the seed.
  /// @return The vertex Z coordinate of the seed.
  float vertexZ() const noexcept { return m_container->vertexZ(m_index); }

  class SpacePointIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ConstSpacePointProxy2;
    using difference_type = std::ptrdiff_t;

    SpacePointIterator() = default;
    SpacePointIterator(const SpacePointContainer2 &spacePointContainer,
                       const SpacePointIndex2 *indexPointer) noexcept
        : m_spacePointContainer{&spacePointContainer},
          m_indexPointer{indexPointer} {}

    SpacePointIterator &operator++() noexcept {
      ++m_indexPointer;
      return *this;
    }

    SpacePointIterator operator++(int) noexcept {
      SpacePointIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    value_type operator*() const noexcept {
      return (*m_spacePointContainer)[*m_indexPointer];
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    const SpacePointIndex2 *m_indexPointer{nullptr};

    friend bool operator==(const SpacePointIterator &a,
                           const SpacePointIterator &b) noexcept {
      return a.m_indexPointer == b.m_indexPointer;
    }
  };

  class SpacePointRange {
   public:
    SpacePointRange(
        const SpacePointContainer2 &spacePointContainer,
        std::span<const SpacePointIndex2> spacePointIndices) noexcept
        : m_spacePointContainer{&spacePointContainer},
          m_spacePointIndices{spacePointIndices} {}

    std::size_t size() const noexcept { return m_spacePointIndices.size(); }
    bool empty() const noexcept { return size() == 0; }

    ConstSpacePointProxy2 operator[](std::size_t index) const noexcept {
      return (*m_spacePointContainer)[m_spacePointIndices[index]];
    }

    SpacePointIterator begin() const noexcept {
      return SpacePointIterator(*m_spacePointContainer,
                                m_spacePointIndices.data());
    }
    SpacePointIterator end() const noexcept {
      return SpacePointIterator(
          *m_spacePointContainer,
          m_spacePointIndices.data() + m_spacePointIndices.size());
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    std::span<const SpacePointIndex2> m_spacePointIndices;
  };

  SpacePointRange spacePoints(
      const SpacePointContainer2 &spacePointContainer) const noexcept {
    return SpacePointRange(spacePointContainer,
                           m_container->spacePointIndices(m_index));
  }

 private:
  ContainerType *m_container{nullptr};
  IndexType m_index{0};
};

}  // namespace Acts::Experimental
