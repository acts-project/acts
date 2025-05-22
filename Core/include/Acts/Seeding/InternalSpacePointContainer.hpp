// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <iterator>
#include <optional>
#include <utility>

namespace Acts {

using SpacePointIndex = std::size_t;

class InternalSpacePointContainer;

template <bool read_only>
class InternalSpacePointProxy {
 public:
  /// Indicates whether this spacepoint proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SpacePointIndex;

  using ContainerType = const_if_t<ReadOnly, InternalSpacePointContainer>;

  InternalSpacePointProxy(ContainerType &container, IndexType index)
      : m_container(&container), m_index(index) {}

  InternalSpacePointProxy(const InternalSpacePointProxy &other) = default;

  InternalSpacePointProxy(const InternalSpacePointProxy<false> &other)
    requires(ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  InternalSpacePointContainer &container() { return *m_container; }

  const InternalSpacePointContainer &container() const { return *m_container; }
  SpacePointIndex index() const { return m_index; }

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
  float &phi()
    requires(!ReadOnly)
  {
    return m_container->phi(m_index);
  }
  float &radius()
    requires(!ReadOnly)
  {
    return m_container->radius(m_index);
  }
  float &varianceR()
    requires(!ReadOnly)
  {
    return m_container->varianceR(m_index);
  }
  float &varianceZ()
    requires(!ReadOnly)
  {
    return m_container->varianceZ(m_index);
  }

  float x() const { return m_container->x(m_index); }
  float y() const { return m_container->y(m_index); }
  float z() const { return m_container->z(m_index); }
  std::optional<float> t() const {
    // TODO
    return std::nullopt;
  }
  float phi() const { return m_container->phi(m_index); }
  float radius() const { return m_container->radius(m_index); }
  float varianceR() const { return m_container->varianceR(m_index); }
  float varianceZ() const { return m_container->varianceZ(m_index); }

  Vector3 topStripVector() const { return Vector3::Zero(); }
  Vector3 topStripCenterPosition() const { return Vector3::Zero(); }
  Vector3 bottomStripVector() const { return Vector3::Zero(); }
  Vector3 stripCenterDistance() const { return Vector3::Zero(); }

  template <bool read_only_it>
  class SourceLinkIterator {
   public:
    static constexpr bool ReadOnly = read_only_it;
    using ContainerType = const_if_t<ReadOnly, InternalSpacePointContainer>;

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
    using ContainerType = const_if_t<ReadOnly, InternalSpacePointContainer>;

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

using MutableInternalSpacePointProxy = InternalSpacePointProxy<false>;
using ConstInternalSpacePointProxy = InternalSpacePointProxy<true>;

class InternalSpacePointContainer {
 public:
  using IndexType = SpacePointIndex;
  using MutableProxyType = MutableInternalSpacePointProxy;
  using ConstProxyType = ConstInternalSpacePointProxy;

  std::size_t size() const { return m_entries.size(); }
  bool empty() const { return size() == 0; }

  void reserve(std::size_t size) {
    m_entries.reserve(size);
    m_xyz.reserve(size * 3);
    m_phiR.reserve(size * 2);
    m_varianceRZ.reserve(size * 2);
    m_sourceLinks.reserve(size);
  }
  void clear() {
    m_entries.clear();
    m_xyz.clear();
    m_phiR.clear();
    m_varianceRZ.clear();
    m_sourceLinks.clear();
  }

  void setBeamPos(const Vector2 &beamPos) { m_beamPos = beamPos; }

  IndexType addSpacePoint(SourceLink sourceLink) {
    return addSpacePoint(std::move(sourceLink), 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
                         0.f);
  }
  IndexType addSpacePoint(SourceLink sourceLink, float x, float y, float z) {
    return addSpacePoint(std::move(sourceLink), x, y, z,
                         VectorHelpers::phi(Eigen::Vector2f{x, y}),
                         VectorHelpers::perp(Eigen::Vector2f{x, y}), 0.f, 0.f);
  }
  IndexType addSpacePoint(SourceLink sourceLink, float x, float y, float z,
                          float varianceR, float varianceZ) {
    return addSpacePoint(std::move(sourceLink), x, y, z,
                         VectorHelpers::phi(Eigen::Vector2f{x, y}),
                         VectorHelpers::perp(Eigen::Vector2f{x, y}), varianceR,
                         varianceZ);
  }
  IndexType addSpacePoint(SourceLink sourceLink, float x, float y, float z,
                          float phi, float radius, float varianceR,
                          float varianceZ) {
    m_entries.emplace_back(m_sourceLinks.size(), static_cast<std::size_t>(1));
    m_xyz.push_back(x);
    m_xyz.push_back(y);
    m_xyz.push_back(z);
    m_phiR.push_back(phi);
    m_phiR.push_back(radius);
    m_varianceRZ.push_back(varianceR);
    m_varianceRZ.push_back(varianceZ);
    m_sourceLinks.emplace_back(std::move(sourceLink));
    return size() - 1;
  }

  MutableProxyType makeSpacePoint(SourceLink sourceLink) {
    return at(addSpacePoint(std::move(sourceLink)));
  }
  MutableProxyType makeSpacePoint(SourceLink sourceLink, float x, float y,
                                  float z) {
    return at(addSpacePoint(std::move(sourceLink), x, y, z));
  }
  MutableProxyType makeSpacePoint(SourceLink sourceLink, float x, float y,
                                  float z, float varianceR, float varianceZ) {
    return at(
        addSpacePoint(std::move(sourceLink), x, y, z, varianceR, varianceZ));
  }
  MutableProxyType makeSpacePoint(SourceLink sourceLink, float x, float y,
                                  float z, float phi, float radius,
                                  float varianceR, float varianceZ) {
    return at(addSpacePoint(std::move(sourceLink), x, y, z, phi, radius,
                            varianceR, varianceZ));
  }

  MutableProxyType at(IndexType index) {
    return MutableProxyType(*this, index);
  }

  SourceLink &sourceLink(IndexType index, std::size_t sourceLinkIndex) {
    return m_sourceLinks[m_entries[index].sourceLinkOffset + sourceLinkIndex];
  }
  float &x(IndexType index) { return m_xyz[index * 3]; }
  float &y(IndexType index) { return m_xyz[index * 3 + 1]; }
  float &z(IndexType index) { return m_xyz[index * 3 + 2]; }
  float &phi(IndexType index) { return m_phiR[index * 2]; }
  float &radius(IndexType index) { return m_phiR[index * 2 + 1]; }
  float &varianceR(IndexType index) { return m_varianceRZ[index * 2]; }
  float &varianceZ(IndexType index) { return m_varianceRZ[index * 2 + 1]; }

  const Vector2 &beamPos() const { return m_beamPos; }

  ConstProxyType at(IndexType index) const {
    if (index >= size()) {
      throw std::out_of_range("Index out of range");
    }
    return ConstProxyType(*this, index);
  }

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
  float phi(IndexType index) const { return m_phiR[index * 2]; }
  float radius(IndexType index) const { return m_phiR[index * 2 + 1]; }
  float varianceR(IndexType index) const { return m_varianceRZ[index * 2]; }
  float varianceZ(IndexType index) const { return m_varianceRZ[index * 2 + 1]; }

  template <bool read_only>
  class SpacePointIterator {
   public:
    static constexpr bool ReadOnly = read_only;

    using ContainerType = const_if_t<ReadOnly, InternalSpacePointContainer>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = InternalSpacePointProxy<ReadOnly>;
    using difference_type = std::ptrdiff_t;

    SpacePointIterator() = default;
    SpacePointIterator(ContainerType &container, IndexType index)
        : m_container(&container), m_index(index) {}

    SpacePointIterator &operator++() {
      ++m_index;
      return *this;
    }
    SpacePointIterator operator++(int) {
      SpacePointIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    bool operator==(const SpacePointIterator &other) const {
      return m_index == other.m_index && m_container == other.m_container;
    }
    bool operator!=(const SpacePointIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
  };
  using iterator = SpacePointIterator<false>;
  using const_iterator = SpacePointIterator<true>;

  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, size()); }

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

 private:
  Vector2 m_beamPos{0, 0};

  struct Entry {
    std::size_t sourceLinkOffset{};
    std::size_t sourceLinkCount{};
  };

  std::vector<Entry> m_entries;
  std::vector<float> m_xyz;
  std::vector<float> m_phiR;
  std::vector<float> m_varianceRZ;
  std::vector<SourceLink> m_sourceLinks;
};

}  // namespace Acts
