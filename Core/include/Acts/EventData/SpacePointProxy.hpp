// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Definitions/Algebra.hpp"

#include <any>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace Acts {

template <typename container_t, bool read_only>
class SpacePointProxy {
 public:
  using ContainerType = container_t;
  using ValueType = typename ContainerType::ValueType;

 public:
  // Never take the ownership of the container
  SpacePointProxy(ContainerType&& container, std::size_t index) = delete;
  // Only get the reference
  SpacePointProxy(ContainerType& container, std::size_t index);
  SpacePointProxy(const ContainerType& container, std::size_t index);
  // copy and move operations are defaults

  ValueType& sp() { return container().sp(m_index); }
  const ValueType& sp() const { return container().sp(m_index); }

  std::size_t index() const;
  float x() const;
  float y() const;
  float z() const;
  float radius() const;
  float varianceR() const;
  float varianceZ() const;

  float topHalfStripLength() const;
  float bottomHalfStripLength() const;
  Acts::Vector3 topStripDirection() const;
  Acts::Vector3 bottomStripDirection() const;
  Acts::Vector3 stripCenterDistance() const;
  Acts::Vector3 topStripCenterPosition() const;

  // component methods for additional quantities
  template <typename T>
  const T& component(HashedString key) const {
    return container().template component<T>(key, m_index);
  }

  template <typename T, bool RO = read_only, std::enable_if_t<!RO, bool> = true>
  T& component(HashedString key) {
    return container().template component<T>(key, m_index);
  }

 private:
  template <bool RO = read_only, std::enable_if_t<!RO, bool> = true>
  ContainerType& container() {
    return *m_container;
  }

  const ContainerType& container() const;

 private:
  Acts::detail::RefHolder<ContainerType> m_container;
  std::size_t m_index;
};

// Implementation
template <typename container_t, bool read_only>
SpacePointProxy<container_t, read_only>::SpacePointProxy(
    typename SpacePointProxy<container_t, read_only>::ContainerType& container,
    std::size_t index)
    : m_container(container), m_index(index) {}

template <typename container_t, bool read_only>
SpacePointProxy<container_t, read_only>::SpacePointProxy(
    const typename SpacePointProxy<container_t, read_only>::ContainerType&
        container,
    std::size_t index)
    : m_container(container), m_index(index) {}

template <typename container_t, bool read_only>
inline std::size_t SpacePointProxy<container_t, read_only>::index() const {
  return m_index;
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::x() const {
  return container().x(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::y() const {
  return container().y(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::z() const {
  return container().z(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::radius() const {
  return container().radius(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::varianceR() const {
  return container().varianceR(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::varianceZ() const {
  return container().varianceZ(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::topHalfStripLength()
    const {
  return container().topHalfStripLength(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::bottomHalfStripLength()
    const {
  return container().bottomHalfStripLength(m_index);
}

template <typename container_t, bool read_only>
inline Acts::Vector3
SpacePointProxy<container_t, read_only>::topStripDirection() const {
  return container().topStripDirection(m_index);
}

template <typename container_t, bool read_only>
inline Acts::Vector3
SpacePointProxy<container_t, read_only>::bottomStripDirection() const {
  return container().bottomStripDirection(m_index);
}

template <typename container_t, bool read_only>
inline Acts::Vector3
SpacePointProxy<container_t, read_only>::stripCenterDistance() const {
  return container().stripCenterDistance(m_index);
}

template <typename container_t, bool read_only>
inline Acts::Vector3
SpacePointProxy<container_t, read_only>::topStripCenterPosition() const {
  return container().topStripCenterPosition(m_index);
}

template <typename container_t, bool read_only>
inline const typename SpacePointProxy<container_t, read_only>::ContainerType&
SpacePointProxy<container_t, read_only>::container() const {
  return *m_container;
}

}  // namespace Acts
