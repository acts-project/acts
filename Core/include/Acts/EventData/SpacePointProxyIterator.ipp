// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

// Implementation
template <typename container_t, bool read_only>
SpacePointProxyIterator<container_t, read_only>::SpacePointProxyIterator(
    typename SpacePointProxyIterator<container_t, read_only>::ContainerType&
        container,
    std::size_t index)
    : m_container(container), m_index(index) {}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>&
SpacePointProxyIterator<container_t, read_only>::operator++() {
  ++m_index;
  return *this;
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>&
SpacePointProxyIterator<container_t, read_only>::operator--() {
  --m_index;
  return *this;
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>
SpacePointProxyIterator<container_t, read_only>::operator++(int) {
  SpacePointProxyIterator other(*this);
  ++m_index;
  return other;
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>
SpacePointProxyIterator<container_t, read_only>::operator--(int) {
  SpacePointProxyIterator other(*this);
  --m_index;
  return other;
}

template <typename container_t, bool read_only>
inline bool SpacePointProxyIterator<container_t, read_only>::operator==(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return m_container.ptr == other.m_container.ptr and m_index == other.m_index;
  ;
}

template <typename container_t, bool read_only>
inline bool SpacePointProxyIterator<container_t, read_only>::operator!=(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return not(*this == other);
}

template <typename container_t, bool read_only>
inline bool SpacePointProxyIterator<container_t, read_only>::operator<(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return m_index < other.m_index;
}

template <typename container_t, bool read_only>
inline bool SpacePointProxyIterator<container_t, read_only>::operator>(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return m_index > other.m_index;
}

template <typename container_t, bool read_only>
inline bool SpacePointProxyIterator<container_t, read_only>::operator<=(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return m_index <= other.m_index;
}

template <typename container_t, bool read_only>
inline bool SpacePointProxyIterator<container_t, read_only>::operator>=(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return m_index >= other.m_index;
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>&
SpacePointProxyIterator<container_t, read_only>::operator+=(
    const std::size_t& offset) {
  m_index += offset;
  return *this;
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>&
SpacePointProxyIterator<container_t, read_only>::operator-=(
    const std::size_t& offset) {
  m_index -= offset;
  return *this;
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>
SpacePointProxyIterator<container_t, read_only>::operator+(
    const std::size_t& offset) const {
  return SpacePointProxyIterator(*m_container, m_index + offset);
}

template <typename container_t, bool read_only>
inline SpacePointProxyIterator<container_t, read_only>
SpacePointProxyIterator<container_t, read_only>::operator-(
    const std::size_t& offset) const {
  return SpacePointProxyIterator(*m_container, m_index - offset);
}

template <typename container_t, bool read_only>
inline typename SpacePointProxyIterator<container_t, read_only>::difference_type
SpacePointProxyIterator<container_t, read_only>::operator-(
    const SpacePointProxyIterator<container_t, read_only>& other) const {
  return m_index - other.m_index;
}

template <typename container_t, bool read_only>
template <bool, typename>
inline typename SpacePointProxyIterator<container_t, read_only>::value_type&
SpacePointProxyIterator<container_t, read_only>::operator*() {
  return m_container->proxy(m_index);
}

template <typename container_t, bool read_only>
inline typename SpacePointProxyIterator<container_t, read_only>::value_type&
SpacePointProxyIterator<container_t, read_only>::operator*() const {
  return m_container->proxy(m_index);
}

}  // namespace Acts
