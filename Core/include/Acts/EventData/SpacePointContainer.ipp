// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    const Acts::SpacePointContainerOptions& options,
    container_t& container)
  : m_options(options.toInternalUnits()),
    m_container(container)
{}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
								const Acts::SpacePointContainerOptions& options,
    container_t&& container)
  : m_options(options.toInternalUnits()),
  m_container(std::move(container))
{}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    SpacePointContainer<container_t, holder_t>& other)
  : m_options(other.m_options),
    m_container(*other.m_container.ptr)
{}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>&
SpacePointContainer<container_t, holder_t>::operator=(
    SpacePointContainer<container_t, holder_t>& other)
{
  m_options = other.m_options;
  m_container.ptr = other.m_container.ptr;
  return *this;
}

template <typename container_t, template <typename> class holder_t>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    SpacePointContainer<container_t, holder_t>&& other) noexcept
  : m_options(other.m_options),
    m_container(std::exchange(other.m_container.ptr, nullptr))
{}
  
template <typename container_t, template <typename> class holder_t>
template <typename T>
T SpacePointContainer<container_t, holder_t>::component(HashedString key,
                                                        std::size_t n) const {
  using namespace Acts::HashedStringLiteral;
  switch (key) {
    case "TopHalfStripLength"_hash:
    case "BottomHalfStripLength"_hash:
    case "TopStripDirection"_hash:
    case "BottomStripDirection"_hash:
    case "StripCenterDistance"_hash:
    case "TopStripCenterPosition"_hash:
      return std::any_cast<T>(container().component_impl(key, n));
    default:
      throw std::runtime_error("no such component " + std::to_string(key));
  }
}

template <typename container_t, template <typename> class holder_t>
SpacePointContainer<container_t, holder_t>&
SpacePointContainer<container_t, holder_t>::operator=(
    SpacePointContainer<container_t, holder_t>&& other) noexcept {
  m_container = std::exchange(other.m_container.ptr, nullptr);
  return *this;
}

template <typename container_t, template <typename> class holder_t>
inline std::size_t SpacePointContainer<container_t, holder_t>::size() const {
  return container().size_impl();
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::iterator
SpacePointContainer<container_t, holder_t>::begin() {
  return {*this, 0};
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::iterator
SpacePointContainer<container_t, holder_t>::end() {
  return {*this, size()};
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::const_iterator
SpacePointContainer<container_t, holder_t>::begin() const {
  return {*this, 0};
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::const_iterator
SpacePointContainer<container_t, holder_t>::end() const {
  return {*this, size()};
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline container_t& SpacePointContainer<container_t, holder_t>::container() {
  return *m_container;
}

template <typename container_t, template <typename> class holder_t>
inline const container_t&
SpacePointContainer<container_t, holder_t>::container() const {
  return *m_container;
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::ValueType&
SpacePointContainer<container_t, holder_t>::sp(std::size_t n) {
  return container().get_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::ValueType&
SpacePointContainer<container_t, holder_t>::sp(std::size_t n) const {
  return container().get_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::x(
    std::size_t n) const {
  return container().x_impl(n) - m_options.beamPos[0];
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::y(
    std::size_t n) const {
  return container().y_impl(n) - m_options.beamPos[1];
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::z(
    std::size_t n) const {
  return container().z_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::radius(
    std::size_t n) const {
  float x = this->x(n);
  float y = this->y(n);
  return std::sqrt(x*x + y*y);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::varianceR(
    std::size_t n) const {
  return container().varianceR_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::varianceZ(
    std::size_t n) const {
  return container().varianceZ_impl(n);
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::ProxyType
SpacePointContainer<container_t, holder_t>::proxy(std::size_t n) {
  return {*this, n};
}

template <typename container_t, template <typename> class holder_t>
inline const typename SpacePointContainer<container_t, holder_t>::ProxyType
SpacePointContainer<container_t, holder_t>::proxy(std::size_t n) const {
  return {*this, n};
}

}  // namespace Acts
