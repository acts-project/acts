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
SpacePointContainer<container_t, holder_t>::SpacePointContainer(container_t& container)
  : m_container(container)
{
  std::size_t n = this->size();
  m_proxies.reserve(n);
  for (std::size_t i(0); i < n; ++i) {
    m_proxies.emplace_back(*this, i);
  }
}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(container_t&& container)
  : m_container(std::move(container))
{
  std::size_t n = this->size();
  m_proxies.reserve(n);
  for (std::size_t i(0); i < n; ++i) {
    m_proxies.emplace_back(*this, i);
  }
}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(SpacePointContainer<container_t, holder_t>& other)
  : m_container(*m_container.ptr),
    m_proxies(other.m_proxies.begin(), other.m_proxies.end())
{}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>&
SpacePointContainer<container_t, holder_t>::operator=(SpacePointContainer<container_t, holder_t>& other)
{
  m_container.ptr = other.m_container.ptr;
  m_proxies.insert(m_proxies.end(), other.m_proxies.begin(),
		   other.m_proxies.end());
  return *this;
}

template <typename container_t, template <typename> class holder_t>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(SpacePointContainer<container_t, holder_t>&& other) noexcept
  : m_container(std::exchange(other.m_container.ptr, nullptr)),
    m_proxies(std::move(other.m_proxies))
{}

template <typename container_t, template <typename> class holder_t>
SpacePointContainer<container_t, holder_t>&
SpacePointContainer<container_t, holder_t>::operator=(SpacePointContainer<container_t, holder_t>&& other) noexcept
{
  m_container = std::exchange(other.m_container.ptr, nullptr);
  m_proxies = std::move(other.m_proxies);
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
inline container_t&
SpacePointContainer<container_t, holder_t>::container() {
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
  return container().storage()[n];
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::ValueType&
SpacePointContainer<container_t, holder_t>::sp(std::size_t n) const {
  return container().storage()[n];
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::x(
    std::size_t n) const {
  return container().x_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::y(
    std::size_t n) const {
  return container().y_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::z(
    std::size_t n) const {
  return container().z_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::radius(
    std::size_t n) const {
  return container().radius_impl(n);
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
inline float SpacePointContainer<container_t, holder_t>::topHalfStripLength(
    std::size_t n) const {
  return container().topHalfStripLength_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::bottomHalfStripLength(
    std::size_t n) const {
  return container().bottomHalfStripLength_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline Acts::Vector3
SpacePointContainer<container_t, holder_t>::topStripDirection(
    std::size_t n) const {
  return container().topStripDirection_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline Acts::Vector3
SpacePointContainer<container_t, holder_t>::bottomStripDirection(
    std::size_t n) const {
  return container().bottomStripDirection_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline Acts::Vector3
SpacePointContainer<container_t, holder_t>::stripCenterDistance(
    std::size_t n) const {
  return container().stripCenterDistance_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline Acts::Vector3
SpacePointContainer<container_t, holder_t>::topStripCenterPosition(
    std::size_t n) const {
  return container().topStripCenterPosition_impl(n);
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::ProxyType&
SpacePointContainer<container_t, holder_t>::get(std::size_t n) {
  return m_proxies[n];
}

template <typename container_t, template <typename> class holder_t>
inline const typename SpacePointContainer<container_t, holder_t>::ProxyType&
SpacePointContainer<container_t, holder_t>::get(std::size_t n) const {
  return m_proxies.at(n);
}

}  // namespace Acts
