// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointProxy.hpp"
#include "Acts/EventData/SpacePointProxyIterator.hpp"
#include "Acts/EventData/Utils.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Definitions/Algebra.hpp"

#include <any>
#include <vector>

namespace Acts {
// add read and write here?
template <typename container_t, template <typename> class holder_t>
class SpacePointContainer {
 public:
  friend class Acts::SpacePointProxy<
      Acts::SpacePointContainer<container_t, holder_t>, false>;
  friend class Acts::SpacePointProxy<
      Acts::SpacePointContainer<container_t, holder_t>, true>;
  friend class Acts::SpacePointProxyIterator<
      Acts::SpacePointContainer<container_t, holder_t>, false>;
  friend class Acts::SpacePointProxyIterator<
      Acts::SpacePointContainer<container_t, holder_t>, true>;

  static constexpr bool ReadOnly = true;

  using SpacePointProxyType =
      Acts::SpacePointProxy<Acts::SpacePointContainer<container_t, holder_t>,
                            false>;
  using ConstSpacePointProxyType =
      Acts::SpacePointProxy<Acts::SpacePointContainer<container_t, holder_t>,
                            true>;
  using iterator = Acts::SpacePointProxyIterator<
      Acts::SpacePointContainer<container_t, holder_t>, false>;
  using const_iterator = Acts::SpacePointProxyIterator<
      Acts::SpacePointContainer<container_t, holder_t>, true>;

  using ValueType = typename container_t::ValueType;

 public:
  // Constructors
  // It makes sense to support both options of
  // taking or not the ownership

  // Do not take ownership
  // Activate only if holder_t is RefHolder
  template <template <typename> class H = holder_t,
            std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::RefHolder>::value>* = nullptr>
  SpacePointContainer(container_t& container) : m_container(container) {
    std::size_t n = size();
    m_proxies.reserve(n);
    for (std::size_t i(0); i < n; i++)
      m_proxies.emplace_back(*this, i);
  }

  // Take the ownership
  // Activate only if holder_t is ValueHolder
  template <template <typename> class H = holder_t,
            std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::ValueHolder>::value>* = nullptr>
  SpacePointContainer(container_t&& container) : m_container(container) {
    std::size_t n = size();
    m_proxies.reserve(n);
    for (std::size_t i(0); i < n; i++)
      m_proxies.emplace_back(*this, i);
  }

  // If we take ownership, forbid copy operations
  // Need to define copy operations only if holder_t is RefHolder !!!
  template <template <typename> class H = holder_t,
            std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::RefHolder>::value>* = nullptr>
  SpacePointContainer(SpacePointContainer& other)
      : m_container(*m_container.ptr),
        m_proxies(other.m_proxies.begin(), other.m_proxies.end()) {}

  template <template <typename> class H = holder_t,
            std::enable_if_t<Acts::detail::is_same_template<
                                 H, Acts::detail::RefHolder>::value,
                             bool> = true>
  SpacePointContainer& operator=(SpacePointContainer& other) {
    m_container.ptr = other.m_container.ptr;
    m_proxies.insert(m_proxies.end(), other.m_proxies.begin(),
                     other.m_proxies.end());
    return *this;
  }

  // move operations
  SpacePointContainer(SpacePointContainer&& other) noexcept
      : m_container(std::exchange(other.m_container.ptr, nullptr)),
        m_proxies(std::move(other.m_proxies)) {}

  SpacePointContainer& operator=(SpacePointContainer&& other) noexcept {
    m_container = std::exchange(other.m_container.ptr, nullptr);
    m_proxies = std::move(other.m_proxies);
    return *this;
  }

  // Destructor
  ~SpacePointContainer() = default;

  std::size_t size() const;

  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  SpacePointProxyType& get(std::size_t n);
  const SpacePointProxyType& get(std::size_t n) const;

  ValueType& sp(std::size_t n) { return container().storage()[n]; }
  const ValueType& sp(std::size_t n) const { return container().storage()[n]; }

  // do these need to be private or public?
  const container_t& container() const;

 private:
  float x(std::size_t n) const;
  float y(std::size_t n) const;
  float z(std::size_t n) const;
  float radius(std::size_t n) const;
  float varianceR(std::size_t n) const;
  float varianceZ(std::size_t n) const;

  // Dynamic variables
  float topHalfStripLength(std::size_t n) const;
  float bottomHalfStripLength(std::size_t n) const;
  Acts::Vector3 topStripDirection(std::size_t n) const;
  Acts::Vector3 bottomStripDirection(std::size_t n) const;
  Acts::Vector3 stripCenterDistance(std::size_t n) const;
  Acts::Vector3 topStripCenterPosition(std::size_t n) const;

  // // component methods for additional quantities
  // template<typename T>
  // const T& component(HashedString key, std::size_t n) const
  // {
  //   using namespace Acts::HashedStringLiteral;
  //   switch (key) {
  //   case "TopHalfStripLength"_hash:
  //   case "BottomHalfStripLength"_hash:
  //   case "TopStripDirection"_hash:
  //   case "BottomStripDirection"_hash:
  //   case "StripCenterDistance"_hash:
  //   case "TopStripCenterPosition"_hash:
  // 	return *std::any_cast<T*>(container().component_impl(key, n));
  //   default:
  // 	throw std::runtime_error("no such component " + std::to_string(key));
  //   }
  // }

 private:
  holder_t<container_t> m_container;
  std::vector<SpacePointProxyType> m_proxies{};  // this will go away ?
};

// Deduction rules
template <typename container_t>
SpacePointContainer(container_t& container)
    -> SpacePointContainer<container_t, Acts::detail::RefHolder>;

template <typename container_t>
SpacePointContainer(container_t&& container)
    -> SpacePointContainer<container_t, Acts::detail::ValueHolder>;

// Implementations
template <typename container_t, template <typename> class holder_t>
inline std::size_t SpacePointContainer<container_t, holder_t>::size() const {
  return container().size_impl();
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::iterator
SpacePointContainer<container_t, holder_t>::begin() {
  return {*this, 0};
}

template <typename container_t, template <typename> class holder_t>
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
inline const container_t&
SpacePointContainer<container_t, holder_t>::container() const {
  return *m_container;
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
inline typename SpacePointContainer<container_t, holder_t>::SpacePointProxyType&
SpacePointContainer<container_t, holder_t>::get(std::size_t n) {
  return m_proxies[n];
}

template <typename container_t, template <typename> class holder_t>
inline const typename SpacePointContainer<container_t,
                                          holder_t>::SpacePointProxyType&
SpacePointContainer<container_t, holder_t>::get(std::size_t n) const {
  return m_proxies.at(n);
}

}  // namespace Acts
