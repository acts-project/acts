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

template <typename container_t,
	  template <typename> class holder_t>
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

public:
  static constexpr bool read_only = true;

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

  using ValueType = typename std::conditional<read_only,
					      typename std::conditional<std::is_const<typename container_t::ValueType>::value,
							       typename container_t::ValueType,
							       const typename container_t::ValueType>::type,
					      typename container_t::ValueType>::type;
  using ProxyType = typename std::conditional<read_only,
					      ConstSpacePointProxyType,
					      SpacePointProxyType>::type;
  
 public:
  // Constructors
  // It makes sense to support both options of
  // taking or not the ownership

  // Do not take ownership
  // Activate only if holder_t is RefHolder
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::RefHolder>::value>>
  SpacePointContainer(container_t& container);

  // Take the ownership
  // Activate only if holder_t is ValueHolder
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
					  H, Acts::detail::ValueHolder>::value>>
  SpacePointContainer(container_t&& container);

  // If we take ownership, forbid copy operations
  // Need to define copy operations only if holder_t is RefHolder !!!
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::RefHolder>::value>>
  SpacePointContainer(SpacePointContainer& other);

  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
					  H, Acts::detail::RefHolder>::value,
					bool>>
  SpacePointContainer& operator=(SpacePointContainer& other);

  // move operations
  SpacePointContainer(SpacePointContainer&& other) noexcept;
  SpacePointContainer& operator=(SpacePointContainer&& other) noexcept;

  // Destructor
  ~SpacePointContainer() = default;

  std::size_t size() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  iterator begin();

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  ValueType& sp(std::size_t n);

  ValueType& sp(std::size_t n) const;

private:
  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  container_t& container();

  const container_t& container() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  ProxyType& proxy(std::size_t n);
  
  const ProxyType& proxy(std::size_t n) const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  std::vector<ProxyType>& proxies();

  const std::vector<ProxyType>& proxies() const;

 private:
  float x(std::size_t n) const;
  float y(std::size_t n) const;
  float z(std::size_t n) const;
  float radius(std::size_t n) const;
  float varianceR(std::size_t n) const;
  float varianceZ(std::size_t n) const;

  // component methods for additional quantities
  template<typename T>
  T component(HashedString key, std::size_t n) const
  {
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

 private:
  holder_t<container_t> m_container;
  std::vector<ProxyType> m_proxies{};
};

}  // namespace Acts

#include "Acts/EventData/SpacePointContainer.ipp"
