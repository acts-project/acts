// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointProxy.hpp"
#include "Acts/Utilities/Holders.hpp"

namespace Acts {

template <typename container_t, bool read_only>
class SpacePointProxyIterator {
 public:
  using ContainerType = typename std::conditional<read_only, const container_t,
                                                  container_t>::type;
  using ProxyType = typename std::conditional<
      read_only, typename container_t::ConstSpacePointProxyType,
      typename container_t::SpacePointProxyType>::type;

  using iterator_category = std::random_access_iterator_tag;
  using value_type =
      typename std::conditional<read_only, const ProxyType, ProxyType>::type;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

 public:
  // Constructors
  SpacePointProxyIterator(ContainerType&& container,
                          std::size_t index) = delete;
  SpacePointProxyIterator(ContainerType& container, std::size_t index);

  SpacePointProxyIterator& operator++();
  SpacePointProxyIterator& operator--();
  SpacePointProxyIterator operator++(int);
  SpacePointProxyIterator operator--(int);

  bool operator==(const SpacePointProxyIterator& other) const;
  bool operator!=(const SpacePointProxyIterator& other) const;
  bool operator<(const SpacePointProxyIterator& other) const;
  bool operator>(const SpacePointProxyIterator& other) const;
  bool operator<=(const SpacePointProxyIterator& other) const;
  bool operator>=(const SpacePointProxyIterator& other) const;

  SpacePointProxyIterator& operator+=(const std::size_t& offset);
  SpacePointProxyIterator& operator-=(const std::size_t& offset);

  SpacePointProxyIterator operator+(const std::size_t& offset) const;
  SpacePointProxyIterator operator-(const std::size_t& offset) const;

  difference_type operator-(const SpacePointProxyIterator& other) const;

  value_type& operator*() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  value_type& operator*();

 private:
  Acts::detail::RefHolder<ContainerType> m_container;
  std::size_t m_index;
};

}  // namespace Acts

#include "Acts/EventData/SpacePointProxyIterator.ipp"
