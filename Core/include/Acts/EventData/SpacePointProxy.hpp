// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <any>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace Acts {

template <typename container_t, bool read_only>
class SpacePointProxy {
 public:
  using ContainerType = typename std::conditional<read_only, const container_t,
                                                  container_t>::type;
  using ValueType = typename std::conditional<
      read_only,
      typename std::conditional<
          std::is_const<typename ContainerType::ValueType>::value,
          typename ContainerType::ValueType,
          const typename ContainerType::ValueType>::type,
      typename ContainerType::ValueType>::type;

 public:
  // Never take the ownership of the container
  SpacePointProxy(ContainerType&& container, std::size_t index) = delete;
  // Only get the reference
  SpacePointProxy(ContainerType& container, std::size_t index);
  // copy and move operations are defaults

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  ValueType& externalSpacePoint();

  ValueType& externalSpacePoint() const;

  std::size_t index() const;

  float x() const;
  float y() const;
  float z() const;
  float phi() const;
  float radius() const;
  float varianceR() const;
  float varianceZ() const;

  const Acts::Vector3& topStripVector() const;
  const Acts::Vector3& bottomStripVector() const;
  const Acts::Vector3& stripCenterDistance() const;
  const Acts::Vector3& topStripCenterPosition() const;

 private:
  ContainerType& container() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  ContainerType& container();

 private:
  Acts::detail::RefHolder<ContainerType> m_container;
  std::size_t m_index;
};

}  // namespace Acts

#include "Acts/EventData/SpacePointProxy.ipp"
