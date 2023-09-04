// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/EventData/SpacePointProxy.hpp"
#include "Acts/EventData/SpacePointProxyIterator.hpp"
#include "Acts/EventData/Utils.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <vector>

#include <math.h>

namespace Acts {
struct SpacePointContainerConfig {
  bool useDetailedDoubleMeasurementInfo = false;
  bool isInInternalUnits = false;

  SpacePointContainerConfig toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for "
          "SpacePointContainerConfig");
    }
    using namespace Acts::UnitLiterals;
    SpacePointContainerConfig config = *this;
    config.isInInternalUnits = true;
    return config;
  };
};

struct SpacePointContainerOptions {
  // location of beam in x,y plane.
  // used as offset for Space Points
  Acts::Vector2 beamPos{0 * Acts::UnitConstants::mm,
                        0 * Acts::UnitConstants::mm};
  bool isInInternalUnits = false;

  SpacePointContainerOptions toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for "
          "SpacePointContainerOptions");
    }
    using namespace Acts::UnitLiterals;
    SpacePointContainerOptions options = *this;
    options.isInInternalUnits = true;
    options.beamPos[0] /= 1_mm;
    options.beamPos[1] /= 1_mm;
    return options;
  }
};

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

  using ValueType = typename std::conditional<
      read_only,
      typename std::conditional<
          std::is_const<typename container_t::ValueType>::value,
          typename container_t::ValueType,
          const typename container_t::ValueType>::type,
      typename container_t::ValueType>::type;
  using ProxyType =
      typename std::conditional<read_only, ConstSpacePointProxyType,
                                SpacePointProxyType>::type;
  using value_type = ProxyType;

 public:
  // Constructors
  // It makes sense to support both options of
  // taking or not the ownership

  // Do not take ownership
  // Activate only if holder_t is RefHolder
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::RefHolder>::value>>
  SpacePointContainer(const Acts::SpacePointContainerConfig& config,
                      const Acts::SpacePointContainerOptions& options,
                      container_t& container);

  // Take the ownership
  // Activate only if holder_t is ValueHolder
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::ValueHolder>::value>>
  SpacePointContainer(const Acts::SpacePointContainerConfig& config,
                      const Acts::SpacePointContainerOptions& options,
                      container_t&& container);

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
  ValueType& sp(const std::size_t n);

  ValueType& sp(const std::size_t n) const;

 private:
  void initialize();

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  container_t& container();

  const container_t& container() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  ProxyType& proxy(const std::size_t n);

  const ProxyType& proxy(const std::size_t n) const;

  const std::vector<ProxyType>& proxies() const;

  template <bool RO = read_only, typename = std::enable_if_t<!RO>>
  std::vector<ProxyType>& proxies();

 private:
  float x(const std::size_t n) const;
  float y(const std::size_t n) const;
  float z(const std::size_t n) const;
  float phi(const std::size_t n) const;
  float radius(const std::size_t n) const;
  float varianceR(const std::size_t n) const;
  float varianceZ(const std::size_t n) const;

  const Acts::Vector3& topStripVector(const std::size_t n) const;
  const Acts::Vector3& bottomStripVector(const std::size_t n) const;
  const Acts::Vector3& stripCenterDistance(const std::size_t n) const;
  const Acts::Vector3& topStripCenterPosition(const std::size_t n) const;

 private:
  Acts::SpacePointContainerConfig m_config;
  Acts::SpacePointContainerOptions m_options;
  Acts::SpacePointData m_data;
  holder_t<container_t> m_container;
  std::vector<ProxyType> m_proxies;
};

}  // namespace Acts

#include "Acts/EventData/SpacePointContainer.ipp"
