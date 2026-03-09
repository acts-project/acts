// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/EventData/SpacePointProxy.hpp"
#include "Acts/EventData/Utils.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"

#include <vector>

namespace Acts {

/// Configuration flags for the space point container.
struct SpacePointContainerConfig {
  /// Whether to use detailed double measurement information
  bool useDetailedDoubleMeasurementInfo = false;

  /// Whether the space point data is in internal units
  bool isInInternalUnits = true;
};

/// Construction options for the space point container.
struct SpacePointContainerOptions {
  /// Location of beam in x,y plane, used as offset for space points
  Acts::Vector2 beamPos{0 * Acts::UnitConstants::mm,
                        0 * Acts::UnitConstants::mm};

  /// Whether the space point data is in internal units
  bool isInInternalUnits = true;
};

/// Container wrapper providing space point proxy access.
template <typename container_t, template <typename> class holder_t>
class SpacePointContainer {
 public:
  friend class Acts::SpacePointProxy<
      Acts::SpacePointContainer<container_t, holder_t>>;

 public:
  /// space point proxy type
  using SpacePointProxyType =
      Acts::SpacePointProxy<Acts::SpacePointContainer<container_t, holder_t>>;

  /// Iterator type
  using iterator = detail::ContainerIterator<
      Acts::SpacePointContainer<container_t, holder_t>, SpacePointProxyType&,
      std::size_t, false>;
  /// Const iterator type
  using const_iterator = detail::ContainerIterator<
      Acts::SpacePointContainer<container_t, holder_t>,
      const SpacePointProxyType&, std::size_t, true>;

  /// Value type of the underlying container
  using ValueType = typename container_t::ValueType;
  /// Proxy type
  using ProxyType = SpacePointProxyType;
  /// Value type
  using value_type = ProxyType;
  /// Size type
  using size_type = std::size_t;

 public:
  // Constructors
  // It makes sense to support both options of
  // taking or not the ownership

  // Do not take ownership
  // Activate only if holder_t is RefHolder
  /// Construct from config, options, and container reference
  /// @param config Configuration for the space point container
  /// @param options Options for the space point container
  /// @param container Reference to the underlying space point container
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::RefHolder>::value>>
  SpacePointContainer(const Acts::SpacePointContainerConfig& config,
                      const Acts::SpacePointContainerOptions& options,
                      const container_t& container);

  // Take the ownership
  // Activate only if holder_t is ValueHolder
  /// Construct from config, options, and container by moving ownership
  /// @param config Configuration for the space point container
  /// @param options Options for the space point container
  /// @param container Rvalue reference to the underlying space point container
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<Acts::detail::is_same_template<
                H, Acts::detail::ValueHolder>::value>>
  SpacePointContainer(const Acts::SpacePointContainerConfig& config,
                      const Acts::SpacePointContainerOptions& options,
                      container_t&& container);

  // No copy operations
  SpacePointContainer(SpacePointContainer& other) = delete;
  SpacePointContainer& operator=(SpacePointContainer& other) = delete;

  // move operations
  /// Move constructor
  /// @param other space point container to move from
  SpacePointContainer(SpacePointContainer&& other) noexcept = default;
  /// Move assignment operator
  /// @param other space point container to move from
  /// @return Reference to this container
  SpacePointContainer& operator=(SpacePointContainer&& other) noexcept =
      default;

  // Destructor
  ~SpacePointContainer() = default;

  /// Return the number of space points
  /// @return Number of space points
  std::size_t size() const;

  /// Return an iterator to the beginning
  /// @return Iterator to the first space point
  iterator begin();
  /// Return an iterator to the end
  /// @return Iterator to one past the last space point
  iterator end();
  /// Return a const iterator to the beginning
  /// @return Const iterator to the first space point
  const_iterator cbegin() const;
  /// Return a const iterator to the end
  /// @return Const iterator to one past the last space point
  const_iterator cend() const;
  /// Return a const iterator to the beginning
  /// @return Const iterator to the first space point
  const_iterator begin() const;
  /// Return a const iterator to the end
  /// @return Const iterator to one past the last space point
  const_iterator end() const;

  /// Access space point proxy at index n
  /// @param n Index of the space point
  /// @return Reference to the space point proxy
  ProxyType& at(const std::size_t n);
  /// Access const space point proxy at index n
  /// @param n Index of the space point
  /// @return Const reference to the space point proxy
  const ProxyType& at(const std::size_t n) const;
  /// Access underlying space point value at index n
  /// @param n Index of the space point
  /// @return Const reference to the underlying space point value
  const ValueType& sp(const std::size_t n) const;

 private:
  void initialize();

  const container_t& container() const;
  const ProxyType& proxy(const std::size_t n) const;
  std::vector<ProxyType>& proxies();
  const std::vector<ProxyType>& proxies() const;

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

  Acts::SpacePointContainerConfig m_config;
  Acts::SpacePointContainerOptions m_options;
  Acts::SpacePointData m_data;
  holder_t<const container_t> m_container;
  std::vector<ProxyType> m_proxies;
};

}  // namespace Acts

#include "Acts/EventData/SpacePointContainer.ipp"
