// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/TrackProxyCommon.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace Acts {

class Surface;

namespace detail_anytrack {

using ParametersMap = detail_tpc::ParametersMap;
using ConstParametersMap = detail_tpc::ConstParametersMap;
using CovarianceMap = detail_tpc::CovarianceMap;
using ConstCovarianceMap = detail_tpc::ConstCovarianceMap;

/// Base class for read-only track handlers
/// Provides accessors that return const references to the underlying data
class TrackHandlerConstBase {
 public:
  virtual ~TrackHandlerConstBase() = default;

  /// Get the reference surface
  virtual const Surface* referenceSurface(const void* container,
                                          TrackIndexType index) const = 0;

  /// Check if track has a reference surface
  virtual bool hasReferenceSurface(const void* container,
                                   TrackIndexType index) const = 0;

  /// Get the particle hypothesis
  virtual ParticleHypothesis particleHypothesis(const void* container,
                                                TrackIndexType index) const = 0;

  /// Get parameter vector
  virtual ConstParametersMap parameters(const void* container,
                                        TrackIndexType index) const = 0;

  /// Get covariance matrix
  virtual ConstCovarianceMap covariance(const void* container,
                                        TrackIndexType index) const = 0;

  /// Get number of track states
  virtual unsigned int nTrackStates(const void* container,
                                    TrackIndexType index) const = 0;

  /// Check if track has a specific dynamic column
  virtual bool hasColumn(const void* container, HashedString key) const = 0;

  /// Get a dynamic column component (type-erased)
  virtual std::any component(const void* container, TrackIndexType index,
                             HashedString key) const = 0;
};

/// Base class for mutable track handlers.
/// Extends the const interface with mutable references to the data.
class TrackHandlerMutableBase : public TrackHandlerConstBase {
 public:
  using TrackHandlerConstBase::component;
  using TrackHandlerConstBase::covariance;
  using TrackHandlerConstBase::parameters;

  /// Get mutable parameter vector
  virtual ParametersMap parameters(void* container,
                                   TrackIndexType index) const = 0;

  /// Get mutable covariance matrix
  virtual CovarianceMap covariance(void* container,
                                   TrackIndexType index) const = 0;

  /// Get mutable dynamic column component (type-erased)
  virtual std::any component(void* container, TrackIndexType index,
                             HashedString key) const = 0;
};

template <typename container_t>
struct TrackHandlerTraits {
  using Container = std::remove_const_t<container_t>;
  static constexpr bool ReadOnly =
      std::is_const_v<container_t> || Container::ReadOnly;
};

/// Concrete handler for a specific track container type
/// This templated class provides static instances that implement the virtual
/// interface for a specific track container type. The static instance approach
/// avoids heap allocation per handle.
template <typename container_t,
          bool read_only = TrackHandlerTraits<container_t>::ReadOnly>
class TrackHandler;

template <typename container_t>
class TrackHandler<container_t, true> final : public TrackHandlerConstBase {
  using ContainerType = typename TrackHandlerTraits<container_t>::Container;

 public:
  /// Get the singleton instance for this track container type
  static const TrackHandler& instance() {
    static const TrackHandler s_instance;
    return s_instance;
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return &tc->getTrack(index).referenceSurface();
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).hasReferenceSurface();
  }

  ParticleHypothesis particleHypothesis(const void* container,
                                        TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).particleHypothesis();
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).parameters();
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).covariance();
  }

  unsigned int nTrackStates(const void* container,
                            TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).nTrackStates();
  }

  bool hasColumn(const void* container, HashedString key) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->hasColumn(key);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->container().component_impl(key, index);
  }

 private:
  TrackHandler() = default;
};

template <typename container_t>
class TrackHandler<container_t, false> final : public TrackHandlerMutableBase {
  using ContainerType = typename TrackHandlerTraits<container_t>::Container;

 public:
  /// Get the singleton instance for this track container type
  static const TrackHandler& instance() {
    static const TrackHandler s_instance;
    return s_instance;
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return &tc->getTrack(index).referenceSurface();
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).hasReferenceSurface();
  }

  ParticleHypothesis particleHypothesis(const void* container,
                                        TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).particleHypothesis();
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).parameters();
  }

  ParametersMap parameters(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<ContainerType*>(container);
    return tc->getTrack(index).parameters();
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).covariance();
  }

  CovarianceMap covariance(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<ContainerType*>(container);
    return tc->getTrack(index).covariance();
  }

  unsigned int nTrackStates(const void* container,
                            TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->getTrack(index).nTrackStates();
  }

  bool hasColumn(const void* container, HashedString key) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->hasColumn(key);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const ContainerType*>(container);
    return tc->container().component_impl(key, index);
  }

  std::any component(void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    auto* tc = static_cast<ContainerType*>(container);
    return tc->container().component_impl(key, index);
  }

 private:
  TrackHandler() = default;
};

}  // namespace detail_anytrack

/// Type-erased track object
/// This class provides a type-erased interface to track proxies without
/// requiring heap allocation per instance. It stores a pointer to the track
/// container and the track index, similar to how TrackProxy works internally.
///
/// The object does not take ownership of the container - the container must
/// outlive the object.
///
/// @tparam read_only If true, provides read-only access to the track
///
/// Usage:
/// @code
/// TrackContainer tracks = ...;
/// auto track = tracks.getTrack(0);
/// AnyMutableTrack AnyTrackProxy(track);  // Mutable track
/// AnyConstTrack constTrack(track);  // Read-only track
/// std::cout << "Chi2: " << AnyTrackProxy.chi2() << std::endl;
/// @endcode
template <bool read_only = true>
class AnyTrackProxy : public TrackProxyCommon<AnyTrackProxy<read_only>,
                                              TrackIndexType, read_only> {
  using Base =
      TrackProxyCommon<AnyTrackProxy<read_only>, TrackIndexType, read_only>;

 public:
  /// Indicates whether this track is read-only
  static constexpr bool ReadOnly = read_only;

  /// Alias for the mutable version
  using MutableTrackHandle = AnyTrackProxy<false>;

  /// Alias for the const version
  using ConstTrackHandle = AnyTrackProxy<true>;
  /// Alias for the const proxy type
  using ConstProxyType = AnyTrackProxy<true>;

  /// Mutable parameters map type
  using ParametersMap = detail_anytrack::ParametersMap;
  /// Const parameters map type
  using ConstParametersMap = detail_anytrack::ConstParametersMap;
  /// Mutable covariance map type
  using CovarianceMap = detail_anytrack::CovarianceMap;
  /// Const covariance map type
  using ConstCovarianceMap = detail_anytrack::ConstCovarianceMap;

  /// Container pointer type
  using ContainerPointer = std::conditional_t<ReadOnly, const void*, void*>;
  /// Handler pointer type
  using HandlerPointer =
      std::conditional_t<ReadOnly,
                         const detail_anytrack::TrackHandlerConstBase*,
                         const detail_anytrack::TrackHandlerMutableBase*>;

  /// Copy constructor: const to const or mutable to mutable
  /// @param other the other track
  AnyTrackProxy(const AnyTrackProxy& other) = default;

  /// Copy assignment operator: const to const or mutable to mutable
  /// @param other the other track
  /// @return reference to this track
  AnyTrackProxy& operator=(const AnyTrackProxy& other) = default;

 private:
  /// Constructor from mutable track
  /// @note Only available if this is read-only
  /// @param other the other track
  explicit AnyTrackProxy(const MutableTrackHandle& other)
    requires ReadOnly
      : m_container(other.m_container),
        m_index(other.m_index),
        m_handler(other.m_handler) {}

  /// Copy assignment operator from mutable track
  /// @note Only available if this is read-only
  /// @param other the other track
  /// @return reference to this track
  AnyTrackProxy& operator=(const MutableTrackHandle& other)
    requires ReadOnly
  {
    m_container = other.m_container;
    m_index = other.m_index;
    m_handler = other.m_handler;
    return *this;
  }

 public:
  /// Construct from a concrete track proxy
  /// @tparam track_proxy_t The concrete track proxy type
  /// @param track The track proxy to wrap
  /// @note Does not take ownership. The track container must outlive this object.
  /// @note AnyConstTrack can be constructed from both mutable and const track proxies.
  ///       AnyMutableTrack can only be constructed from mutable track proxies.
  template <TrackProxyConcept track_proxy_t>
    requires(ReadOnly || !track_proxy_t::ReadOnly)
  explicit AnyTrackProxy(track_proxy_t track)
      : m_container(nullptr), m_index(track.m_index), m_handler(nullptr) {
    using container_t = std::remove_reference_t<decltype(*track.m_container)>;
    auto* containerPtr = &(*track.m_container);
    if constexpr (ReadOnly) {
      m_container = static_cast<const void*>(containerPtr);
    } else {
      m_container = static_cast<void*>(containerPtr);
    }
    m_handler = &detail_anytrack::TrackHandler<container_t>::instance();
  }

  /// Get the index of this track
  /// @return The track index
  TrackIndexType index() const { return m_index; }

  /// Get the reference surface
  /// @return Reference to the reference surface
  const Surface& referenceSurface() const {
    const Surface* surface =
        constHandler()->referenceSurface(containerPtr(), m_index);
    assert(surface != nullptr);
    return *surface;
  }

  /// Check if track has a reference surface
  /// @return true if a reference surface exists
  bool hasReferenceSurface() const {
    return constHandler()->hasReferenceSurface(containerPtr(), m_index);
  }

  /// Get the particle hypothesis
  /// @return The particle hypothesis
  ParticleHypothesis particleHypothesis() const {
    return constHandler()->particleHypothesis(containerPtr(), m_index);
  }

  /// Get the bound parameters map
  /// @return The const parameters map
  ConstParametersMap parameters() const {
    return constHandler()->parameters(containerPtr(), m_index);
  }

  /// Get the mutable bound parameters map
  /// @return The mutable parameters map
  ParametersMap parameters()
    requires(!ReadOnly)
  {
    return mutableHandler()->parameters(mutableContainerPtr(), m_index);
  }

  /// Get a parameter value by index
  /// @param i The parameter index
  /// @return The parameter value
  double parameter(std::size_t i) const { return parameters()[i]; }

  /// Get a mutable reference to a parameter component
  /// @param i The parameter index
  /// @return Mutable reference to the value
  double& parameter(std::size_t i)
    requires(!ReadOnly)
  {
    auto params = parameters();
    return params[i];
  }

  /// Get the covariance map
  /// @return The const covariance map
  ConstCovarianceMap covariance() const {
    return constHandler()->covariance(containerPtr(), m_index);
  }

  /// Get the mutable covariance map
  /// @return The mutable covariance map
  CovarianceMap covariance()
    requires(!ReadOnly)
  {
    return mutableHandler()->covariance(mutableContainerPtr(), m_index);
  }

  /// Get a covariance matrix element
  /// @param i Row index
  /// @param j Column index
  /// @return The covariance element
  double covariance(std::size_t i, std::size_t j) const {
    return covariance()(i, j);
  }

  /// Get a mutable reference to a covariance element
  /// @param i Row index
  /// @param j Column index
  /// @return Mutable reference to the covariance value
  double& covariance(std::size_t i, std::size_t j)
    requires(!ReadOnly)
  {
    auto cov = covariance();
    return cov(i, j);
  }

  using Base::absoluteMomentum;
  using Base::charge;
  using Base::fourMomentum;
  using Base::loc0;
  using Base::loc1;
  using Base::phi;
  using Base::qOverP;
  using Base::theta;
  using Base::time;
  using Base::transverseMomentum;

  /// Get the number of track states
  /// @return The number of track states
  unsigned int nTrackStates() const {
    return constHandler()->nTrackStates(containerPtr(), m_index);
  }

  /// Check if the track has a specific dynamic column
  /// @param key The hashed column key
  /// @return true if the column exists
  bool hasColumn(HashedString key) const {
    return constHandler()->hasColumn(containerPtr(), key);
  }

  /// Get a mutable reference to a dynamic column component
  /// @tparam T The type of the component
  /// @tparam key String key for the component to access
  /// @return Mutable reference to the component
  /// @note Only available if ReadOnly is false
  template <typename T, HashedString key>
  T& component()
    requires(!ReadOnly)
  {
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  /// Get a mutable reference to a dynamic column component
  /// @tparam T The type of the component
  /// @param key String key for the component to access
  /// @return Mutable reference to the component
  /// @note Only available if ReadOnly is false
  template <typename T>
  T& component(HashedString key)
    requires(!ReadOnly)
  {
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  /// Get a const reference to a dynamic column component
  /// @tparam T The type of the component
  /// @tparam key String key for the component to access
  /// @return Const reference to the component
  template <typename T, HashedString key>
  const T& component() const {
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  /// Get a const reference to a dynamic column component
  /// @tparam T The type of the component
  /// @param key String key for the component to access
  /// @return Const reference to the component
  template <typename T>
  const T& component(HashedString key) const {
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

 private:
  template <bool R>
  friend class AnyTrackProxy;

  const detail_anytrack::TrackHandlerConstBase* constHandler() const {
    return m_handler;
  }

  const detail_anytrack::TrackHandlerMutableBase* mutableHandler() const
    requires(!ReadOnly)
  {
    return m_handler;
  }

  const void* containerPtr() const {
    if constexpr (ReadOnly) {
      return m_container;
    } else {
      return static_cast<const void*>(m_container);
    }
  }

  void* mutableContainerPtr() const
    requires(!ReadOnly)
  {
    return m_container;
  }

  /// Pointer to the track container (non-owning)
  ContainerPointer m_container;

  /// Track index within the container
  TrackIndexType m_index;

  /// Pointer to the static handler instance
  HandlerPointer m_handler;
};

/// Alias for read-only type-erased track
using AnyConstTrackProxy = AnyTrackProxy<true>;

/// Alias for mutable type-erased track (though currently all operations are
/// const)
using AnyMutableTrackProxy = AnyTrackProxy<false>;

static_assert(ConstTrackProxyConcept<AnyConstTrackProxy>);
static_assert(MutableTrackProxyConcept<AnyMutableTrackProxy>);

}  // namespace Acts
