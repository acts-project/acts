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
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <cassert>
#include <memory>
#include <type_traits>

namespace Acts {

class Surface;

namespace detail {

using ParametersMap = detail_tpc::ParametersMap;
using ConstParametersMap = detail_tpc::ConstParametersMap;
using CovarianceMap = detail_tpc::CovarianceMap;
using ConstCovarianceMap = detail_tpc::ConstCovarianceMap;

/// Base class for read-only track handlers
/// Provides accessors that return const references to the underlying data
class TrackHandlerConstBase {
 public:
  virtual ~TrackHandlerConstBase() = default;

  /// Get the tip index of the track
  virtual const TrackIndexType& tipIndex(const void* container,
                                         TrackIndexType index) const = 0;

  /// Get the stem index of the track
  virtual const TrackIndexType& stemIndex(const void* container,
                                          TrackIndexType index) const = 0;

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

  /// Get theta parameter
  virtual double theta(const void* container, TrackIndexType index) const = 0;

  /// Get phi parameter
  virtual double phi(const void* container, TrackIndexType index) const = 0;

  /// Get qOverP parameter
  virtual double qOverP(const void* container, TrackIndexType index) const = 0;

  /// Get charge
  virtual double charge(const void* container, TrackIndexType index) const = 0;

  /// Get absolute momentum
  virtual double absoluteMomentum(const void* container,
                                  TrackIndexType index) const = 0;

  /// Get transverse momentum
  virtual double transverseMomentum(const void* container,
                                    TrackIndexType index) const = 0;

  /// Get number of measurements
  virtual const unsigned int& nMeasurements(const void* container,
                                            TrackIndexType index) const = 0;

  /// Get number of holes
  virtual const unsigned int& nHoles(const void* container,
                                     TrackIndexType index) const = 0;

  /// Get number of outliers
  virtual const unsigned int& nOutliers(const void* container,
                                        TrackIndexType index) const = 0;

  /// Get number of shared hits
  virtual const unsigned int& nSharedHits(const void* container,
                                          TrackIndexType index) const = 0;

  /// Get chi2
  virtual const float& chi2(const void* container,
                            TrackIndexType index) const = 0;

  /// Get number of degrees of freedom
  virtual const unsigned int& nDoF(const void* container,
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
  using TrackHandlerConstBase::chi2;
  using TrackHandlerConstBase::component;
  using TrackHandlerConstBase::covariance;
  using TrackHandlerConstBase::nDoF;
  using TrackHandlerConstBase::nHoles;
  using TrackHandlerConstBase::nMeasurements;
  using TrackHandlerConstBase::nOutliers;
  using TrackHandlerConstBase::nSharedHits;
  using TrackHandlerConstBase::parameters;
  using TrackHandlerConstBase::stemIndex;
  using TrackHandlerConstBase::tipIndex;

  /// Get mutable tip index
  virtual TrackIndexType& tipIndex(void* container,
                                   TrackIndexType index) const = 0;

  /// Get mutable stem index
  virtual TrackIndexType& stemIndex(void* container,
                                    TrackIndexType index) const = 0;

  /// Get mutable parameter vector
  virtual ParametersMap parameters(void* container,
                                   TrackIndexType index) const = 0;

  /// Get mutable covariance matrix
  virtual CovarianceMap covariance(void* container,
                                   TrackIndexType index) const = 0;

  /// Get mutable number of measurements
  virtual unsigned int& nMeasurements(void* container,
                                      TrackIndexType index) const = 0;

  /// Get mutable number of holes
  virtual unsigned int& nHoles(void* container, TrackIndexType index) const = 0;

  /// Get mutable number of outliers
  virtual unsigned int& nOutliers(void* container,
                                  TrackIndexType index) const = 0;

  /// Get mutable number of shared hits
  virtual unsigned int& nSharedHits(void* container,
                                    TrackIndexType index) const = 0;

  /// Get mutable chi2
  virtual float& chi2(void* container, TrackIndexType index) const = 0;

  /// Get mutable number of degrees of freedom
  virtual unsigned int& nDoF(void* container, TrackIndexType index) const = 0;

  /// Get mutable dynamic column component (type-erased)
  virtual std::any component(void* container, TrackIndexType index,
                             HashedString key) const = 0;
};

/// Concrete handler for a specific track container type
/// This templated class provides static instances that implement the virtual
/// interface for a specific track container type. The static instance approach
/// avoids heap allocation per handle.
template <typename container_t>
class TrackHandler final : public TrackHandlerMutableBase {
 public:
  /// Get the singleton instance for this track container type
  static const TrackHandler& instance() {
    static const TrackHandler s_instance;
    return s_instance;
  }

  const TrackIndexType& tipIndex(const void* container,
                                 TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<TrackIndexType>(tc, index, detail_tp::kTipIndexKey);
  }

  TrackIndexType& tipIndex(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<TrackIndexType>(tc, index, detail_tp::kTipIndexKey);
  }

  const TrackIndexType& stemIndex(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<TrackIndexType>(tc, index, detail_tp::kStemIndexKey);
  }

  TrackIndexType& stemIndex(void* container,
                            TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<TrackIndexType>(tc, index,
                                            detail_tp::kStemIndexKey);
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return &tc->getTrack(index).referenceSurface();
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).hasReferenceSurface();
  }

  ParticleHypothesis particleHypothesis(const void* container,
                                        TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).particleHypothesis();
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).parameters();
  }

  ParametersMap parameters(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return tc->getTrack(index).parameters();
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).covariance();
  }

  CovarianceMap covariance(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return tc->getTrack(index).covariance();
  }

  double theta(const void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).theta();
  }

  double phi(const void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).phi();
  }

  double qOverP(const void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).qOverP();
  }

  double charge(const void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).charge();
  }

  double absoluteMomentum(const void* container,
                          TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).absoluteMomentum();
  }

  double transverseMomentum(const void* container,
                            TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).transverseMomentum();
  }

  const unsigned int& nMeasurements(const void* container,
                                    TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<unsigned int>(tc, index, detail_tp::kMeasurementsKey);
  }

  unsigned int& nMeasurements(void* container,
                              TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<unsigned int>(tc, index,
                                          detail_tp::kMeasurementsKey);
  }

  const unsigned int& nHoles(const void* container,
                             TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<unsigned int>(tc, index, detail_tp::kHolesKey);
  }

  unsigned int& nHoles(void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<unsigned int>(tc, index, detail_tp::kHolesKey);
  }

  const unsigned int& nOutliers(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<unsigned int>(tc, index, detail_tp::kOutliersKey);
  }

  unsigned int& nOutliers(void* container,
                          TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<unsigned int>(tc, index, detail_tp::kOutliersKey);
  }

  const unsigned int& nSharedHits(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<unsigned int>(tc, index, detail_tp::kSharedHitsKey);
  }

  unsigned int& nSharedHits(void* container,
                            TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<unsigned int>(tc, index, detail_tp::kSharedHitsKey);
  }

  const float& chi2(const void* container,
                    TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<float>(tc, index, detail_tp::kChi2Key);
  }

  float& chi2(void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<float>(tc, index, detail_tp::kChi2Key);
  }

  const unsigned int& nDoF(const void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return getConstColumn<unsigned int>(tc, index, detail_tp::kNdfKey);
  }

  unsigned int& nDoF(void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return getMutableColumn<unsigned int>(tc, index, detail_tp::kNdfKey);
  }

  unsigned int nTrackStates(const void* container,
                            TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->getTrack(index).nTrackStates();
  }

  bool hasColumn(const void* container, HashedString key) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->hasColumn(key);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    const auto* tc = static_cast<const container_t*>(container);
    return tc->container().component_impl(key, index);
  }

  std::any component(void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    auto* tc = static_cast<container_t*>(container);
    return tc->container().component_impl(key, index);
  }

 private:
  template <typename value_t>
  static const value_t& getConstColumn(const container_t* tc,
                                       TrackIndexType index, HashedString key) {
    std::any result = tc->container().component_impl(key, index);
    return *std::any_cast<const value_t*>(result);
  }

  template <typename value_t>
  static value_t& getMutableColumn(container_t* tc, TrackIndexType index,
                                   HashedString key) {
    std::any result = tc->container().component_impl(key, index);
    return *std::any_cast<value_t*>(result);
  }

  TrackHandler() = default;
};

}  // namespace detail

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
/// AnyMutableTrack anyTrack(track);  // Mutable track
/// AnyConstTrack constTrack(track);  // Read-only track
/// std::cout << "Chi2: " << anyTrack.chi2() << std::endl;
/// @endcode
template <bool read_only = true>
class AnyTrack {
 public:
  /// Indicates whether this track is read-only
  static constexpr bool ReadOnly = read_only;

  /// Alias for the mutable version
  using MutableTrackHandle = AnyTrack<false>;

  /// Alias for the const version
  using ConstTrackHandle = AnyTrack<true>;
  using ConstProxyType = AnyTrack<true>;

  using ParametersMap = detail::ParametersMap;
  using ConstParametersMap = detail::ConstParametersMap;
  using CovarianceMap = detail::CovarianceMap;
  using ConstCovarianceMap = detail::ConstCovarianceMap;

  using ContainerPointer = std::conditional_t<ReadOnly, const void*, void*>;
  using HandlerPointer =
      std::conditional_t<ReadOnly, const detail::TrackHandlerConstBase*,
                         const detail::TrackHandlerMutableBase*>;

  class ContainerView {
   public:
    ContainerView() = default;
    ContainerView(const detail::TrackHandlerConstBase* handler,
                  const void* container)
        : m_handler(handler), m_container(container) {}

    bool hasColumn(HashedString key) const {
      assert(m_handler != nullptr && m_container != nullptr);
      return m_handler->hasColumn(m_container, key);
    }

   private:
    const detail::TrackHandlerConstBase* m_handler = nullptr;
    const void* m_container = nullptr;
  };

  /// Default constructor creates an invalid track
  AnyTrack()
      : m_container(nullptr), m_index(kTrackIndexInvalid), m_handler(nullptr) {}

  /// Copy constructor: const to const or mutable to mutable
  /// @param other the other track
  AnyTrack(const AnyTrack& other) = default;

  /// Copy assignment operator: const to const or mutable to mutable
  /// @param other the other track
  /// @return reference to this track
  AnyTrack& operator=(const AnyTrack& other) = default;

  /// Constructor from mutable track
  /// @note Only available if this is read-only
  /// @param other the other track
  explicit AnyTrack(const MutableTrackHandle& other)
    requires ReadOnly
      : m_container(other.m_container),
        m_index(other.m_index),
        m_handler(other.m_handler) {}

  /// Copy assignment operator from mutable track
  /// @note Only available if this is read-only
  /// @param other the other track
  /// @return reference to this track
  AnyTrack& operator=(const MutableTrackHandle& other)
    requires ReadOnly
  {
    m_container = other.m_container;
    m_index = other.m_index;
    m_handler = other.m_handler;
    return *this;
  }

  /// Construct from a concrete track proxy
  /// @tparam track_proxy_t The concrete track proxy type
  /// @param track The track proxy to wrap
  /// @note Does not take ownership. The track container must outlive this object.
  /// @note AnyConstTrack can be constructed from both mutable and const track proxies.
  ///       AnyMutableTrack can only be constructed from mutable track proxies.
  template <TrackProxyConcept track_proxy_t>
    requires(ReadOnly || !track_proxy_t::ReadOnly)
  explicit AnyTrack(track_proxy_t track)
      : m_container(nullptr), m_index(track.m_index), m_handler(nullptr) {
    using container_t = std::remove_const_t<
        std::remove_reference_t<decltype(*track.m_container)>>;
    auto* containerPtr = &(*track.m_container);
    if constexpr (ReadOnly) {
      m_container = static_cast<const void*>(containerPtr);
    } else {
      m_container = static_cast<void*>(containerPtr);
    }
    m_handler = &detail::TrackHandler<container_t>::instance();
  }

  /// Check if the track is valid
  /// @return true if this points to a valid track
  bool isValid() const {
    return m_container != nullptr && m_handler != nullptr &&
           m_index != kTrackIndexInvalid;
  }

  /// Explicit conversion to bool for validity checking
  explicit operator bool() const { return isValid(); }

  /// Get the tip index of the track
  /// @return The tip index
  const TrackIndexType& tipIndex() const {
    assert(isValid());
    return constHandler()->tipIndex(containerPtr(), m_index);
  }

  /// Get a mutable reference to the tip index
  /// @return The tip index reference
  TrackIndexType& tipIndex()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->tipIndex(mutableContainerPtr(), m_index);
  }

  /// Get the stem index of the track
  /// @return The stem index
  const TrackIndexType& stemIndex() const {
    assert(isValid());
    return constHandler()->stemIndex(containerPtr(), m_index);
  }

  /// Get a mutable reference to the stem index
  /// @return The stem index reference
  TrackIndexType& stemIndex()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->stemIndex(mutableContainerPtr(), m_index);
  }

  /// Get the index of this track
  /// @return The track index
  TrackIndexType index() const {
    assert(isValid());
    return m_index;
  }

  /// Get the reference surface
  /// @return Reference to the reference surface
  const Surface& referenceSurface() const {
    assert(isValid());
    const Surface* surface =
        constHandler()->referenceSurface(containerPtr(), m_index);
    assert(surface != nullptr);
    return *surface;
  }

  /// Check if track has a reference surface
  /// @return true if a reference surface exists
  bool hasReferenceSurface() const {
    assert(isValid());
    return constHandler()->hasReferenceSurface(containerPtr(), m_index);
  }

  /// Get the particle hypothesis
  /// @return The particle hypothesis
  ParticleHypothesis particleHypothesis() const {
    assert(isValid());
    return constHandler()->particleHypothesis(containerPtr(), m_index);
  }

  /// Get the bound parameters map
  ConstParametersMap parameters() const {
    assert(isValid());
    return constHandler()->parameters(containerPtr(), m_index);
  }

  /// Get the mutable bound parameters map
  ParametersMap parameters()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->parameters(mutableContainerPtr(), m_index);
  }

  /// Get a parameter value by index
  /// @param i The parameter index
  /// @return The parameter value
  const double& parameter(std::size_t i) const { return parameters()[i]; }

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
  ConstCovarianceMap covariance() const {
    assert(isValid());
    return constHandler()->covariance(containerPtr(), m_index);
  }

  /// Get the mutable covariance map
  CovarianceMap covariance()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->covariance(mutableContainerPtr(), m_index);
  }

  /// Get a covariance matrix element
  /// @param i Row index
  /// @param j Column index
  /// @return The covariance element
  const double& covariance(std::size_t i, std::size_t j) const {
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

  /// Get the theta parameter
  /// @return The theta value
  double theta() const {
    assert(isValid());
    return constHandler()->theta(containerPtr(), m_index);
  }

  /// Get the phi parameter
  /// @return The phi value
  double phi() const {
    assert(isValid());
    return constHandler()->phi(containerPtr(), m_index);
  }

  /// Get the q/p parameter
  /// @return The q/p value
  double qOverP() const {
    assert(isValid());
    return constHandler()->qOverP(containerPtr(), m_index);
  }

  /// Get the charge
  /// @return The charge value
  double charge() const {
    assert(isValid());
    return constHandler()->charge(containerPtr(), m_index);
  }

  /// Get the absolute momentum
  /// @return The absolute momentum value
  double absoluteMomentum() const {
    assert(isValid());
    return constHandler()->absoluteMomentum(containerPtr(), m_index);
  }

  /// Get the transverse momentum
  /// @return The transverse momentum value
  double transverseMomentum() const {
    assert(isValid());
    return constHandler()->transverseMomentum(containerPtr(), m_index);
  }

  /// Get the number of measurements
  /// @return The number of measurements
  const unsigned int& nMeasurements() const {
    assert(isValid());
    return constHandler()->nMeasurements(containerPtr(), m_index);
  }

  unsigned int& nMeasurements()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->nMeasurements(mutableContainerPtr(), m_index);
  }

  /// Get the number of holes
  /// @return The number of holes
  const unsigned int& nHoles() const {
    assert(isValid());
    return constHandler()->nHoles(containerPtr(), m_index);
  }

  unsigned int& nHoles()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->nHoles(mutableContainerPtr(), m_index);
  }

  /// Get the number of outliers
  /// @return The number of outliers
  const unsigned int& nOutliers() const {
    assert(isValid());
    return constHandler()->nOutliers(containerPtr(), m_index);
  }

  unsigned int& nOutliers()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->nOutliers(mutableContainerPtr(), m_index);
  }

  /// Get the number of shared hits
  /// @return The number of shared hits
  const unsigned int& nSharedHits() const {
    assert(isValid());
    return constHandler()->nSharedHits(containerPtr(), m_index);
  }

  unsigned int& nSharedHits()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->nSharedHits(mutableContainerPtr(), m_index);
  }

  /// Get the chi2 value
  /// @return The chi2 value
  const float& chi2() const {
    assert(isValid());
    return constHandler()->chi2(containerPtr(), m_index);
  }

  float& chi2()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->chi2(mutableContainerPtr(), m_index);
  }

  /// Get the number of degrees of freedom
  /// @return The number of degrees of freedom
  const unsigned int& nDoF() const {
    assert(isValid());
    return constHandler()->nDoF(containerPtr(), m_index);
  }

  unsigned int& nDoF()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableHandler()->nDoF(mutableContainerPtr(), m_index);
  }

  /// Get the number of track states
  /// @return The number of track states
  unsigned int nTrackStates() const {
    assert(isValid());
    return constHandler()->nTrackStates(containerPtr(), m_index);
  }

  /// Check if the track has a specific dynamic column
  /// @param key The hashed column key
  /// @return true if the column exists
  bool hasColumn(HashedString key) const {
    assert(isValid());
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
    assert(isValid());
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
    assert(isValid());
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
    assert(isValid());
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  /// Get a const reference to a dynamic column component
  /// @tparam T The type of the component
  /// @param key String key for the component to access
  /// @return Const reference to the component
  template <typename T>
  const T& component(HashedString key) const {
    assert(isValid());
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  // ContainerView container() const {
  //   assert(isValid());
  //   return ContainerView{constHandler(), containerPtr()};
  // }

 private:
  template <bool R>
  friend class AnyTrack;

  const detail::TrackHandlerConstBase* constHandler() const {
    return m_handler;
  }

  const detail::TrackHandlerMutableBase* mutableHandler() const
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
using AnyConstTrack = AnyTrack<true>;

/// Alias for mutable type-erased track (though currently all operations are
/// const)
using AnyMutableTrack = AnyTrack<false>;

static_assert(ConstTrackProxyConcept<AnyConstTrack>);
static_assert(MutableTrackProxyConcept<AnyMutableTrack>);

}  // namespace Acts
