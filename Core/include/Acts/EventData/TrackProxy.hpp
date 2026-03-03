// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryBackendConcept.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackContainerBackendConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxyCommon.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <iterator>

namespace Acts {

template <TrackContainerBackend track_container_t,
          CommonMultiTrajectoryBackend traj_t,
          template <typename> class holder_t>
  requires HolderFor<holder_t, track_container_t> && HolderFor<holder_t, traj_t>
class TrackContainer;

template <bool read_only>
class AnyTrack;

/// Proxy class representing a single track.
/// This class provides a **view** into an associated @ref TrackContainer, and
/// has **reference semantics**. You can think of it as a pointer to a vector
/// of tracks, which exposes an object-oriented interface for accessing the
/// track properties.
///
/// @tparam track_container_t the container backend
/// @tparam trajectory_t the track state container backend
/// @tparam holder_t ownership management class for the backend
/// @tparam read_only true if this track container is not mutable
template <typename track_container_t, typename trajectory_t,
          template <typename> class holder_t, bool read_only = true>
class TrackProxy
    : public TrackProxyCommon<
          TrackProxy<track_container_t, trajectory_t, holder_t, read_only>,
          typename track_container_t::IndexType, read_only> {
  using Base = TrackProxyCommon<
      TrackProxy<track_container_t, trajectory_t, holder_t, read_only>,
      typename track_container_t::IndexType, read_only>;

 public:
  /// Indicates whether this track proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// The track container backend given as a template parameter
  using Container = track_container_t;

  /// The track state container backend given as a template parameter
  using Trajectory = trajectory_t;

  /// The index type of the track container
  using IndexType = TrackIndexType;

  /// Sentinel value that indicates an invalid index
  static constexpr IndexType kInvalid = kTrackIndexInvalid;

  /// Alias for the mutable version of this track proxy, with the same backends
  using MutableTrackProxy =
      TrackProxy<track_container_t, trajectory_t, holder_t, false>;

  /// Alias for the const version of this track proxy, with the same backends
  using ConstTrackProxy =
      TrackProxy<track_container_t, trajectory_t, holder_t, true>;

  /// Alias for an associated const track proxy, with the same backends
  using ConstProxyType = ConstTrackProxy;

  /// Alias for an associated mutable track state proxy, with the same backends
  using TrackStateProxy = typename Trajectory::TrackStateProxy;

  /// Alias for an associated const track state proxy, with the same backends
  using ConstTrackStateProxy = typename Trajectory::ConstTrackStateProxy;

  /// Map-type for a bound parameter vector. This has reference semantics, i.e.
  /// points at a matrix by an internal pointer.
  using Parameters =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CoefficientsMap;

  /// Same as @ref Parameters, but with const semantics
  using ConstParameters =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CoefficientsMap;

  /// Map-type for a bound covariance. This has reference semantics, i.e.
  /// points at a matrix by an internal pointer.
  using Covariance =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CovarianceMap;

  /// Same as @ref Covariance, but with const semantics
  using ConstCovariance =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CovarianceMap;

#ifndef DOXYGEN
  friend TrackContainer<Container, Trajectory, holder_t>;
  friend MutableTrackProxy;
  friend ConstTrackProxy;
  // Track proxies are friends, not food!
  template <typename A, typename B, template <typename> class H, bool R>
  friend class TrackProxy;
  template <bool R>
  friend class AnyTrackProxy;
#endif

  /// @anchor track_proxy_construct
  /// @name Constructors and assignment operator
  ///
  /// Public constructors and assignment operators for @c TrackProxy only
  /// allow construction from another @c TrackProxy. You should generally
  /// not have to construct @c TrackProxy manually.
  ///
  /// @{

  /// Copy constructor: const to const or mutable to mutable
  /// @param other the other track proxy
  TrackProxy(const TrackProxy& other) = default;

  /// Copy assignment operator: const to const or mutable to mutable
  /// @param other the other track proxy
  /// @return reference to this track proxy
  TrackProxy& operator=(const TrackProxy& other) = default;

  /// Constructor from mutable track proxy
  /// @note Only available if the track proxy is read-only
  /// @param other the other track proxy
  explicit TrackProxy(const MutableTrackProxy& other)
    requires ReadOnly
      : m_container{other.m_container}, m_index{other.m_index} {}

  /// Copy assignment operator from mutable track proxy
  /// @note Only available if the track proxy is read-only
  /// @param other the other track proxy
  /// @return reference to this track proxy
  TrackProxy& operator=(const MutableTrackProxy& other)
    requires ReadOnly
  {
    m_container = other.m_container;
    m_index = other.m_index;
    return *this;
  }

  /// @}

  /// Equality operator with another track proxy
  /// Checks the container identity and the track index
  /// @param other Other track proxy to compare with
  /// @return True if the track proxies refer to the same track
  bool operator==(const TrackProxy& other) const {
    return &(*m_container) == &(*other.m_container) && m_index == other.m_index;
  }

  /// @anchor track_proxy_props
  /// @name TrackProxy properties
  /// Methods that give access to the properties of a track represented by
  /// @c TrackProxy.
  ///
  /// Many of these methods come in a @c const and a non-@c const version. The
  /// non-@c const version is only available if you have an instance of
  /// @c TrackProxy that does not have the @c read_only template parameter set to
  /// @c true, even if you hold it as an lvalue.
  ///
  /// @{

  /// Get the reference surface of the track (e.g. the perigee)
  /// @return the reference surface
  const Surface& referenceSurface() const {
    return *m_container->container().referenceSurface_impl(m_index);
  }

  // NOLINTBEGIN(performance-unnecessary-value-param)
  // looks like a false-positive. clang-tidy believes `srf` is not movable.
  /// Set a new reference surface for this track
  /// @param srf The surface to set
  void setReferenceSurface(std::shared_ptr<const Surface> srf)
    requires(!ReadOnly)
  {
    m_container->container().setReferenceSurface_impl(m_index, std::move(srf));
  }
  // NOLINTEND(performance-unnecessary-value-param)

  /// Returns whether the track has a reference surface or not
  /// @return whether a surface exists or not
  bool hasReferenceSurface() const {
    // @TODO: This could be more efficient
    return m_container->container().referenceSurface_impl(m_index) != nullptr;
  }

  /// Get the parameters of the track at the reference surface (e.g. perigee).
  /// Const version
  /// @return Proxy vector for the parameters
  ConstParameters parameters() const {
    return m_container->parameters(m_index);
  }

  /// Get the covariance of the track at the reference surface (e.g. perigee).
  /// Const version
  /// @return Proxy matrix for the covariance
  ConstCovariance covariance() const {
    return m_container->covariance(m_index);
  }

  /// Get the parameters of the track at the reference surface (e.g. perigee).
  /// Mutable version
  /// @note Only available if the track proxy is not read-only
  /// @return Proxy vector for the parameters
  Parameters parameters()
    requires(!ReadOnly)
  {
    return m_container->parameters(m_index);
  }

  /// Get the covariance of the track at the reference surface (e.g. perigee).
  /// Mutable version
  /// @note Only available if the track proxy is not read-only
  /// @return Proxy matrix for the covariance
  Covariance covariance()
    requires(!ReadOnly)
  {
    return m_container->covariance(m_index);
  }

  /// Get the particle hypothesis
  /// @return the particle hypothesis
  ParticleHypothesis particleHypothesis() const {
    return m_container->container().particleHypothesis_impl(m_index);
  }

  /// Set a new particle hypothesis for this track
  /// @note Only available if the track proxy is not read-only
  /// @param particleHypothesis The particle hypothesis to set
  void setParticleHypothesis(const ParticleHypothesis& particleHypothesis)
    requires(!ReadOnly)
  {
    m_container->container().setParticleHypothesis_impl(m_index,
                                                        particleHypothesis);
  }

  using Base::absoluteMomentum;
  using Base::charge;
  using Base::chi2;
  using Base::direction;
  using Base::fourMomentum;
  using Base::loc0;
  using Base::loc1;
  using Base::momentum;
  using Base::nDoF;
  using Base::nHoles;
  using Base::nMeasurements;
  using Base::nOutliers;
  using Base::nSharedHits;
  using Base::phi;
  using Base::qOverP;
  using Base::stemIndex;
  using Base::theta;
  using Base::time;
  using Base::tipIndex;
  using Base::transverseMomentum;

  /// Return the number of track states associated to this track
  /// @note This is calculated by iterating over the track states which is
  ///       somewhat expensive. Consider caching this value if you need It
  ///       more than once.
  /// @return The number of track states
  unsigned int nTrackStates() const {
    // @TODO: This should probably be cached, distance is expensive
    //        without random access
    if (tipIndex() == kInvalid) {
      // no tip index -> no track states
      return 0;
    }
    auto tsRange = trackStatesReversed();
    return std::distance(tsRange.begin(), tsRange.end());
  }

  /// Return the index of this track in the track container
  /// @note This is separate from the tip index
  /// @return the track index
  IndexType index() const { return m_index; }

  /// @}

  /// @anchor track_proxy_track_states
  /// @name TrackProxy track state access
  /// Methods that give access to the track states of a track represented by @c TrackProxy.
  /// @{

  /// Return a const track state proxy to the outermost track state
  /// @return The outermost track state proxy
  ConstTrackStateProxy outermostTrackState() const {
    return m_container->trackStateContainer().getTrackState(tipIndex());
  }

  /// Return a mutable track state proxy to the outermost track state
  /// @return The outermost track state proxy
  TrackStateProxy outermostTrackState()
    requires(!ReadOnly)
  {
    return m_container->trackStateContainer().getTrackState(tipIndex());
  }

  /// Return a const track state proxy to the innermost track state
  /// @note This is only available, if the track is forward linked
  /// @return The innermost track state proxy
  auto innermostTrackState() const {
    using proxy_t = decltype(m_container->trackStateContainer().getTrackState(
        std::declval<IndexType>()));

    IndexType stem = component<IndexType, detail_tp::kStemIndexKey>();
    if (stem == kInvalid) {
      return std::optional<proxy_t>{};
    } else {
      return std::optional<proxy_t>{
          m_container->trackStateContainer().getTrackState(stem)};
    }
  }

  /// Return a mutable track state proxy to the innermost track state
  /// @note This is only available, if the track is forward linked
  /// @note Only available if the track proxy is not read-only
  /// @return The innermost track state proxy
  auto innermostTrackState()
    requires(!ReadOnly)
  {
    using proxy_t = decltype(m_container->trackStateContainer().getTrackState(
        std::declval<IndexType>()));

    IndexType stem = component<IndexType>(detail_tp::kStemIndexKey);
    if (stem == kInvalid) {
      return std::optional<proxy_t>{};
    } else {
      return std::optional<proxy_t>{
          m_container->trackStateContainer().getTrackState(stem)};
    }
  }

  /// Get a range over the track states of this track. Return value is
  /// compatible with range based for loop. Const version
  /// @note This range is from the outside inwards!
  /// @return Track state range to iterate over
  auto trackStatesReversed() const {
    return m_container->reverseTrackStateRange(m_index);
  }

  /// Get a range over the track states of this track. Return value is
  /// compatible with range based for loop. Mutable version
  /// @note Only available if the track proxy is not read-only
  /// @note This range is from the outside inwards!
  /// @return Track state range to iterate over
  auto trackStatesReversed()
    requires(!ReadOnly)
  {
    return m_container->reverseTrackStateRange(m_index);
  }

  /// Get a range over the track states of this track. Return value is
  /// compatible with range based for loop. This overload returns a const-only
  /// track state range, which means you cannot modify the track states obtained
  /// in the iteration.
  /// @note This range is from the inside out!
  /// @warning This access direction is only possible if the track states are
  ///          **forward-linked**.
  /// @return Track state range to iterate over
  auto trackStates() const {
    return m_container->forwardTrackStateRange(m_index);
  }

  /// Get a range over the track states of this track.
  /// Return value is compatible with range based for loop.
  /// This overload returns a mutable track state range, which means you
  /// can modify the track states obtained in the iteration.
  /// @note Only available if the track proxy is not read-only
  /// @note This range is from the inside out!
  /// @warning This access direction is only possible if the track states are
  ///          **forward-linked**.
  /// @return Track state range to iterate over
  auto trackStates()
    requires(!ReadOnly)
  {
    return m_container->forwardTrackStateRange(m_index);
  }

  /// @}

  /// @anchor track_proxy_track_state_manipulation
  /// @name TrackProxy track state manipulation
  /// Methods that manipulate the track states of a track represented by @c
  /// TrackProxy.
  ///
  /// **Copy Methods Overview:**
  ///
  /// Three main copy methods are available with different behaviors:
  /// - @c copyFrom(): Deep copy including all track states (creates new track
  ///                 states)
  /// - @c copyFromWithoutStates(): Copy only track properties, invalidate
  ///                               track state indices
  /// - @c copyFromShallow(): Shallow copy sharing the same track states (copy
  ///                         indices only)
  ///
  /// Choose based on your needs:
  /// - Use @c copyFrom() for independent track copies with separate track
  ///   states
  /// - Use @c copyFromWithoutStates() to update track metadata without
  ///   affecting trajectories
  /// - Use @c copyFromShallow() for lightweight copies when track states can
  ///   be shared
  /// @{

  /// Forward connect a track.
  /// This means setting indices from the inside out on all track states.
  /// @note Only available if the track proxy is not read-only
  void linkForward()
    requires(!ReadOnly)
  {
    IndexType last = kInvalid;
    for (auto ts : trackStatesReversed()) {
      ts.template component<IndexType>(detail_tp::kNextKey) = last;
      last = ts.index();
    }
    stemIndex() = last;
  }

  /// Append a track state to this track.
  /// This will modify the tip index to point at the newly created track state,
  /// which will be directly after the previous track state at tip index.
  /// @note Only available if the track proxy is not read-only
  /// @param mask The allocation prop mask for the new track state
  /// @return The newly added track state
  auto appendTrackState(TrackStatePropMask mask = TrackStatePropMask::All)
    requires(!ReadOnly)
  {
    auto& tsc = m_container->trackStateContainer();
    auto ts = tsc.makeTrackState(mask, tipIndex());
    tipIndex() = ts.index();
    return ts;
  }

  /// Create a complete deep copy of another track, including all track states.
  /// This creates new track states in the destination trajectory and copies
  /// all data from the source track states. The track state sequence order
  /// is preserved.
  ///
  /// **Implementation details:**
  /// - Track states are initially copied in reversed order for efficiency
  /// - The track state links are then updated using reverseTrackStates()
  /// - As a consequence, the resulting track is forward-linked
  ///
  /// **What gets copied:**
  /// - All track-level properties (parameters, covariance, particle hypothesis,
  ///   etc.)
  /// - Reference surface (shared pointer is copied)
  /// - Track summary data (nMeasurements, nHoles, chi2, etc.)
  /// - All dynamic track columns
  /// - Complete sequence of track states with all their data
  /// - All dynamic track state columns
  ///
  /// **Result:**
  /// - The destination track will have newly created track states
  /// - tipIndex() and stemIndex() will point to the new track
  /// states
  /// - Track state indices will be different from the source
  /// - All track state data will be identical to the source
  /// - The track will be forward-linked (stemIndex() will be valid)
  ///
  /// @note Only available if the track proxy is not read-only
  /// @note Both track containers must have compatible dynamic columns
  /// @tparam track_proxy_t the other track proxy's type
  /// @param other The source track proxy to copy from
  template <TrackProxyConcept track_proxy_t>
  void copyFrom(const track_proxy_t& other)
    requires(!ReadOnly)
  {
    copyFromWithoutStates(other);

    // append track states (cheap), but they're in the wrong order
    for (const auto& srcTrackState : other.trackStatesReversed()) {
      auto destTrackState = appendTrackState(srcTrackState.getMask());
      destTrackState.copyFrom(srcTrackState, Acts::TrackStatePropMask::All,
                              true);
    }

    // reverse using standard linked list reversal algorithm
    reverseTrackStates();
  }

  /// Copy track-level properties from another track, but not the track states.
  /// This copies all track metadata and properties but leaves the track state
  /// sequence unchanged. Useful when you want to copy track properties to an
  /// existing track that may already have track states.
  ///
  /// **What gets copied:**
  /// - Track parameters at reference surface
  /// - Covariance matrix at reference surface
  /// - Particle hypothesis
  /// - Reference surface (shared pointer is copied)
  /// - Track summary data (nMeasurements, nHoles, nOutliers, nSharedHits, chi2,
  /// nDoF)
  /// - All dynamic track columns
  ///
  /// **What does NOT get copied:**
  /// - Track states (existing track states remain unchanged in the container)
  ///
  /// **Result:**
  /// - All track-level properties are updated to match the source
  /// - tipIndex() and stemIndex() are set to kInvalid (track states
  /// become inaccessible)
  /// - Existing track states remain in the container but are no longer linked
  /// to this track
  /// - nTrackStates() will return 0 due to invalid indices
  ///
  /// @note Only available if the track proxy is not read-only
  /// @note Both track containers must have compatible dynamic columns
  /// @tparam track_proxy_t the other track proxy's type
  /// @param other The source track proxy to copy properties from
  template <TrackProxyConcept track_proxy_t>
  void copyFromWithoutStates(const track_proxy_t& other)
    requires(!ReadOnly)
  {
    setParticleHypothesis(other.particleHypothesis());

    if (other.hasReferenceSurface()) {
      setReferenceSurface(other.referenceSurface().getSharedPtr());
      parameters() = other.parameters();
      covariance() = other.covariance();
    } else {
      setReferenceSurface(nullptr);
    }

    nMeasurements() = other.nMeasurements();
    nHoles() = other.nHoles();
    nOutliers() = other.nOutliers();
    nSharedHits() = other.nSharedHits();
    chi2() = other.chi2();
    nDoF() = other.nDoF();

    m_container->copyDynamicFrom(m_index, other.m_container->container(),
                                 other.m_index);

    tipIndex() = kInvalid;
    stemIndex() = kInvalid;
  }

  /// Create a shallow copy from another track, sharing the same track states.
  /// This copies all track-level properties and makes the destination track
  /// point to the same track state sequence as the source. The track states
  /// themselves are not duplicated - both tracks will reference the same
  /// track state objects in memory.
  ///
  /// **What gets copied:**
  /// - All track-level properties (parameters, covariance, particle hypothesis,
  /// etc.)
  /// - Reference surface (shared pointer is copied)
  /// - Track summary data (nMeasurements, nHoles, chi2, etc.)
  /// - All dynamic track columns
  /// - tipIndex() and stemIndex() (track state linking information)
  ///
  /// **What gets shared (not duplicated):**
  /// - Track states (both tracks reference the same track state objects)
  ///
  /// **Result:**
  /// - The destination track will have the same nTrackStates() as the source
  /// - Both tracks will iterate over the same track state sequence
  /// - Modifications to track states will be visible in both tracks
  /// - Track state indices will be identical between tracks
  /// - The destination track will have a different track index than the source
  ///
  /// @warning Modifying track states through either track will affect both tracks
  ///          since they share the same track state objects
  /// @warning It is the user's responsibility to ensure that the tip and stem
  ///          indices from the source track are valid in the destination
  ///          track's track state container. No validation is performed -
  ///          invalid indices will lead to undefined behavior when accessing
  ///          track states
  /// @note Only available if the track proxy is not read-only
  /// @note Both track containers must have compatible dynamic columns
  /// @tparam track_proxy_t the other track proxy's type
  /// @param other The source track proxy to create a shallow copy from
  template <TrackProxyConcept track_proxy_t>
  void copyFromShallow(const track_proxy_t& other)
    requires(!ReadOnly)
  {
    copyFromWithoutStates(other);
    tipIndex() = other.tipIndex();
    stemIndex() = other.stemIndex();
  }

  /// Reverse the ordering of track states for this track
  /// Afterwards, the previous endpoint of the track state sequence will be
  /// the "innermost" track state
  /// @note Only available if the track proxy is not read-only
  /// @note This is dangerous with branching track state sequences, as it will break them
  /// @note This also automatically forward-links the track!
  /// @param invertJacobians Whether to invert the Jacobians of the track states
  void reverseTrackStates(bool invertJacobians = false)
    requires(!ReadOnly)
  {
    IndexType current = tipIndex();
    IndexType next = kInvalid;
    IndexType prev = kInvalid;

    stemIndex() = tipIndex();

    // @TODO: Maybe refactor to not need this variable if invertJacobians == false
    BoundMatrix nextJacobian;

    while (current != kInvalid) {
      auto ts = m_container->trackStateContainer().getTrackState(current);
      prev = ts.previous();
      ts.template component<IndexType>(detail_tp::kNextKey) = prev;
      ts.previous() = next;
      if (invertJacobians) {
        if (next != kInvalid) {
          BoundMatrix curJacobian = ts.jacobian();
          ts.jacobian() = nextJacobian.inverse();
          nextJacobian = curJacobian;
        } else {
          nextJacobian = ts.jacobian();
          ts.jacobian().setZero();
        }
      }
      next = current;
      tipIndex() = current;
      current = prev;
    }
  }

  /// @}

  /// @anchor track_proxy_generic_component
  /// @name TrackProxy generic component access
  /// Methods that give access to generic components of a track represented by
  /// @c TrackProxy.  Internally, a compile-time hash of the component name is
  /// used to identify which component is being requested. Most of the named
  /// methods in @ref track_proxy_props "TrackProxy properties" use these
  /// methods to retrieve the actual data.
  ///
  /// A number of overloads exist, where you can either supply the
  /// @ref HashedString @c key as a template parameter or a runtime argument.  The
  /// former has the advantage of being guaranteed to be evaluated at
  /// compile-time.
  ///
  /// @{

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr T& component()
    requires(!ReadOnly)
  {
    return m_container->template component<T, key>(m_index);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T>
  constexpr T& component(HashedString key)
    requires(!ReadOnly)
  {
    return m_container->template component<T>(key, m_index);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return Mutable reference to the component given by @p key
  template <typename T>
  constexpr T& component(std::string_view key)
    requires(!ReadOnly)
  {
    return m_container->template component<T>(hashStringDynamic(key), m_index);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Const reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr const T& component() const {
    return m_container->template component<T, key>(m_index);
  }

  /// Check whether a dynamic column exists
  /// @param key String key for the component to check
  /// @return whether the column exists
  bool hasColumn(HashedString key) const { return m_container->hasColumn(key); }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(HashedString key) const {
    return m_container->template component<T>(key, m_index);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(std::string_view key) const {
    return m_container->template component<T>(hashStringDynamic(key), m_index);
  }

  /// @}

  /// Return the track parameters at the reference surface
  /// @note The parameters are created on the fly
  /// @return the track parameters
  BoundTrackParameters createParametersAtReference() const {
    return BoundTrackParameters(referenceSurface().getSharedPtr(), parameters(),
                                covariance(), particleHypothesis());
  }

  /// Convert a track state into track parameters
  /// @note The parameters are created on the fly
  /// @param trackState Track state to convert to parameters
  /// @return the track parameters
  BoundTrackParameters createParametersFromState(
      const ConstTrackStateProxy& trackState) const {
    return BoundTrackParameters(trackState.referenceSurface().getSharedPtr(),
                                trackState.parameters(),
                                trackState.covariance(), particleHypothesis());
  }

  /// Return a reference to the track container backend, mutable version.
  /// @note Only available if the track proxy is not read-only
  /// @return reference to the track container backend
  auto& container()
    requires(!ReadOnly)
  {
    return *m_container;
  }

  /// Return a reference to the track container backend, const version.
  /// @return reference to the track container backend
  const auto& container() const { return *m_container; }

 private:
  TrackProxy(
      const_if_t<ReadOnly, TrackContainer<Container, Trajectory, holder_t>>&
          container,
      IndexType itrack)
      : m_container{&container}, m_index{itrack} {}

  detail_lt::TransitiveConstPointer<
      const_if_t<ReadOnly, TrackContainer<Container, Trajectory, holder_t>>>
      m_container;
  IndexType m_index;
};

}  // namespace Acts
