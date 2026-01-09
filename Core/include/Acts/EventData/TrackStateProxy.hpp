// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxyCommon.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/EigenConcepts.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cstddef>
#include <ranges>
#include <span>

#include <Eigen/Core>

namespace Acts {

template <typename derived_t>
class MultiTrajectory;

template <bool read_only>
class AnyTrackStateProxy;

namespace detail_lt {

/// Helper type to make a member pointers constness transitive.
template <typename T>
class TransitiveConstPointer {
 public:
  using element_type = T;
  TransitiveConstPointer() = default;
  explicit TransitiveConstPointer(T* ptr) : m_ptr{ptr} {}

  template <typename U>
  explicit TransitiveConstPointer(const TransitiveConstPointer<U>& other)
      : m_ptr{other.ptr()} {}

  template <typename U>
  TransitiveConstPointer& operator=(const TransitiveConstPointer<U>& other) {
    m_ptr = other.m_ptr;
    return *this;
  }

  template <typename U>
  bool operator==(const TransitiveConstPointer<U>& other) const {
    return m_ptr == other.m_ptr;
  }

  const T* operator->() const { return m_ptr; }

  T* operator->() { return m_ptr; }

  template <typename U>
  friend class TransitiveConstPointer;
  template <bool R>
  friend class AnyTrackStateProxy;

  const T& operator*() const { return *m_ptr; }

  T& operator*() { return *m_ptr; }

  explicit operator bool() const { return m_ptr != nullptr; }

  const T* ptr() const { return m_ptr; }
  T* ptr() { return m_ptr; }

 private:
  T* m_ptr{nullptr};
};

}  // namespace detail_lt

/// Proxy object to access a single point on the trajectory.
///
/// @tparam SourceLink Type to link back to an original measurement
/// @tparam M          Maximum number of measurement dimensions
/// @tparam read_only  true for read-only access to underlying storage
template <typename trajectory_t, std::size_t M, bool read_only = true>
class TrackStateProxy
    : public TrackStateProxyCommon<TrackStateProxy<trajectory_t, M, read_only>,
                                   read_only> {
  using Base =
      TrackStateProxyCommon<TrackStateProxy<trajectory_t, M, read_only>,
                            read_only>;

  friend class TrackStateProxyCommon<
      TrackStateProxy<trajectory_t, M, read_only>, read_only>;

 public:
  /// Indicates whether this track state proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Alias for an associated const track state proxy, with the same backends
  using ConstProxyType = TrackStateProxy<trajectory_t, M, true>;

  /// Map-type for a bound parameter vector. This has reference semantics, i.e.
  /// points at a matrix by an internal pointer.
  using Parameters = typename TrackStateTraits<M, false>::Parameters;

  /// Same as @ref Parameters, but with const semantics
  using ConstParameters = typename TrackStateTraits<M, true>::Parameters;

  /// Map-type for a bound covariance. This has reference semantics, i.e.
  /// points at a matrix by an internal pointer.
  using Covariance = typename TrackStateTraits<M, false>::Covariance;

  /// Jacobian shape is identical to @ref Covariance
  using Jacobian = typename TrackStateTraits<M, false>::Covariance;

  /// Same as @ref Covariance, but with const semantics
  using ConstCovariance = typename TrackStateTraits<M, true>::Covariance;

  /// Jacobian shape is identical to @ref Covariance
  using ConstJacobian = typename TrackStateTraits<M, true>::Covariance;

  /// Map-type for a calibrated measurement vector, where the local measurement
  /// dimension is variable.
  template <std::size_t N>
  using Calibrated = typename TrackStateTraits<N, false>::Calibrated;

  /// Same as @c Calibrated, but with const semantics
  template <std::size_t N>
  using ConstCalibrated = typename TrackStateTraits<N, true>::Calibrated;

  /// Map-type for a calibrated measurement covariance matrix, where the local
  /// measurement dimension is variable.
  template <std::size_t N>
  using CalibratedCovariance =
      typename TrackStateTraits<N, false>::CalibratedCovariance;

  /// Same as @ref CalibratedCovariance, but with const semantics
  template <std::size_t N>
  using ConstCalibratedCovariance =
      typename TrackStateTraits<N, true>::CalibratedCovariance;

  /// Map-type for a measurement vector, where the local measurement dimension
  /// is variable.
  using EffectiveCalibrated =
      typename TrackStateTraits<M, false>::EffectiveCalibrated;

  /// Same as @c EffectiveCalibrated, but with const semantics
  using ConstEffectiveCalibrated =
      typename TrackStateTraits<M, true>::EffectiveCalibrated;

  /// Map-type for a measurement covariance matrix, where the local measurement
  /// dimension is variable.
  using EffectiveCalibratedCovariance =
      typename TrackStateTraits<M, false>::EffectiveCalibratedCovariance;

  /// Same as @ref EffectiveCalibratedCovariance, but with const semantics
  using ConstEffectiveCalibratedCovariance =
      typename TrackStateTraits<M, true>::EffectiveCalibratedCovariance;

  /// The index type of the track state container
  using IndexType = TrackIndexType;

  /// Sentinel value that indicates an invalid index
  static constexpr IndexType kInvalid = kTrackIndexInvalid;

  /// The track state container backend given as a template parameter
  using Trajectory = trajectory_t;

  using Base::allocateCalibrated;
  using Base::calibrated;
  using Base::calibratedCovariance;
  using Base::chi2;
  using Base::covariance;
  using Base::effectiveCalibrated;
  using Base::effectiveCalibratedCovariance;
  using Base::filtered;
  using Base::filteredCovariance;
  using Base::getMask;
  using Base::hasCalibrated;
  using Base::hasFiltered;
  using Base::hasJacobian;
  using Base::hasPredicted;
  using Base::hasPrevious;
  using Base::hasProjector;
  using Base::hasSmoothed;
  using Base::parameters;
  using Base::pathLength;
  using Base::predicted;
  using Base::predictedCovariance;
  using Base::previous;
  using Base::projectorSubspaceHelper;
  using Base::projectorSubspaceIndices;
  using Base::setProjectorSubspaceIndices;
  using Base::smoothed;
  using Base::smoothedCovariance;
  using Base::typeFlags;

  /// @anchor track_state_proxy_construct
  /// @name Constructors and assignment operator
  ///
  /// Public constructors and assignment operators for @c TrackStateProxy only
  /// allow construction from another @c TrackStateProxy. You should generally
  /// not have to construct @c TrackStateProxy manually.
  ///
  /// @{

  /// Copy constructor: const to const or mutable to mutable
  /// @param other The other TrackStateProxy to construct from
  TrackStateProxy(const TrackStateProxy& other) = default;

  /// Copy assignment operator: const to const or mutable to mutable
  /// @param other The other TrackStateProxy to assign from
  /// @return Reference to this TrackStateProxy
  TrackStateProxy& operator=(const TrackStateProxy& other) = default;

  /// Constructor from mutable TrackStateProxy
  /// @note Only available if the track state proxy is read-only
  /// @param other The other TrackStateProxy to construct from
  explicit TrackStateProxy(const TrackStateProxy<Trajectory, M, false>& other)
    requires ReadOnly
      : m_traj{other.m_traj}, m_istate{other.m_istate} {}

  /// Assignment operator to from mutable @c TrackStateProxy
  /// @param other The other TrackStateProxy to assign from
  /// @note Only available if the track state proxy is read-only
  /// @return Reference to this TrackStateProxy
  TrackStateProxy& operator=(const TrackStateProxy<Trajectory, M, false>& other)
    requires ReadOnly
  {
    m_traj = other.m_traj;
    m_istate = other.m_istate;

    return *this;
  }

  /// @}

  /// @anchor track_state_proxy_props
  /// @name Track state properties
  ///
  /// Properties of the track state represented by @c TrackStateProxy.
  ///
  /// Many of these methods come in a @c const and a non-@c const version. The
  /// non-@c const version is only available if you have an instance of
  /// @c TrackStateProxy that does not have the @c read_only template parameter set to
  /// @c true, even if you hold it as an lvalue.
  ///
  /// The track states each have an index in the track state container. The
  /// sequence of track states is implemented as a one or two-way linked list,
  /// which uses indices into the same container.
  ///
  /// Each track state has a @c previous index, which points at the track state
  /// immediately preceding. A track state with a @c previous index of @c
  /// kInvalid is the first (innermost) track state in a track or track
  /// candidate. This is also referred to as a *stem* at the track level.
  ///
  /// During track finding and fitting, track states are usually appended to the
  /// sequence, populating the @c previous index of the new track state. Combinatorial
  /// track finding can produce track states which fork in this way, by having
  /// more than one track state with the same @c previous index.
  ///
  /// The track states have static, optional and dynamic properties. Static
  /// properties are always present, and can always be retrieved. Optional
  /// components use an extra indirection mechanism that coordinates with the
  /// backend to allow both not having the component set, or sharing it with
  /// other track states. An example is a branching trajectory from track
  /// finding which shares the same predicted parameter vector and associated
  /// covariance.
  ///
  /// Optional components are
  /// - predicted parameters and covariance
  /// - filtered parameters and covariance
  /// - smoothed parameters and covariance
  /// - jacobian
  /// - calibrated measurement info including projector
  ///
  /// They can be unset via @ref unset, @ref getMask can be used to check which
  /// components are present. The first four are shareable between track
  /// states via @ref shareFrom.
  ///
  /// @{

  /// Index within the trajectory.
  /// @return the index
  IndexType index() const { return m_istate; }

  /// Unset an optional track state component
  /// @note Only available if the track state proxy is not read-only
  /// @param target The component to unset
  void unset(TrackStatePropMask target)
    requires(!ReadOnly)
  {
    m_traj->self().unset(target, m_istate);
  }

  /// Add additional components to the track state
  /// @note Only available if the track state proxy is not read-only
  /// @param mask The bitmask that instructs which components to allocate
  void addComponents(TrackStatePropMask mask)
    requires(!ReadOnly)
  {
    m_traj->self().addTrackStateComponents_impl(m_istate, mask);
  }

  /// Reference surface.
  /// @return the reference surface
  const Surface& referenceSurface() const {
    assert(hasReferenceSurface() &&
           "TrackState does not have reference surface");
    return *m_traj->referenceSurface(m_istate);
  }

  /// Returns if the track state has a non nullptr surface associated
  /// @return whether a surface exists or not
  bool hasReferenceSurface() const {
    return m_traj->referenceSurface(m_istate) != nullptr;
  }

  // NOLINTBEGIN(performance-unnecessary-value-param)
  // looks like a false-positive. clang-tidy believes `srf` is not movable.

  /// Set the reference surface to a given value
  /// @param srf Shared pointer to the surface to set
  /// @note This overload is only present in case @c ReadOnly is false.
  void setReferenceSurface(std::shared_ptr<const Surface> srf)
    requires(!ReadOnly)
  {
    m_traj->setReferenceSurface(m_istate, std::move(srf));
  }
  // NOLINTEND(performance-unnecessary-value-param)

  /// @}

  /// Returns the jacobian from the previous trackstate to this one
  /// @return The jacobian matrix
  /// @note Const version
  ConstCovariance jacobian() const {
    assert(has<detail_tsp::kJacobianKey>());
    return m_traj->self().jacobian(m_istate);
  }

  /// Returns the jacobian from the previous trackstate to this one
  /// @return The jacobian matrix
  /// @note Mutable version
  Covariance jacobian()
    requires(!ReadOnly)
  {
    assert(has<detail_tsp::kJacobianKey>());
    return m_traj->self().jacobian(m_istate);
  }

  /// @anchor track_state_proxy_meas
  /// @name Track state measurement properties
  ///
  /// Properties of the measurement associated with the track state represented.
  /// This consists of a vector and an associated square matrix of a measurement
  /// dimension which is between one and the size of the track parametrization.
  /// The measurement coordinate frame is required to be a strict subset of the
  /// bound track parametrization on the local geometry coordinate frame, i.e.
  /// using a pure projector matrix to convert from the bound parametrization to
  /// the measurement frame is possible.
  ///
  /// The track state stores the parameter vector and covariance, and the
  /// backend is given the possibility to do so in a jagged way, i.e. only
  /// storing the number of values needed. This requires calling
  /// @ref allocateCalibrated before storing the measurements
  /// (even if it might be a no-op).
  ///
  /// The projector matrix is packed as a bitset, which is converted to a matrix
  /// on-demand (and therefore returned by value).
  ///
  /// The track state also includes a @ref SourceLink which acts as a proxy
  /// to the original uncalibrated measurement that the calibrated measurement
  /// was derived from. It is set and returned by value, to allow unpacking /
  /// repacking by the backend, if needed.
  ///
  /// @{

  /// Uncalibrated measurement in the form of a source link. Const version
  /// @return The uncalibrated measurement source link
  SourceLink getUncalibratedSourceLink() const;

  /// Set an uncalibrated source link
  /// @param sourceLink The uncalibrated source link to set
  void setUncalibratedSourceLink(SourceLink&& sourceLink)
    requires(!ReadOnly)
  {
    m_traj->setUncalibratedSourceLink(m_istate, std::move(sourceLink));
  }

  /// Check if the point has an associated uncalibrated measurement.
  /// @return Whether it is set
  bool hasUncalibratedSourceLink() const {
    return has<detail_tsp::kUncalibratedKey>();
  }

  /// Return the (dynamic) number of dimensions stored for this measurement.
  /// @note Depending on the backend, this size is used to determine the
  ///       memory range of the measurement vector and covariance.
  /// @return The number of dimensions
  IndexType calibratedSize() const { return m_traj->calibratedSize(m_istate); }

  /// Allocate storage to be able to store a measurement of size @p measdim.
  /// This must be called **before** setting the measurement content.
  /// @param measdim Number of measurement dimensions to allocate
  /// @note This does not allocate if an allocation of the same size already exists
  /// @note This will zero-initialize the allocated storage
  /// @note This is an error if an existing allocation has different size
  void allocateCalibrated(std::size_t measdim)
    requires(!ReadOnly)
  {
    m_traj->allocateCalibrated(m_istate, measdim);
  }

  /// @}

  /// @anchor track_state_share_copy
  /// @name Sharing and copying
  ///
  /// Methods to share and copy track state components. Sharing means setting up
  /// more than one track state to point to the same component.
  ///
  /// Shareable components are
  /// - predicted parameters and covariance
  /// - filtered parameters and covariance
  /// - smoothed parameters and covariance
  /// - jacobian
  ///
  /// See @ref TrackStatePropMask.
  ///
  /// @{

  /// Share a shareable component **within** this track state
  /// @param shareSource Which component to share from
  /// @param shareTarget Which component to share as. This should be different from
  ///                    as @p shareSource, e.g. predicted can be shared as filtered.
  void shareFrom(TrackStatePropMask shareSource, TrackStatePropMask shareTarget)
    requires(!ReadOnly)
  {
    shareFrom(*this, shareSource, shareTarget);
  }

  /// Share a shareable component from another track state.
  /// @param other Track state proxy to share component from
  /// @param component Which component to share.
  /// @note The track states both need to be stored in the
  ///       same @c MultiTrajectory instance
  template <bool ReadOnlyOther>
  void shareFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                 TrackStatePropMask component)
    requires(!ReadOnly)
  {
    shareFrom(other, component, component);
  }

  /// Share a shareable component from another track state
  /// @param other Track state proxy to share component(s) from
  /// @param shareSource Which component to share from
  /// @param shareTarget Which component to share as. This can be be different from
  ///                    as @p shareSource, e.g. predicted can be shared as filtered.
  /// @note Shareable components are predicted, filtered, smoothed, calibrated, jacobian,
  ///       or projector. See @c TrackStatePropMask.
  template <bool ReadOnlyOther>
  void shareFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                 TrackStatePropMask shareSource, TrackStatePropMask shareTarget)
    requires(!ReadOnly)
  {
    assert(m_traj == other.m_traj &&
           "Cannot share components across MultiTrajectories");

    assert(ACTS_CHECK_BIT(other.getMask(), shareSource) &&
           "Source has incompatible allocation");

    m_traj->self().shareFrom(m_istate, other.m_istate, shareSource,
                             shareTarget);
  }

  /// Copy the contents of another track state proxy into this one
  /// @param other The other track state to copy from
  /// @param mask An optional mask to determine what to copy from
  /// @param onlyAllocated Whether to only copy allocated components
  /// @note If the this track state proxy does not have compatible allocations
  ///       with the source track state proxy, and @p onlyAllocated is false,
  ///       an exception is thrown.
  /// @note The mask parameter will not cause a copy of components that are
  ///       not allocated in the source track state proxy.
  template <TrackStateProxyConcept track_state_proxy_t>
  void copyFrom(const track_state_proxy_t& other,
                TrackStatePropMask mask = TrackStatePropMask::All,
                bool onlyAllocated = true)
    requires(!ReadOnly)
  {
    using PM = TrackStatePropMask;

    if (onlyAllocated) {
      auto dest = getMask();
      auto src = other.getMask() &
                 mask;  // combine what we have with what we want to copy

      if (ACTS_CHECK_BIT(src, PM::Calibrated)) {
        // on-demand allocate calibrated
        dest |= PM::Calibrated;
      }

      if ((static_cast<std::underlying_type_t<TrackStatePropMask>>(
               (src ^ dest) & src) != 0 ||
           dest == TrackStatePropMask::None ||
           src == TrackStatePropMask::None) &&
          mask != TrackStatePropMask::None) {
        throw std::runtime_error(
            "Attempt track state copy with incompatible allocations");
      }

      // we're sure now this has correct allocations, so just copy
      if (ACTS_CHECK_BIT(src, PM::Predicted)) {
        predicted() = other.predicted();
        predictedCovariance() = other.predictedCovariance();
      }

      if (ACTS_CHECK_BIT(src, PM::Filtered)) {
        filtered() = other.filtered();
        filteredCovariance() = other.filteredCovariance();
      }

      if (ACTS_CHECK_BIT(src, PM::Smoothed)) {
        smoothed() = other.smoothed();
        smoothedCovariance() = other.smoothedCovariance();
      }

      if (other.hasUncalibratedSourceLink()) {
        setUncalibratedSourceLink(other.getUncalibratedSourceLink());
      }

      if (ACTS_CHECK_BIT(src, PM::Jacobian)) {
        jacobian() = other.jacobian();
      }

      if (ACTS_CHECK_BIT(src, PM::Calibrated)) {
        visit_measurement(other.calibratedSize(), [&](auto N) {
          constexpr int measdim = decltype(N)::value;
          allocateCalibrated(
              other.template calibrated<measdim>().eval(),
              other.template calibratedCovariance<measdim>().eval());
        });

        setProjectorSubspaceIndices(other.projectorSubspaceIndices());
      }
    } else {
      if (ACTS_CHECK_BIT(mask, PM::Predicted) &&
          has<detail_tsp::kPredictedKey>() &&
          other.template has<detail_tsp::kPredictedKey>()) {
        predicted() = other.predicted();
        predictedCovariance() = other.predictedCovariance();
      }

      if (ACTS_CHECK_BIT(mask, PM::Filtered) &&
          has<detail_tsp::kFilteredKey>() &&
          other.template has<detail_tsp::kFilteredKey>()) {
        filtered() = other.filtered();
        filteredCovariance() = other.filteredCovariance();
      }

      if (ACTS_CHECK_BIT(mask, PM::Smoothed) &&
          has<detail_tsp::kSmoothedKey>() &&
          other.template has<detail_tsp::kSmoothedKey>()) {
        smoothed() = other.smoothed();
        smoothedCovariance() = other.smoothedCovariance();
      }

      if (other.hasUncalibratedSourceLink()) {
        setUncalibratedSourceLink(other.getUncalibratedSourceLink());
      }

      if (ACTS_CHECK_BIT(mask, PM::Jacobian) &&
          has<detail_tsp::kJacobianKey>() &&
          other.template has<detail_tsp::kJacobianKey>()) {
        jacobian() = other.jacobian();
      }

      // NOTE: we should not check hasCalibrated on this, since it
      // may be not yet allocated
      if (ACTS_CHECK_BIT(mask, PM::Calibrated) &&
          other.template has<detail_tsp::kCalibratedKey>()) {
        visit_measurement(other.calibratedSize(), [&](auto N) {
          constexpr int measdim = decltype(N)::value;
          allocateCalibrated(
              other.template calibrated<measdim>().eval(),
              other.template calibratedCovariance<measdim>().eval());
        });

        setProjectorSubspaceIndices(other.projectorSubspaceIndices());
      }
    }

    chi2() = other.chi2();
    pathLength() = other.pathLength();
    typeFlags() = other.typeFlags();

    if (other.hasReferenceSurface()) {
      setReferenceSurface(other.referenceSurface().getSharedPtr());
    }

    m_traj->copyDynamicFrom(m_istate, other.container(), other.index());
  }

  /// @}

  /// @anchor track_state_proxy_generic_component
  /// @name Track state proxy Generic component access
  /// @{

  /// Check if a component is set
  /// @tparam key Hashed string key to check for
  /// @return true if the component exists, false if not
  template <HashedString key>
  constexpr bool has() const {
    return m_traj->template has<key>(m_istate);
  }

  /// Check if a component is set
  /// @param key Hashed string key to check for
  /// @return true if the component exists, false if not
  constexpr bool has(HashedString key) const {
    return m_traj->has(key, m_istate);
  }

  /// Check if a component is set
  /// @param key String key to check for
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return true if the component exists, false if not
  constexpr bool has(std::string_view key) const {
    return has(hashStringDynamic(key));
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr T& component()
    requires(!ReadOnly)
  {
    return m_traj->template component<T, key>(m_istate);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T>
  constexpr T& component(HashedString key)
    requires(!ReadOnly)
  {
    return m_traj->template component<T>(key, m_istate);
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
    return m_traj->template component<T>(hashStringDynamic(key), m_istate);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Const reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr const T& component() const {
    return m_traj->template component<T, key>(m_istate);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(HashedString key) const {
    return m_traj->template component<T>(key, m_istate);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(std::string_view key) const {
    return m_traj->template component<T>(hashStringDynamic(key), m_istate);
  }

  /// @}

  /// Return a mutable reference to the underlying backend container
  /// @return A reference to the backend container
  MultiTrajectory<Trajectory>& trajectory()
    requires(!ReadOnly)
  {
    return *m_traj;
  }

  /// Return a const reference to the underlying backend container
  /// @return A const reference to the backend container
  const MultiTrajectory<Trajectory>& trajectory() const { return *m_traj; }

  /// Get a mutable reference to the track state container backend
  /// @return a mutable reference to the backend
  auto& container()
    requires(!ReadOnly)
  {
    return *m_traj;
  }

  /// Get a const reference to the track state container backend
  /// @return a const reference to the backend
  const auto& container() const { return *m_traj; }

  /// Check if the track state has a specific dynamic column
  /// @param key The hashed column key
  /// @return true if the column exists
  bool hasColumn(HashedString key) const { return container().hasColumn(key); }

 protected:
  ConstParameters parametersAtIndex(IndexType parIndex) const {
    return m_traj->parameters(parIndex);
  }

  Parameters parametersAtIndexMutable(IndexType parIndex)
    requires(!ReadOnly)
  {
    return m_traj->parameters(parIndex);
  }

  ConstCovariance covarianceAtIndex(IndexType covIndex) const {
    return m_traj->covariance(covIndex);
  }

  Covariance covarianceAtIndexMutable(IndexType covIndex)
    requires(!ReadOnly)
  {
    return m_traj->covariance(covIndex);
  }

  double* calibratedDataMutable()
    requires(!ReadOnly)
  {
    return m_traj->template calibrated<M>(m_istate).data();
  }

  const double* calibratedData() const {
    return m_traj->template calibrated<M>(m_istate).data();
  }

  double* calibratedCovarianceDataMutable()
    requires(!ReadOnly)
  {
    return m_traj->template calibratedCovariance<M>(m_istate).data();
  }

  const double* calibratedCovarianceData() const {
    return m_traj->template calibratedCovariance<M>(m_istate).data();
  }

 private:
  // Private since it can only be created by the trajectory.
  TrackStateProxy(const_if_t<ReadOnly, MultiTrajectory<Trajectory>>& trajectory,
                  IndexType istate);

  detail_lt::TransitiveConstPointer<
      const_if_t<ReadOnly, MultiTrajectory<Trajectory>>>
      m_traj;
  IndexType m_istate;

  friend class Acts::MultiTrajectory<Trajectory>;
  friend class TrackStateProxy<Trajectory, M, true>;
  friend class TrackStateProxy<Trajectory, M, false>;
  template <bool R>
  friend class AnyTrackStateProxy;

  const auto* rawTrajectoryPtr() const { return m_traj.ptr(); }
  auto* rawTrajectoryPtr() { return m_traj.ptr(); }
};
}  // namespace Acts

#include "Acts/EventData/TrackStateProxy.ipp"
