// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <bitset>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

namespace Acts {

/// @enum TrackStateFlag
///
/// This enum describes the type of TrackState
enum TrackStateFlag {
  MeasurementFlag = 0,
  ParameterFlag = 1,
  OutlierFlag = 2,
  HoleFlag = 3,
  MaterialFlag = 4,
  SharedHitFlag = 5,
  NumTrackStateFlags = 6
};

class ConstTrackStateType;

/// View type over a bitset stored in a 64 bit integer
/// This view allows modifications.
class TrackStateType {
 public:
  using raw_type = std::uint64_t;
  static constexpr std::size_t kRawBits =
      std::numeric_limits<std::make_unsigned<raw_type>::type>::digits;
  /// Constructor from a reference to the underlying value container
  /// @param raw the value container
  TrackStateType(raw_type& raw) : m_raw{&raw} {}

  /// Assign the value from another set of flags
  /// @param other the other set of flags to assign
  /// @return this object
  TrackStateType& operator=(const TrackStateType& other) {
    assert(other.m_raw != nullptr);
    *m_raw = *other.m_raw;
    return *this;
  }

  /// Assign the value from another set of flags
  /// @param other the other set of flags to assign
  /// @return this object
  TrackStateType& operator=(const ConstTrackStateType& other);

  /// Automatically convert to const track state type
  operator ConstTrackStateType();

  /// Return if the bit at position @p pos is 1
  /// @param pos the bit position
  /// @return if the bit at @p pos is one or not
  bool test(std::size_t pos) const {
    assert(m_raw != nullptr);
    std::bitset<kRawBits> bs{*m_raw};
    return bs.test(pos);
  }

  /// Change the value of the bit at position @p pos to @p value.
  /// @param pos the position of the bit to change
  /// @param value the value to change the bit to
  void set(std::size_t pos, bool value = true) {
    assert(m_raw != nullptr);
    std::bitset<kRawBits> bs{*m_raw};
    bs.set(pos, value);
    *m_raw = bs.to_ullong();
  }

  /// Change the value of the bit at position at @p pos to @c false
  /// @param pos the position of the bit to change
  void reset(std::size_t pos) { set(pos, false); }

 private:
  raw_type* m_raw{nullptr};
};

/// View type over a bitset stored in a 64 bit integer
/// This view does not allow modifications
class ConstTrackStateType {
 public:
  using raw_type = std::uint64_t;
  static constexpr std::size_t kRawBits =
      std::numeric_limits<std::make_unsigned<raw_type>::type>::digits;

  /// Constructor from a reference to the underlying value container
  /// @param raw the value container
  ConstTrackStateType(const raw_type& raw) : m_raw{&raw} {}

  /// Return if the bit at position @p pos is 1
  /// @param pos the bit position
  /// @return if the bit at @p pos is one or not
  bool test(std::size_t pos) const {
    assert(m_raw != nullptr);
    std::bitset<kRawBits> bs{*m_raw};
    return bs.test(pos);
  }

  friend std::ostream& operator<<(std::ostream& os, ConstTrackStateType t) {
    assert(t.m_raw != nullptr);
    std::bitset<kRawBits> bs{*t.m_raw};
    std::bitset<TrackStateFlag::NumTrackStateFlags> trunc;
    for (size_t i = 0; i < TrackStateFlag::NumTrackStateFlags; i++) {
      trunc[i] = bs[i];
    }
    os << "MPOHMS=" << trunc;
    return os;
  }

 private:
  friend class TrackStateType;
  const raw_type* m_raw{nullptr};
};

inline TrackStateType& TrackStateType::operator=(
    const ConstTrackStateType& other) {
  assert(other.m_raw != nullptr);
  *m_raw = *other.m_raw;
  return *this;
}
inline TrackStateType::operator ConstTrackStateType() {
  return {*m_raw};
}

// using TrackStateType = std::bitset<TrackStateFlag::NumTrackStateFlags>;

// forward declarations
template <typename derived_t>
class MultiTrajectory;
class Surface;

using ProjectorBitset = uint64_t;

namespace detail_lt {
/// Either type T or const T depending on the boolean.
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

/// Helper type to make a member pointers constness transitive.
template <typename T>
class TransitiveConstPointer {
 public:
  TransitiveConstPointer() = default;
  TransitiveConstPointer(T* ptr) : m_ptr{ptr} {}

  template <typename U>
  TransitiveConstPointer(const TransitiveConstPointer<U>& other)
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

  const T& operator*() const { return *m_ptr; }

  T& operator*() { return *m_ptr; }

 private:
  T* ptr() const { return m_ptr; }

  T* m_ptr;
};

/// Type construction helper for coefficients and associated covariances.
template <size_t Size, bool ReadOnlyMaps = true>
struct Types {
  enum {
    Flags = Eigen::ColMajor | Eigen::AutoAlign,
  };

  using Scalar = ActsScalar;
  // single items
  using Coefficients = Eigen::Matrix<Scalar, Size, 1, Flags>;
  using Covariance = Eigen::Matrix<Scalar, Size, Size, Flags>;
  using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
  using CovarianceMap = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>>;

  using DynamicCoefficients = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Flags>;
  using DynamicCovariance =
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Flags>;
  using DynamicCoefficientsMap =
      Eigen::Map<ConstIf<DynamicCoefficients, ReadOnlyMaps>>;
  using DynamicCovarianceMap =
      Eigen::Map<ConstIf<DynamicCovariance, ReadOnlyMaps>>;
};
}  // namespace detail_lt

// This is public
template <size_t M, bool ReadOnly = true>
struct TrackStateTraits {
  using Parameters =
      typename detail_lt::Types<eBoundSize, ReadOnly>::CoefficientsMap;
  using Covariance =
      typename detail_lt::Types<eBoundSize, ReadOnly>::CovarianceMap;
  using Measurement = typename detail_lt::Types<M, ReadOnly>::CoefficientsMap;
  using MeasurementCovariance =
      typename detail_lt::Types<M, ReadOnly>::CovarianceMap;

  constexpr static auto ProjectorFlags = Eigen::RowMajor | Eigen::AutoAlign;
  using Projector =
      Eigen::Matrix<typename Covariance::Scalar, M, eBoundSize, ProjectorFlags>;
};

namespace detail_lt {

/// Proxy object to access a single point on the trajectory.
///
/// @tparam SourceLink Type to link back to an original measurement
/// @tparam M         Maximum number of measurement dimensions
/// @tparam ReadOnly  true for read-only access to underlying storage
template <typename trajectory_t, size_t M, bool ReadOnly = true>
class TrackStateProxy {
 public:
  using Parameters = typename TrackStateTraits<M, ReadOnly>::Parameters;
  using Covariance = typename TrackStateTraits<M, ReadOnly>::Covariance;
  using ConstParameters = typename TrackStateTraits<M, true>::Parameters;
  using ConstCovariance = typename TrackStateTraits<M, true>::Covariance;

  template <size_t N>
  using Measurement = typename TrackStateTraits<N, ReadOnly>::Measurement;
  template <size_t N>
  using MeasurementCovariance =
      typename TrackStateTraits<N, ReadOnly>::MeasurementCovariance;
  template <size_t N>
  using ConstMeasurement = typename TrackStateTraits<N, true>::Measurement;
  template <size_t N>
  using ConstMeasurementCovariance =
      typename TrackStateTraits<N, true>::MeasurementCovariance;

  using IndexType = TrackIndexType;
  static constexpr IndexType kInvalid = kTrackIndexInvalid;

  // as opposed to the types above, this is an actual Matrix (rather than a
  // map)
  // @TODO: Does not copy flags, because this fails: can't have col major row
  // vector, but that's required for 1xN projection matrices below.
  using Projector = typename TrackStateTraits<M, ReadOnly>::Projector;
  using EffectiveProjector =
      Eigen::Matrix<typename Projector::Scalar, Eigen::Dynamic, Eigen::Dynamic,
                    TrackStateTraits<M, ReadOnly>::ProjectorFlags, M,
                    eBoundSize>;

  using Trajectory = trajectory_t;

  // Constructor and assignment operator to construct ReadOnly TrackStateProxy
  // from ReadWrite (mutable -> const)
  TrackStateProxy(const TrackStateProxy<Trajectory, M, false>& other)
      : m_traj{other.m_traj}, m_istate{other.m_istate} {}

  TrackStateProxy& operator=(
      const TrackStateProxy<Trajectory, M, false>& other) {
    m_traj = other.m_traj;
    m_istate = other.m_istate;

    return *this;
  }

  /// Index within the trajectory.
  /// @return the index
  IndexType index() const { return m_istate; }

  /// Return the index of the track state 'previous' in the track sequence
  /// @return The index of the previous track state.
  IndexType previous() const {
    return component<IndexType, hashString("previous")>();
  }

  /// Return a mutable reference to the index of the track state 'previous' in
  /// the track sequence
  /// @return The index of the previous track state.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType& previous() {
    return component<IndexType, hashString("previous")>();
  }

  /// Return whather this track state has a previous (parent) track state.
  /// @return Boolean indicating whether a previous track state exists
  bool hasPrevious() const {
    return component<IndexType, hashString("previous")>() != kInvalid;
  }

  /// Build a mask that represents all the allocated components of this track
  /// state proxy
  /// @return The generated mask
  TrackStatePropMask getMask() const;

  /// Share a shareable component within this track state
  /// @param shareSource Which component to share from
  /// @param shareTarget Which component to share as. This should be different from
  ///                    as @p shareSource, e.g. predicted can be shared as filtered.
  /// @note Shareable components are predicted, filtered, smoothed, calibrated, jacobian,
  ///       or projector. See @c TrackStatePropMask.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void shareFrom(TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) {
    shareFrom(*this, shareSource, shareTarget);
  }

  /// Share a shareable component from another track state.
  /// @param other Track state proxy to share component from
  /// @param component Which component to share.
  /// @note Shareable components are predicted, filtered, smoothed, calibrated, jacobian,
  ///       or projector. See @c TrackStatePropMask.
  /// @note The track states both need to be stored in the
  ///       same @c MultiTrajectory instance
  template <bool RO = ReadOnly, bool ReadOnlyOther,
            typename = std::enable_if_t<!RO>>
  void shareFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                 TrackStatePropMask component) {
    shareFrom(other, component, component);
  }

  /// Share a shareable component from anothe track state
  /// @param shareSource Which component to share from
  /// @param shareTarget Which component to share as. This can be be different from
  ///                    as @p shareSource, e.g. predicted can be shared as filtered.
  /// @note Shareable components are predicted, filtered, smoothed, calibrated, jacobian,
  ///       or projector. See @c TrackStatePropMask.
  template <bool RO = ReadOnly, bool ReadOnlyOther,
            typename = std::enable_if_t<!RO>>
  void shareFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                 TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) {
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
  template <typename track_state_proxy_t, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void copyFrom(const track_state_proxy_t& other,
                TrackStatePropMask mask = TrackStatePropMask::All,
                bool onlyAllocated = true) {
    using PM = TrackStatePropMask;

    // @TODO: How to support arbitrary columns here?

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
        allocateCalibrated(other.calibratedSize());

        // workaround for gcc8 bug:
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
        auto* self = this;
        visit_measurement(other.calibratedSize(), [&](auto N) {
          constexpr int measdim = decltype(N)::value;
          self->template calibrated<measdim>() =
              other.template calibrated<measdim>();
          self->template calibratedCovariance<measdim>() =
              other.template calibratedCovariance<measdim>();
        });

        setProjectorBitset(other.projectorBitset());
      }
    } else {
      if (ACTS_CHECK_BIT(mask, PM::Predicted) &&
          has<hashString("predicted")>() &&
          other.template has<hashString("predicted")>()) {
        predicted() = other.predicted();
        predictedCovariance() = other.predictedCovariance();
      }

      if (ACTS_CHECK_BIT(mask, PM::Filtered) && has<hashString("filtered")>() &&
          other.template has<hashString("filtered")>()) {
        filtered() = other.filtered();
        filteredCovariance() = other.filteredCovariance();
      }

      if (ACTS_CHECK_BIT(mask, PM::Smoothed) && has<hashString("smoothed")>() &&
          other.template has<hashString("smoothed")>()) {
        smoothed() = other.smoothed();
        smoothedCovariance() = other.smoothedCovariance();
      }

      if (other.hasUncalibratedSourceLink()) {
        setUncalibratedSourceLink(other.getUncalibratedSourceLink());
      }

      if (ACTS_CHECK_BIT(mask, PM::Jacobian) && has<hashString("jacobian")>() &&
          other.template has<hashString("jacobian")>()) {
        jacobian() = other.jacobian();
      }

      if (ACTS_CHECK_BIT(mask, PM::Calibrated) &&
          has<hashString("calibrated")>() &&
          other.template has<hashString("calibrated")>()) {
        allocateCalibrated(other.calibratedSize());

        // workaround for gcc8 bug:
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
        auto* self = this;
        visit_measurement(other.calibratedSize(), [&](auto N) {
          constexpr int measdim = decltype(N)::value;
          self->template calibrated<measdim>() =
              other.template calibrated<measdim>();
          self->template calibratedCovariance<measdim>() =
              other.template calibratedCovariance<measdim>();
        });

        setProjectorBitset(other.projectorBitset());
      }
    }

    chi2() = other.chi2();
    pathLength() = other.pathLength();
    typeFlags() = other.typeFlags();

    if (other.hasReferenceSurface()) {
      setReferenceSurface(other.referenceSurface().getSharedPtr());
    }
  }

  /// Unset an optional track state component
  /// @param target The component to unset
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void unset(TrackStatePropMask target) {
    m_traj->self().unset(target, m_istate);
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
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(std::shared_ptr<const Surface> srf) {
    m_traj->setReferenceSurface(m_istate, std::move(srf));
  }
  // NOLINTEND(performance-unnecessary-value-param)

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
    return has(hashString(key));
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component() {
    return m_traj->template component<T, key>(m_istate);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key) {
    return m_traj->template component<T>(key, m_istate);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return Mutable reference to the component given by @p key
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(std::string_view key) {
    return m_traj->template component<T>(hashString(key), m_istate);
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
    return m_traj->template component<T>(hashString(key), m_istate);
  }

  /// Track parameters vector. This tries to be somewhat smart and return the
  /// first parameters that are set in this order: predicted -> filtered ->
  /// smoothed
  /// @return one of predicted, filtered or smoothed parameters
  ConstParameters parameters() const;

  /// Track parameters covariance matrix. This tries to be somewhat smart and
  /// return the
  /// first parameters that are set in this order: predicted -> filtered ->
  /// smoothed
  /// @return one of predicted, filtered or smoothed covariances
  ConstCovariance covariance() const;

  /// Predicted track parameters vector
  /// @return The predicted parameters
  ConstParameters predicted() const {
    assert(has<hashString("predicted")>());
    return m_traj->self().parameters(
        component<IndexType, hashString("predicted")>());
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Parameters predicted() {
    assert(has<hashString("predicted")>());
    return m_traj->self().parameters(
        component<IndexType, hashString("predicted")>());
  }

  /// Predicted track parameters covariance matrix.
  /// @return The predicted track parameter covariance
  ConstCovariance predictedCovariance() const {
    assert(has<hashString("predicted")>());
    return m_traj->self().covariance(
        component<IndexType, hashString("predicted")>());
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Covariance predictedCovariance() {
    assert(has<hashString("predicted")>());
    return m_traj->self().covariance(
        component<IndexType, hashString("predicted")>());
  }

  /// Check whether the predicted parameters+covariance is set
  /// @return Whether it is set or not
  bool hasPredicted() const { return has<hashString("predicted")>(); }

  /// Filtered track parameters vector
  /// @return The filtered parameters
  /// @note Const version
  ConstParameters filtered() const {
    assert(has<hashString("filtered")>());
    return m_traj->self().parameters(
        component<IndexType, hashString("filtered")>());
  }

  /// Filtered track parameters vector
  /// @return The filtered parameters
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Parameters filtered() {
    assert(has<hashString("filtered")>());
    return m_traj->self().parameters(
        component<IndexType, hashString("filtered")>());
  }

  /// Filtered track parameters covariance matrix
  /// @return The filtered parameters covariance
  /// @note Const version
  ConstCovariance filteredCovariance() const {
    assert(has<hashString("filtered")>());
    return m_traj->self().covariance(
        component<IndexType, hashString("filtered")>());
  }

  /// Filtered track parameters covariance matrix
  /// @return The filtered parameters covariance
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Covariance filteredCovariance() {
    assert(has<hashString("filtered")>());
    return m_traj->self().covariance(
        component<IndexType, hashString("filtered")>());
  }

  /// Return whether filtered parameters+covariance is set
  /// @return Whether it is set
  bool hasFiltered() const { return has<hashString("filtered")>(); }

  /// Smoothed track parameters vector
  /// @return The smoothed parameters
  /// @note Const version
  ConstParameters smoothed() const {
    assert(has<hashString("smoothed")>());
    return m_traj->self().parameters(
        component<IndexType, hashString("smoothed")>());
  }

  /// Smoothed track parameters vector
  /// @return The smoothed parameters
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Parameters smoothed() {
    assert(has<hashString("smoothed")>());
    return m_traj->self().parameters(
        component<IndexType, hashString("smoothed")>());
  }

  /// Smoothed track parameters covariance matrix
  /// @return the parameter covariance matrix
  /// @note Const version
  ConstCovariance smoothedCovariance() const {
    assert(has<hashString("smoothed")>());
    return m_traj->self().covariance(
        component<IndexType, hashString("smoothed")>());
  }

  /// Smoothed track parameters covariance matrix
  /// @return the parameter covariance matrix
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Covariance smoothedCovariance() {
    assert(has<hashString("smoothed")>());
    return m_traj->self().covariance(
        component<IndexType, hashString("smoothed")>());
  }

  /// Return whether smoothed parameters+covariance is set
  /// @return Whether it is set
  bool hasSmoothed() const { return has<hashString("smoothed")>(); }

  /// Returns the jacobian from the previous trackstate to this one
  /// @return The jacobian matrix
  /// @note Const version
  ConstCovariance jacobian() const {
    assert(has<hashString("jacobian")>());
    return m_traj->self().jacobian(m_istate);
  }

  /// Returns the jacobian from the previous trackstate to this one
  /// @return The jacobian matrix
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Covariance jacobian() {
    assert(has<hashString("jacobian")>());
    return m_traj->self().jacobian(m_istate);
  }

  /// Returns whether a jacobian is set for this trackstate
  /// @return Whether it is set
  bool hasJacobian() const { return has<hashString("jacobian")>(); }

  /// Returns the projector (measurement mapping function) for this track
  /// state. It is derived from the uncalibrated measurement
  /// @note This function returns the overallocated projector. This means it
  /// is of dimension MxM, where M is the maximum number of measurement
  /// dimensions. The NxM submatrix, where N is the actual dimension of the
  /// measurement, is located in the top left corner, everything else is zero.
  /// @return The overallocated projector
  Projector projector() const;

  /// Returns whether a projector is set
  /// @return Whether it is set
  bool hasProjector() const { return has<hashString("projector")>(); }

  /// Returns the projector (measurement mapping function) for this track
  /// state. It is derived from the uncalibrated measurement
  /// @note This function returns the effective projector. This means it
  /// is of dimension NxM, where N is the actual dimension of the
  /// measurement.
  /// @return The effective projector
  EffectiveProjector effectiveProjector() const {
    return projector().topLeftCorner(calibratedSize(), M);
  }

  /// Set the projector on this track state
  /// This will convert the projector to a more compact bitset representation
  /// and store it.
  /// @param projector The projector in the form of a dense matrix
  /// @note @p projector is assumed to only have 0s or 1s as components.
  template <typename Derived, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void setProjector(const Eigen::MatrixBase<Derived>& projector) {
    constexpr int rows = Eigen::MatrixBase<Derived>::RowsAtCompileTime;
    constexpr int cols = Eigen::MatrixBase<Derived>::ColsAtCompileTime;

    static_assert(rows != -1 && cols != -1,
                  "Assignment of dynamic matrices is currently not supported.");

    assert(has<hashString("projector")>());

    static_assert(rows <= M, "Given projector has too many rows");
    static_assert(cols <= eBoundSize, "Given projector has too many columns");

    // set up full size projector with only zeros
    typename TrackStateProxy::Projector fullProjector =
        decltype(fullProjector)::Zero();

    // assign (potentially) smaller actual projector to matrix, preserving
    // zeroes outside of smaller matrix block.
    fullProjector.template topLeftCorner<rows, cols>() = projector;

    // convert to bitset before storing
    auto projectorBitset = matrixToBitset(fullProjector);
    component<ProjectorBitset, hashString("projector")>() =
        projectorBitset.to_ullong();
  }

  /// Get the projector bitset, a compressed form of a projection matrix
  /// @note This is mainly to copy explicitly a projector from one state
  /// to another. Use the `projector` or `effectiveProjector` method if
  /// you want to access the matrix.
  /// @return The projector bitset
  ProjectorBitset projectorBitset() const {
    assert(has<hashString("projector")>());
    return component<ProjectorBitset, hashString("projector")>();
  }

  /// Set the projector bitset, a compressed form of a projection matrix
  /// @param proj The projector bitset
  ///
  /// @note This is mainly to copy explicitly a projector from one state
  /// to another. If you have a projection matrix, set it with `setProjector`.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setProjectorBitset(ProjectorBitset proj) {
    assert(has<hashString("projector")>());
    component<ProjectorBitset, hashString("projector")>() = proj;
  }

  /// Uncalibrated measurement in the form of a source link. Const version
  /// @return The uncalibrated measurement source link
  SourceLink getUncalibratedSourceLink() const;

  /// Set an uncalibrated source link
  /// @param sourceLink The uncalibrated source link to set
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setUncalibratedSourceLink(SourceLink sourceLink) {
    m_traj->setUncalibratedSourceLink(m_istate, std::move(sourceLink));
  }

  /// Check if the point has an associated uncalibrated measurement.
  /// @return Whether it is set
  bool hasUncalibratedSourceLink() const {
    return has<hashString("uncalibratedSourceLink")>();
  }

  /// Check if the point has an associated calibrated measurement.
  /// @return Whether it is set
  bool hasCalibrated() const { return has<hashString("calibrated")>(); }

  /// Full calibrated measurement vector. Might contain additional zeroed
  /// dimensions.
  /// @return The measurement vector
  /// @note Const version
  template <size_t measdim>
  ConstMeasurement<measdim> calibrated() const {
    assert(has<hashString("calibrated")>());
    return m_traj->self().template measurement<measdim>(m_istate);
  }

  /// Full calibrated measurement vector. Might contain additional zeroed
  /// dimensions.
  /// @return The measurement vector
  /// @note Mutable version
  template <size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  Measurement<measdim> calibrated() {
    assert(has<hashString("calibrated")>());
    return m_traj->self().template measurement<measdim>(m_istate);
  }

  /// Full calibrated measurement covariance matrix. The effective covariance
  /// is located in the top left corner, everything else is zeroed.
  /// @return The measurement covariance matrix
  /// @note Const version
  template <size_t measdim>
  ConstMeasurementCovariance<measdim> calibratedCovariance() const {
    assert(has<hashString("calibratedCov")>());
    return m_traj->self().template measurementCovariance<measdim>(m_istate);
  }

  /// Full calibrated measurement covariance matrix. The effective covariance
  /// is located in the top left corner, everything else is zeroed.
  /// @return The measurement covariance matrix
  /// @note Mutable version
  template <size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  MeasurementCovariance<measdim> calibratedCovariance() {
    assert(has<hashString("calibratedCov")>());
    return m_traj->self().template measurementCovariance<measdim>(m_istate);
  }

  /// Dynamic measurement vector with only the valid dimensions.
  /// @return The effective calibrated measurement vector
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto effectiveCalibrated() {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename Types<M, ReadOnly>::DynamicCoefficientsMap{
          self->template calibrated<kMeasurementSize>().data(),
          kMeasurementSize};
    });
  }

  /// Dynamic measurement vector with only the valid dimensions.
  /// @return The effective calibrated measurement vector
  /// @note Const version
  auto effectiveCalibrated() const {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename Types<M, true>::DynamicCoefficientsMap{
          self->template calibrated<kMeasurementSize>().data(),
          kMeasurementSize};
    });
  }

  /// Dynamic measurement covariance matrix with only the valid dimensions.
  /// @return The effective calibrated covariance matrix
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto effectiveCalibratedCovariance() {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename Types<M, ReadOnly>::DynamicCovarianceMap{
          self->template calibratedCovariance<kMeasurementSize>().data(),
          kMeasurementSize, kMeasurementSize};
    });
  }

  /// Dynamic measurement covariance matrix with only the valid dimensions.
  /// @return The effective calibrated covariance matrix
  /// @note Const version
  auto effectiveCalibratedCovariance() const {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename Types<M, true>::DynamicCovarianceMap{
          self->template calibratedCovariance<kMeasurementSize>().data(),
          kMeasurementSize, kMeasurementSize};
    });
  }

  /// Return the (dynamic) number of dimensions stored for this measurement.
  /// @note The underlying storage is overallocated to MeasurementSizeMax
  /// regardless of this value
  /// @return The number of dimensions
  IndexType calibratedSize() const { return m_traj->calibratedSize(m_istate); }

  /// Overwrite existing measurement data.
  ///
  /// @tparam kMeasurementSize Size of the calibrated measurement
  /// @param meas The measurement object to set
  ///
  /// @note This assumes this TrackState stores it's own calibrated
  ///   measurement. **If storage is shared with another TrackState, both will
  ///   be overwritten!**. Also assumes none of the calibrated components is
  ///   *invalid* (i.e. unset) for this TrackState.
  /// @note This does not set the reference surface.
  template <size_t kMeasurementSize, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void setCalibrated(
      const Acts::Measurement<BoundIndices, kMeasurementSize>& meas) {
    static_assert(kMeasurementSize <= M,
                  "Input measurement must be within the allowed size");

    allocateCalibrated(kMeasurementSize);
    assert(hasCalibrated());

    calibrated<kMeasurementSize>().setZero();
    calibrated<kMeasurementSize>().template head<kMeasurementSize>() =
        meas.parameters();
    calibratedCovariance<kMeasurementSize>().setZero();
    calibratedCovariance<kMeasurementSize>()
        .template topLeftCorner<kMeasurementSize, kMeasurementSize>() =
        meas.covariance();
    setProjector(meas.projector());
  }

  void allocateCalibrated(size_t measdim) {
    m_traj->allocateCalibrated(m_istate, measdim);
  }

  /// Getter/setter for chi2 value associated with the track state
  /// This overload returns a mutable reference, which allows setting a new
  /// value directly into the backing store.
  /// @note this overload is only enabled in case the proxy is not read-only
  /// @return Mutable reference to the chi2 value
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  double& chi2() {
    return component<double, hashString("chi2")>();
  }

  /// Getter for the chi2 value associated with the track state.
  /// This overload returns a copy of the chi2 value, and thus does not allow
  /// modification of the value in the backing storage.
  /// @return the chi2 value of the track state
  double chi2() const { return component<double, hashString("chi2")>(); }

  /// Getter for the path length associated with the track state.
  /// This overloaded is only enabled if not read-only, and returns a mutable
  /// reference.
  /// @return Mutable reference to the pathlength.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  double& pathLength() {
    return component<double, hashString("pathLength")>();
  }

  /// Getter for the path length. Returns a copy of the path length value.
  /// @return The path length of this track state
  double pathLength() const {
    return component<double, hashString("pathLength")>();
  }

  /// Getter for the type flags associated with the track state.
  /// This overloaded is only enabled if not read-only, and returns a mutable
  /// reference.
  /// @return reference to the type flags.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackStateType typeFlags() {
    return TrackStateType{
        component<TrackStateType::raw_type, hashString("typeFlags")>()};
  }

  /// Getter for the type flags. Returns a copy of the type flags value.
  /// @return The type flags of this track state
  ConstTrackStateType typeFlags() const {
    return ConstTrackStateType{
        component<TrackStateType::raw_type, hashString("typeFlags")>()};
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  MultiTrajectory<Trajectory>& trajectory() {
    return *m_traj;
  }

  const MultiTrajectory<Trajectory>& trajectory() const { return *m_traj; }

 private:
  // Private since it can only be created by the trajectory.
  TrackStateProxy(ConstIf<MultiTrajectory<Trajectory>, ReadOnly>& trajectory,
                  IndexType istate);

  TransitiveConstPointer<ConstIf<MultiTrajectory<Trajectory>, ReadOnly>> m_traj;
  IndexType m_istate;

  friend class Acts::MultiTrajectory<Trajectory>;
  friend class TrackStateProxy<Trajectory, M, true>;
  friend class TrackStateProxy<Trajectory, M, false>;
};

/// Helper type that wraps two iterators
template <typename trajectory_t, size_t M, bool ReadOnly>
class TrackStateRange {
  using ProxyType = TrackStateProxy<trajectory_t, M, ReadOnly>;

 public:
  /// Iterator that wraps a track state proxy. The nullopt case signifies the
  /// end of the range, i.e. the "past-the-end" iterator
  struct Iterator {
    std::optional<ProxyType> proxy;

    using iterator_category = std::forward_iterator_tag;
    using value_type = ProxyType;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = void;

    Iterator& operator++() {
      if (!proxy) {
        return *this;
      }
      if (proxy->hasPrevious()) {
        proxy = proxy->trajectory().getTrackState(proxy->previous());
        return *this;
      } else {
        proxy = std::nullopt;
        return *this;
      }
    }

    bool operator==(const Iterator& other) const {
      if (!proxy && !other.proxy) {
        return true;
      }
      if (proxy && other.proxy) {
        return proxy->index() == other.proxy->index();
      }
      return false;
    }

    bool operator!=(const Iterator& other) const { return !(*this == other); }

    ProxyType operator*() const { return *proxy; }
    ProxyType operator*() { return *proxy; }
  };

  TrackStateRange(ProxyType _begin) : m_begin{_begin} {}
  TrackStateRange() : m_begin{std::nullopt} {}

  Iterator begin() { return m_begin; }
  Iterator end() { return Iterator{std::nullopt}; }

 private:
  Iterator m_begin;
};

// implement track state visitor concept
template <typename T, typename TS>
using call_operator_t = decltype(std::declval<T>()(std::declval<TS>()));

template <typename T, typename TS>
constexpr bool VisitorConcept = Concepts ::require<
    Concepts ::either<Concepts ::identical_to<bool, call_operator_t, T, TS>,
                      Concepts ::identical_to<void, call_operator_t, T, TS>>>;

}  // namespace detail_lt

/// This namespace contains typedefs and constant values that are used by
/// other parts of the @c MultiTrajectory implementation. It extracts these
/// from @c TrackStateTraits using the default maximum measurement dimension.
namespace MultiTrajectoryTraits {
constexpr unsigned int MeasurementSizeMax = eBoundSize;
using IndexType = TrackIndexType;
constexpr IndexType kInvalid = kTrackIndexInvalid;
}  // namespace MultiTrajectoryTraits

template <typename T>
struct IsReadOnlyMultiTrajectory;

/// Store a trajectory of track states with multiple components.
///
/// This container supports both simple, sequential trajectories as well
/// as combinatorial or multi-component trajectories. Each point can store
/// a parent point such that the trajectory forms a directed, acyclic graph
/// of sub-trajectories. From a set of endpoints, all possible sub-components
/// can be easily identified. Some functionality is provided to simplify
/// iterating over specific sub-components.
template <typename derived_t>
class MultiTrajectory {
 public:
  using Derived = derived_t;

  static constexpr bool ReadOnly = IsReadOnlyMultiTrajectory<Derived>::value;

  // Pull out type alias and re-expose them for ease of use.
  //
  static constexpr unsigned int MeasurementSizeMax =
      MultiTrajectoryTraits::MeasurementSizeMax;

  using ConstTrackStateProxy =
      detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, true>;
  using TrackStateProxy =
      detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, false>;

  using IndexType = typename TrackStateProxy::IndexType;
  static constexpr IndexType kInvalid = TrackStateProxy::kInvalid;

 protected:
  MultiTrajectory() = default;  // pseudo abstract base class

 private:
  /// Helper to static cast this to the Derived class for CRTP
  constexpr Derived& self() { return static_cast<Derived&>(*this); }
  /// Helper to static cast this to the Derived class for CRTP. Const version.
  constexpr const Derived& self() const {
    return static_cast<const Derived&>(*this);
  }

  /// Helper function to check if a component exists IF it is an optional one.
  /// Used in assertions
  constexpr bool checkOptional(HashedString key, IndexType istate) const {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "calibrated"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
        return self().has_impl(key, istate);
      default:
        return true;
    }
  }

 public:
  /// Access a read-only point on the trajectory by index.
  /// @param istate The index to access
  /// @return Read only proxy to the stored track state
  ConstTrackStateProxy getTrackState(IndexType istate) const {
    return {*this, istate};
  }

  /// Access a writable point on the trajectory by index.
  /// @param istate The index to access
  /// @return Read-write proxy to the stored track state
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackStateProxy getTrackState(IndexType istate) {
    return {*this, istate};
  }

  /// Visit all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   non-modifying functor to be called with each point
  template <typename F>
  void visitBackwards(IndexType iendpoint, F&& callable) const;

  /// Range for the track states from @p iendpoint to the trajectory start
  /// @param iendpoint Trajectory entry point to start from
  /// @return Iterator pair to iterate over
  /// @note Const version
  auto trackStateRange(IndexType iendpoint) const {
    using range_t =
        decltype(detail_lt::TrackStateRange{getTrackState(iendpoint)});
    if (iendpoint == kInvalid) {
      return range_t{};
    }

    return range_t{getTrackState(iendpoint)};
  }

  /// Range for the track states from @p iendpoint to the trajectory start
  /// @param iendpoint Trajectory entry point to start from
  /// @return Iterator pair to iterate over
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto trackStateRange(IndexType iendpoint) {
    using range_t =
        decltype(detail_lt::TrackStateRange{getTrackState(iendpoint)});
    if (iendpoint == kInvalid) {
      return range_t{};
    }

    return range_t{getTrackState(iendpoint)};
  }

  /// Apply a function to all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   modifying functor to be called with each point
  ///
  /// @warning If the trajectory contains multiple components with common
  ///          points, this can have an impact on the other components.
  template <typename F, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void applyBackwards(IndexType iendpoint, F&& callable) {
    static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                  "Callable needs to satisfy VisitorConcept");

    if (iendpoint == MultiTrajectoryTraits::kInvalid) {
      throw std::runtime_error(
          "Cannot apply backwards with kInvalid as endpoint");
    }

    while (true) {
      auto ts = getTrackState(iendpoint);
      if constexpr (std::is_same_v<std::invoke_result_t<F, TrackStateProxy>,
                                   bool>) {
        bool proceed = callable(ts);
        // this point has no parent and ends the trajectory, or a break was
        // requested
        if (!proceed || !ts.hasPrevious()) {
          break;
        }
      } else {
        callable(ts);
        // this point has no parent and ends the trajectory
        if (!ts.hasPrevious()) {
          break;
        }
      }
      iendpoint = ts.previous();
    }
  }

  auto&& convertToReadOnly() const {
    auto&& cv = self().convertToReadOnly_impl();
    static_assert(
        IsReadOnlyMultiTrajectory<decltype(cv)>::value,
        "convertToReadOnly_impl does not return something that reports "
        "being ReadOnly.");
    return cv;
  }

  /// Clear the @c MultiTrajectory. Leaves the underlying storage untouched
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void clear() {
    self().clear_impl();
  }

  /// Returns the number of track states contained
  constexpr IndexType size() const { return self().size_impl(); }

  /// Add a track state without providing explicit information. Which components
  /// of the track state are initialized/allocated can be controlled via @p mask
  /// @param mask The bitmask that instructs which components to allocate and
  /// which to leave invalid
  /// @param iprevious index of the previous state, kInvalid if first
  /// @return Index of the newly added track state
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr IndexType addTrackState(
      TrackStatePropMask mask = TrackStatePropMask::All,
      IndexType iprevious = kInvalid) {
    return self().addTrackState_impl(mask, iprevious);
  }

  /// Add a column to the @c MultiTrajectory
  /// @tparam T Type of the column values to add
  /// @note This takes a string argument rather than a hashed string to maintain
  ///       compatibility with backends.
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(const std::string& key) {
    self().template addColumn_impl<T>(key);
  }

  /// Check if a column with a key @p key exists.
  /// @param key Key to check for a column with
  /// @return True if the column exists, false if not.
  constexpr bool hasColumn(HashedString key) const {
    return self().hasColumn_impl(key);
  }

 protected:
  // These are internal helper functions which the @c TrackStateProxy class talks to

  /// Check for component existence of @p key in track satet @p istate
  /// @param key The key for which to check
  /// @param istate The track state index to check
  /// @return True if the component exists, false if not
  constexpr bool has(HashedString key, IndexType istate) const {
    return self().has_impl(key, istate);
  }

  /// Check for component existence of @p key in track satet @p istate
  /// @tparam key The key for which to check
  /// @param istate The track state index to check
  /// @return True if the component exists, false if not
  template <HashedString key>
  constexpr bool has(IndexType istate) const {
    return self().has_impl(key, istate);
  }

  /// Retrieve a parameter proxy instance for parameters at a given index
  /// @param parIdx Index into the parameter column
  /// @return Mutable proxy
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::Parameters parameters(IndexType parIdx) {
    return self().parameters_impl(parIdx);
  }

  /// Retrieve a parameter proxy instance for parameters at a given index
  /// @param parIdx Index into the parameter column
  /// @return Const proxy
  constexpr typename ConstTrackStateProxy::Parameters parameters(
      IndexType parIdx) const {
    return self().parameters_impl(parIdx);
  }

  /// Retrieve a covariance proxy instance for a covariance at a given index
  /// @param covIdx Index into the covariance column
  /// @return Mutable proxy
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::Covariance covariance(IndexType covIdx) {
    return self().covariance_impl(covIdx);
  }

  /// Retrieve a covariance proxy instance for a covariance at a given index
  /// @param covIdx Index into the covariance column
  /// @return Const proxy
  constexpr typename ConstTrackStateProxy::Covariance covariance(
      IndexType covIdx) const {
    return self().covariance_impl(covIdx);
  }

  /// Retrieve a jacobian proxy instance for a jacobian at a given index
  /// @param jacIdx Index into the jacobian column
  /// @return Mutable proxy
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::Covariance jacobian(IndexType jacIdx) {
    return self().jacobian_impl(jacIdx);
  }

  /// Retrieve a jacobian proxy instance for a jacobian at a given index
  /// @param jacIdx Index into the jacobian column
  /// @return Const proxy
  constexpr typename ConstTrackStateProxy::Covariance jacobian(
      IndexType jacIdx) const {
    return self().jacobian_impl(jacIdx);
  }

  /// Retrieve a measurement proxy instance for a measurement at a given index
  /// @tparam measdim the measurement dimension
  /// @param measIdx Index into the measurement column
  /// @return Mutable proxy
  template <size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::template Measurement<measdim> measurement(
      IndexType measIdx) {
    return self().template measurement_impl<measdim>(measIdx);
  }

  /// Retrieve a measurement proxy instance for a measurement at a given index
  /// @tparam measdim the measurement dimension
  /// @param measIdx Index into the measurement column
  /// @return Const proxy
  template <size_t measdim>
  constexpr typename ConstTrackStateProxy::template Measurement<measdim>
  measurement(IndexType measIdx) const {
    return self().template measurement_impl<measdim>(measIdx);
  }

  /// Retrieve a measurement covariance proxy instance for a measurement at a
  /// given index
  /// @tparam measdim the measurement dimension
  /// @param covIdx Index into the measurement covariance column
  /// @return Mutable proxy
  template <size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::template MeasurementCovariance<measdim>
  measurementCovariance(IndexType covIdx) {
    return self().template measurementCovariance_impl<measdim>(covIdx);
  }

  /// Retrieve a measurement covariance proxy instance for a measurement at a
  /// given index
  /// @param covIdx Index into the measurement covariance column
  /// @return Const proxy
  template <size_t measdim>
  constexpr
      typename ConstTrackStateProxy::template MeasurementCovariance<measdim>
      measurementCovariance(IndexType covIdx) const {
    return self().template measurementCovariance_impl<measdim>(covIdx);
  }

  /// Get the calibrated measurement size for a track state
  /// @param istate The track state
  /// @return the calibrated size
  IndexType calibratedSize(IndexType istate) const {
    return self().calibratedSize_impl(istate);
  }

  /// Share a shareable component from between track state.
  /// @param iself The track state index to share "into"
  /// @param iother The track state index to share from
  /// @param shareSource Which component to share from
  /// @param shareTarget Which component to share as. This doesn't have to be the same
  ///                    as @p shareSource, e.g. predicted can be shared as filtered.
  /// @note Shareable components are predicted, filtered, smoothed, calibrated, jacobian,
  ///       or projector. See @c TrackStatePropMask.
  /// @note The track states both need to be stored in the
  ///       same @c MultiTrajectory instance
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void shareFrom(IndexType iself, IndexType iother,
                           TrackStatePropMask shareSource,
                           TrackStatePropMask shareTarget) {
    self().shareFrom_impl(iself, iother, shareSource, shareTarget);
  }

  /// Unset an optional track state component
  /// @param target The component to unset
  /// @param istate The track state index to operate on
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void unset(TrackStatePropMask target, IndexType istate) {
    self().unset_impl(target, istate);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Mutable reference to the component given by @p key
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType istate) {
    assert(checkOptional(key, istate));
    return *std::any_cast<T*>(self().component_impl(key, istate));
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Mutable reference to the component given by @p key
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType istate) {
    assert(checkOptional(key, istate));
    return *std::any_cast<T*>(self().component_impl(key, istate));
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Const reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr const T& component(IndexType istate) const {
    assert(checkOptional(key, istate));
    return *std::any_cast<const T*>(self().component_impl(key, istate));
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(HashedString key, IndexType istate) const {
    assert(checkOptional(key, istate));
    return *std::any_cast<const T*>(self().component_impl(key, istate));
  }

  /// Allocate storage for a calibrated measurement of specified dimension
  /// @param istate The track state to store for
  /// @param measdim the dimension of the measurement to store
  /// @note Is a noop if the track state already has an allocation
  ///       an the dimension is the same.
  void allocateCalibrated(IndexType istate, size_t measdim) {
    self().allocateCalibrated_impl(istate, measdim);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setUncalibratedSourceLink(IndexType istate, SourceLink sourceLink) {
    self().setUncalibratedSourceLink_impl(istate, std::move(sourceLink));
  }

  SourceLink getUncalibratedSourceLink(IndexType istate) const {
    return self().getUncalibratedSourceLink_impl(istate);
  }

  const Surface* referenceSurface(IndexType istate) const {
    return self().referenceSurface_impl(istate);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(IndexType istate,
                           std::shared_ptr<const Surface> surface) {
    self().setReferenceSurface_impl(istate, std::move(surface));
  }

 private:
  friend class detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, true>;
  friend class detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, false>;
};

}  // namespace Acts

#include "Acts/EventData/MultiTrajectory.ipp"
