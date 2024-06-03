// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cstddef>

#include <Eigen/Core>

namespace Acts {

template <typename derived_t>
class MultiTrajectory;

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
template <std::size_t Size, bool ReadOnlyMaps = true>
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
template <std::size_t M, bool ReadOnly = true>
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

/// Proxy object to access a single point on the trajectory.
///
/// @tparam SourceLink Type to link back to an original measurement
/// @tparam M         Maximum number of measurement dimensions
/// @tparam read_only  true for read-only access to underlying storage
template <typename trajectory_t, std::size_t M, bool read_only = true>
class TrackStateProxy {
 public:
  /// Indicates whether this track state proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Alias for an associated const track state proxy, with the same backends
  using ConstProxyType = TrackStateProxy<trajectory_t, M, true>;

  /// Map-type for a bound parameter vector. This has reference semantics, i.e.
  /// points at a matrix by an internal pointer.
  using Parameters = typename TrackStateTraits<M, ReadOnly>::Parameters;

  /// Same as @ref Parameters, but with const semantics
  using ConstParameters = typename TrackStateTraits<M, true>::Parameters;

  /// Map-type for a bound covariance. This has reference semantics, i.e.
  /// points at a matrix by an internal pointer.
  using Covariance = typename TrackStateTraits<M, ReadOnly>::Covariance;

  /// Same as @ref Covariance, but with const semantics
  using ConstCovariance = typename TrackStateTraits<M, true>::Covariance;

  /// Map-type for a measurement vector, where the local measurement dimension
  /// is variable.
  template <std::size_t N>
  using Measurement = typename TrackStateTraits<N, ReadOnly>::Measurement;

  /// Same as @c Measurement, but with const semantics
  template <std::size_t N>
  using ConstMeasurement = typename TrackStateTraits<N, true>::Measurement;

  /// Map-type for a measurement covariance matrix, where the local measurement
  /// dimension is variable.
  template <std::size_t N>
  using MeasurementCovariance =
      typename TrackStateTraits<N, ReadOnly>::MeasurementCovariance;

  /// Same as @ref MeasurementCovariance, but with const semantics
  template <std::size_t N>
  using ConstMeasurementCovariance =
      typename TrackStateTraits<N, true>::MeasurementCovariance;

  /// The index type of the track state container
  using IndexType = TrackIndexType;

  /// Sentinel value that indicates an invalid index
  static constexpr IndexType kInvalid = kTrackIndexInvalid;

  /// Matrix representing the projector (measurement mapping function) for a
  /// measurement.  This is not a map type, but an actual matrix. This matrix
  /// is always \f$M \times M\f$, even if the local measurement dimension is lower.
  /// The actual \f$N\times M\f$ projector is given by the top \f$N\f$ rows.
  using Projector = typename TrackStateTraits<M, ReadOnly>::Projector;

  /// Dynamic variant of the projector matrix
  /// @warning Using this type is discouraged, as it has a runtime overhead
  using EffectiveProjector =
      Eigen::Matrix<typename Projector::Scalar, Eigen::Dynamic, Eigen::Dynamic,
                    TrackStateTraits<M, ReadOnly>::ProjectorFlags, M,
                    eBoundSize>;

  /// The track state container backend given as a template parameter
  using Trajectory = trajectory_t;

  /// @anchor track_state_proxy_construct
  /// @name Constructors and assignment operator
  ///
  /// Public constructors and assignment operators for @c TrackStateProxy only
  /// allow construction from another @c TrackStateProxy. You should generally
  /// not have to construct @c TrackStateProxy manually.
  ///
  /// @{

  /// Constructor and assignment operator to construct TrackStateProxy
  /// from mutable
  /// @param other The other TrackStateProxy to construct from
  TrackStateProxy(const TrackStateProxy<Trajectory, M, false>& other)
      : m_traj{other.m_traj}, m_istate{other.m_istate} {}

  /// Assignment operator to from mutable @c TrackStateProxy
  /// @param other The other TrackStateProxy to construct from
  /// @return Reference to this TrackStateProxy
  TrackStateProxy& operator=(
      const TrackStateProxy<Trajectory, M, false>& other) {
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

  /// Return the index of the track state `previous` in the track sequence
  /// @return The index of the previous track state.
  IndexType previous() const {
    return component<IndexType, hashString("previous")>();
  }

  /// Return a mutable reference to the index of the track state 'previous' in
  /// the track sequence
  /// @note Only available if the track state proxy is not read-only
  /// @return The index of the previous track state.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType& previous() {
    return component<IndexType, hashString("previous")>();
  }

  /// Return whether this track state has a previous (parent) track state.
  /// @return Boolean indicating whether a previous track state exists
  bool hasPrevious() const {
    return component<IndexType, hashString("previous")>() != kInvalid;
  }

  /// Build a mask that represents all the allocated components of this track
  /// state proxy
  /// @return The generated mask
  TrackStatePropMask getMask() const;

  /// Unset an optional track state component
  /// @note Only available if the track state proxy is not read-only
  /// @param target The component to unset
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void unset(TrackStatePropMask target) {
    m_traj->self().unset(target, m_istate);
  }

  /// Add additional components to the track state
  /// @note Only available if the track state proxy is not read-only
  /// @param mask The bitmask that instructs which components to allocate
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void addComponents(TrackStatePropMask mask) {
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
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(std::shared_ptr<const Surface> srf) {
    m_traj->setReferenceSurface(m_istate, std::move(srf));
  }
  // NOLINTEND(performance-unnecessary-value-param)

  /// Getter/setter for chi2 value associated with the track state
  /// This overload returns a mutable reference, which allows setting a new
  /// value directly into the backing store.
  /// @note this overload is only enabled in case the proxy is not read-only
  /// @return Mutable reference to the chi2 value
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  float& chi2() {
    return component<float, hashString("chi2")>();
  }

  /// Getter for the chi2 value associated with the track state.
  /// This overload returns a copy of the chi2 value, and thus does not allow
  /// modification of the value in the backing storage.
  /// @return the chi2 value of the track state
  float chi2() const { return component<float, hashString("chi2")>(); }

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

  /// @}

  /// @anchor track_state_proxy_params
  /// @name Track state parameters
  /// @{

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

  /// @}

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
  /// A convenience function to assign this from the @ref Measurement class
  /// is provided, although it's use is discouraged.
  ///
  /// The track state also includes a @ref SourceLink which acts as a proxy
  /// to the original uncalibrated measurement that the calibrated measurement
  /// was derived from. It is set and returned by value, to allow unpacking /
  /// repacking by the backend, if needed.
  ///
  /// @{

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
  /// @warning This function returns the effective projector. This means it
  /// is of dimension \f$N\times M\f$, where \f$N\f$ is the actual dimension of the
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

  /// Directly get the projector bitset, a compressed form of a projection
  /// matrix
  /// @note This is mainly to copy explicitly a projector from one state
  ///       to another. Use the `projector` or `effectiveProjector` method if
  ///       you want to access the matrix.
  /// @return The projector bitset
  ProjectorBitset projectorBitset() const {
    assert(has<hashString("projector")>());
    return component<ProjectorBitset, hashString("projector")>();
  }

  /// Set the projector bitset, a compressed form of a projection matrix
  /// @param proj The projector bitset
  ///
  /// @note This is mainly to copy explicitly a projector from one state
  ///       to another. If you have a projection matrix, set it with
  ///       `setProjector`.
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
  template <std::size_t measdim>
  ConstMeasurement<measdim> calibrated() const {
    assert(has<hashString("calibrated")>());
    return m_traj->self().template measurement<measdim>(m_istate);
  }

  /// Full calibrated measurement vector. Might contain additional zeroed
  /// dimensions.
  /// @return The measurement vector
  /// @note Mutable version
  template <std::size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  Measurement<measdim> calibrated() {
    assert(has<hashString("calibrated")>());
    return m_traj->self().template measurement<measdim>(m_istate);
  }

  /// Const full calibrated measurement covariance matrix. The effective
  /// covariance is located in the top left corner, everything else is zeroed.
  /// @return The measurement covariance matrix
  template <std::size_t measdim>
  ConstMeasurementCovariance<measdim> calibratedCovariance() const {
    assert(has<hashString("calibratedCov")>());
    return m_traj->self().template measurementCovariance<measdim>(m_istate);
  }

  /// Mutable full calibrated measurement covariance matrix. The effective
  /// covariance is located in the top left corner, everything else is zeroed.
  /// @return The measurement covariance matrix
  template <std::size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  MeasurementCovariance<measdim> calibratedCovariance() {
    assert(has<hashString("calibratedCov")>());
    return m_traj->self().template measurementCovariance<measdim>(m_istate);
  }

  /// Mutable dynamic measurement vector with only the valid dimensions.
  /// @warning The dynamic vector has a runtime overhead!
  /// @return The effective calibrated measurement vector
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto effectiveCalibrated() {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename detail_lt::Types<M, ReadOnly>::DynamicCoefficientsMap{
          self->template calibrated<kMeasurementSize>().data(),
          kMeasurementSize};
    });
  }

  /// Const dynamic measurement vector with only the valid dimensions.
  /// @warning The dynamic matrix has a runtime overhead!
  /// @return The effective calibrated measurement vector
  auto effectiveCalibrated() const {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename detail_lt::Types<M, true>::DynamicCoefficientsMap{
          self->template calibrated<kMeasurementSize>().data(),
          kMeasurementSize};
    });
  }

  /// Mutable dynamic measurement covariance matrix with only the valid
  /// dimensions.
  /// @warning The dynamic matrix has a runtime overhead!
  /// @return The effective calibrated covariance matrix
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto effectiveCalibratedCovariance() {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename detail_lt::Types<M, ReadOnly>::DynamicCovarianceMap{
          self->template calibratedCovariance<kMeasurementSize>().data(),
          kMeasurementSize, kMeasurementSize};
    });
  }

  /// Const dynamic measurement covariance matrix with only the valid
  /// dimensions.
  /// @warning The dynamic matrix has a runtime overhead!
  /// @return The effective calibrated covariance matrix
  auto effectiveCalibratedCovariance() const {
    // repackage the data pointer to a dynamic map type
    // workaround for gcc8 bug:
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86594
    auto* self = this;
    return visit_measurement(calibratedSize(), [&](auto N) {
      constexpr int kMeasurementSize = decltype(N)::value;
      return typename detail_lt::Types<M, true>::DynamicCovarianceMap{
          self->template calibratedCovariance<kMeasurementSize>().data(),
          kMeasurementSize, kMeasurementSize};
    });
  }

  /// Return the (dynamic) number of dimensions stored for this measurement.
  /// @note Depending on the backend, this size is used to determine the
  ///       memory range of the measurement vector and covariance.
  /// @return The number of dimensions
  IndexType calibratedSize() const { return m_traj->calibratedSize(m_istate); }

  /// Set the calibrated measurement of this track state from a measurement
  /// object. This is a convenience function to set the calibrated measurement
  /// from the @c Acts::Measurement object. In general, you should implement this
  /// functionality specifically for an (experiment-specific) uncalibrated
  /// measurement EDM.
  ///
  /// @tparam kMeasurementSize Size of the calibrated measurement
  /// @param meas The measurement object to set
  ///
  /// @warning This assumes this TrackState stores it's own calibrated
  ///       measurement. **If storage is shared with another TrackState, both
  ///       will be overwritten!**. Also assumes none of the calibrated
  ///       components is *invalid* (i.e. unset) for this TrackState.
  /// @note This does not set the reference surface.
  template <std::size_t kMeasurementSize, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void setCalibrated(
      const Acts::Measurement<BoundIndices, kMeasurementSize>& meas) {
    static_assert(kMeasurementSize <= M,
                  "Input measurement must be within the allowed size");

    allocateCalibrated(kMeasurementSize);
    assert(hasCalibrated());

    calibrated<kMeasurementSize>() = meas.parameters();
    calibratedCovariance<kMeasurementSize>() = meas.covariance();
    setProjector(meas.projector());
  }

  /// Allocate storage to be able to store a measurement of size @p measdim.
  /// This must be called **before** setting the measurement content.
  void allocateCalibrated(std::size_t measdim) {
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
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void shareFrom(TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) {
    shareFrom(*this, shareSource, shareTarget);
  }

  /// Share a shareable component from another track state.
  /// @param other Track state proxy to share component from
  /// @param component Which component to share.
  /// @note The track states both need to be stored in the
  ///       same @c MultiTrajectory instance
  template <bool RO = ReadOnly, bool ReadOnlyOther,
            typename = std::enable_if_t<!RO>>
  void shareFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                 TrackStatePropMask component) {
    shareFrom(other, component, component);
  }

  /// Share a shareable component from another track state
  /// @param other Track state proxy to share component(s) from
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
  template <ACTS_CONCEPT(TrackStateProxyConcept) track_state_proxy_t,
            bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void copyFrom(const track_state_proxy_t& other,
                TrackStatePropMask mask = TrackStatePropMask::All,
                bool onlyAllocated = true) {
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

  /// @}

  /// Return a mutable reference to the underlying backend container
  /// @return A reference to the backend container
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  MultiTrajectory<Trajectory>& trajectory() {
    return *m_traj;
  }

  /// Return a const reference to the underlying backend container
  /// @return A const reference to the backend container
  const MultiTrajectory<Trajectory>& trajectory() const { return *m_traj; }

  /// Get a mutable reference to the track state container backend
  /// @return a mutable reference to the backend
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& container() {
    return *m_traj;
  }

  /// Get a const reference to the track state container backend
  /// @return a const reference to the backend
  const auto& container() const { return *m_traj; }

 private:
  // Private since it can only be created by the trajectory.
  TrackStateProxy(
      detail_lt::ConstIf<MultiTrajectory<Trajectory>, ReadOnly>& trajectory,
      IndexType istate);

  detail_lt::TransitiveConstPointer<
      detail_lt::ConstIf<MultiTrajectory<Trajectory>, ReadOnly>>
      m_traj;
  IndexType m_istate;

  friend class Acts::MultiTrajectory<Trajectory>;
  friend class TrackStateProxy<Trajectory, M, true>;
  friend class TrackStateProxy<Trajectory, M, false>;
};
}  // namespace Acts

#include "Acts/EventData/TrackStateProxy.ipp"
