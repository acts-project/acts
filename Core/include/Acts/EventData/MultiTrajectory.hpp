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
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <bitset>
#include <cstdint>
#include <memory>
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

using TrackStateType = std::bitset<TrackStateFlag::NumTrackStateFlags>;

// forward declarations
template <typename derived_t>
class MultiTrajectory;
class Surface;

using ProjectorBitset = std::bitset<eBoundSize * eBoundSize>;

namespace detail_lt {
/// Either type T or const T depending on the boolean.
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

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

  using IndexType = std::uint16_t;
  static constexpr IndexType kInvalid = std::numeric_limits<IndexType>::max();
  static constexpr IndexType kNoPrevious = kInvalid - 1;

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
  using Measurement = typename TrackStateTraits<M, ReadOnly>::Parameters;
  using MeasurementCovariance =
      typename TrackStateTraits<M, ReadOnly>::Covariance;

  using IndexType = typename TrackStateTraits<M, ReadOnly>::IndexType;
  static constexpr IndexType kInvalid = TrackStateTraits<M, ReadOnly>::kInvalid;
  static constexpr IndexType kNoPrevious =
      TrackStateTraits<M, ReadOnly>::kInvalid;

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
  size_t index() const { return m_istate; }

  /// Return the index of the track state 'previous' in the track sequence
  /// @return The index of the previous track state.
  size_t previous() const {
    return component<size_t, hashString("previous")>();
  }

  bool hasPrevious() const {
    return component<size_t, hashString("previous")>() < kNoPrevious;
  }

  /// Build a mask that represents all the allocated components of this track
  /// state proxy
  /// @return The generated mask
  TrackStatePropMask getMask() const;

  template <bool RO = ReadOnly, typename = std::enable_if<!RO>>
  void shareFrom(TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) {
    shareFrom(*this, shareSource, shareTarget);
  }

  template <bool RO = ReadOnly, bool ReadOnlyOther,
            typename = std::enable_if<!RO>>
  void shareFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                 TrackStatePropMask component) {
    shareFrom(other, component, component);
  }

  template <bool RO = ReadOnly, bool ReadOnlyOther,
            typename = std::enable_if<!RO>>
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
  template <bool RO = ReadOnly, bool ReadOnlyOther,
            typename = std::enable_if<!RO>>
  void copyFrom(const TrackStateProxy<Trajectory, M, ReadOnlyOther>& other,
                TrackStatePropMask mask = TrackStatePropMask::All,
                bool onlyAllocated = true) {
    using PM = TrackStatePropMask;

    // @TODO: How to support arbitrary columns here?

    if (onlyAllocated) {
      auto dest = getMask();
      auto src = other.getMask() &
                 mask;  // combine what we have with what we want to copy
      if (static_cast<std::underlying_type_t<TrackStatePropMask>>((src ^ dest) &
                                                                  src) != 0) {
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

      // need to do it this way since other might be nullptr
      component<const SourceLink*, hashString("sourceLink")>() =
          other.template component<const SourceLink*,
                                   hashString("sourceLink")>();

      if (ACTS_CHECK_BIT(src, PM::Jacobian)) {
        jacobian() = other.jacobian();
      }

      if (ACTS_CHECK_BIT(src, PM::Calibrated)) {
        // need to do it this way since other might be nullptr
        component<const SourceLink*, hashString("calibratedSourceLink")>() =
            other.template component<const SourceLink*,
                                     hashString("calibratedSourceLink")>();
        calibrated() = other.calibrated();
        calibratedCovariance() = other.calibratedCovariance();
        calibratedSize() = other.calibratedSize();
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

      // need to do it this way since other might be nullptr
      component<const SourceLink*, hashString("sourceLink")>() =
          other.template component<const SourceLink*,
                                   hashString("sourceLink")>();

      if (ACTS_CHECK_BIT(mask, PM::Jacobian) && has<hashString("jacobian")>() &&
          other.template has<hashString("jacobian")>()) {
        jacobian() = other.jacobian();
      }

      if (ACTS_CHECK_BIT(mask, PM::Calibrated) &&
          has<hashString("calibrated")>() &&
          other.template has<hashString("calibrated")>() &&
          has<hashString("calibratedSourceLink")>() &&
          other.template has<hashString("calibratedSourceLink")>()) {
        // need to do it this way since other might be nullptr
        component<const SourceLink*, hashString("calibratedSourceLink")>() =
            other.template component<const SourceLink*,
                                     hashString("calibratedSourceLink")>();
        calibrated() = other.calibrated();
        calibratedCovariance() = other.calibratedCovariance();
        calibratedSize() = other.calibratedSize();
        setProjectorBitset(other.projectorBitset());
      }
    }

    chi2() = other.chi2();
    pathLength() = other.pathLength();
    typeFlags() = other.typeFlags();

    // can be nullptr, but we just take that
    setReferenceSurface(other.referenceSurfacePointer());
  }

  /// Unset an optional track state component
  /// @param target The component to unset
  template <bool RO = ReadOnly, typename = std::enable_if<!RO>>
  void unset(TrackStatePropMask target) {
    m_traj->self().unset(target, m_istate);
  }

  /// Reference surface.
  /// @return the reference surface
  const Surface& referenceSurface() const {
    assert(has<hashString("referenceSurface")>());
    return *component<std::shared_ptr<const Surface>,
                      hashString("referenceSurface")>();
  }

  /// Set the reference surface to a given value
  /// @param srf Shared pointer to the surface to set
  /// @note This overload is only present in case @c ReadOnly is false.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(std::shared_ptr<const Surface> srf) {
    assert(has<hashString("referenceSurface")>());
    component<std::shared_ptr<const Surface>,
              hashString("referenceSurface")>() = std::move(srf);
  }

  template <HashedString key>
  constexpr bool has() const {
    return has(key);
  }

  constexpr bool has(HashedString key) const {
    return m_traj->self().has(key, m_istate);
  }

  constexpr bool has(std::string_view key) const {
    return has(hashString(key));
  }

  template <typename T>
  constexpr bool has(T key) const {
    return has(hashString(key));
  }

  template <typename T, HashedString key>
  constexpr T& component() {
    return m_traj->self().template component<T, key>(m_istate);
  }

  template <typename T>
  constexpr T& component(HashedString key) {
    return m_traj->self().template component<T>(key, m_istate);
  }

  template <typename T>
  constexpr T& component(std::string_view key) {
    return m_traj->self().template component<T>(hashString(key), m_istate);
  }

  template <typename T, typename K>
  constexpr T& component(K key) {
    return m_traj->self().template component<T>(hashString(key), m_istate);
  }

  template <typename T, HashedString key>
  constexpr const T& component() const {
    return m_traj->self().template component<T, key>(m_istate);
  }

  template <typename T>
  constexpr const T& component(HashedString key) const {
    return m_traj->self().template component<T>(key, m_istate);
  }

  template <typename T, typename K>
  constexpr const T& component(K key) const {
    return m_traj->self().template component<T>(hashString(key), m_istate);
  }

  template <typename T>
  constexpr const T& component(std::string_view key) const {
    return m_traj->self().template component<T>(hashString(key), m_istate);
  }

  /// Track parameters vector. This tries to be somewhat smart and return the
  /// first parameters that are set in this order: predicted -> filtered ->
  /// smoothed
  /// @return one of predicted, filtered or smoothed parameters
  Parameters parameters() const;

  /// Track parameters covariance matrix. This tries to be somewhat smart and
  /// return the
  /// first parameters that are set in this order: predicted -> filtered ->
  /// smoothed
  /// @return one of predicted, filtered or smoothed covariances
  Covariance covariance() const;

  /// Predicted track parameters vector
  /// @return The predicted parameters
  Parameters predicted() const;

  /// Predicted track parameters covariance matrix.
  /// @return The predicted track parameter covariance
  Covariance predictedCovariance() const;

  /// Check whether the predicted parameters+covariance is set
  /// @return Whether it is set or not
  bool hasPredicted() const { return has<hashString("predicted")>(); }

  /// Filtered track parameters vector
  /// @return The filtered parameters
  Parameters filtered() const;

  /// Filtered track parameters covariance matrix
  /// @return The filtered parameters covariance
  Covariance filteredCovariance() const;

  /// Return whether filtered parameters+covariance is set
  /// @return Whether it is set
  bool hasFiltered() const { return has<hashString("filtered")>(); }

  /// Smoothed track parameters vector
  /// @return The smoothed parameters
  Parameters smoothed() const;

  /// Smoothed track parameters covariance matrix
  /// @return the parameter covariance matrix
  Covariance smoothedCovariance() const;

  /// Return whether smoothed parameters+covariance is set
  /// @return Whether it is set
  bool hasSmoothed() const { return has<hashString("smoothed")>(); }

  /// Returns the jacobian from the previous trackstate to this one
  /// @return The jacobian matrix
  Covariance jacobian() const;

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
    component<ProjectorBitset, hashString("projector")>() =
        matrixToBitset(fullProjector);
  }

  /// Uncalibrated measurement in the form of a source link. Const version
  /// @return The uncalibrated measurement source link
  const SourceLink& uncalibrated() const;

  /// Set an uncalibrated source link
  /// @param sourceLink The uncalibrated source link to set
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setUncalibrated(const SourceLink& sourceLink) {
    assert(has<hashString("sourceLink")>());
    using T = const SourceLink*;
    T& sl = component<const SourceLink*, hashString("sourceLink")>();
    sl = &sourceLink;

    assert(
        (component<const SourceLink*, hashString("sourceLink")>() != nullptr));
  }

  /// Check if the point has an associated calibrated measurement.
  /// @return Whether it is set
  bool hasCalibrated() const { return has<hashString("calibrated")>(); }

  /// The source link of the calibrated measurement. Const version
  /// @note This does not necessarily have to be the uncalibrated source link.
  /// @return The source link
  const SourceLink& calibratedSourceLink() const;

  /// Set a calibrated source link
  /// @param sourceLink The source link to set
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setCalibratedSourceLink(const SourceLink& sourceLink) {
    assert(has<hashString("sourceLink")>());
    m_traj->self().template component<const SourceLink*>(
        "calibratedSourceLink", m_istate) = &sourceLink;

    assert(m_traj->self().template component<const SourceLink*>(
               "calibratedSourceLink", m_istate) != nullptr);
  }

  /// Full calibrated measurement vector. Might contain additional zeroed
  /// dimensions.
  /// @return The measurement vector
  Measurement calibrated() const;

  /// Full calibrated measurement covariance matrix. The effective covariance
  /// is located in the top left corner, everything else is zeroed.
  /// @return The measurement covariance matrix
  MeasurementCovariance calibratedCovariance() const;

  /// Dynamic measurement vector with only the valid dimensions.
  /// @return The effective calibrated measurement vector
  auto effectiveCalibrated() const {
    return calibrated().head(calibratedSize());
  }

  /// Dynamic measurement covariance matrix with only the valid dimensions.
  /// @return The effective calibrated covariance matrix
  auto effectiveCalibratedCovariance() const {
    const size_t measdim = calibratedSize();
    return calibratedCovariance().topLeftCorner(measdim, measdim);
  }

  /// Return the (dynamic) number of dimensions stored for this measurement.
  /// @note The underlying storage is overallocated to MeasurementSizeMax
  /// regardless of this value
  /// @return The number of dimensions
  IndexType calibratedSize() const {
    return component<IndexType, hashString("measdim")>();
  }

  /// Return reference to the (dynamic) number of dimensions stored for this
  /// measurement.
  /// @note The underlying storage is overallocated to MeasurementSizeMax
  /// regardless of this value
  /// @return The number of dimensions
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType& calibratedSize() {
    return component<IndexType, hashString("measdim")>();
  }

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

    calibratedSize() = kMeasurementSize;

    assert(has<hashString("calibratedSourceLink")>());
    component<const SourceLink*, hashString("calibratedSourceLink")>() =
        &meas.sourceLink();
    assert(
        (component<const SourceLink*, hashString("calibratedSourceLink")>() !=
         nullptr));

    assert(hasCalibrated());
    calibrated().setZero();
    calibrated().template head<kMeasurementSize>() = meas.parameters();
    calibratedCovariance().setZero();
    calibratedCovariance()
        .template topLeftCorner<kMeasurementSize, kMeasurementSize>() =
        meas.covariance();
    setProjector(meas.projector());
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
  TrackStateType& typeFlags() {
    return component<TrackStateType, hashString("typeFlags")>();
  }

  /// Getter for the type flags. Returns a copy of the type flags value.
  /// @return The type flags of this track state
  TrackStateType typeFlags() const {
    return component<TrackStateType, hashString("typeFlags")>();
  }

 private:
  // Private since it can only be created by the trajectory.
  TrackStateProxy(ConstIf<MultiTrajectory<Trajectory>, ReadOnly>& trajectory,
                  size_t istate);

  const std::shared_ptr<const Surface>& referenceSurfacePointer() const {
    assert(has<hashString("referenceSurface")>());
    return component<std::shared_ptr<const Surface>,
                     hashString("referenceSurface")>();
  }

  ProjectorBitset projectorBitset() const {
    assert(has<hashString("projector")>());
    return component<ProjectorBitset, hashString("projector")>();
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setProjectorBitset(ProjectorBitset proj) {
    assert(has<hashString("projector")>());
    component<ProjectorBitset, hashString("projector")>() = proj;
  }

  ConstIf<MultiTrajectory<Trajectory>, ReadOnly>* m_traj;
  size_t m_istate;

  friend class Acts::MultiTrajectory<Trajectory>;
  friend class TrackStateProxy<Trajectory, M, true>;
  friend class TrackStateProxy<Trajectory, M, false>;
};

// implement track state visitor concept
template <typename T, typename TS>
using call_operator_t = decltype(std::declval<T>()(std::declval<TS>()));

template <typename T, typename TS>
constexpr bool VisitorConcept = Concepts ::require<
    Concepts ::either<Concepts ::identical_to<bool, call_operator_t, T, TS>,
                      Concepts ::identical_to<void, call_operator_t, T, TS>>>;

}  // namespace detail_lt

// template <typename derived_t>
// class MultiTrajectoryBackend {
// public:
// using Derived = derived_t;

// template <typename T>
// constexpr void addColumn(HashedString key) {
// Derived& self = static_cast<Derived&>(*this);
// assert(self.size() == 0 &&
// "Adding columns not supported after track states have been added");
// self.template addColumnImpl<T>(key);
// }

// template <HashedString key, typename T>
// constexpr void addColumn() {
// return addColumn<T>(key);
// }

// template <typename T>
// constexpr void addColumn(std::string_view key) {
// return addColumn<T>(hashString(key));
// }

// template <typename T, typename K>
// constexpr void addColumn(K key) {
// return addColumn<T>(hashString(key));
// }

// template <HashedString key>
// constexpr bool hasColumn() {
// return hasColumn(key);
// }

// constexpr bool hasColumn(HashedString key) {
// Derived& self = static_cast<Derived&>(*this);
// return self.hasColumnImpl(key);
// }

// constexpr bool hasColumn(std::string_view key) {
// Derived& self = static_cast<Derived&>(*this);
// return self.hasColumnImpl(hashString(key));
// }

// template <typename K>
// constexpr bool hasColumn(K key) {
// Derived& self = static_cast<Derived&>(*this);
// return self.hasColumnImpl(hashString(key));
// }
// };

namespace MultiTrajectoryTraits {
constexpr unsigned int MeasurementSizeMax = eBoundSize;
}

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

  static constexpr unsigned int MeasurementSizeMax =
      MultiTrajectoryTraits::MeasurementSizeMax;

  using ConstTrackStateProxy =
      detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, true>;
  using TrackStateProxy =
      detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, false>;

  using IndexType = typename TrackStateProxy::IndexType;
  static constexpr IndexType kInvalid = TrackStateProxy::kInvalid;
  static constexpr IndexType kNoPrevious = kInvalid - 1;

  template <HashedString K, typename T>
  struct Column {
    constexpr static HashedString key = K;
    using type = T;
  };

  // This is just for convenience, maybe remove
  template <HashedString K, typename T>
  using C = MultiTrajectory::Column<K, T>;

  // template <typename T, typename... Args>
  // std::unique_ptr<MultiTrajectory> static createWithBackend(
  // std::unique_ptr<T> backend, Args&&... args) {
  // MultiTrajectoryBackend<typename T::Derived>& impl = *backend;
  // auto addColumns = [&impl](auto&& column) {
  // using column_t = std::decay_t<decltype(column)>;
  // impl.template addColumn<column_t::key, typename column_t::type>();
  // };
  // (addColumns(std::forward<Args>(args)), ...);
  // return backend;
  // }

  // virtual ~MultiTrajectory() = 0;
 protected:
  MultiTrajectory() = default;  // pseudo abstract base class

 private:
  Derived& self() { return static_cast<Derived&>(*this); }
  const Derived& self() const { return static_cast<const Derived&>(*this); }

 public:
  /// Access a read-only point on the trajectory by index.
  /// @param istate The index to access
  /// @return Read only proxy to the stored track state
  ConstTrackStateProxy getTrackState(size_t istate) const {
    return {*this, istate};
  }

  /// Access a writable point on the trajectory by index.
  /// @param istate The index to access
  /// @return Read-write proxy to the stored track state
  TrackStateProxy getTrackState(size_t istate) { return {*this, istate}; }

  /// Visit all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   non-modifying functor to be called with each point
  template <typename F>
  void visitBackwards(size_t iendpoint, F&& callable) const;

  /// Apply a function to all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   modifying functor to be called with each point
  ///
  /// @warning If the trajectory contains multiple components with common
  ///          points, this can have an impact on the other components.
  template <typename F>
  void applyBackwards(size_t iendpoint, F&& callable);

  /// Clear the @c MultiTrajectory. Leaves the underlying storage untouched
  constexpr void clear() { self().clear_impl(); }

  /// Returns the number of track states contained
  constexpr size_t size() const { return self().size_impl(); }

  /// Add a track state without providing explicit information. Which components
  /// of the track state are initialized/allocated can be controlled via @p mask
  /// @param mask The bitmask that instructs which components to allocate and
  /// which to leave invalid
  /// @param iprevious index of the previous state, kInvalid if first
  /// @return Index of the newly added track state
  constexpr size_t addTrackState(
      TrackStatePropMask mask = TrackStatePropMask::All,
      size_t iprevious = kNoPrevious) {
    return self().addTrackState_impl(mask, iprevious);
  }

 protected:
  constexpr bool has(HashedString key, IndexType istate) const {
    return self().has_impl(key, istate);
  }

  constexpr typename TrackStateProxy::Parameters parameters(IndexType parIdx) {
    return self().parameters_impl(parIdx);
  };
  constexpr typename ConstTrackStateProxy::Parameters parameters(
      IndexType parIdx) const {
    return self().parameters_impl(parIdx);
  }

  constexpr typename TrackStateProxy::Covariance covariance(IndexType covIdx) {
    return self().covariance_impl(covIdx);
  }
  constexpr typename ConstTrackStateProxy::Covariance covariance(
      IndexType covIdx) const {
    return self().covariance_impl(covIdx);
  }

  constexpr typename TrackStateProxy::Covariance jacobian(IndexType covIdx) {
    return self().jacobian_impl(covIdx);
  }
  constexpr typename ConstTrackStateProxy::Covariance jacobian(
      IndexType covIdx) const {
    return self().jacobian_impl(covIdx);
  }

  constexpr typename TrackStateProxy::Measurement measurement(
      IndexType measIdx) {
    return self().measurement_impl(measIdx);
  }
  constexpr typename ConstTrackStateProxy::Measurement measurement(
      IndexType measIdx) const {
    return self().measurement_impl(measIdx);
  }

  constexpr typename TrackStateProxy::MeasurementCovariance
  measurementCovariance(IndexType covIdx) {
    return self().measurementCovariance_impl(covIdx);
  }
  constexpr typename ConstTrackStateProxy::MeasurementCovariance
  measurementCovariance(IndexType covIdx) const {
    return self().measurementCovariance_impl(covIdx);
  }

  constexpr void shareFrom(IndexType iself, IndexType iother,
                           TrackStatePropMask shareSource,
                           TrackStatePropMask shareTarget) {
    self().shareFrom_impl(iself, iother, shareSource, shareTarget);
  }

  constexpr void unset(TrackStatePropMask target, IndexType istate) {
    self().unset_impl(target, istate);
  }

  template <typename T, HashedString key>
  constexpr T& component(IndexType istate) {
    assert(self().has(key, istate));
    return *static_cast<T*>(self().component_impl(key, istate));
  }

  template <typename T>
  constexpr T& component(HashedString key, IndexType istate) {
    assert(self().has(key, istate));
    return *static_cast<T*>(self().component_impl(key, istate));
  }

  template <typename T, HashedString key>
  constexpr const T& component(IndexType istate) const {
    assert(self().has(key, istate));
    return *static_cast<const T*>(self().component_impl(key, istate));
  }

  template <typename T>
  constexpr const T& component(HashedString key, IndexType istate) const {
    assert(self().has(key, istate));
    return *static_cast<const T*>(self().component_impl(key, istate));
  }

 private:
  friend class detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, true>;
  friend class detail_lt::TrackStateProxy<Derived, MeasurementSizeMax, false>;
};  // namespace Acts

}  // namespace Acts

#include "Acts/EventData/MultiTrajectory.ipp"
