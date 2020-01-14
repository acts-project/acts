// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "Acts/EventData/TrackState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

// forward declarations
template <typename source_link_t>
class MultiTrajectory;
class Surface;

namespace detail_lt {
/// Either type T or const T depending on the boolean.
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;
/// wrapper for a dynamic Eigen type that adds support for automatic growth
///
/// \warning Assumes the underlying storage has a fixed number of rows
template <typename Storage, size_t kSizeIncrement>
struct GrowableColumns {
  /// Make sure storage for @p n additional columns is allocated. Will update
  /// the size of the container accordingly. The indices added by this call
  /// can safely be written to.
  /// @param n Number of columns to add, defaults to 1.
  /// @return View into the last allocated column
  auto addCol(size_t n = 1) {
    size_t index = m_size + (n - 1);
    while (capacity() <= index) {
      data.conservativeResize(Eigen::NoChange, data.cols() + kSizeIncrement);
    }
    m_size = index + 1;

    // @TODO: do this or not? If we assume this happens only when something is
    // written, the expectation is that everything is zero
    data.col(index).setZero();

    return data.col(index);
  }

  /// Writable access to a column w/o checking its existence first.
  auto col(size_t index) { return data.col(index); }

  /// Read-only access to a column w/o checking its existence first.
  auto col(size_t index) const { return data.col(index); }

  /// Return the current allocated storage capacity
  size_t capacity() const { return static_cast<size_t>(data.cols()); }

  size_t size() const { return m_size; }

 private:
  Storage data;
  size_t m_size{0};
};

/// Type construction helper for coefficients and associated covariances.
template <size_t Size, bool ReadOnlyMaps = true>
struct Types {
  enum {
    Flags = Eigen::ColMajor | Eigen::AutoAlign,
    SizeIncrement = 8,
  };
  using Scalar = double;
  // single items
  using Coefficients = Eigen::Matrix<Scalar, Size, 1, Flags>;
  using Covariance = Eigen::Matrix<Scalar, Size, Size, Flags>;
  using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
  using CovarianceMap = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>>;
  // storage of multiple items in flat arrays
  using StorageCoefficients =
      GrowableColumns<Eigen::Array<Scalar, Size, Eigen::Dynamic, Flags>,
                      SizeIncrement>;
  using StorageCovariance =
      GrowableColumns<Eigen::Array<Scalar, Size * Size, Eigen::Dynamic, Flags>,
                      SizeIncrement>;
};

struct IndexData {
  using IndexType = uint16_t;

  static constexpr IndexType kInvalid = UINT16_MAX;

  IndexType irefsurface = kInvalid;
  IndexType iprevious = kInvalid;
  IndexType ipredicted = kInvalid;
  IndexType ifiltered = kInvalid;
  IndexType ismoothed = kInvalid;
  IndexType ijacobian = kInvalid;
  IndexType iprojector = kInvalid;

  double chi2;
  double pathLength;
  TrackStateType typeFlags;

  IndexType iuncalibrated = kInvalid;
  IndexType icalibrated = kInvalid;
  IndexType icalibratedsourcelink = kInvalid;
  IndexType measdim = 0;
};

/// Proxy object to access a single point on the trajectory.
///
/// @tparam source_link_t Type to link back to an original measurement
/// @tparam N         Number of track parameters
/// @tparam M         Maximum number of measurement dimensions
/// @tparam ReadOnly  true for read-only access to underlying storage
template <typename source_link_t, size_t N, size_t M, bool ReadOnly = true>
class TrackStateProxy {
 public:
  using SourceLink = source_link_t;
  using Parameters = typename Types<N, ReadOnly>::CoefficientsMap;
  using Covariance = typename Types<N, ReadOnly>::CovarianceMap;
  using Measurement = typename Types<M, ReadOnly>::CoefficientsMap;
  using MeasurementCovariance = typename Types<M, ReadOnly>::CovarianceMap;

  // as opposed to the types above, this is an actual Matrix (rather than a
  // map)
  // @TODO: Does not copy flags, because this fails: can't have col major row
  // vector, but that's required for 1xN projection matrices below.
  constexpr static auto ProjectorFlags = Eigen::RowMajor | Eigen::AutoAlign;
  using Projector =
      Eigen::Matrix<typename Covariance::Scalar, M, N, ProjectorFlags>;
  using EffectiveProjector =
      Eigen::Matrix<typename Projector::Scalar, Eigen::Dynamic, Eigen::Dynamic,
                    ProjectorFlags, M, N>;

  /// Index within the trajectory.
  /// @return the index
  size_t index() const { return m_istate; }

  /// Return the index of the track state 'previous' in the track sequence
  /// @return The index of the previous track state.
  size_t previous() const { return data().iprevious; }

  /// Return the index tuple that makes up this track state
  /// @return Mutable ref to index tuple from the parent @c MultiTrajectory
  /// @note This overload is only present in case @c ReadOnly is false.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexData& data() {
    return m_traj->m_index[m_istate];
  }

  /// Return the index tuple that makes up this track state
  /// @return Immutable ref to index tuple from the parent @c MultiTrajectory
  const IndexData& data() const { return m_traj->m_index[m_istate]; }

  /// Reference surface.
  /// @return the reference surface
  const Surface& referenceSurface() const {
    assert(data().irefsurface != IndexData::kInvalid);
    return *m_traj->m_referenceSurfaces[data().irefsurface];
  }

  /// Set the reference surface to a given value
  /// @param srf Shared pointer to the surface to set
  /// @note This overload is only present in case @c ReadOnly is false.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(std::shared_ptr<const Surface> srf) {
    m_traj->m_referenceSurfaces[data().irefsurface] = std::move(srf);
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

  /// Predicted parameters in the form of BoundParameters
  /// @return BoundParameters instance of the predicted parameters
  BoundParameters predictedParameters(const Acts::GeometryContext& gctx) const;

  /// Predicted track parameters covariance matrix.
  /// @return The predicted track parameter covariance
  Covariance predictedCovariance() const;

  /// Check whether the predicted parameters+covariance is set
  /// @return Whether it is set or not
  bool hasPredicted() const { return data().ipredicted != IndexData::kInvalid; }

  /// Filtered track parameters vector
  /// @return The filtered parameters
  Parameters filtered() const;

  /// Filtered track parameters covariance matrix
  /// @return The filtered parameters covariance
  Covariance filteredCovariance() const;

  /// Filtered parameters in the form of BoundParameters
  /// @return BoundParameters instance of the filtered parameters
  BoundParameters filteredParameters(const Acts::GeometryContext& gctx) const;

  /// Return whether filtered parameters+covariance is set
  /// @return Whether it is set
  bool hasFiltered() const { return data().ifiltered != IndexData::kInvalid; }

  /// Smoothed track parameters vector
  /// @return the parameter vector
  Parameters smoothed() const;

  /// Smoothed track parameters covariance matrix
  /// @return the parameter covariance matrix
  Covariance smoothedCovariance() const;

  /// Smoothed parameters in the form of BoundParameters
  /// @return BoundParameters instance of the smoothed parameters
  BoundParameters smoothedParameters(const Acts::GeometryContext& gctx) const;

  /// Return whether smoothed parameters+covariance is set
  /// @return Whether it is set
  bool hasSmoothed() const { return data().ismoothed != IndexData::kInvalid; }

  /// Returns the jacobian from the previous trackstate to this one
  /// @return The jacobian matrix
  Covariance jacobian() const;

  /// Returns whether a jacobian is set for this trackstate
  /// @return Whether it is set
  bool hasJacobian() const { return data().ijacobian != IndexData::kInvalid; }

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
  bool hasProjector() const { return data().iprojector != IndexData::kInvalid; }

  /// Returns the projector (measurement mapping function) for this track
  /// state. It is derived from the uncalibrated measurement
  /// @note This function returns the effective projector. This means it
  /// is of dimension NxM, where N is the actual dimension of the
  /// measurement.
  /// @return The effective projector
  EffectiveProjector effectiveProjector() const {
    return projector().topLeftCorner(data().measdim, M);
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

    IndexData& dataref = data();
    assert(dataref.iprojector != IndexData::kInvalid);

    constexpr int max_measdim = MultiTrajectory<SourceLink>::MeasurementSizeMax;

    static_assert(rows <= max_measdim, "Given projector has too many rows");
    static_assert(cols <= max_measdim, "Given projector has too many columns");

    // set up full size projector with only zeros
    typename TrackStateProxy::Projector fullProjector =
        decltype(fullProjector)::Zero();

    // assign (potentially) smaller actual projector to matrix, preserving
    // zeroes outside of smaller matrix block.
    fullProjector.template topLeftCorner<rows, max_measdim>() = projector;

    // convert to bitset before storing
    m_traj->m_projectors[dataref.iprojector] = matrixToBitset(fullProjector);
  }

  /// Return whether an uncalibrated measurement (source link) is set
  /// @return Whether it is set
  bool hasUncalibrated() const {
    return data().iuncalibrated != IndexData::kInvalid;
  }

  /// Uncalibrated measurement in the form of a source link. Const version
  /// @return The uncalibrated measurement source link
  const SourceLink& uncalibrated() const;

  /// Uncalibrated measurement in the form of a source link. Mutable version
  /// @note This overload is only available if @c ReadOnly is false
  /// @return The uncalibrated measurement source link
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  SourceLink& uncalibrated() {
    assert(data().iuncalibrated != IndexData::kInvalid);
    return m_traj->m_sourceLinks[data().iuncalibrated];
  }

  /// Check if the point has an associated calibrated measurement.
  /// @return Whether it is set
  bool hasCalibrated() const {
    return data().icalibrated != IndexData::kInvalid;
  }

  /// The source link of the calibrated measurement. Const version
  /// @note This does not necessarily have to be the uncalibrated source link.
  /// @return The source link
  const SourceLink& calibratedSourceLink() const;

  /// The source link of the calibrated measurement. Mutable version
  /// @note This does not necessarily have to be the uncalibrated source link.
  /// @note This overload is only available if @c ReadOnly is false
  /// @return The source link
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  SourceLink& calibratedSourceLink() {
    assert(data().icalibratedsourcelink != IndexData::kInvalid);
    return m_traj->m_sourceLinks[data().icalibratedsourcelink];
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
  auto effectiveCalibrated() const { return calibrated().head(data().measdim); }

  /// Dynamic measurement covariance matrix with only the valid dimensions.
  /// @return The effective calibrated covariance matrix
  auto effectiveCalibratedCovariance() const {
    return calibratedCovariance().topLeftCorner(data().measdim, data().measdim);
  }

  /// Return the (dynamic) number of dimensions stored for this measurement.
  /// @note The underlying storage is overallocated to MeasurementSizeMax
  /// regardless of this value
  /// @return The number of dimensions
  size_t calibratedSize() const { return data().measdim; }

  /// Setter for a full measurement object
  /// @note This assumes this TrackState stores it's own calibrated
  /// measurement. **If storage is shared with another TrackState, both will
  /// be overwritten!**. Also assumes none of the calibrated components is
  /// *invalid* (i.e. unset) for this TrackState..
  /// @tparam params The parameter tags of the measurement
  /// @param meas The measurement object to set
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>,
            ParID_t... params>
  void setCalibrated(const Acts::Measurement<SourceLink, params...>& meas) {
    IndexData& dataref = data();
    constexpr size_t measdim = Acts::Measurement<SourceLink, params...>::size();

    dataref.measdim = measdim;

    assert(hasCalibrated());
    calibrated().setZero();
    calibrated().template head<measdim>() = meas.parameters();

    calibratedCovariance().setZero();
    calibratedCovariance().template topLeftCorner<measdim, measdim>() =
        meas.covariance();

    setProjector(meas.projector());

    // this shouldn't change
    assert(data().irefsurface != IndexData::kInvalid);
    std::shared_ptr<const Surface>& refSrf =
        m_traj->m_referenceSurfaces[dataref.irefsurface];
    // either unset, or the same, otherwise this is inconsistent assignment
    assert(!refSrf || refSrf.get() == &meas.referenceSurface());
    if (!refSrf) {
      // ref surface is not set, set it now
      refSrf = meas.referenceSurface().getSharedPtr();
    }

    assert(dataref.icalibratedsourcelink != IndexData::kInvalid);
    calibratedSourceLink() = meas.sourceLink();
  }

  /// Setter for a full measurement object.
  /// @note This allocates new storage for the calibrated measurement. If this
  /// TrackState previously already had unique storage for these components,
  /// they will **not be removed**, but may become unaccessible.
  /// @tparam params The parameter tags of the measurement
  /// @param meas The measurement object to set
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>,
            ParID_t... params>
  void resetCalibrated(const Acts::Measurement<SourceLink, params...>& meas) {
    IndexData& dataref = data();
    auto& traj = *m_traj;
    // force reallocate, whether currently invalid or shared index
    traj.m_meas.addCol();
    traj.m_measCov.addCol();
    // shared index between meas par
    // and cov
    dataref.icalibrated = traj.m_meas.size() - 1;

    traj.m_sourceLinks.emplace_back();
    dataref.icalibratedsourcelink = traj.m_sourceLinks.size() - 1;

    traj.m_projectors.emplace_back();
    dataref.iprojector = traj.m_projectors.size() - 1;

    // now actually assign to the allocated entries
    setCalibrated(meas);
  }

  /// Getter/setter for chi2 value associated with the track state
  /// This overload returns a mutable reference, which allows setting a new
  /// value directly into the backing store.
  /// @note this overload is only enabled in case the proxy is not read-only
  /// @return Mutable reference to the chi2 value
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  double& chi2() {
    return data().chi2;
  }

  /// Getter for the chi2 value associated with the track state.
  /// This overload returns a copy of the chi2 value, and thus does not allow
  /// modification of the value in the backing storage.
  /// @return the chi2 value of the track state
  double chi2() const { return data().chi2; }

  /// Getter for the path length associated with the track state.
  /// This overloaded is only enabled if not read-only, and returns a mutable
  /// reference.
  /// @return Mutable reference to the pathlength.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  double& pathLength() {
    return data().pathLength;
  }

  /// Getter for the path length. Returns a copy of the path length value.
  /// @return The path length of this track state
  double pathLength() const { return data().pathLength; }

  /// Getter for the type flags associated with the track state.
  /// This overloaded is only enabled if not read-only, and returns a mutable
  /// reference.
  /// @return reference to the type flags.
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackStateType& typeFlags() {
    return data().typeFlags;
  }

  /// Getter for the type flags. Returns a copy of the type flags value.
  /// @return The type flags of this track state
  TrackStateType typeFlags() const { return data().typeFlags; }

 private:
  // Private since it can only be created by the trajectory.
  TrackStateProxy(ConstIf<MultiTrajectory<SourceLink>, ReadOnly>& trajectory,
                  size_t istate);

  ConstIf<MultiTrajectory<SourceLink>, ReadOnly>* m_traj;
  size_t m_istate;

  friend class Acts::MultiTrajectory<SourceLink>;
};

// implement track state visitor concept
template <typename T, typename TS>
using call_operator_t = decltype(std::declval<T>()(std::declval<TS>()));

template <typename T, typename TS>
constexpr bool VisitorConcept = concept ::require<
    concept ::either<concept ::identical_to<bool, call_operator_t, T, TS>,
                     concept ::identical_to<void, call_operator_t, T, TS>>>;

}  // namespace detail_lt

/// Collection of bit masks to enable steering which components of a track state
/// should be initialized, and which should be left invalid.
/// These mask values can be combined using binary operators, so
/// (TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian) will instruct
/// allocating storage for both predicted parameters (including covariance) and
/// a jacobian.
namespace TrackStatePropMask {
/// Type of the bitmasks
using Type = std::bitset<8>;

constexpr Type Predicted{1 << 0};
constexpr Type Filtered{1 << 1};
constexpr Type Smoothed{1 << 2};
constexpr Type Jacobian{1 << 3};

constexpr Type Uncalibrated{1 << 4};
constexpr Type Calibrated{1 << 5};

// Initialize to all 1s. This only works as long as number of bits <= 64
constexpr Type All{static_cast<unsigned long long>(-1)};
constexpr Type None{0};
}  // namespace TrackStatePropMask

/// Store a trajectory of track states with multiple components.
///
/// This container supports both simple, sequential trajectories as well
/// as combinatorial or multi-component trajectories. Each point can store
/// a parent point such that the trajectory forms a directed, acyclic graph
/// of sub-trajectories. From a set of endpoints, all possible sub-components
/// can be easily identified. Some functionality is provided to simplify
/// iterating over specific sub-components.
/// @tparam source_link_t Type to link back to an original measurement
template <typename source_link_t>
class MultiTrajectory {
 public:
  enum {
    ParametersSize = BoundParsDim,
    MeasurementSizeMax = BoundParsDim,
  };
  using SourceLink = source_link_t;
  using ConstTrackStateProxy =
      detail_lt::TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax,
                                 true>;
  using TrackStateProxy = detail_lt::TrackStateProxy<SourceLink, ParametersSize,
                                                     MeasurementSizeMax, false>;

  using ProjectorBitset = std::bitset<ParametersSize * MeasurementSizeMax>;

  /// Create an empty trajectory.
  MultiTrajectory() = default;

  /// Add a track state using information from a separate track state object.
  ///
  /// @tparam parameters_t The parameter type used for the trackstate
  /// @param trackParameters  at the local point
  /// @param iprevious        index of the previous state, SIZE_MAX if first
  /// @note The parameter type from @p parameters_t is not currently stored in
  /// MultiTrajectory.
  /// @return Index of the newly added track state
  template <typename parameters_t>
  size_t addTrackState(const TrackState<SourceLink, parameters_t>& ts,
                       size_t iprevious = SIZE_MAX);

  /// Add a track state without providing explicit information. Which components
  /// of the track state are initialized/allocated can be controlled via @p mask
  /// @param mask The bitmask that instructs which components to allocate and
  /// which to leave invalid
  /// @param iprevious index of the previous state, SIZE_MAX if first
  /// @return Index of the newly added track state
  size_t addTrackState(
      const TrackStatePropMask::Type& mask = TrackStatePropMask::All,
      size_t iprevious = SIZE_MAX);

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

 private:
  /// index to map track states to the corresponding
  std::vector<detail_lt::IndexData> m_index;
  typename detail_lt::Types<ParametersSize>::StorageCoefficients m_params;
  typename detail_lt::Types<ParametersSize>::StorageCovariance m_cov;
  typename detail_lt::Types<MeasurementSizeMax>::StorageCoefficients m_meas;
  typename detail_lt::Types<MeasurementSizeMax>::StorageCovariance m_measCov;
  typename detail_lt::Types<ParametersSize>::StorageCovariance m_jac;
  std::vector<SourceLink> m_sourceLinks;
  std::vector<ProjectorBitset> m_projectors;

  // owning vector of shared pointers to surfaces
  // @TODO: This might be problematic when appending a large number of surfaces
  // trackstates, because vector has to reallocated and thus copy. This might
  // be handled in a smart way by moving but not sure.
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  friend class detail_lt::TrackStateProxy<SourceLink, ParametersSize,
                                          MeasurementSizeMax, true>;
  friend class detail_lt::TrackStateProxy<SourceLink, ParametersSize,
                                          MeasurementSizeMax, false>;
};

}  // namespace Acts

#include "Acts/EventData/MultiTrajectory.ipp"
