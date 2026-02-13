// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/EigenConcepts.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <algorithm>
#include <ranges>
#include <string_view>

namespace Acts {

namespace detail_tsp {
inline constexpr HashedString kPreviousKey = hashString("previous");
inline constexpr HashedString kChi2Key = hashString("chi2");
inline constexpr HashedString kPathLengthKey = hashString("pathLength");
inline constexpr HashedString kTypeFlagsKey = hashString("typeFlags");
inline constexpr HashedString kPredictedKey = hashString("predicted");
inline constexpr HashedString kFilteredKey = hashString("filtered");
inline constexpr HashedString kSmoothedKey = hashString("smoothed");
inline constexpr HashedString kJacobianKey = hashString("jacobian");
inline constexpr HashedString kProjectorKey = hashString("projector");
inline constexpr HashedString kUncalibratedKey =
    hashString("uncalibratedSourceLink");
inline constexpr HashedString kCalibratedKey = hashString("calibrated");
inline constexpr HashedString kCalibratedCovKey = hashString("calibratedCov");
inline constexpr HashedString kMeasDimKey = hashString("measdim");
inline constexpr HashedString kNextKey = hashString("next");

}  // namespace detail_tsp

/// Common CRTP implementation shared by track state proxy front-ends.
/// This base class provides access to track state components including
/// predicted, filtered, and smoothed parameters and covariances, as well as
/// calibrated measurement data and various metadata. The derived proxy must
/// expose the underlying storage access methods.
///
/// @tparam Derived The proxy implementation inheriting from this base
/// @tparam read_only Whether the proxy provides mutable access
template <typename Derived, bool read_only>
class TrackStateProxyCommon {
 protected:
  /// Index type for track states
  using IndexType = Acts::TrackIndexType;

  /// Cast to derived class
  /// @return Reference to derived proxy
  constexpr Derived& derived() { return static_cast<Derived&>(*this); }
  /// Cast to derived class (const)
  /// @return Const reference to derived proxy
  constexpr const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

 public:
  /// Retrieve the previous track state index in the linked trajectory.
  /// @return Index of the previous state or `kTrackIndexInvalid`.
  TrackIndexType previous() const {
    return derived()
        .template component<TrackIndexType, detail_tsp::kPreviousKey>();
  }

  /// Retrieve a mutable reference to the previous track state index.
  /// @return Mutable index of the previous state.
  TrackIndexType& previous()
    requires(!read_only)
  {
    return derived()
        .template component<TrackIndexType, detail_tsp::kPreviousKey>();
  }

  /// Check whether this state links to a previous state.
  /// @return True if the previous index is valid.
  bool hasPrevious() const { return previous() != kTrackIndexInvalid; }

  /// Check for presence of predicted track parameters.
  /// @return True if the predicted component exists.
  bool hasPredicted() const { return derived().has(detail_tsp::kPredictedKey); }

  /// Check for presence of filtered track parameters.
  /// @return True if the filtered component exists.
  bool hasFiltered() const { return derived().has(detail_tsp::kFilteredKey); }

  /// Check for presence of smoothed track parameters.
  /// @return True if the smoothed component exists.
  bool hasSmoothed() const { return derived().has(detail_tsp::kSmoothedKey); }

  /// Check for presence of a transport Jacobian.
  /// @return True if a Jacobian is stored.
  bool hasJacobian() const { return derived().has(detail_tsp::kJacobianKey); }

  /// Check for presence of a measurement projector.
  /// @return True if projector indices are stored.
  bool hasProjector() const { return derived().has(detail_tsp::kProjectorKey); }

  /// Check for presence of an uncalibrated source link.
  /// @return True if an uncalibrated source link is stored.
  bool hasUncalibratedSourceLink() const {
    return derived().has(detail_tsp::kUncalibratedKey);
  }

  /// Check for presence of calibrated measurement data.
  /// @return True if calibrated measurements exist.
  bool hasCalibrated() const {
    return derived().has(detail_tsp::kCalibratedKey);
  }

  /// @name Type aliases for parameter and covariance access
  /// @{

  /// Mutable Eigen map type for bound track parameters.
  using ParametersMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CoefficientsMap;
  /// Const Eigen map type for bound track parameters.
  using ConstParametersMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CoefficientsMap;
  /// Mutable Eigen map type for bound track parameter covariance.
  using CovarianceMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CovarianceMap;
  /// Const Eigen map type for bound track parameter covariance.
  using ConstCovarianceMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CovarianceMap;

  /// Mutable Eigen map type for calibrated measurements (dynamic size).
  using EffectiveCalibratedMap =
      typename detail_tsp::DynamicSizeTypes<false>::CoefficientsMap;
  /// Const Eigen map type for calibrated measurements (dynamic size).
  using ConstEffectiveCalibratedMap =
      typename detail_tsp::DynamicSizeTypes<true>::CoefficientsMap;
  /// Mutable Eigen map type for calibrated measurement covariance (dynamic
  /// size).
  using EffectiveCalibratedCovarianceMap =
      typename detail_tsp::DynamicSizeTypes<false>::CovarianceMap;
  /// Const Eigen map type for calibrated measurement covariance (dynamic size).
  using ConstEffectiveCalibratedCovarianceMap =
      typename detail_tsp::DynamicSizeTypes<true>::CovarianceMap;

  /// @}

  /// Access the predicted parameter vector.
  /// @return Bound parameter map for the predicted state.
  ConstParametersMap predicted() const {
    assert(hasPredicted());
    const auto parIndex =
        derived().template component<IndexType>(detail_tsp::kPredictedKey);
    return derived().parametersAtIndex(parIndex);
  }

  /// Access the predicted parameter vector.
  /// @return Mutable bound parameter map for the predicted state.
  ParametersMap predicted()
    requires(!read_only)
  {
    assert(hasPredicted());
    const auto parIndex =
        derived().template component<IndexType>(detail_tsp::kPredictedKey);
    return derived().parametersAtIndexMutable(parIndex);
  }

  /// Access the predicted covariance matrix.
  /// @return Bound covariance map for the predicted state.
  ConstCovarianceMap predictedCovariance() const {
    assert(hasPredicted());
    const auto covIndex =
        derived().template component<IndexType>(detail_tsp::kPredictedKey);
    return derived().covarianceAtIndex(covIndex);
  }

  /// Access the predicted covariance matrix.
  /// @return Mutable bound covariance map for the predicted state.
  CovarianceMap predictedCovariance()
    requires(!read_only)
  {
    assert(hasPredicted());
    const auto covIndex =
        derived().template component<IndexType>(detail_tsp::kPredictedKey);
    return derived().covarianceAtIndexMutable(covIndex);
  }

  /// Access the filtered parameter vector.
  /// @return Bound parameter map for the filtered state.
  ConstParametersMap filtered() const {
    assert(hasFiltered());
    const auto parIndex =
        derived().template component<IndexType>(detail_tsp::kFilteredKey);
    return derived().parametersAtIndex(parIndex);
  }

  /// Access the filtered parameter vector.
  /// @return Mutable bound parameter map for the filtered state.
  ParametersMap filtered()
    requires(!read_only)
  {
    const auto parIndex =
        derived().template component<IndexType>(detail_tsp::kFilteredKey);
    return derived().parametersAtIndexMutable(parIndex);
  }

  /// Access the filtered covariance matrix.
  /// @return Bound covariance map for the filtered state.
  ConstCovarianceMap filteredCovariance() const {
    assert(hasFiltered());
    const auto covIndex =
        derived().template component<IndexType>(detail_tsp::kFilteredKey);
    return derived().covarianceAtIndex(covIndex);
  }

  /// Access the filtered covariance matrix.
  /// @return Mutable bound covariance map for the filtered state.
  CovarianceMap filteredCovariance()
    requires(!read_only)
  {
    const auto covIndex =
        derived().template component<IndexType>(detail_tsp::kFilteredKey);
    return derived().covarianceAtIndexMutable(covIndex);
  }

  /// Access the smoothed parameter vector.
  /// @return Bound parameter map for the smoothed state.
  ConstParametersMap smoothed() const {
    assert(hasSmoothed());
    const auto parIndex =
        derived().template component<IndexType>(detail_tsp::kSmoothedKey);
    return derived().parametersAtIndex(parIndex);
  }

  /// Access the smoothed parameter vector.
  /// @return Mutable bound parameter map for the smoothed state.
  ParametersMap smoothed()
    requires(!read_only)
  {
    assert(hasSmoothed());
    const auto parIndex =
        derived().template component<IndexType>(detail_tsp::kSmoothedKey);
    return derived().parametersAtIndexMutable(parIndex);
  }

  /// Access the smoothed covariance matrix.
  /// @return Bound covariance map for the smoothed state.
  ConstCovarianceMap smoothedCovariance() const {
    assert(hasSmoothed());
    const auto covIndex =
        derived().template component<IndexType>(detail_tsp::kSmoothedKey);
    return derived().covarianceAtIndex(covIndex);
  }

  /// Access the smoothed covariance matrix.
  /// @return Mutable bound covariance map for the smoothed state.
  CovarianceMap smoothedCovariance()
    requires(!read_only)
  {
    assert(hasSmoothed());
    const auto covIndex =
        derived().template component<IndexType>(detail_tsp::kSmoothedKey);
    return derived().covarianceAtIndexMutable(covIndex);
  }

  /// Retrieve the accumulated path length.
  /// @return Path length stored on the state.
  double pathLength() const {
    return derived().template component<double, detail_tsp::kPathLengthKey>();
  }

  /// Retrieve a mutable reference to the accumulated path length.
  /// @return Mutable path length.
  double& pathLength()
    requires(!read_only)
  {
    return derived().template component<double, detail_tsp::kPathLengthKey>();
  }

  /// Retrieve the track-state type flags.
  /// @return Bit mask describing the state type.
  ConstTrackStateTypeMap typeFlags() const {
    const auto& raw = derived()
                          .template component<ConstTrackStateTypeMap::raw_type,
                                              detail_tsp::kTypeFlagsKey>();
    return ConstTrackStateTypeMap{raw};
  }

  /// Retrieve mutable track-state type flags.
  /// @return Mutable bit mask describing the state type.
  MutableTrackStateTypeMap typeFlags()
    requires(!read_only)
  {
    auto& raw = derived()
                    .template component<MutableTrackStateTypeMap::raw_type,
                                        detail_tsp::kTypeFlagsKey>();
    return MutableTrackStateTypeMap{raw};
  }

  /// Retrieve the local chi2 contribution.
  /// @return Chi2 value associated with this state.
  float chi2() const {
    return derived().template component<float, detail_tsp::kChi2Key>();
  }

  /// Retrieve a mutable reference to the local chi2 contribution.
  /// @return Mutable chi2 value.
  float& chi2()
    requires(!read_only)
  {
    return derived().template component<float, detail_tsp::kChi2Key>();
  }

  /// Decode the measurement projector indices.
  /// @return Bound parameter indices used for projection.
  BoundSubspaceIndices projectorSubspaceIndices() const {
    assert(hasProjector());
    const auto& serialized =
        derived()
            .template component<SerializedSubspaceIndices,
                                detail_tsp::kProjectorKey>();
    return deserializeSubspaceIndices<eBoundSize>(serialized);
  }

  /// Returns the projector subspace indices
  /// @return The projector subspace indices
  template <std::size_t measdim>
  SubspaceIndices<measdim> projectorSubspaceIndices() const {
    BoundSubspaceIndices boundSubspace = projectorSubspaceIndices();
    SubspaceIndices<measdim> subspace;
    std::copy(boundSubspace.begin(), boundSubspace.begin() + measdim,
              subspace.begin());
    return subspace;
  }

  /// Creates a variable size subspace helper
  /// @return The subspace helper
  VariableBoundSubspaceHelper projectorSubspaceHelper() const {
    BoundSubspaceIndices boundSubspace = projectorSubspaceIndices();
    std::span<std::uint8_t> validSubspaceIndices(
        boundSubspace.begin(),
        boundSubspace.begin() + derived().calibratedSize());
    return VariableBoundSubspaceHelper(validSubspaceIndices);
  }

  /// Creates a fixed size subspace helper
  /// @return The subspace helper
  template <std::size_t measdim>
  FixedBoundSubspaceHelper<measdim> projectorSubspaceHelper() const {
    SubspaceIndices<measdim> subspace = projectorSubspaceIndices<measdim>();
    return FixedBoundSubspaceHelper<measdim>(subspace);
  }

  /// Store subspace indices describing the measurement projector.
  /// @tparam index_range_t Range of indices to encode.
  /// @param subspaceIndices Collection of bound indices forming the projector rows.
  template <std::ranges::sized_range index_range_t>
  void setProjectorSubspaceIndices(const index_range_t& subspaceIndices)
    requires(!read_only &&
             std::convertible_to<std::ranges::range_value_t<index_range_t>,
                                 std::uint8_t>)
  {
    assert(derived().template has<detail_tsp::kProjectorKey>());
    assert(subspaceIndices.size() <= eBoundSize);
    BoundSubspaceIndices boundSubspace{};
    std::transform(subspaceIndices.begin(), subspaceIndices.end(),
                   boundSubspace.begin(),
                   [](auto i) { return static_cast<std::uint8_t>(i); });
    derived()
        .template component<SerializedSubspaceIndices,
                            detail_tsp::kProjectorKey>() =
        serializeSubspaceIndices(boundSubspace);
  }

  /// Compute the property mask describing which components are present.
  /// @return Bit mask of available properties.
  TrackStatePropMask getMask() const {
    using PM = TrackStatePropMask;
    PM mask = PM::None;
    if (hasPredicted()) {
      mask |= PM::Predicted;
    }
    if (hasFiltered()) {
      mask |= PM::Filtered;
    }
    if (hasSmoothed()) {
      mask |= PM::Smoothed;
    }
    if (hasJacobian()) {
      mask |= PM::Jacobian;
    }
    if (hasCalibrated()) {
      mask |= PM::Calibrated;
    }
    return mask;
  }

  /// Access the best available parameters (smoothed, filtered, or predicted).
  /// @return Bound parameter map for the state.
  ConstParametersMap parameters() const {
    if (hasSmoothed()) {
      return smoothed();
    } else if (hasFiltered()) {
      return filtered();
    }
    return predicted();
  }

  /// Access the best available covariance (smoothed, filtered, or predicted).
  /// @return Bound covariance map for the state.
  ConstCovarianceMap covariance() const {
    if (hasSmoothed()) {
      return smoothedCovariance();
    } else if (hasFiltered()) {
      return filteredCovariance();
    }
    return predictedCovariance();
  }

  /// Access the calibrated measurement values with runtime dimension.
  /// @return Eigen map referencing the calibrated measurement vector.
  ConstEffectiveCalibratedMap effectiveCalibrated() const {
    const double* data = derived().calibratedData();
    const auto size = derived().calibratedSize();
    return ConstEffectiveCalibratedMap{data, size};
  }

  /// Access mutable calibrated measurement values with runtime dimension.
  /// @return Eigen map referencing the calibrated measurement vector.
  EffectiveCalibratedMap effectiveCalibrated()
    requires(!read_only)
  {
    double* data = derived().calibratedDataMutable();
    const auto size = derived().calibratedSize();
    return EffectiveCalibratedMap{data, size};
  }

  /// Access the calibrated covariance with runtime dimension.
  /// @return Eigen map referencing the measurement covariance matrix.
  ConstEffectiveCalibratedCovarianceMap effectiveCalibratedCovariance() const {
    const double* data = derived().calibratedCovarianceData();
    const auto size = derived().calibratedSize();
    return ConstEffectiveCalibratedCovarianceMap{data, size, size};
  }

  /// Access mutable calibrated covariance with runtime dimension.
  /// @return Eigen map referencing the measurement covariance matrix.
  EffectiveCalibratedCovarianceMap effectiveCalibratedCovariance()
    requires(!read_only)
  {
    double* data = derived().calibratedCovarianceDataMutable();
    const auto size = derived().calibratedSize();
    return EffectiveCalibratedCovarianceMap{data, size, size};
  }

  /// Access calibrated measurement data with compile-time dimension.
  /// @tparam measdim Measurement dimension.
  /// @return Eigen map referencing the calibrated measurement vector.
  template <std::size_t measdim>
  typename TrackStateTraits<measdim, true>::Calibrated calibrated() const {
    assert(derived().calibratedSize() == static_cast<TrackIndexType>(measdim));
    const double* data = derived().calibratedData();
    return typename TrackStateTraits<measdim, true>::Calibrated(data);
  }

  /// Access calibrated measurement data with compile-time dimension.
  /// @tparam measdim Measurement dimension.
  /// @return Mutable Eigen map referencing the calibrated measurement vector.
  template <std::size_t measdim>
  typename TrackStateTraits<measdim, false>::Calibrated calibrated()
    requires(!read_only)
  {
    assert(derived().calibratedSize() == static_cast<TrackIndexType>(measdim));
    double* data = derived().calibratedDataMutable();
    return typename TrackStateTraits<measdim, false>::Calibrated(data);
  }

  /// Access calibrated covariance data with compile-time dimension.
  /// @tparam measdim Measurement dimension.
  /// @return Eigen map referencing the covariance matrix.
  template <std::size_t measdim>
  typename TrackStateTraits<measdim, true>::CalibratedCovariance
  calibratedCovariance() const {
    assert(derived().calibratedSize() == static_cast<TrackIndexType>(measdim));
    const double* data = derived().calibratedCovarianceData();
    return typename TrackStateTraits<measdim, true>::CalibratedCovariance(data);
  }

  /// Access calibrated covariance data with compile-time dimension.
  /// @tparam measdim Measurement dimension.
  /// @return Mutable Eigen map referencing the covariance matrix.
  template <std::size_t measdim>
  typename TrackStateTraits<measdim, false>::CalibratedCovariance
  calibratedCovariance()
    requires(!read_only)
  {
    assert(derived().calibratedSize() == static_cast<TrackIndexType>(measdim));
    double* data = derived().calibratedCovarianceDataMutable();
    return
        typename TrackStateTraits<measdim, false>::CalibratedCovariance(data);
  }

  /// Allocate and initialize calibrated data from static-size Eigen objects.
  /// @tparam val_t Eigen vector type holding calibrated values.
  /// @tparam cov_t Eigen matrix type holding the covariance.
  /// @param val Vector to copy into the calibrated storage.
  /// @param cov Covariance matrix to copy into the calibrated storage.
  template <typename val_t, typename cov_t>
  void allocateCalibrated(const Eigen::DenseBase<val_t>& val,
                          const Eigen::DenseBase<cov_t>& cov)
    requires(!read_only && Concepts::eigen_base_is_fixed_size<val_t> &&
             Concepts::eigen_bases_have_same_num_rows<val_t, cov_t> &&
             Concepts::eigen_base_is_square<cov_t> &&
             Eigen::PlainObjectBase<val_t>::RowsAtCompileTime <=
                 static_cast<std::underlying_type_t<BoundIndices>>(eBoundSize))
  {
    constexpr std::size_t measdim =
        static_cast<std::size_t>(val_t::RowsAtCompileTime);
    derived().allocateCalibrated(measdim);
    calibrated<measdim>() = val;
    calibratedCovariance<measdim>() = cov;
  }
};

}  // namespace Acts
