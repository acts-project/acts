// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/HashedString.hpp"

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

namespace detail_tsp {

template <typename Derived, bool read_only>
class TrackStateProxyCommon {
 protected:
  using IndexType = Acts::TrackIndexType;
  // using ConstParameters = const_parameters_t;
  // using Parameters = parameters_t;
  // using ConstCovariance = const_covariance_t;
  // using Covariance = covariance_t;
  //
  constexpr Derived& derived() { return static_cast<Derived&>(*this); }
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

  /// Check for presence of calibrated measurement data.
  /// @return True if calibrated measurements exist.
  bool hasCalibrated() const {
    return derived().has(detail_tsp::kCalibratedKey);
  }

  using ParametersMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CoefficientsMap;
  using ConstParametersMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CoefficientsMap;
  using CovarianceMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CovarianceMap;
  using ConstCovarianceMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CovarianceMap;
  using EffectiveCalibratedMap =
      typename detail_tsp::DynamicSizeTypes<false>::CoefficientsMap;
  using ConstEffectiveCalibratedMap =
      typename detail_tsp::DynamicSizeTypes<true>::CoefficientsMap;
  using EffectiveCalibratedCovarianceMap =
      typename detail_tsp::DynamicSizeTypes<false>::CovarianceMap;
  using ConstEffectiveCalibratedCovarianceMap =
      typename detail_tsp::DynamicSizeTypes<true>::CovarianceMap;

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
};

}  // namespace detail_tsp
}  // namespace Acts
