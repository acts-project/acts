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

/// Either type T or const T depending on the boolean.
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

/// Type construction helper for fixed size coefficients and associated
/// covariances.
template <std::size_t Size, bool ReadOnlyMaps = true>
struct FixedSizeTypes {
  constexpr static auto Flags = Eigen::ColMajor | Eigen::AutoAlign;

  // single items
  using Coefficients = Eigen::Matrix<double, Size, 1, Flags>;
  using Covariance = Eigen::Matrix<double, Size, Size, Flags>;
  using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
  using CovarianceMap = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>>;

  using DynamicCoefficients = Eigen::Matrix<double, Eigen::Dynamic, 1, Flags>;
  using DynamicCovariance =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Flags>;
  using DynamicCoefficientsMap =
      Eigen::Map<ConstIf<DynamicCoefficients, ReadOnlyMaps>>;
  using DynamicCovarianceMap =
      Eigen::Map<ConstIf<DynamicCovariance, ReadOnlyMaps>>;
};

// Type construction helper for dynamic sized coefficients and associated
/// covariances.
template <bool ReadOnlyMaps = true>
struct DynamicSizeTypes {
  constexpr static auto Flags = Eigen::ColMajor | Eigen::AutoAlign;

  using Coefficients = Eigen::Matrix<double, Eigen::Dynamic, 1, Flags>;
  using Covariance =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Flags>;
  using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
  using CovarianceMap = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>>;
};

}  // namespace detail_tsp

template <std::size_t M, bool ReadOnly = true>
struct TrackStateTraits {
  using Parameters =
      typename detail_tsp::FixedSizeTypes<eBoundSize,
                                          ReadOnly>::CoefficientsMap;
  using Covariance =
      typename detail_tsp::FixedSizeTypes<eBoundSize, ReadOnly>::CovarianceMap;
  using Calibrated =
      typename detail_tsp::FixedSizeTypes<M, ReadOnly>::CoefficientsMap;
  using CalibratedCovariance =
      typename detail_tsp::FixedSizeTypes<M, ReadOnly>::CovarianceMap;
  using EffectiveCalibrated =
      typename detail_tsp::DynamicSizeTypes<ReadOnly>::CoefficientsMap;
  using EffectiveCalibratedCovariance =
      typename detail_tsp::DynamicSizeTypes<ReadOnly>::CovarianceMap;
};

}  // namespace Acts
