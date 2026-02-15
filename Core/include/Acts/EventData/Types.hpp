// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cstdint>
#include <limits>
#include <span>

namespace Acts {

/// Type alias for track index values
using TrackIndexType = std::uint32_t;

/// Sentinel value for an invalid / unset track EDM related index
static constexpr TrackIndexType kTrackIndexInvalid =
    std::numeric_limits<TrackIndexType>::max();

/// Maximum size of measurement dimension supported
static constexpr std::size_t kMeasurementSizeMax = eBoundSize;

/// @brief Type alias for bitset representing parameter projections
/// @details Used to store parameter projection information in a compact bit format
using ProjectorBitset = std::uint64_t;

/// Type alias for subspace index
using SubspaceIndex = std::uint8_t;
/// Template alias for subspace indices array
template <std::size_t measdim>
using SubspaceIndices = std::array<SubspaceIndex, measdim>;
/// Type alias for bound parameter subspace indices
using BoundSubspaceIndices = SubspaceIndices<eBoundSize>;

// template <std::size_t measdim>
/// @brief type alias for indices of bound track parameters subspace
/// @details used to specify which components of the bound track parameters are being referenced
// using BoundSubspaceIndices = std::array<std::size_t, eBoundSize>;
static constexpr BoundSubspaceIndices kBoundSubspaceIndicesInvalid = {
    eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize};

/// @brief Type alias for serialized subspace indices
/// @details Compact representation of subspace indices as a 64-bit unsigned integer
using SerializedSubspaceIndices = std::uint64_t;

/// Index type for spacepoints
using SpacePointIndex2 = std::uint32_t;
/// Range of spacepoint indices defined by a pair of start and end indices
using SpacePointIndexRange2 = std::pair<SpacePointIndex2, SpacePointIndex2>;
/// Subset of spacepoint indices represented as a span
using SpacePointIndexSubset2 = std::span<const SpacePointIndex2>;

/// Index type for seeds
using SeedIndex2 = std::uint32_t;

namespace detail_tsp {

template <std::size_t Size, bool ReadOnlyMaps = true>
struct FixedSizeTypes {
  constexpr static auto Flags = Eigen::ColMajor | Eigen::AutoAlign;

  // single items
  using Coefficients = Eigen::Matrix<double, Size, 1, Flags>;
  using Covariance = Eigen::Matrix<double, Size, Size, Flags>;
  using CoefficientsMap = Eigen::Map<const_if_t<ReadOnlyMaps, Coefficients>>;
  using CovarianceMap = Eigen::Map<const_if_t<ReadOnlyMaps, Covariance>>;

  using DynamicCoefficients = Eigen::Matrix<double, Eigen::Dynamic, 1, Flags>;
  using DynamicCovariance =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Flags>;
  using DynamicCoefficientsMap =
      Eigen::Map<const_if_t<ReadOnlyMaps, DynamicCoefficients>>;
  using DynamicCovarianceMap =
      Eigen::Map<const_if_t<ReadOnlyMaps, DynamicCovariance>>;
};

// Type construction helper for dynamic sized coefficients and associated
/// covariances.
template <bool ReadOnlyMaps = true>
struct DynamicSizeTypes {
  constexpr static auto Flags = Eigen::ColMajor | Eigen::AutoAlign;

  using Coefficients = Eigen::Matrix<double, Eigen::Dynamic, 1, Flags>;
  using Covariance =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Flags>;
  using CoefficientsMap = Eigen::Map<const_if_t<ReadOnlyMaps, Coefficients>>;
  using CovarianceMap = Eigen::Map<const_if_t<ReadOnlyMaps, Covariance>>;
};

}  // namespace detail_tsp

/// Type aliases describing track state component mappings.
template <std::size_t M, bool ReadOnly = true>
struct TrackStateTraits {
  /// Track parameters vector type
  using Parameters =
      typename detail_tsp::FixedSizeTypes<eBoundSize,
                                          ReadOnly>::CoefficientsMap;
  /// Track parameters covariance matrix type
  using Covariance =
      typename detail_tsp::FixedSizeTypes<eBoundSize, ReadOnly>::CovarianceMap;
  /// Calibrated measurement vector type
  using Calibrated =
      typename detail_tsp::FixedSizeTypes<M, ReadOnly>::CoefficientsMap;
  /// Calibrated measurement covariance matrix type
  using CalibratedCovariance =
      typename detail_tsp::FixedSizeTypes<M, ReadOnly>::CovarianceMap;
  /// Effective calibrated measurement vector type
  using EffectiveCalibrated =
      typename detail_tsp::DynamicSizeTypes<ReadOnly>::CoefficientsMap;
  /// Effective calibrated measurement covariance matrix type
  using EffectiveCalibratedCovariance =
      typename detail_tsp::DynamicSizeTypes<ReadOnly>::CovarianceMap;
};

}  // namespace Acts
