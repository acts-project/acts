// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <concepts>
#include <cstdint>
#include <ranges>
#include <type_traits>
#include <utility>

namespace Acts {

namespace detail {

/// The parameter vector type of a measurement like type
/// @tparam T the measurement like type
template <typename T>
using MeasurementParametersType =
    std::decay_t<decltype(std::declval<const T&>().parameters())>;

/// The covariance matrix type of a measurement like type
/// @tparam T the measurement like type
template <typename T>
using MeasurementCovarianceType =
    std::decay_t<decltype(std::declval<const T&>().covariance())>;

}  // namespace detail

/// Requirements on a calibrated measurement expressed in a subspace of a
/// parameter space.
///
/// This is the contract between a track finder and whatever the client uses to
/// represent a calibrated measurement. It is deliberately minimal so that
/// existing measurement proxies satisfy it without a copy, and so that clients
/// are free to return their own type.
///
/// The dimension may be known at compile time or only at runtime; see
/// @ref MeasurementSizeAtCompileTime and @ref StaticMeasurementConcept. Note
/// that @ref subspaceIndices must yield exactly @c size() entries, not a padded
/// range.
///
/// @tparam T the measurement like type
template <typename T>
concept MeasurementConcept = requires(const T& measurement) {
  /// The measurement dimension
  { measurement.size() } -> std::integral;
  /// The parameters of the parameter space the measurement constrains,
  /// exactly `size()` of them
  { measurement.subspaceIndices() } -> std::ranges::sized_range;
  /// The measured values
  { measurement.parameters() };
  /// The covariance of the measured values
  { measurement.covariance() };

  requires std::convertible_to<
      std::ranges::range_value_t<decltype(measurement.subspaceIndices())>,
      std::uint8_t>;
  requires detail::MeasurementParametersType<T>::ColsAtCompileTime == 1;
  requires detail::MeasurementParametersType<T>::RowsAtCompileTime ==
               detail::MeasurementCovarianceType<T>::RowsAtCompileTime;
};

/// The dimension of a measurement if it is known at compile time,
/// `Eigen::Dynamic` otherwise.
///
/// @tparam measurement_t the measurement type
template <MeasurementConcept measurement_t>
constexpr int MeasurementSizeAtCompileTime =
    detail::MeasurementParametersType<measurement_t>::RowsAtCompileTime;

/// A @ref MeasurementConcept whose dimension is known at compile time.
///
/// Algorithms can use this to skip the runtime dispatch over the measurement
/// dimension, which is what makes it worthwhile for a client to hand out a
/// fixed size measurement type in the first place.
///
/// @tparam measurement_t the measurement type
template <typename measurement_t>
concept StaticMeasurementConcept =
    MeasurementConcept<measurement_t> &&
    (MeasurementSizeAtCompileTime<measurement_t> != Eigen::Dynamic);

}  // namespace Acts
