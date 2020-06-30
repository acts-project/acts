// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLinkConcept.hpp"

#include <variant>

namespace Acts {

class Surface;

namespace MeasurementHelpers {

/// @brief Extract surface from a type erased measurement object
/// @tparam T The FittableMeasurement type
/// @return const pointer to the extracted surface
template <typename T>
const Surface* getSurface(const T& fittable_measurement) {
  return std::visit([](const auto& meas) { return &meas.referenceSurface(); },
                    fittable_measurement);
}

template <typename T>
size_t getSize(const T& fittable_measurement) {
  return std::visit([](const auto& meas) { return meas.size(); },
                    fittable_measurement);
}
}  // namespace MeasurementHelpers

struct MinimalSourceLink {
  const FittableMeasurement<MinimalSourceLink>* meas{nullptr};

  bool operator==(const MinimalSourceLink& rhs) const;

  const Surface& referenceSurface() const;

  const FittableMeasurement<MinimalSourceLink>& operator*() const;
};

inline std::ostream& operator<<(std::ostream& os, const MinimalSourceLink& sl) {
  os << "SourceLink(" << sl.meas << ")";
  return os;
}

static_assert(SourceLinkConcept<MinimalSourceLink>,
              "MinimalSourceLink does not fulfill SourceLinkConcept");

namespace detail {

/// Helper functor for @c visit_measurement. This is the actual functor given
/// to @c template_switch.
/// @tparam I Compile time int value
template <size_t I>
struct visit_measurement_callable {
  /// The invoked function. It will perform the head/top-left corner
  /// extraction, and pass thee results to the given lambda.
  /// @tparam L The lambda type
  /// @tparam A The parameter vector type
  /// @tparam B The covariance matrix type
  /// @note No requirements on @c A and @c B are made, to enable a single
  /// overload for both const and non-const matrices/vectors.
  /// @param param The parameter vector
  /// @param cov The covariance matrix
  /// @param lambda The lambda to call with the statically sized subsets
  template <typename L, typename A, typename B>
  auto static constexpr invoke(A&& param, B&& cov, L&& lambda) {
    return lambda(param.template head<I>(), cov.template topLeftCorner<I, I>());
  }
};
}  // namespace detail

/// Dispatch a lambda call on an overallocated parameter vector and covariance
/// matrix, based on a runtime dimension value. Inside the lambda call, the
/// vector and matrix will have fixed dimensions, but will still point back to
/// the originally given overallocated values.
/// @tparam L The lambda type
/// @tparam A The parameter vector type
/// @tparam B The covariance matrix type
/// @note No requirements on @c A and @c B are made, to enable a single
/// overload for both const and non-const matrices/vectors.
/// @param param The parameter vector
/// @param cov The covariance matrix
/// @param dim The actual dimension as a runtime value
/// @param lambda The lambda to call with the statically sized subsets
template <typename L, typename A, typename B>
auto visit_measurement(A&& param, B&& cov, size_t dim, L&& lambda) {
  return template_switch<detail::visit_measurement_callable, 1,
                         eBoundParametersSize>(dim, param, cov, lambda);
}

}  // namespace Acts
