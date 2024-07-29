// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/detail/CalculateResiduals.hpp"
#include "Acts/EventData/detail/ParameterTraits.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <variant>
#include <vector>

namespace ActsExamples {

/// A measurement of a fixed-size subspace of the full parameters.
///
/// @tparam indices_t Parameter index type, determines the full parameter space
/// @tparam kSize Size of the parameter subset.
///
/// The measurement intentionally does not store a pointer/reference to the
/// reference object in the geometry hierarchy, i.e. the surface or volume. The
/// reference object can already be identified via the geometry identifier
/// provided by the source link. Since a measurement **must** be anchored within
/// the geometry hierarchy, all measurement surfaces and volumes **must**
/// provide valid geometry identifiers. In all use-cases, e.g. Kalman filtering,
/// a pointer/reference to the reference object is available before the
/// measurement is accessed; e.g. the propagator provides the surface pointer
/// during navigation, which is then used to lookup possible measurements.
///
/// The pointed-to geometry object would differ depending on the parameter type.
/// This means either, that there needs to be an additional variable type or
/// that a pointer to a base object is stored (requiring a `dynamic_cast` later
/// on). Both variants add additional complications. Since the geometry object
/// is not required anyway (as discussed above), not storing it removes all
/// these complications altogether.
template <typename indices_t, std::size_t kSize>
class FixedSizeMeasurement {
  static constexpr std::size_t kFullSize =
      Acts::detail::kParametersSize<indices_t>;

 public:
  using Scalar = Acts::ActsScalar;
  /// Vector type containing for measured parameter values.
  using ParametersVector = Acts::ActsVector<kSize>;
  /// Matrix type for the measurement covariance.
  using CovarianceMatrix = Acts::ActsSquareMatrix<kSize>;
  /// Vector type containing all parameters in the same space.
  using FullParametersVector = Acts::ActsVector<kFullSize>;
  using ProjectionMatrix = Acts::ActsMatrix<kSize, kFullSize>;
  using ExpansionMatrix = Acts::ActsMatrix<kFullSize, kSize>;

  /// Construct from source link, subset indices, and measured data.
  ///
  /// @tparam parameters_t Input parameters vector type
  /// @tparam covariance_t Input covariance matrix type
  /// @param source The link that connects to the underlying detector readout
  /// @param indices Which parameters are measured
  /// @param params Measured parameters values
  /// @param cov Measured parameters covariance
  ///
  /// @note The indices must be ordered and must describe/match the content
  ///   of parameters and covariance.
  template <typename parameters_t, typename covariance_t>
  FixedSizeMeasurement(Acts::SourceLink source,
                       const std::array<indices_t, kSize>& indices,
                       const Eigen::MatrixBase<parameters_t>& params,
                       const Eigen::MatrixBase<covariance_t>& cov)
      : m_source(std::move(source)), m_params(params), m_cov(cov) {
    for (std::size_t i = 0u; i < kSize; ++i) {
      assert(indices[i] >= 0 && indices[i] < kFullSize &&
             "Index out of bounds");
      m_indices[i] = static_cast<std::uint8_t>(indices[i]);
    }
  }
  /// A measurement can only be constructed with valid parameters.
  FixedSizeMeasurement() = delete;
  FixedSizeMeasurement(const FixedSizeMeasurement&) = default;
  FixedSizeMeasurement(FixedSizeMeasurement&&) = default;
  ~FixedSizeMeasurement() = default;
  FixedSizeMeasurement& operator=(const FixedSizeMeasurement&) = default;
  FixedSizeMeasurement& operator=(FixedSizeMeasurement&&) = default;

  /// Source link that connects to the underlying detector readout.
  const Acts::SourceLink& sourceLink() const { return m_source; }

  /// Number of measured parameters.
  static constexpr std::size_t size() { return kSize; }

  /// Check if a specific parameter is part of this measurement.
  bool contains(indices_t i) const {
    return std::find(m_indices.begin(), m_indices.end(), i) != m_indices.end();
  }

  /// The measurement indices
  constexpr std::array<Acts::BoundIndices, kSize> indices() const {
    std::array<Acts::BoundIndices, kSize> result;
    for (std::size_t i = 0u; i < kSize; ++i) {
      result[i] = static_cast<Acts::BoundIndices>(m_indices[i]);
    }
    return result;
  }

  /// Measured parameters values.
  const ParametersVector& parameters() const { return m_params; }

  /// Measured parameters covariance.
  const CovarianceMatrix& covariance() const { return m_cov; }

  /// Projection matrix from the full space into the measured subspace.
  ProjectionMatrix projector() const {
    ProjectionMatrix proj = ProjectionMatrix::Zero();
    for (std::size_t i = 0u; i < kSize; ++i) {
      proj(i, m_indices[i]) = 1;
    }
    return proj;
  }

  /// Expansion matrix from the measured subspace into the full space.
  ///
  /// This is equivalent to the transpose of the projection matrix only in the
  /// case of a trivial projection matrix. While this is the case here, it is
  /// still recommended to use the expansion matrix directly in cases where it
  /// is explicitly used.
  ExpansionMatrix expander() const {
    ExpansionMatrix expn = ExpansionMatrix::Zero();
    for (std::size_t i = 0u; i < kSize; ++i) {
      expn(m_indices[i], i) = 1;
    }
    return expn;
  }

  /// Compute residuals in the measured subspace.
  ///
  /// @param reference Reference parameters in the full space.
  ///
  /// This computes the difference `measured - reference` taking into account
  /// the allowed parameter ranges. Only the reference values in the measured
  /// subspace are used for the computation.
  ParametersVector residuals(const FullParametersVector& reference) const {
    ParametersVector res = ParametersVector::Zero();
    Acts::detail::calculateResiduals(static_cast<indices_t>(kSize), m_indices,
                                     reference, m_params, res);
    return res;
  }

  std::ostream& operator<<(std::ostream& os) const {
    Acts::detail::printMeasurement(os, static_cast<indices_t>(kSize),
                                   m_indices.data(), m_params.data(),
                                   m_cov.data());
    return os;
  }

 private:
  Acts::SourceLink m_source;
  std::array<std::uint8_t, kSize> m_indices;
  ParametersVector m_params;
  CovarianceMatrix m_cov;
};

/// Construct a fixed-size measurement for the given indices.
///
/// @tparam parameters_t Input parameters vector type
/// @tparam covariance_t Input covariance matrix type
/// @tparam indices_t Parameter index type, determines the full parameter space
/// @tparam tail_indices_t Helper types required to support variadic arguments;
///   all types must be convertibale to `indices_t`.
/// @param source The link that connects to the underlying detector readout
/// @param params Measured parameters values
/// @param cov Measured parameters covariance
/// @param index0 Required parameter index, measurement must be at least 1d
/// @param tailIndices Additional parameter indices for larger measurements
/// @return Fixed-size measurement w/ the correct type and given inputs
///
/// This helper function can be used to create a fixed-size measurement using an
/// explicit set of indices, e.g.
///
///     auto m = makeMeasurement(s, p, c, eBoundLoc0, eBoundTime);
///
/// for a 2d measurement w/ one position and time.
///
/// @note The indices must be ordered and must be consistent with the content of
///   parameters and covariance.
template <typename parameters_t, typename covariance_t, typename indices_t,
          typename... tail_indices_t>
auto makeFixedSizeMeasurement(Acts::SourceLink source,
                              const Eigen::MatrixBase<parameters_t>& params,
                              const Eigen::MatrixBase<covariance_t>& cov,
                              indices_t index0, tail_indices_t... tailIndices)
    -> FixedSizeMeasurement<indices_t, 1u + sizeof...(tail_indices_t)> {
  using IndexContainer = std::array<indices_t, 1u + sizeof...(tail_indices_t)>;
  return {std::move(source), IndexContainer{index0, tailIndices...}, params,
          cov};
}

namespace detail {
/// @cond

// Recursive construction of the measurement variant. `kN` is counted down until
// zero while the sizes are accumulated in the parameter pack.
//
// Example:
//
//        VariantMeasurementGenerator<..., 4>
//     -> VariantMeasurementGenerator<..., 3, 4>
//     -> VariantMeasurementGenerator<..., 2, 3, 4>
//     -> VariantMeasurementGenerator<..., 1, 2, 3, 4>
//     -> VariantMeasurementGenerator<..., 0, 1, 2, 3, 4>
//
template <typename indices_t, std::size_t kN, std::size_t... kSizes>
struct VariantMeasurementGenerator
    : VariantMeasurementGenerator<indices_t, kN - 1u, kN, kSizes...> {};
template <typename indices_t, std::size_t... kSizes>
struct VariantMeasurementGenerator<indices_t, 0u, kSizes...> {
  using Type = std::variant<FixedSizeMeasurement<indices_t, kSizes>...>;
};

/// @endcond
}  // namespace detail

/// Variant that can contain all possible measurements in a parameter space.
///
/// @tparam indices_t Parameter index type, determines the full parameter space
template <typename indices_t>
using VariantMeasurement = typename detail::VariantMeasurementGenerator<
    indices_t, Acts::detail::kParametersSize<indices_t>>::Type;

/// Variant that can hold all possible bound measurements.
///
using BoundVariantMeasurement = VariantMeasurement<Acts::BoundIndices>;

/// Variant that can hold all possible free measurements.
///
using FreeVariantMeasurement = VariantMeasurement<Acts::FreeIndices>;

template <typename indices_t>
std::ostream& operator<<(std::ostream& os,
                         const VariantMeasurement<indices_t>& vm) {
  return std::visit([&](const auto& m) { return (os << m); }, vm);
}

/// Variable measurement type that can contain all possible combinations.
using Measurement = BoundVariantMeasurement;

/// Container of measurements.
///
/// In contrast to the source links, the measurements themself must not be
/// orderable. The source links stored in the measurements are treated
/// as opaque here and no ordering is enforced on the stored measurements.
using MeasurementContainer = std::vector<Measurement>;

}  // namespace ActsExamples
