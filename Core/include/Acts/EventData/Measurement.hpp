// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
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
#include "Acts/Utilities/detail/Subspace.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <variant>

namespace Acts {

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
class Measurement {
  static constexpr std::size_t kFullSize = detail::kParametersSize<indices_t>;

  using Subspace = detail::FixedSizeSubspace<kFullSize, kSize>;

 public:
  using Scalar = ActsScalar;
  /// Vector type containing for measured parameter values.
  using ParametersVector = ActsVector<kSize>;
  /// Matrix type for the measurement covariance.
  using CovarianceMatrix = ActsSquareMatrix<kSize>;
  /// Vector type containing all parameters in the same space.
  using FullParametersVector = ActsVector<kFullSize>;
  using ProjectionMatrix = ActsMatrix<kSize, kFullSize>;
  using ExpansionMatrix = ActsMatrix<kFullSize, kSize>;

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
  Measurement(SourceLink source, const std::array<indices_t, kSize>& indices,
              const Eigen::MatrixBase<parameters_t>& params,
              const Eigen::MatrixBase<covariance_t>& cov)
      : m_source(std::move(source)),
        m_subspace(indices),
        m_params(params),
        m_cov(cov) {
    // TODO we should be able to support arbitrary ordering, by sorting the
    //   indices and reordering parameters/covariance. since the parameter order
    //   can be modified by the user, the user can not always know what the
    //   right order is. another option is too remove the requirement for index
    //   ordering from the subspace types, but that will make it harder to
    //   refactor their implementation later on.
  }
  /// A measurement can only be constructed with valid parameters.
  Measurement() = delete;
  Measurement(const Measurement&) = default;
  Measurement(Measurement&&) = default;
  ~Measurement() = default;
  Measurement& operator=(const Measurement&) = default;
  Measurement& operator=(Measurement&&) = default;

  /// Source link that connects to the underlying detector readout.
  const SourceLink& sourceLink() const { return m_source; }

  /// Number of measured parameters.
  static constexpr std::size_t size() { return kSize; }

  /// Check if a specific parameter is part of this measurement.
  bool contains(indices_t i) const { return m_subspace.contains(i); }

  /// The measurement indices
  constexpr std::array<indices_t, kSize> indices() const {
    std::array<uint8_t, kSize> subInds = m_subspace.indices();
    std::array<indices_t, kSize> inds{};
    for (std::size_t i = 0; i < kSize; i++) {
      inds[i] = static_cast<indices_t>(subInds[i]);
    }
    return inds;
  }

  /// Measured parameters values.
  const ParametersVector& parameters() const { return m_params; }

  /// Measured parameters covariance.
  const CovarianceMatrix& covariance() const { return m_cov; }

  /// Projection matrix from the full space into the measured subspace.
  ProjectionMatrix projector() const {
    return m_subspace.template projector<Scalar>();
  }

  /// Expansion matrix from the measured subspace into the full space.
  ///
  /// This is equivalent to the transpose of the projection matrix only in the
  /// case of a trivial projection matrix. While this is the case here, it is
  /// still recommended to use the expansion matrix directly in cases where it
  /// is explicitly used.
  ExpansionMatrix expander() const {
    return m_subspace.template expander<Scalar>();
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
    detail::calculateResiduals(static_cast<indices_t>(kSize),
                               m_subspace.indices(), reference, m_params, res);
    return res;
  }

  std::ostream& operator<<(std::ostream& os) const {
    detail::printMeasurement(os, static_cast<indices_t>(kSize),
                             m_subspace.indices().data(), m_params.data(),
                             m_cov.data());
    return os;
  }

 private:
  SourceLink m_source;
  Subspace m_subspace;
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
auto makeMeasurement(SourceLink source,
                     const Eigen::MatrixBase<parameters_t>& params,
                     const Eigen::MatrixBase<covariance_t>& cov,
                     indices_t index0, tail_indices_t... tailIndices)
    -> Measurement<indices_t, 1u + sizeof...(tail_indices_t)> {
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
  using Type = std::variant<Measurement<indices_t, kSizes>...>;
};

/// @endcond
}  // namespace detail

/// Variant that can contain all possible measurements in a parameter space.
///
/// @tparam indices_t Parameter index type, determines the full parameter space
template <typename indices_t>
using VariantMeasurement = typename detail::VariantMeasurementGenerator<
    indices_t, detail::kParametersSize<indices_t>>::Type;

/// Variant that can hold all possible bound measurements.
///
using BoundVariantMeasurement = VariantMeasurement<BoundIndices>;

/// Variant that can hold all possible free measurements.
///
using FreeVariantMeasurement = VariantMeasurement<FreeIndices>;

template <typename indices_t>
std::ostream& operator<<(std::ostream& os,
                         const VariantMeasurement<indices_t>& vm) {
  return std::visit([&](const auto& m) { return (os << m); }, vm);
}

}  // namespace Acts
