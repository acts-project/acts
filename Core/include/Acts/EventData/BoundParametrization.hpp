// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

namespace Acts {

/// Interpretation of the coordinates of a bound state.
///
/// This value owns its interpretation and borrows neither a surface nor a
/// GeometryContext. Construction and copying do not allocate. Vectors and
/// covariance matrices retain their existing Eigen storage. The private storage
/// representation has no stable binary layout contract.
///
/// Currently only canonical global-spherical coordinates are supported. Values
/// can represent the poles, but the inverse direction derivative is singular
/// there. No alternative chart can be selected yet.
class BoundParametrization {
 private:
  struct GlobalSpherical {
    template <typename Derived>
    Vector3 direction(const Eigen::MatrixBase<Derived>& parameters) const {
      const double phi = parameters[eBoundDir0];
      const double theta = parameters[eBoundDir1];
      requireFiniteAngles(phi, theta);
      return makeDirectionFromPhiTheta(phi, theta);
    }

    Vector2 directionCoordinates(const Vector3& direction) const {
      const Vector3 unit = unitDirection(direction);
      return {(unit.x() == 0. && unit.y() == 0.)
                  ? 0.
                  : detail::radian_sym(VectorHelpers::phi(unit)),
              VectorHelpers::theta(unit)};
    }

    template <typename Derived>
    Matrix<3, 2> directionJacobian(
        const Eigen::MatrixBase<Derived>& parameters) const {
      const double phi = parameters[eBoundDir0];
      const double theta = parameters[eBoundDir1];
      requireFiniteAngles(phi, theta);
      const double cp = std::cos(phi), sp = std::sin(phi);
      const double ct = std::cos(theta), st = std::sin(theta);
      Matrix<3, 2> jacobian;
      jacobian << -sp * st, cp * ct, cp * st, sp * ct, 0., -st;
      return jacobian;
    }

    Matrix<2, 3> directionCoordinatesJacobian(
        const Vector3& direction, double minTransverse = 1e-12) const {
      checkThreshold(minTransverse);
      const Vector3 unit = unitDirection(direction);
      const double transverse = std::hypot(unit.x(), unit.y());
      if (transverse <= minTransverse) {
        throw std::domain_error("Singular global-spherical direction Jacobian");
      }
      const double cp = unit.x() / transverse;
      const double sp = unit.y() / transverse;
      Matrix<2, 3> jacobian;
      jacobian << -sp / transverse, cp / transverse, 0., cp * unit.z(),
          sp * unit.z(), -transverse;
      return jacobian / direction.stableNorm();
    }

    template <typename Derived>
    bool isInDomain(const Eigen::MatrixBase<Derived>& parameters,
                    double minTransverse = 1e-12) const {
      checkThreshold(minTransverse);
      return isValid(parameters, true, 0.) && parameters[eBoundDir1] > 0. &&
             parameters[eBoundDir1] < std::numbers::pi &&
             std::sin(parameters[eBoundDir1]) > minTransverse;
    }

    template <typename Derived>
    bool isValid(
        const Eigen::MatrixBase<Derived>& parameters,
        bool validateCoordinateRange = true, double epsilon = 1e-6,
        double maxAbsEta = std::numeric_limits<double>::infinity()) const {
      if (!parameters.allFinite()) {
        return false;
      }
      const double phi = parameters[eBoundDir0];
      const double theta = parameters[eBoundDir1];
      constexpr double pi = std::numbers::pi;
      if (validateCoordinateRange &&
          !(phi + epsilon >= -pi && phi - epsilon < pi &&
            theta + epsilon >= 0. && theta - epsilon <= pi)) {
        return false;
      }
      if (std::isfinite(maxAbsEta)) {
        const double eta = AngleHelpers::etaFromTheta(theta);
        return std::isfinite(eta) && std::abs(eta) - epsilon <= maxAbsEta;
      }
      return true;
    }

    template <typename Derived>
    BoundVector canonicalize(
        const Eigen::MatrixBase<Derived>& parameters) const {
      requireFinite(parameters);
      BoundVector result = parameters;
      std::tie(result[eBoundDir0], result[eBoundDir1]) =
          detail::normalizePhiTheta(result[eBoundDir0], result[eBoundDir1]);
      return result;
    }

    template <typename Derived>
    BoundMatrix canonicalizationJacobian(
        const Eigen::MatrixBase<Derived>& parameters) const {
      requireFinite(parameters);
      const double theta = detail::radian_pos(parameters[eBoundDir1]);
      if (theta == 0. || theta == std::numbers::pi) {
        throw std::domain_error(
            "Canonicalization derivative is undefined at a pole");
      }
      BoundMatrix jacobian = BoundMatrix::Identity();
      if (theta > std::numbers::pi) {
        jacobian(eBoundDir1, eBoundDir1) = -1.;
      }
      return jacobian;
    }

    template <typename Derived, typename Covariance>
    std::pair<BoundVector, BoundMatrix> canonicalize(
        const Eigen::MatrixBase<Derived>& parameters,
        const Eigen::MatrixBase<Covariance>& covariance) const {
      requireFinite(covariance);
      const BoundMatrix jacobian = canonicalizationJacobian(parameters);
      BoundMatrix transformed = covariance;
      if (jacobian(eBoundDir1, eBoundDir1) < 0.) {
        transformed.row(eBoundDir1) *= -1.;
        transformed.col(eBoundDir1) *= -1.;
      }
      return {canonicalize(parameters), transformed};
    }

    template <typename Left, typename Right>
    BoundVector difference(const Eigen::MatrixBase<Left>& lhs,
                           const Eigen::MatrixBase<Right>& rhs) const {
      requireFinite(lhs);
      requireFinite(rhs);
      BoundVector result = lhs - rhs;
      requireFinite(result);
      result[eBoundDir0] = detail::radian_sym(result[eBoundDir0]);
      result[eBoundDir1] = detail::radian_sym(result[eBoundDir1]);
      return result;
    }

    template <typename Parameters, typename Delta>
    BoundVector update(const Eigen::MatrixBase<Parameters>& parameters,
                       const Eigen::MatrixBase<Delta>& delta) const {
      return canonicalize(parameters + delta);
    }

   private:
    static void requireFiniteAngles(double phi, double theta) {
      if (!std::isfinite(phi) || !std::isfinite(theta)) {
        throw std::domain_error("Nonfinite bound direction coordinates");
      }
    }
    template <typename Derived>
    static void requireFinite(const Eigen::MatrixBase<Derived>& values) {
      if (!values.allFinite()) {
        throw std::domain_error("Nonfinite bound coordinates or covariance");
      }
    }
    static void checkThreshold(double threshold) {
      if (!std::isfinite(threshold) || threshold < 0. || threshold >= 1.) {
        throw std::invalid_argument("Invalid direction conditioning threshold");
      }
    }
    static Vector3 unitDirection(const Vector3& direction) {
      const double norm = direction.stableNorm();
      if (!direction.allFinite() || !std::isfinite(norm) || norm == 0.) {
        throw std::domain_error("Direction must be finite and nonzero");
      }
      return direction / norm;
    }
  };
  using Variant = std::variant<GlobalSpherical>;
  Variant m_chart;

  explicit constexpr BoundParametrization(GlobalSpherical chart)
      : m_chart(chart) {}

 public:
  /// Construct the canonical global-spherical interpretation.
  /// @return Allocation-free parametrization value
  static constexpr BoundParametrization globalSpherical() noexcept {
    return BoundParametrization(GlobalSpherical{});
  }

  /// Whether algorithms requiring global spherical coordinates are supported.
  /// @return True for the canonical convention
  constexpr bool isGlobalSpherical() const noexcept {
    return std::holds_alternative<GlobalSpherical>(m_chart);
  }

  /// Decode the global direction without copying the bound vector.
  /// @tparam Derived Eigen vector, map or expression type
  /// @param parameters Bound coordinates
  /// @return Global unit direction, including at the poles
  /// @throws std::domain_error if the direction coordinates are not finite
  template <typename Derived>
  Vector3 direction(const Eigen::MatrixBase<Derived>& parameters) const {
    return std::get<GlobalSpherical>(m_chart).direction(parameters);
  }

  /// Encode a global direction into the two direction coordinates.
  /// @param direction Nonzero global direction; normalization is ignored
  /// @return Canonical direction coordinates; phi is zero at an exact pole
  /// @throws std::domain_error for a zero or nonfinite direction
  Vector2 directionCoordinates(const Vector3& direction) const {
    return std::get<GlobalSpherical>(m_chart).directionCoordinates(direction);
  }

  /// Derivative of global unit direction with respect to direction coordinates.
  /// @tparam Derived Eigen vector, map or expression type
  /// @param parameters Bound coordinates
  /// @return A 3x2 Jacobian (rank deficient at the poles)
  /// @throws std::domain_error for nonfinite direction coordinates
  template <typename Derived>
  Matrix<3, 2> directionJacobian(
      const Eigen::MatrixBase<Derived>& parameters) const {
    return std::get<GlobalSpherical>(m_chart).directionJacobian(parameters);
  }

  /// Derivative of encoding with respect to a global Cartesian direction.
  /// @param direction Nonzero global direction; need not be normalized
  /// @param minTransverse Minimum transverse fraction for acceptable conditioning
  /// @return A 2x3 Jacobian including the derivative of normalization
  /// @throws std::domain_error at a pole or below the conditioning threshold
  /// @throws std::invalid_argument for an invalid threshold
  Matrix<2, 3> directionCoordinatesJacobian(
      const Vector3& direction, double minTransverse = 1e-12) const {
    return std::get<GlobalSpherical>(m_chart).directionCoordinatesJacobian(
        direction, minTransverse);
  }

  /// Whether a bound state admits a sufficiently conditioned inverse
  /// derivative.
  /// @tparam Derived Eigen vector, map or expression type
  /// @param parameters Bound coordinates, required to be canonical and finite
  /// @param minTransverse Minimum transverse fraction
  /// @return True inside the differentiable chart domain
  /// @throws std::invalid_argument for an invalid threshold
  template <typename Derived>
  bool isInDomain(const Eigen::MatrixBase<Derived>& parameters,
                  double minTransverse = 1e-12) const {
    return std::get<GlobalSpherical>(m_chart).isInDomain(parameters,
                                                         minTransverse);
  }

  /// Validate storage values and optional canonical ranges and eta acceptance.
  /// This is separate from derivative conditioning: valid values include poles.
  /// @tparam Derived Eigen vector, map or expression type
  /// @param parameters Bound coordinates
  /// @param validateCoordinateRange Whether to check canonical coordinate ranges
  /// @param epsilon Numerical tolerance for ranges and eta acceptance
  /// @param maxAbsEta Maximum absolute global pseudorapidity
  /// @return Whether all values satisfy the requested checks
  template <typename Derived>
  bool isValid(
      const Eigen::MatrixBase<Derived>& parameters,
      bool validateCoordinateRange = true, double epsilon = 1e-6,
      double maxAbsEta = std::numeric_limits<double>::infinity()) const {
    return std::get<GlobalSpherical>(m_chart).isValid(
        parameters, validateCoordinateRange, epsilon, maxAbsEta);
  }

  /// Canonicalize coordinates without changing the represented direction.
  /// @tparam Derived Eigen vector, map or expression type
  /// @param parameters Bound coordinates
  /// @return Canonical bound coordinates
  /// @throws std::domain_error for nonfinite input
  template <typename Derived>
  BoundVector canonicalize(const Eigen::MatrixBase<Derived>& parameters) const {
    return std::get<GlobalSpherical>(m_chart).canonicalize(parameters);
  }

  /// Local covariance-basis transformation associated with canonicalization.
  /// The theta row changes sign when folding across a pole. The phi wrapping
  /// seam is treated using local periodic coordinates, with derivative +1.
  /// @tparam Derived Eigen vector, map or expression type
  /// @param parameters Original bound coordinates
  /// @return Jacobian from the original to canonical coordinate basis
  /// @throws std::domain_error for nonfinite input or exactly at a pole
  template <typename Derived>
  BoundMatrix canonicalizationJacobian(
      const Eigen::MatrixBase<Derived>& parameters) const {
    return std::get<GlobalSpherical>(m_chart).canonicalizationJacobian(
        parameters);
  }

  /// Canonicalize a state and transform its covariance into the same basis.
  /// @tparam Derived Eigen vector, map or expression type
  /// @tparam Covariance Eigen covariance matrix or map type
  /// @param parameters Original bound coordinates
  /// @param covariance Covariance in the original coordinate basis
  /// @return Canonical coordinates and transformed covariance
  /// @throws std::domain_error for nonfinite input or an undefined derivative
  template <typename Derived, typename Covariance>
  std::pair<BoundVector, BoundMatrix> canonicalize(
      const Eigen::MatrixBase<Derived>& parameters,
      const Eigen::MatrixBase<Covariance>& covariance) const {
    return std::get<GlobalSpherical>(m_chart).canonicalize(parameters,
                                                           covariance);
  }

  /// Coordinate residual lhs - rhs, respecting this parametrization's topology.
  /// @tparam Left Left Eigen vector or map type
  /// @tparam Right Right Eigen vector or map type
  /// @param lhs Left bound coordinates
  /// @param rhs Right bound coordinates
  /// @return Residual suitable for a small filtering/smoothing update
  /// @throws std::domain_error for nonfinite input
  template <typename Left, typename Right>
  BoundVector difference(const Eigen::MatrixBase<Left>& lhs,
                         const Eigen::MatrixBase<Right>& rhs) const {
    return std::get<GlobalSpherical>(m_chart).difference(lhs, rhs);
  }

  /// Apply a small coordinate update and canonicalize the result.
  /// @tparam Parameters Eigen vector or map type
  /// @tparam Delta Eigen update vector or expression type
  /// @param parameters Original coordinates
  /// @param delta Small coordinate increment
  /// @return Canonical updated coordinates
  /// @throws std::domain_error for nonfinite input
  template <typename Parameters, typename Delta>
  BoundVector update(const Eigen::MatrixBase<Parameters>& parameters,
                     const Eigen::MatrixBase<Delta>& delta) const {
    return std::get<GlobalSpherical>(m_chart).update(parameters, delta);
  }
};

static_assert(std::is_trivially_copyable_v<BoundParametrization>);

}  // namespace Acts
