// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <optional>
#include <variant>

namespace Acts {

/// @brief Variant-like type to capture different types of boundary tolerances
///
/// Since our track hypothesis comes with uncertainties, we sometimes need to
/// check if the track is not just within the boundary of the surface but also
/// within a certain tolerance. This class captures different parameterizations
/// of such tolerances. The surface class will then use these tolerances to
/// check if a ray is within the boundary+tolerance of the surface.
///
/// Different types of boundary tolerances implemented:
/// - Infinite: Infinite tolerance i.e. no boundary check will be performed.
/// - None: No tolerance i.e. exact boundary check will be performed.
/// - AbsoluteBound: Absolute tolerance in bound coordinates.
///   The tolerance is defined as a pair of absolute values for the bound
///   coordinates. Only if both coordinates are within the tolerance, the
///   boundary check is considered as passed.
/// - AbsoluteCartesian: Absolute tolerance in Cartesian coordinates.
///   The tolerance is defined as a pair of absolute values for the Cartesian
///   coordinates. The transformation to Cartesian coordinates can be done via
///   the Jacobian for small distances. Only if both coordinates are within
///   the tolerance, the boundary check is considered as passed.
/// - AbsoluteEuclidean: Absolute tolerance in Euclidean distance.
///   The tolerance is defined as a single absolute value for the Euclidean
///   distance. The Euclidean distance can be calculated via the local bound
///   Jacobian and the bound coordinate residual. If the distance is within
///   the tolerance, the boundary check is considered as passed.
/// - Chi2Bound: Chi2 tolerance in bound coordinates.
///   The tolerance is defined as a maximum chi2 value and a weight matrix,
///   which is the inverse of the bound covariance matrix. The chi2 value is
///   calculated from the bound coordinates residual and the weight matrix.
///   If the chi2 value is below the maximum chi2 value, the boundary check
///   is considered as passed.
///
/// The bound coordinates residual is defined as the difference between the
/// point checked and the closest point on the boundary. The Jacobian is the
/// derivative of the bound coordinates with respect to the Cartesian
/// coordinates.
///
class BoundaryTolerance {
 public:
  struct InfiniteParams {};

  struct NoneParams {};

  struct AbsoluteBoundParams {
    long double tolerance0{};
    long double tolerance1{};

    AbsoluteBoundParams() = default;
    AbsoluteBoundParams(long double tolerance0_, long double tolerance1_)
        : tolerance0(tolerance0_), tolerance1(tolerance1_) {
      if (tolerance0 < 0 || tolerance1 < 0) {
        throw std::invalid_argument(
            "AbsoluteBound: Tolerance must be non-negative");
      }
    }
  };

  struct AbsoluteCartesianParams {
    long double tolerance0{};
    long double tolerance1{};

    AbsoluteCartesianParams() = default;
    AbsoluteCartesianParams(long double tolerance0_, long double tolerance1_)
        : tolerance0(tolerance0_), tolerance1(tolerance1_) {
      if (tolerance0 < 0 || tolerance1 < 0) {
        throw std::invalid_argument(
            "AbsoluteCartesian: Tolerance must be non-negative");
      }
      if ((tolerance0 == 0) != (tolerance1 == 0)) {
        throw std::invalid_argument(
            "AbsoluteCartesian: Both tolerances must be zero or non-zero");
      }
    }
  };

  struct AbsoluteEuclideanParams {
    long double tolerance{};

    AbsoluteEuclideanParams() = default;
    explicit AbsoluteEuclideanParams(long double tolerance_)
        : tolerance(tolerance_) {}
  };

  struct Chi2BoundParams {
    long double maxChi2{};
    SquareMatrix2 weight = SquareMatrix2::Identity();

    Chi2BoundParams() = default;
    Chi2BoundParams(const SquareMatrix2& weight_, long double maxChi2_)
        : maxChi2(maxChi2_), weight(weight_) {}
  };

 private:
  /// Underlying variant type
  using Variant = std::variant<InfiniteParams, NoneParams, AbsoluteBoundParams,
                               AbsoluteCartesianParams, AbsoluteEuclideanParams,
                               Chi2BoundParams>;

  /// Construct from variant
  explicit BoundaryTolerance(Variant variant);

 public:
  /// Infinite tolerance i.e. no boundary check
  static auto Infinite() { return BoundaryTolerance{InfiniteParams{}}; }

  /// No tolerance i.e. exact boundary check
  static auto None() { return BoundaryTolerance{NoneParams{}}; }

  /// Absolute tolerance in bound coordinates
  static auto AbsoluteBound(long double tolerance0, long double tolerance1) {
    return BoundaryTolerance{AbsoluteBoundParams{tolerance0, tolerance1}};
  }

  /// Absolute tolerance in Cartesian coordinates
  static auto AbsoluteCartesian(long double tolerance0,
                                long double tolerance1) {
    return BoundaryTolerance{AbsoluteCartesianParams{tolerance0, tolerance1}};
  }

  /// Absolute tolerance in Euclidean distance
  static auto AbsoluteEuclidean(long double tolerance) {
    return BoundaryTolerance{AbsoluteEuclideanParams{tolerance}};
  }

  /// Chi2 tolerance in bound coordinates
  static auto Chi2Bound(const SquareMatrix2& weight, long double maxChi2) {
    return BoundaryTolerance{Chi2BoundParams{weight, maxChi2}};
  }

  enum class ToleranceMode {
    Extend,  // Extend the boundary
    None,    // No tolerance
    Shrink   // Shrink the boundary
  };

  /// Check if the tolerance is infinite.
  bool isInfinite() const;
  /// Check if the is no tolerance.
  bool isNone() const;
  /// Check if the tolerance is absolute with bound coordinates.
  bool hasAbsoluteBound(bool isCartesian = false) const;
  /// Check if the tolerance is absolute with Cartesian coordinates.
  bool hasAbsoluteCartesian() const;
  /// Check if the tolerance is absolute with Euclidean distance.
  bool hasAbsoluteEuclidean() const;
  /// Check if the tolerance is chi2 with bound coordinates.
  bool hasChi2Bound() const;

  /// Check if any tolerance is set.
  ToleranceMode toleranceMode() const;

  /// Get the tolerance as absolute bound.
  AbsoluteBoundParams asAbsoluteBound(bool isCartesian = false) const;
  /// Get the tolerance as absolute Cartesian.
  const AbsoluteCartesianParams& asAbsoluteCartesian() const;
  /// Get the tolerance as absolute Euclidean.
  const AbsoluteEuclideanParams& asAbsoluteEuclidean() const;
  /// Get the tolerance as chi2 bound.
  const Chi2BoundParams& asChi2Bound() const;

  /// Get the tolerance as absolute bound if possible.
  std::optional<AbsoluteBoundParams> asAbsoluteBoundOpt(
      bool isCartesian = false) const;

  /// Check if the distance is tolerated.
  bool isTolerated(const Vector2& distance,
                   const std::optional<SquareMatrix2>& jacobianOpt) const;

  /// Check if there is a metric assigned with this tolerance.
  bool hasMetric(bool hasJacobian) const;

  /// Get the metric for the tolerance.
  SquareMatrix2 getMetric(const std::optional<SquareMatrix2>& jacobian) const;

 private:
  Variant m_variant;

  /// Check if the boundary check is of a specific type.
  template <typename T>
  bool holdsVariant() const {
    return std::holds_alternative<T>(m_variant);
  }

  /// Get the specific underlying type.
  template <typename T>
  const T& getVariant() const {
    return std::get<T>(m_variant);
  }

  template <typename T>
  const T* getVariantPtr() const {
    return holdsVariant<T>() ? &getVariant<T>() : nullptr;
  }
};

}  // namespace Acts
