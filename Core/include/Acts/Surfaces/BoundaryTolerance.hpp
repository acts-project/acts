// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <type_traits>
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
/// - Chi2Cartesian: Chi2 tolerance in Cartesian coordinates.
///   Similar to Chi2Bound, but the chi2 value is calculated in Cartesian
///   coordinates.
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

  /// Parameters for absolute Euclidean boundary tolerance.
  struct AbsoluteEuclideanParams {
    /// Tolerance value
    double tolerance{};
  };

  /// Parameters for chi2 boundary tolerance in bound coordinates.
  struct Chi2BoundParams {
    /// Maximum chi2 value
    double maxChi2{};
    /// Weight matrix stored as flat array
    std::array<double, 4> weight{};

    /// Get weight matrix as Eigen matrix
    /// @return Mapped weight matrix
    Eigen::Map<SquareMatrix2> weightMatrix() {
      return Eigen::Map<SquareMatrix2>(weight.data());
    }

    /// Get weight matrix as const Eigen matrix
    /// @return Mapped weight matrix
    Eigen::Map<const SquareMatrix2> weightMatrix() const {
      return Eigen::Map<const SquareMatrix2>(weight.data());
    }
  };

  /// Parameters for chi2 boundary tolerance in Cartesian coordinates.
  struct Chi2CartesianParams {
    /// Maximum chi2 value
    double maxChi2{};
    /// Weight matrix stored as flat array
    std::array<double, 4> weight{};

    /// Get weight matrix as Eigen matrix
    /// @return Mapped weight matrix
    Eigen::Map<SquareMatrix2> weightMatrix() {
      return Eigen::Map<SquareMatrix2>(weight.data());
    }

    /// Get weight matrix as const Eigen matrix
    /// @return Mapped weight matrix
    Eigen::Map<const SquareMatrix2> weightMatrix() const {
      return Eigen::Map<const SquareMatrix2>(weight.data());
    }
  };

  static_assert(std::is_trivially_copyable_v<Chi2BoundParams>);

 private:
  /// Underlying variant type
  using Variant =
      std::variant<InfiniteParams, NoneParams, AbsoluteEuclideanParams,
                   Chi2BoundParams, Chi2CartesianParams>;
  static_assert(std::is_trivially_copyable_v<Variant>);

  /// Construct from variant
  constexpr explicit BoundaryTolerance(Variant variant) : m_variant{variant} {}

 public:
  /// Infinite tolerance i.e. no boundary check
  /// @return BoundaryTolerance configured for infinite tolerance
  constexpr static auto Infinite() noexcept {
    return BoundaryTolerance{InfiniteParams{}};
  }

  /// No tolerance i.e. exact boundary check
  /// @return BoundaryTolerance configured for no tolerance (exact checking)
  constexpr static auto None() noexcept {
    return BoundaryTolerance{NoneParams{}};
  }

  /// Absolute tolerance in Euclidean distance
  /// @param tolerance The tolerance value in Euclidean distance
  /// @return BoundaryTolerance configured for absolute Euclidean tolerance
  constexpr static auto AbsoluteEuclidean(double tolerance) noexcept {
    return BoundaryTolerance{AbsoluteEuclideanParams{tolerance}};
  }

  /// Chi2 tolerance in bound coordinates
  /// @param weight The weight matrix for the chi2 calculation
  /// @param maxChi2 The maximum chi2 value allowed
  /// @return BoundaryTolerance configured for chi2 tolerance in bound coordinates
  static auto Chi2Bound(const SquareMatrix2& weight, double maxChi2) noexcept {
    Chi2BoundParams tolerance{maxChi2, {}};
    tolerance.weightMatrix() = weight;
    return BoundaryTolerance{tolerance};
  }

  /// Chi2 tolerance in Cartesian coordinates
  /// @param weight The weight matrix for the chi2 calculation
  /// @param maxChi2 The maximum chi2 value allowed
  /// @return BoundaryTolerance configured for chi2 tolerance in Cartesian coordinates
  static auto Chi2Cartesian(const SquareMatrix2& weight,
                            double maxChi2) noexcept {
    Chi2CartesianParams tolerance{maxChi2, {}};
    tolerance.weightMatrix() = weight;
    return BoundaryTolerance{tolerance};
  }

  /// Copy constructor
  /// @param other The BoundaryTolerance object to copy
  BoundaryTolerance(const BoundaryTolerance& other) noexcept = default;
  /// Copy assignment operator
  /// @param other The BoundaryTolerance object to copy
  /// @return Reference to this BoundaryTolerance after copying
  BoundaryTolerance& operator=(const BoundaryTolerance& other) noexcept =
      default;
  /// Move constructor
  /// @param other The BoundaryTolerance object to move
  BoundaryTolerance(BoundaryTolerance&& other) noexcept = default;
  /// Move assignment operator
  /// @param other The BoundaryTolerance object to move
  /// @return Reference to this BoundaryTolerance after moving
  BoundaryTolerance& operator=(BoundaryTolerance&& other) noexcept = default;

  /// @enum ToleranceMode
  /// Enumeration defining how tolerance should be applied to boundaries
  enum class ToleranceMode {
    Extend,  //< Extend the boundary
    None,    //< No tolerance
    Shrink   //< Shrink the boundary
  };

  /// Check if the tolerance is infinite.
  /// @return True if configured for infinite tolerance, false otherwise
  constexpr bool isInfinite() const { return holdsVariant<InfiniteParams>(); }
  /// Check if the is no tolerance.
  /// @return True if configured for no tolerance (exact checking), false otherwise
  constexpr bool isNone() const { return holdsVariant<NoneParams>(); }
  /// Check if the tolerance is absolute with Euclidean distance.
  /// @return True if configured for absolute Euclidean tolerance, false otherwise
  constexpr bool hasAbsoluteEuclidean() const {
    return holdsVariant<AbsoluteEuclideanParams>();
  }
  /// Check if the tolerance is chi2 with bound coordinates.
  /// @return True if configured for chi2 tolerance in bound coordinates, false otherwise
  constexpr bool hasChi2Bound() const {
    return holdsVariant<Chi2BoundParams>();
  }
  /// Check if the tolerance is chi2 with Cartesian coordinates.
  /// @return True if configured for chi2 tolerance in Cartesian coordinates, false otherwise
  constexpr bool hasChi2Cartesian() const {
    return holdsVariant<Chi2CartesianParams>();
  }

  /// Get the tolerance mode.
  /// @return ToleranceMode indicating how the tolerance should be applied
  ToleranceMode toleranceMode() const;

  /// Get the tolerance as absolute Euclidean.
  /// @return Reference to the absolute Euclidean tolerance parameters
  constexpr const AbsoluteEuclideanParams& asAbsoluteEuclidean() const {
    return getVariant<AbsoluteEuclideanParams>();
  }
  /// Get the tolerance as chi2 bound.
  /// @return Reference to the chi2 bound tolerance parameters
  constexpr const Chi2BoundParams& asChi2Bound() const {
    return getVariant<Chi2BoundParams>();
  }
  /// Get the tolerance as chi2 Cartesian.
  /// @return Reference to the chi2 Cartesian tolerance parameters
  constexpr const Chi2CartesianParams& asChi2Cartesian() const {
    return getVariant<Chi2CartesianParams>();
  }

  /// Check if the bound position delta is tolerated.
  /// @param boundDelta The delta in bound coordinates
  /// @param boundToCartesian The transformation matrix from bound to Cartesian
  /// @return True if the delta is within tolerance, false otherwise
  bool isTolerated(const Vector2& boundDelta,
                   const SquareMatrix2& boundToCartesian) const;

 private:
  Variant m_variant;

  /// Check if the boundary check is of a specific type.
  template <typename T>
  constexpr bool holdsVariant() const {
    return std::holds_alternative<T>(m_variant);
  }

  /// Get the specific underlying type.
  template <typename T>
  constexpr const T& getVariant() const {
    return std::get<T>(m_variant);
  }

  template <typename T>
  constexpr const T* getVariantPtr() const {
    return holdsVariant<T>() ? &getVariant<T>() : nullptr;
  }
};

static_assert(std::is_trivially_copyable_v<BoundaryTolerance>);
static_assert(std::is_trivially_move_constructible_v<BoundaryTolerance>);

}  // namespace Acts
