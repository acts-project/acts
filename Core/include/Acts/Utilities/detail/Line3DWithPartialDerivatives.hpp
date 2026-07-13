// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"

#include <cstdint>

namespace Acts::detail {
/// Helper to describe a line in 3D space. The line is described by a point on
/// the line and a corresponding direction vector with unit length.
/// Additionally, the `Line3DWithPartialDerivatives` holds the first and second
/// derivative of the line with respect to the parameters `x0`, `y0`, `theta`,
/// and `phi` (the parameters are defined in the enum `ParIndex`).
template <std::floating_point T>
class Line3DWithPartialDerivatives {
 public:
  /// Abrivation of the Vector
  using Vector = Eigen::Matrix<T, 3, 1>;
  /// Enum to map the indices of the parameter vector
  enum class ParIndex : std::uint8_t {
    y0 = 0,     /// y-axis intercept at z=0
    theta = 1,  /// polar angle w.r.t the z-axis
    x0 = 2,     /// x-axis intercept at z=0
    phi = 3,    /// Azimuthal angle in the x-y plane
    nPars = 4,  /// Number of line parameters
  };
  static constexpr auto s_nPars = static_cast<std::size_t>(ParIndex::nPars);
  /// Abrivation of the parameter vector type
  using ParamVector = std::array<T, s_nPars>;
  /// Default constructor
  Line3DWithPartialDerivatives() = default;
  /// Default copy constructor
  Line3DWithPartialDerivatives(
      const Line3DWithPartialDerivatives& other) noexcept = default;
  /// Default move constructor
  Line3DWithPartialDerivatives(Line3DWithPartialDerivatives&& other) noexcept =
      default;
  /// Default assignment operator
  Line3DWithPartialDerivatives& operator=(
      const Line3DWithPartialDerivatives& other) noexcept = default;
  /// Default move assignemt operator
  Line3DWithPartialDerivatives& operator=(
      Line3DWithPartialDerivatives&& other) noexcept = default;
  /// Constructor taking a parameter array which is at least as big as the
  /// ParamVector. The first 4 indices are interpreted as the corresponding line
  /// parameters
  template <std::size_t N>
  explicit Line3DWithPartialDerivatives(
      const std::array<T, N>& initPars) noexcept;

  /// Update the line & derivatives with the new parameters
  /// @param newPars The new parameters to update the line with
  template <std::size_t N>
  void updateParameters(const std::array<T, N>& newPars) noexcept
    requires(N >= s_nPars);
  /// Returns a point on the line
  const Vector& position() const;
  /// Returns the direction of the line
  const Vector& direction() const;
  /// Returns the first derivative of the line with respect to the passed
  /// @param param: Index of the parameter to get the derivative for
  const Vector& gradient(const ParIndex param) const;
  /// Returns the second derivative of the line with respect to the passed
  /// @param param1: Index of the first parameter to get the derivative for
  /// @param param2: Index of the second parameter to get the derivative for
  const Vector& hessian(const ParIndex param1, const ParIndex param2) const;
  /// Returns a point along the line which is by lambda apart from the line
  /// reference
  /// @param lambda: Separation from the reference
  Vector point(const double lambda) const;
  /// Returns the currently set parameters
  ParamVector parameters() const;

 private:
  /// Current line position
  Vector m_pos{Vector::Zero()};
  /// Current line direction
  Vector m_dir{Vector::UnitZ()};
  /// Storage of the current partial derivatives w.r.t. line parameters
  std::array<Vector, s_nPars> m_gradient{
      filledArray<Vector, s_nPars>(Vector::Zero())};
  /// Storage of the  current secodn order partial derivatives w.r.t. line
  /// parameters
  std::array<Vector, sumUpToN(s_nPars)> m_hessian{
      filledArray<Vector, sumUpToN(s_nPars)>(Vector::Zero())};
};
}  // namespace Acts::detail
#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.ipp"
