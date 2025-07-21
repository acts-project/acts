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

namespace Acts::detail {
/// @brief Helper to describe a line in 3D space. The line is described
///        by a point on the line and a corresponding direction vector with unit
///        length. Additionally, the `Line3DWithPartialDerivatives` holds the
///        first and second derivative of the line with respect to the
///        parameters `x0`, `y0`, `theta`, and `phi` (the parameters are defined
///        in the enum `ParIndices`).
template <std::floating_point T>
class Line3DWithPartialDerivatives {
 public:
  /// @brief Abrivation of the Vector
  using Vector = Eigen::Matrix<T, 3, 1>;
  /// @brief Enum to map the indices of the parameter vector
  enum ParIndices : std::size_t {
    x0 = 0,
    y0 = 1,
    theta = 2,
    phi = 3,
    nPars = 4
  };
  /// @brief Abrivation of the parameter vector type
  using ParamVector = std::array<T, nPars>;

  /// @brief Update the line & derivatives with the new parameters
  /// @param newPars The new parameters to update the line with
  void updateParameters(const ParamVector& newPars);
  /// @brief Returns a point on the line
  const Vector& position() const;
  /// @brief Returns the direction of the line
  const Vector& direction() const;
  /// @brief Returns the first derivative of the line with respect to the passed
  /// @param param: Index of the parameter to get the derivative for
  const Vector& gradient(const std::size_t param) const;
  /// @brief Returns the second derivative of the line with respect to the passed
  /// @param param1: Index of the first parameter to get the derivative for
  /// @param param2: Index of the second parameter to get the derivative for
  const Vector& hessian(const std::size_t param1,
                        const std::size_t param2) const;
  /// @brief Returns a point along the line which is by lambda apart from
  ///        the line reference
  /// @param lambda: Separation from the reference
  Vector point(const double lambda) const;

 private:
  Vector m_pos{Vector::Zero()};
  Vector m_dir{Vector::Zero()};
  static constexpr std::size_t s_nPars = ParIndices::nPars;
  std::array<Vector, s_nPars> m_gradient{
      filledArray<Vector, s_nPars>(Vector::Zero())};

  std::array<Vector, sumUpToN(s_nPars)> m_hessian{
      filledArray<Vector, sumUpToN(s_nPars)>(Vector::Zero())};
};
}  // namespace Acts::detail
#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.ipp"
