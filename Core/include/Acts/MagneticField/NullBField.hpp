// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace Acts {

/// @ingroup MagneticField
/// @brief Null bfield which returns 0 always
class NullBField final {
 public:
  struct Cache {
    /// @brief constructor with context
    Cache(std::reference_wrapper<const MagneticFieldContext> /*mcfg*/) {}
  };

  /// @brief Default constructor
  NullBField() = default;

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global position
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Vector3D getField(const Vector3D& /*position*/) const { return m_BField; }

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global position
  /// @param [in] cache Cache object (is ignored)
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  Vector3D getField(const Vector3D& /*position*/, Cache& /*cache*/) const {
    return m_BField;
  }

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D getFieldGradient(const Vector3D& /*position*/,
                            ActsMatrixD<3, 3>& /*derivative*/) const {
    return m_BField;
  }

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @param [in] cache Cache object (is ignored)
  /// @return magnetic field vector
  ///
  /// @note The @p position is ignored and only kept as argument to provide
  ///       a consistent interface with other magnetic field services.
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D getFieldGradient(const Vector3D& /*position*/,
                            ActsMatrixD<3, 3>& /*derivative*/,
                            Cache& /*cache*/) const {
    return m_BField;
  }

  /// @brief check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position global 3D position
  /// @return @c true if position is inside the defined look-up grid,
  ///         otherwise @c false
  /// @note The method will always return true for the null B-Field
  bool isInside(const Vector3D& /*position*/) const { return true; }

 private:
  /// magnetic field vector
  const Vector3D m_BField = Vector3D::Zero();
};
}  // namespace Acts
