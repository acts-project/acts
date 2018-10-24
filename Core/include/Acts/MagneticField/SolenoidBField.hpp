// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @ingroup MagneticField
//
/// @class SolenoidBField
/// Implements a multi-coil solenoid magnetic field. On every call, the field
/// is evaluated at that exact position. The field has radially symmetry, the
/// field vectors point in +z direction.
/// The config exposes a target field value in the center. This value is used
/// to empirically determine a scale factor which reproduces this field value
/// in the center.
class SolenoidBField
{
public:
  struct Cache
  {
    // empty, we don't need one
  };

  struct Config
  {
    double R;
    double L;
    size_t nCoils;
    double bMagCenter;
  };

  /// @brief the constructur with a shared pointer
  /// @note since it is a shared field, we enforce it to be const
  /// @tparam bField is the shared BField to be stored
  SolenoidBField(Config config);

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  ///
  /// @return magnetic field vector at given position
  Vector3D
  getField(const Vector3D& position) const;

  /// @brief Retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  /// @param [in] cache Cache object, passed through to wrapped BField
  Vector3D
  getField(const Vector3D& position, Cache& /*cache*/) const;

  /// @brief Retrieve magnetic field value in local (r,z) coordinates
  ///
  /// @param [in] position local 2D position
  Vector2D
  getField(const Vector2D& position) const;

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @return magnetic field vector
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D
  getFieldGradient(const Vector3D& position,
                   ActsMatrixD<3, 3>& /*derivative*/) const;

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @param [in] cache Cache object, passed through to wrapped BField
  /// @return magnetic field vector
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D
  getFieldGradient(const Vector3D& position,
                   ActsMatrixD<3, 3>& /*derivative*/,
                   Cache& /*cache*/) const;

private:
  Config m_cfg;
  double m_scale;
  double m_dz;
  double m_R2;

  Vector2D
  multiCoilField(const Vector2D& pos, double scale) const;

  Vector2D
  singleCoilField(const Vector2D& pos, double scale) const;

  double
  B_r(const Vector2D& pos, double scale) const;

  double
  B_z(const Vector2D& pos, double scale) const;

  double
  k2(double r, double z) const;
};

}  // namespace Acts
