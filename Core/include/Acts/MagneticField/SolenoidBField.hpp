// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldContext.hpp"
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
///
/// E_1(k^2) = complete elliptic integral of the 1st kind
/// E_2(k^2) = complete elliptic integral of the 2nd kind
///
/// E_1(k^2) and E_2(k^2) are usually indicated as K(k^2) and E(k^2) in
/// literature,
/// respectively
///              _
///     2       /  pi / 2          2    2          - 1 / 2
/// E (k )  =   |         ( 1  -  k  sin {theta} )         dtheta
///  1         _/  0
///
///              _          ____________________
///     2       /  pi / 2| /       2    2
/// E (k )  =   |        |/ 1  -  k  sin {theta} dtheta
///  2         _/  0
///
/// k^2 = is a function of the point (r, z) and of the radius of the coil R
///
///  2           4Rr
/// k   =  ---------------
///               2      2
///        (R + r)   +  z
/// Using these, you can evaluate the two components B_r and B_z of the
/// magnetic field:
///                            _                             _
///              mu  I        |  /     2 \                    |
///                0     kz   |  |2 - k  |    2          2    |
/// B (r, z)  =  ----- ------ |  |-------|E (k )  -  E (k )   |
///  r            4pi     ___ |  |      2| 2          1       |
///                    | /  3 |_ \2 - 2k /                   _|
///                    |/ Rr
///
///                         _                                       _
///             mu  I      |  /         2      \                     |
///               0     k  |  | (R + r)k  - 2r |     2          2    |
/// B (r,z)  =  ----- ---- |  | -------------- | E (k )  +  E (k )   |
///  z           4pi    __ |  |           2    |  2          1       |
///                   |/Rr |_ \   2r(1 - k )   /                    _|
///
class SolenoidBField {
 public:
  struct Cache {
    /// @brief Constructor with magnetic field context
    ///
    /// @param mcfg the magnetic field context
    Cache(std::reference_wrapper<const MagneticFieldContext> /*mcfg*/) {}
  };

  /// Config struct for the SolenoidBfield.
  struct Config {
    /// Radius at which the coils are located.
    double radius;
    /// Extent of the solenoid in z. It goes from
    /// -length/2 to +length/2 by convention
    double length;
    /// The number of coils that make up the solenoid
    size_t nCoils;
    /// The target magnetic field strength at the center.
    /// This will be used to scale coefficients
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
  Vector3D getField(const Vector3D& position) const;

  /// @brief Retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  /// @param [in] cache Cache object, passed through to wrapped BField
  Vector3D getField(const Vector3D& position, Cache& /*cache*/) const;

  /// @brief Retrieve magnetic field value in local (r,z) coordinates
  ///
  /// @param [in] position local 2D position
  Vector2D getField(const Vector2D& position) const;

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @return magnetic field vector
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D getFieldGradient(const Vector3D& position,
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
  Vector3D getFieldGradient(const Vector3D& position,
                            ActsMatrixD<3, 3>& /*derivative*/,
                            Cache& /*cache*/) const;

 private:
  Config m_cfg;
  double m_scale;
  double m_dz;
  double m_R2;

  Vector2D multiCoilField(const Vector2D& pos, double scale) const;

  Vector2D singleCoilField(const Vector2D& pos, double scale) const;

  double B_r(const Vector2D& pos, double scale) const;

  double B_z(const Vector2D& pos, double scale) const;

  double k2(double r, double z) const;
};

}  // namespace Acts
