// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstddef>

namespace Acts {

/// @class SolenoidBField
/// @brief Analytical solenoid magnetic field implementation
///
/// @ingroup magnetic_field
///
/// @section Concept
///
/// This class implements a multi-coil solenoid magnetic field. On every call,
/// the field is evaluated at that exact position. The field has radially
/// symmetry, the field vectors point in +z direction. The config exposes a
/// target field value in the center. This value is used to empirically
/// determine a scale factor which reproduces this field value in the center.
///
/// @image html bfield/quiver.png "Picture of a solenoid field in rz, with arrows indicating the direction of the field, and their size denoting the strength. The field is almost homogeneous in the center." width=100%
///
/// @note
/// A configuration of
/// ```cpp
/// SolenoidBField::Config cfg;
/// cfg.length = 5.8_m;
/// cfg.radius = (2.56 + 2.46) * 0.5 * 0.5_m;
/// cfg.nCoils = 1154;
/// cfg.bMagCenter = 2_T;
/// SolenoidBField bField(cfg);
/// ```
/// roughly corresponds to the solenoid wrapping the Inner Detector in ATLAS!
///
/// @warning
/// As the evaluation of @f$E_1(k^2)@f$ and @f$E_2(k^2)@f$ is **slow**. The
/// @ref Acts::InterpolatedBFieldMap easily outperforms
/// @ref Acts::SolenoidBField. A helper is provided that builds a map from the
/// analytical implementation and is much faster to lookup:
/// @ref Acts::solenoidFieldMap.
///
/// @section Implementation
///
/// - @f$E_1(k^2)@f$ = complete elliptic integral of the 1st kind
/// - @f$E_2(k^2)@f$ = complete elliptic integral of the 2nd kind
///
/// @f$E_1(k^2)@f$ and @f$E_2(k^2)@f$ are usually indicated as @f$K(k^2)@f$ and @f$E(k^2)@f$ in
/// literature,
/// respectively
///
/// @f[
/// E_1(k^2) = \int_0^{\pi/2} \left( 1 - k^2 \sin^2{\theta} \right )^{-1/2} \mathrm{d}\theta
/// @f]
///
/// @f[
/// E_2(k^2) = \int_0^{\pi/2}\sqrt{1 - k^2 \sin^2{\theta}} \mathrm{d}\theta
/// @f]
///
/// @f$k^2@f$ is a function of the point @f$(r, z)@f$ and of the radius of the coil @f$R@f$
///
/// @f[
/// k^2 = \frac{4Rr}{(R+r)^2 + z^2}
/// @f]
///
/// Using these, you can evaluate the two components @f$B_r@f$ and @f$B_z@f$ of the
/// magnetic field:
///
/// @f[
/// B_r(r, z) = \frac{\mu_0 I}{4\pi} \frac{kz}{\sqrt{Rr^3}} \left[ \left(\frac{2-k^2}{2-2k^2}\right)E_2(k^2) - E_1(k^2) \right ]
/// @f]
///
/// @f[
/// B_z(r,z) = \frac{\mu_0 I}{4\pi} \frac{k}{\sqrt{Rr}} \left[ \left( \frac{(R+r)k^2-2r}{2r(1-k^2)} \right ) E_2(k^2) + E_1(k^2) \right ]
/// @f]
///
///
/// In the implementation the factor of @f$(\mu_0\cdot I)@f$ is defined to be a scaling
/// factor. It is evaluated and defined as the magnetic field in the center of
/// the
/// coil, i.e. the scale set in @ref Acts::SolenoidBField::Config::bMagCenter.
///
class SolenoidBField final : public MagneticFieldProvider {
 public:
  struct Cache {
    /// @brief Constructor with magnetic field context
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  /// Config struct for the SolenoidBfield.
  struct Config {
    /// Radius at which the coils are located.
    double radius;
    /// Extent of the solenoid in z. It goes from
    /// -length/2 to +length/2 by convention
    double length;
    /// The number of coils that make up the solenoid
    std::size_t nCoils;
    /// The target magnetic field strength at the center.
    /// This will be used to scale coefficients
    double bMagCenter;
  };

  /// @brief the constructor with a shared pointer
  /// @note since it is a shared field, we enforce it to be const
  /// @tparam bField is the shared BField to be stored
  /// @param config Configuration struct containing solenoid parameters
  explicit SolenoidBField(Config config);

  /// @brief Retrieve magnetic field value in local (r,z) coordinates
  ///
  /// @param [in] position local 2D position
  /// @return Magnetic field vector in local (r,z) coordinates
  Vector2 getField(const Vector2& position) const;

  /// @copydoc MagneticFieldProvider::makeCache(const MagneticFieldContext&) const
  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override;

  /// @brief Get the B field at a position
  ///
  /// @param position The position to query at
  /// @return Magnetic field vector in global coordinates
  Vector3 getField(const Vector3& position) const;

  /// @copydoc MagneticFieldProvider::getField(const Vector3&,MagneticFieldProvider::Cache&) const
  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override;

 private:
  Config m_cfg;
  double m_scale;
  double m_dz;
  double m_R2;

  Vector2 multiCoilField(const Vector2& pos, double scale) const;

  Vector2 singleCoilField(const Vector2& pos, double scale) const;

  double B_r(const Vector2& pos, double scale) const;

  double B_z(const Vector2& pos, double scale) const;

  double k2(double r, double z) const;
};

}  // namespace Acts
