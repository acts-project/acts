// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

Acts::SolenoidBField::SolenoidBField(Config config) : m_cfg(std::move(config))
{
  m_dz = m_cfg.L / m_cfg.nCoils;
  m_R2 = m_cfg.R * m_cfg.R;
  // we need to scale so we reproduce the expected B field strength
  // at the center of the solenoid
  Vector2D field = multiCoilField({0, 0}, 1.);  // scale = 1
  m_scale        = m_cfg.bMagCenter / field.norm();
}

Acts::Vector3D
Acts::SolenoidBField::getField(const Vector3D& position) const
{
  using VectorHelpers::perp;
  Vector2D rzPos(perp(position), position.z());
  Vector2D rzField = multiCoilField(rzPos, m_scale);
  Vector3D xyzField(0, 0, rzField[1]);

  if (rzPos[0] != 0.) {
    // add xy field component, radially symmetric
    Vector3D rDir = Vector3D(position.x(), position.y(), 0).normalized();
    xyzField += rDir * rzField[0];
  }

  return xyzField;
}

Acts::Vector3D
Acts::SolenoidBField::getField(const Vector3D& position, Cache& /*cache*/) const
{
  return getField(position);
}

Acts::Vector2D
Acts::SolenoidBField::getField(const Vector2D& position) const
{
  return multiCoilField(position, m_scale);
}

Acts::Vector3D
Acts::SolenoidBField::getFieldGradient(const Vector3D& position,
                                       ActsMatrixD<3, 3>& /*derivative*/) const
{
  return getField(position);
}

Acts::Vector2D
Acts::SolenoidBField::multiCoilField(const Vector2D& pos, double scale) const
{
  // iterate over all coils
  Vector2D resultField(0, 0);
  for (size_t coil = 0; coil < m_cfg.nCoils; coil++) {
    Vector2D shiftedPos
        = Vector2D(pos[0], pos[1] + m_cfg.L * 0.5 - m_dz * (coil + 0.5));
    resultField += singleCoilField(shiftedPos, scale);
  }

  return resultField;
}

Acts::Vector2D
Acts::SolenoidBField::singleCoilField(const Vector2D& pos, double scale) const
{
  return {B_r(pos, scale), B_z(pos, scale)};
}

double
Acts::SolenoidBField::B_r(const Vector2D& pos, double scale) const
{
  using boost::math::ellint_1;
  using boost::math::ellint_2;

  double r = std::abs(pos[0]);
  double z = pos[1];

  if (r == 0) return 0.;

  double k_2      = k2(r, z);
  double k        = std::sqrt(k_2);
  double constant = scale * k * z / (4 * M_PI * std::sqrt(m_cfg.R * r * r * r));

  double B = (2. - k_2) / (2. - 2. * k_2) * ellint_2(k_2) - ellint_1(k_2);

  // pos[0] is still signed!
  return r / pos[0] * constant * B;
}

double
Acts::SolenoidBField::B_z(const Vector2D& pos, double scale) const
{
  using boost::math::ellint_1;
  using boost::math::ellint_2;

  double r = std::abs(pos[0]);
  double z = pos[1];

  if (r == 0) {
    double res = scale / 2. * m_R2 / (std::sqrt(m_R2 + z * z) * (m_R2 + z * z));
    return res;
  }

  double k_2      = k2(r, z);
  double k        = std::sqrt(k_2);
  double constant = scale * k / (4 * M_PI * std::sqrt(m_cfg.R * r));
  double B
      = ((m_cfg.R + r) * k_2 - 2. * r) / (2. * r * (1. - k_2)) * ellint_2(k_2)
      + ellint_1(k_2);

  return constant * B;
}

double
Acts::SolenoidBField::k2(double r, double z) const
{
  return 4 * m_cfg.R * r / ((m_cfg.R + r) * (m_cfg.R + r) + z * z);
}
