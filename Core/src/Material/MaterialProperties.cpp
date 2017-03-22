Core/include/ACTS/Material/MaterialProperties.hpp// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialProperties.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Material/MaterialProperties.hpp"
#include <climits>

Acts::MaterialProperties::MaterialProperties()
  : m_material()
  , m_dInX0(0.)
  , m_dInL0(0.)
{
}

Acts::MaterialProperties::MaterialProperties(float Xo,
                                             float Lo,
                                             float averageA,
                                             float averageZ,
                                             float averageRho,
                                             float thickness)
  : m_material(Xo, Lo, averageA, averageZ, averageRho)
  , m_dInX0(Xo * Xo > 10e-10 ? thickness / Xo : 0.)
  , m_dInL0(Lo * Lo > 10e-10 ? thickness / Lo : 0.)
{
}

Acts::MaterialProperties::MaterialProperties(const Material& material,
                                             float thickness)
  : m_material(material)
  , m_dInX0(material.X0() * material.X0() > 10e-10 ? thickness / material.X0() : 0.)
  , m_dInL0(material.L0() * material.L0() > 10e-10 ? thickness / material.L0() : 0.)
{
}

Acts::MaterialProperties::MaterialProperties(
    const Acts::MaterialProperties& mprop)
  : m_material(mprop.m_material)
  , m_dInX0(mprop.m_dInX0)
  , m_dInL0(mprop.m_dInL0)
{
}

Acts::MaterialProperties*
Acts::MaterialProperties::clone() const
{
  return new Acts::MaterialProperties(*this);
}

Acts::MaterialProperties&
Acts::MaterialProperties::operator=(const Acts::MaterialProperties& mprop)
{
  if (this != &mprop) {
    m_material = mprop.m_material;
    m_dInX0    = mprop.m_dInX0;
    m_dInL0    = mprop.m_dInL0;
  }
  return (*this);
}

Acts::MaterialProperties&
Acts::MaterialProperties::operator*=(float scale)
{
  // assuming rescaling of the material thickness
  m_dInX0 *= scale;
  m_dInL0 *= scale;
  return (*this);
}

std::ostream&
Acts::operator<<(std::ostream& sl, const MaterialProperties& mprop)
{
  sl << "Acts::MaterialProperties: " << std::endl;
  sl << "   - thickness/X0                          = " << mprop.thicknessInX0()
     << std::endl;
  sl << "   - thickness                       [mm]  = " << mprop.thickness()
     << std::endl;
  sl << "   - radiation length X0             [mm]  = " << mprop.averageX0()
     << std::endl;
  sl << "   - nuclear interaction length L0   [mm]  = " << mprop.averageL0()
     << std::endl;
  sl << "   - average material Z/A*rho [gram/mm^3]  = "
     << mprop.zOverAtimesRho() << std::endl;
  /*  interface not finalized
  if (mprop.material().composition){
      sl << "   - material composition from " <<
  mprop.material().composition->size() << " elements " << std::endl;
      sl << "       listing them (prob. ordereded ) : " << std::endl;
      for ( auto& eIter : (*mprop.material().composition) )
          sl << "         -> Z : " << eIter.element() << "( fraction : "  <<
  eIter.fraction() << " )" << std::endl;
  }
  */
  return sl;
}
