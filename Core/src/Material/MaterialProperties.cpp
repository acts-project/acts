// This file is part of the ACTS project.
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
  : m_material(), m_dInX0(0.), m_dInL0(0.)
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
                                             float           thickness)
  : m_material(material)
  , m_dInX0(material.X0() * material.X0() > 10e-10 ? thickness / material.X0()
                                                   : 0.)
  , m_dInL0(material.L0() * material.L0() > 10e-10 ? thickness / material.L0()
                                                   : 0.)
{
}

Acts::MaterialProperties::MaterialProperties(
    const Acts::MaterialProperties& mprop)
  : m_material(mprop.m_material), m_dInX0(mprop.m_dInX0), m_dInL0(mprop.m_dInL0)
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

void
Acts::MaterialProperties::add(const Acts::MaterialProperties& mprop)
{
  if (!(*this)) {
    m_material = mprop.m_material;
    m_dInX0    = mprop.m_dInX0;
    m_dInL0    = mprop.m_dInL0;
  }
  if (mprop && mprop.thickness() != 0) {
    // create a new average material, which scales with the thickness
    // scale the density with the thickness

    float oldScaledRho = this->thickness() * this->averageRho();
    float newScaledRho = mprop.thickness() * mprop.averageRho();
    float updateRho    = oldScaledRho + newScaledRho;
    // scale A and Z with the scaled density
    float updateA
        = oldScaledRho * this->averageA() + newScaledRho * mprop.averageA();
    float updateZ
        = oldScaledRho * this->averageZ() + newScaledRho * mprop.averageZ();
    // Sum of thicknesses
    float sumThickness = this->thickness() + mprop.thickness();
    // divide sumA and sumZ by the sum of scaled rho
    updateA /= updateRho;
    updateZ /= updateRho;
    // divide scaled rho by sum of thicknesses
    updateRho /= sumThickness;
    // dInX0 & and dInL0 are already scaled - just add the new properties
    m_dInX0 += mprop.m_dInX0;
    m_dInL0 += mprop.m_dInL0;
    // calculate X0 and L0 to create new material
    float updateX0 = sumThickness / m_dInX0;
    float updateL0 = sumThickness / m_dInL0;
    // create new material with the updated parameters
    m_material = Material(updateX0, updateL0, updateA, updateZ, updateRho);
  }
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
