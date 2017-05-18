// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// HomogeneousSurfaceMaterial.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"

Acts::HomogeneousSurfaceMaterial::HomogeneousSurfaceMaterial(
    const MaterialProperties& full,
    double                    splitFactor)
  : SurfaceMaterial(splitFactor), m_fullMaterial(full)
{
}

Acts::HomogeneousSurfaceMaterial::HomogeneousSurfaceMaterial(
    const Acts::HomogeneousSurfaceMaterial& lmp)
  : SurfaceMaterial(lmp.m_splitFactor), m_fullMaterial(lmp.m_fullMaterial)
{
}

Acts::HomogeneousSurfaceMaterial::~HomogeneousSurfaceMaterial()
{
}

Acts::HomogeneousSurfaceMaterial&
Acts::HomogeneousSurfaceMaterial::
operator=(const HomogeneousSurfaceMaterial& lmp)
{
  if (this != &lmp) {
    // now refill evertything
    m_fullMaterial                       = lmp.m_fullMaterial;
    Acts::SurfaceMaterial::m_splitFactor = lmp.m_splitFactor;
  }
  return (*this);
}

Acts::HomogeneousSurfaceMaterial&
Acts::HomogeneousSurfaceMaterial::operator*=(double scale)
{
  // scale the sub properties
  m_fullMaterial *= scale;
  return (*this);
}

std::ostream&
Acts::HomogeneousSurfaceMaterial::dump(std::ostream& sl) const
{
  sl << "Acts::HomogeneousSurfaceMaterial : " << std::endl;
  sl << "   - fullMaterial : " << m_fullMaterial << std::endl;
  sl << "   - split factor : " << m_splitFactor << std::endl;
  return sl;
}
