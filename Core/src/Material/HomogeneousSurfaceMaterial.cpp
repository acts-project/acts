// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// HomogeneousSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"

Acts::HomogeneousSurfaceMaterial::HomogeneousSurfaceMaterial(
    const MaterialProperties& full,
    double                    splitFactor)
  : SurfaceMaterial(splitFactor), m_fullMaterial(full)
{
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
