// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"

#include "Acts/Material/MaterialSlab.hpp"

#include <ostream>

namespace Acts {

HomogeneousSurfaceMaterial::HomogeneousSurfaceMaterial(const MaterialSlab& full,
                                                       double splitFactor,
                                                       MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType), m_fullMaterial(full) {}

HomogeneousSurfaceMaterial& HomogeneousSurfaceMaterial::scale(double factor) {
  m_fullMaterial.scaleThickness(factor);
  return *this;
}

const MaterialSlab& HomogeneousSurfaceMaterial::materialSlab(
    const Vector2& /*lp*/) const {
  return m_fullMaterial;
}

const MaterialSlab& HomogeneousSurfaceMaterial::materialSlab(
    const Vector3& /*gp*/) const {
  return m_fullMaterial;
}

std::ostream& HomogeneousSurfaceMaterial::toStream(std::ostream& sl) const {
  sl << "HomogeneousSurfaceMaterial : " << std::endl;
  sl << "   - fullMaterial : " << m_fullMaterial << std::endl;
  sl << "   - split factor : " << m_splitFactor << std::endl;
  return sl;
}

}  // namespace Acts
