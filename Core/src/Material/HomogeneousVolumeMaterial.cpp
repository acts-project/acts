// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/HomogeneousVolumeMaterial.hpp"

#include "Acts/Material/Material.hpp"

#include <ostream>

namespace Acts {

HomogeneousVolumeMaterial::HomogeneousVolumeMaterial(const Material& material)
    : m_material(material) {}

const Material HomogeneousVolumeMaterial::material(
    const Vector3& /*position*/) const {
  return m_material;
}

std::ostream& HomogeneousVolumeMaterial::toStream(std::ostream& sl) const {
  sl << "HomogeneousVolumeMaterial : " << std::endl;
  sl << "   - material : " << m_material << std::endl;
  return sl;
}

}  // namespace Acts
