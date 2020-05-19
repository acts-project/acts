// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"

Acts::HomogeneousVolumeMaterial::HomogeneousVolumeMaterial(
    const Material& material)
    : m_material(material) {}

std::ostream& Acts::HomogeneousVolumeMaterial::toStream(
    std::ostream& sl) const {
  sl << "Acts::HomogeneousVolumeMaterial : " << std::endl;
  sl << "   - material : " << m_material << std::endl;
  return sl;
}
