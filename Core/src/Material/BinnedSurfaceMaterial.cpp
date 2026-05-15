// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/BinnedSurfaceMaterial.hpp"

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <ostream>
#include <utility>
#include <vector>

namespace Acts {

BinnedSurfaceMaterial::BinnedSurfaceMaterial(const BinUtility& binUtility,
                                             MaterialSlabVector fullProperties,
                                             double splitFactor,
                                             MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType), m_binUtility(binUtility) {
  // Catch the cases where the bin Utility is only-1D and add either first or
  // second dummy dimension
  if (binUtility.dimensions() == 1u) {
    auto aDir = binUtility.binningData()[0u].binvalue;
    // Only act if the axis direction is not already the first one in an allowed
    // local frame
    if (aDir != AxisDirection::AxisX && aDir != AxisDirection::AxisR &&
        aDir != AxisDirection::AxisRPhi) {
      m_axesSwapped = true;
    }
  }
  m_fullMaterial.push_back(std::move(fullProperties));
}

BinnedSurfaceMaterial::BinnedSurfaceMaterial(const BinUtility& binUtility,
                                             MaterialSlabMatrix fullProperties,
                                             double splitFactor,
                                             MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_binUtility(binUtility),
      m_fullMaterial(std::move(fullProperties)) {}

BinnedSurfaceMaterial& BinnedSurfaceMaterial::scale(double factor) {
  for (auto& materialVector : m_fullMaterial) {
    for (auto& materialBin : materialVector) {
      materialBin.scaleThickness(factor);
    }
  }
  return (*this);
}

const MaterialSlab& BinnedSurfaceMaterial::materialSlab(
    const Vector2& lp) const {
  // the first bin
  std::size_t b0lu = m_axesSwapped ? 1 : 0;
  std::size_t b1lu = m_axesSwapped ? 0 : 1;
  std::size_t ibin0 = m_binUtility.bin(lp, b0lu);
  std::size_t ibin1 =
      m_binUtility.max(1) != 0u ? m_binUtility.bin(lp, b1lu) : 0;
  return m_fullMaterial[ibin1][ibin0];
}

std::vector<AxisDirection> BinnedSurfaceMaterial::materialAxisDirections()
    const {
  std::vector<AxisDirection> axisDirs;
  for (const auto& bd : m_binUtility.binningData()) {
    axisDirs.push_back(bd.binvalue);
  }
  return axisDirs;
}

const MaterialSlab& BinnedSurfaceMaterial::materialSlab(
    const Vector3& gp) const {
  // the first bin
  std::size_t ibin0 = m_binUtility.bin(gp, 0);
  std::size_t ibin1 = m_binUtility.max(1) != 0u ? m_binUtility.bin(gp, 1) : 0;
  return m_fullMaterial[ibin1][ibin0];
}

std::ostream& BinnedSurfaceMaterial::toStream(std::ostream& sl) const {
  sl << "BinnedSurfaceMaterial : " << std::endl;
  sl << "   - Number of Material bins [0,1] : " << m_binUtility.max(0) + 1
     << " / " << m_binUtility.max(1) + 1 << std::endl;
  sl << "   - Parse full update material    : " << std::endl;  //
  // output  the full material
  unsigned int imat1 = 0;
  for (auto& materialVector : m_fullMaterial) {
    unsigned int imat0 = 0;
    // the vector iterator
    for (auto& materialBin : materialVector) {
      sl << " Bin [" << imat1 << "][" << imat0 << "] - " << (materialBin);
      ++imat0;
    }
    ++imat1;
  }
  sl << "  - BinUtility: " << m_binUtility << std::endl;
  return sl;
}

}  // namespace Acts
