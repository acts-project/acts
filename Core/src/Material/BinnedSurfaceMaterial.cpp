// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/BinnedSurfaceMaterial.hpp"

#include "Acts/Material/MaterialSlab.hpp"

#include <ostream>
#include <utility>
#include <vector>

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const DirectedProtoAxis& dProtoAxis, MaterialSlabVector fullProperties,
    double splitFactor, Acts::MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_directedProtoAxes({dProtoAxis}) {
  // fill the material with deep copy
  m_fullMaterial.push_back(std::move(fullProperties));
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const DirectedProtoAxis& dProtoAxis0, const DirectedProtoAxis& dProtoAxis1,
    MaterialSlabMatrix fullProperties, double splitFactor,
    Acts::MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_directedProtoAxes({dProtoAxis0, dProtoAxis1}),
      m_fullMaterial(std::move(fullProperties)) {}

Acts::BinnedSurfaceMaterial& Acts::BinnedSurfaceMaterial::scale(double factor) {
  for (auto& materialVector : m_fullMaterial) {
    for (auto& materialBin : materialVector) {
      materialBin.scaleThickness(factor);
    }
  }
  return (*this);
}

const Acts::MaterialSlab& Acts::BinnedSurfaceMaterial::materialSlab(
    const Vector2& lp) const {
  const auto& dpa0 = m_directedProtoAxes[0];
  std::size_t ibin0 = dpa0.getAxis().getBin(lp[0u]);
  std::size_t ibin1 = 0u;
  if (m_directedProtoAxes.size() > 1) {
    const auto& dpa1 = m_directedProtoAxes[1];
    ibin1 = dpa1.getAxis().getBin(lp[1u]);
  }
  return m_fullMaterial[ibin1][ibin0];
}

const Acts::MaterialSlab& Acts::BinnedSurfaceMaterial::materialSlab(
    const Acts::Vector3& gp) const {
  const auto& dpa0 = m_directedProtoAxes[0];
  double v0 = VectorHelpers::cast(gp, dpa0.getAxisDirection());
  std::size_t ibin0 = dpa0.getAxis().getBin(v0);
  std::size_t ibin1 = 0u;
  if (m_directedProtoAxes.size() > 1) {
    const auto& dpa1 = m_directedProtoAxes[1];
    double v1 = VectorHelpers::cast(gp, dpa1.getAxisDirection());
    ibin1 = dpa1.getAxis().getBin(v1);
  }
  return m_fullMaterial[ibin1][ibin0];
}

std::ostream& Acts::BinnedSurfaceMaterial::toStream(std::ostream& sl) const {
  sl << "Acts::BinnedSurfaceMaterial : " << std::endl;
  std::size_t nbins0 = m_directedProtoAxes[0].getAxis().getNBins();
  std::size_t nbins1 = 0u;
  if (m_directedProtoAxes.size() > 1) {
    nbins1 = m_directedProtoAxes[1].getAxis().getNBins();
  }
  sl << "   - Number of Material bins [0,1] : " << nbins0 << " / " << nbins1
     << std::endl;
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
  sl << "  - ProtoAxes: " << m_directedProtoAxes << std::endl;
  return sl;
}
