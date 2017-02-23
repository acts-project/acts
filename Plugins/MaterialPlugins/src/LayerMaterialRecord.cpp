// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerMaterialRecord.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/LayerMaterialRecord.hpp"

Acts::LayerMaterialRecord::LayerMaterialRecord()
  : m_binUtility(nullptr), m_materialMatrix()
{
}

Acts::LayerMaterialRecord::LayerMaterialRecord(const BinUtility* binutility)
  : m_binUtility(binutility), m_materialMatrix()
{
  // reserve
  m_materialMatrix.reserve(m_binUtility->max(1) + 1);
  for (unsigned int ibin1 = 0; ibin1 < (unsigned int)m_binUtility->max(1) + 1;
       ++ibin1) {
    // create the vector for the push_back
    Acts::MaterialPropertiesVector matVec;
    matVec.reserve(m_binUtility->max(0) + 1);
    for (unsigned int ibin = 0; ibin < (unsigned int)m_binUtility->max(0) + 1;
         ++ibin)
      matVec.push_back(new MaterialProperties(0., 0., 0., 0., 0., 0.));
    m_materialMatrix.push_back(matVec);
  }
}

Acts::LayerMaterialRecord::LayerMaterialRecord(
    const LayerMaterialRecord& lmrecord)
  : m_binUtility(lmrecord.m_binUtility)
  , m_materialMatrix(lmrecord.m_materialMatrix)
{
}

Acts::LayerMaterialRecord*
Acts::LayerMaterialRecord::clone() const
{
  return (new LayerMaterialRecord(*this));
}

Acts::LayerMaterialRecord&
Acts::LayerMaterialRecord::operator=(const LayerMaterialRecord& lmrecord)
{
  if (this != &lmrecord) {
    m_binUtility     = lmrecord.m_binUtility;
    m_materialMatrix = lmrecord.m_materialMatrix;
  }
  return (*this);
}

void
Acts::LayerMaterialRecord::addLayerMaterialProperties(
    const Acts::Vector3D&           pos,
    const Acts::MaterialProperties* newMaterial)
{
  // get the bins corresponding to the position
  size_t bin0 = m_binUtility->bin(pos, 0);
  size_t bin1 = m_binUtility->bin(pos, 1);
  // get the material which might be there already, add new material and
  // weigh it
  const Acts::MaterialProperties* material = m_materialMatrix.at(bin1).at(bin0);
  float                           newThickness = newMaterial->thickness();
  float                           newRho       = newMaterial->averageRho();
  float                           thickness    = 0.;
  float                           rho          = 0.;
  float                           x0           = 0.;
  float                           l0           = 0.;
  float                           A            = 0.;
  float                           Z            = 0.;
  // access the new material properties
  if (material) {
    thickness += material->thickness();
    rho += material->averageRho();
    x0 += material->x0();
    l0 += material->l0();
    A += material->averageA();
    Z += material->averageZ();
  }
  // add the new material properties and weigh them
  thickness += newThickness;
  rho += newRho * newThickness;
  x0 += newMaterial->x0() * newThickness;
  l0 += newMaterial->l0() * newThickness;
  A += newMaterial->averageA() * newRho;
  Z += newMaterial->averageZ() * newRho;
  // set the new current material (not averaged yet)
  const Acts::Material updatedMaterial(x0, l0, A, Z, rho);
  m_materialMatrix.at(bin1).at(bin0)->setMaterial(updatedMaterial, thickness);
}

void
Acts::LayerMaterialRecord::averageMaterial()
{
  // access the bins
  size_t bins0 = m_binUtility->bins(0);
  size_t bins1 = m_binUtility->bins(1);
  // loop through the material properties matrix and average
  for (size_t bin0 = 0; bin0 < bins0; bin0++) {
    for (size_t bin1 = 0; bin1 < bins1; bin1++) {
      const Acts::MaterialProperties* material
          = m_materialMatrix.at(bin1).at(bin0);

      float thickness = material->thickness();
      float rho       = material->averageRho();
      float x0        = material->x0();
      float l0        = material->l0();
      float A         = material->averageA();
      float Z         = material->averageZ();
      // divide
      if (x0 != 0.) x0 /= thickness;
      if (l0 != 0.) l0 /= thickness;
      if (A != 0.) A /= rho;
      if (Z != 0.) Z /= rho;
      if (rho != 0.) rho /= thickness;
      if (thickness != 0.) thickness /= material->entries();
      // set the new current material
      const Acts::Material updatedMaterial(x0, l0, A, Z, rho);
      m_materialMatrix.at(bin1).at(bin0)->setMaterial(updatedMaterial,
                                                      thickness);
    }  // b2
  }    // b1
}

std::shared_ptr<const Acts::BinnedSurfaceMaterial>
Acts::LayerMaterialRecord::layerMaterial() const
{
  // return the binned surface material
  return (std::make_shared<const Acts::BinnedSurfaceMaterial>(
      *m_binUtility, m_materialMatrix));
}
