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
  : m_binUtility(nullptr), m_materialMatrix(), m_matStepsAndAssignedPos()
{
}

Acts::LayerMaterialRecord::LayerMaterialRecord(const BinUtility* binutility)
  : m_binUtility(binutility), m_materialMatrix(), m_matStepsAndAssignedPos()
{
  // reserve
  m_materialMatrix.reserve(m_binUtility->max(1) + 1);
  for (unsigned int ibin1 = 0; ibin1 < (unsigned int)m_binUtility->max(1) + 1;
       ++ibin1) {
    // create the vector for the push_back
    Acts::MaterialPropertiesVector matVec;
    matVec.reserve(m_binUtility->max(0) + 1);
    for (unsigned int ibin0 = 0; ibin0 < (unsigned int)m_binUtility->max(0) + 1;
         ++ibin0)
      matVec.push_back(new MaterialProperties(0., 0., 0., 0., 0., 0., 0., 0));
    m_materialMatrix.push_back(matVec);
  }
}

Acts::LayerMaterialRecord::LayerMaterialRecord(
    const LayerMaterialRecord& lmrecord)
  : m_binUtility(lmrecord.m_binUtility)
  , m_materialMatrix(lmrecord.m_materialMatrix)
  , m_matStepsAndAssignedPos(lmrecord.m_matStepsAndAssignedPos)
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
    m_matStepsAndAssignedPos.clear();
  }
  return (*this);
}

void
Acts::LayerMaterialRecord::addLayerMaterialProperties(
    const Acts::Vector3D&                 pos,
    const std::vector<Acts::MaterialStep> layerMaterialSteps)
{
  m_matStepsAndAssignedPos.push_back(std::make_pair(layerMaterialSteps, pos));
  // sum up all material at this point for this layer
  float newThickness = 0.;
  float newRho       = 0.;
  float newA         = 0.;
  float newZ         = 0.;
  float newtInX0     = 0.;
  float newtInL0     = 0.;

  for (auto& layerStep : layerMaterialSteps) {
    float t       = layerStep.material().thickness();
    float density = layerStep.material().averageRho();
    newThickness += t;

    newRho += density * t;
    newA += layerStep.material().averageA() * density * t;
    newZ += layerStep.material().averageZ() * density * t;
    newtInX0 += (layerStep.material().x0() != 0.)
        ? t / layerStep.material().x0()
        : 0.;

    newtInL0 += (layerStep.material().l0() != 0.)
        ? t / layerStep.material().l0()
        : 0.;
  }

  newA /= (newRho != 0.) ? newRho : 1.;
  newZ /= (newRho != 0.) ? newRho : 1.;

  newRho /= (newThickness != 0.) ? newThickness : 1.;

  // now add it at the corresponding assigned position
  // get the bins corresponding to the position
  size_t bin0 = m_binUtility->bin(pos, 0);
  size_t bin1 = m_binUtility->bin(pos, 1);
  // get the material which might be there already, add new material and
  // weigh it
  const Acts::MaterialProperties* material = m_materialMatrix.at(bin1).at(bin0);
  float                           thickness = 0.;
  float                           rho       = 0.;
  float                           tInX0     = 0.;
  float                           tInL0     = 0.;
  float                           A         = 0.;
  float                           Z         = 0.;
  // access the old material properties
  if (material) {
    thickness += material->thickness();
    rho += material->averageRho();
    A += material->averageA();
    Z += material->averageZ();
    tInX0 += material->thicknessInX0();
    tInL0 += material->thicknessInL0();
  }
  // sum up material properties (different MaterialTracks)
  thickness += newThickness;
  rho += newRho;
  A += newA;
  Z += newZ;
  tInX0 += newtInX0;
  tInL0 += newtInL0;

  float x0 = (thickness != 0. && tInX0 != 0.) ? thickness / tInX0 : 0.;
  float l0 = (thickness != 0. && tInL0 != 0.) ? thickness / tInL0 : 0.;

  // set the new current material (not averaged yet)
  const Acts::Material updatedMaterial(x0, l0, A, Z, rho);
  // pick the number of entries for the next material entry
  size_t entries = m_materialMatrix.at(bin1).at(bin0)->entries();
  // set the material with the right number of entries
  m_materialMatrix.at(bin1).at(bin0)->setMaterial(
      updatedMaterial, thickness, entries);
  // increase the number of entries for this material
  m_materialMatrix.at(bin1).at(bin0)->addEntry();
}

void
Acts::LayerMaterialRecord::averageMaterial()
{
  // access the bins
  size_t bins0 = m_binUtility->bins(0);
  size_t bins1 = m_binUtility->bins(1);
  // loop through the material properties matrix and average
  for (size_t bin1 = 0; bin1 < bins1; bin1++) {
    for (size_t bin0 = 0; bin0 < bins0; bin0++) {
      const Acts::MaterialProperties* material
          = m_materialMatrix.at(bin1).at(bin0);

      float thickness = material->thickness();
      float rho       = material->averageRho();
      float tInX0     = material->thicknessInX0();
      float tInL0     = material->thicknessInL0();
      float A         = material->averageA();
      float Z         = material->averageZ();
      // caclulate mean value by dividing summed up material of all track
      // records through the number of track record for each bin
      size_t n = material->entries();
      if (n != 0) {
        tInX0 /= n;
        tInL0 /= n;
        A /= n;
        Z /= n;
        rho /= n;
        thickness /= n;
      }

      float x0 = (thickness != 0. && tInX0 != 0.) ? thickness / tInX0 : 0.;
      float l0 = (thickness != 0. && tInL0 != 0.) ? thickness / tInL0 : 0.;
      // set the new current material (resetting number of entries)
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
