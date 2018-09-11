// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/AccumulatedSurfaceMaterial.hpp"

// Default Constructor - for homogeneous material
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(double splitFactor)
  : m_splitFactor(splitFactor)
{
  AccumulatedVector accMat = {{AccumulatedMaterialProperties()}};
  m_accumulatedMaterial    = {{accMat}};
}

// Binned Material constructor with split factor
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    const BinUtility& binUtility,
    double            splitFactor)
  : m_binUtility(binUtility), m_splitFactor(splitFactor)
{
  size_t            bins0 = m_binUtility.bins(0);
  size_t            bins1 = m_binUtility.bins(1);
  AccumulatedVector accVec(bins0, AccumulatedMaterialProperties());
  m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
}

// Assign a material properites object
void
Acts::AccumulatedSurfaceMaterial::assign(const Vector2D&           lp,
                                         const MaterialProperties& mp)
{
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0] += mp;
  } else {
    size_t bin0 = m_binUtility.bin(lp, 0);
    size_t bin1 = m_binUtility.bin(lp, 1);
    m_accumulatedMaterial[bin1][bin0] += mp;
  }
}

// Assign a material properites object
void
Acts::AccumulatedSurfaceMaterial::assign(const Vector3D&           gp,
                                         const MaterialProperties& mp)
{
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0] += mp;
  } else {
    std::array<size_t, 3> bTriple = m_binUtility.binTriple(gp);
    m_accumulatedMaterial[bTriple[1]][bTriple[0]] += mp;
  }
}

// Assign a vector of material properites object
void
Acts::AccumulatedSurfaceMaterial::assign(
    const Vector3D& gp,
    const std::vector<std::pair<MaterialProperties, Vector3D>>& mps)
{
  if (m_binUtility.dimensions() == 0) {
    for (auto mp : mps) {
      m_accumulatedMaterial[0][0] += mp.first;
    }
  } else {
    std::array<size_t, 3> bTriple = m_binUtility.binTriple(gp);
    for (auto mp : mps) {
      m_accumulatedMaterial[bTriple[1]][bTriple[0]] += mp.first;
    }
  }
}

// Average the information accumulated during one event
void
Acts::AccumulatedSurfaceMaterial::eventAverage()
{
  for (auto& matVec : m_accumulatedMaterial) {
    for (auto& mat : matVec) {
      mat.eventAverage();
    }
  }
}
