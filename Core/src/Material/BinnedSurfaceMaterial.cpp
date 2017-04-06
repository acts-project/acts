// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedSurfaceMaterial.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinUtility&               binUtility,
    const MaterialPropertiesVector& fullProperties,
    double                          splitFactor,
    size_t                          entries)
  : SurfaceMaterial(splitFactor)
  , m_binUtility(binUtility)
  , m_entries(entries)
{
  // fill the material with deep copy
  m_fullMaterial.push_back(fullProperties);
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinUtility&               binUtility,
    const MaterialPropertiesMatrix& fullProperties,
    double                          splitFactor,
    size_t                          entries)
  : Acts::SurfaceMaterial(splitFactor)
  , m_binUtility(binUtility)
  , m_fullMaterial(fullProperties)
  , m_entries(entries)
{
}

Acts::BinnedSurfaceMaterial::~BinnedSurfaceMaterial()
{
  clearMaterial();
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinnedSurfaceMaterial& lmp)
  : SurfaceMaterial(lmp)
  , m_binUtility(lmp.m_binUtility)
  , m_entries(lmp.m_entries)
{
  // clear the material
  clearMaterial();
  // and fill the material
  fillMaterial(lmp.m_fullMaterial);
}

Acts::BinnedSurfaceMaterial&
Acts::BinnedSurfaceMaterial::operator=(const BinnedSurfaceMaterial& lmp)
{
  if (this != &lmp) {
    // reassign
    SurfaceMaterial::operator=(lmp);
    m_binUtility = lmp.binUtility();
    m_entries = lmp.m_entries;
    // clear
    clearMaterial();
    // reassign the material
    fillMaterial(lmp.m_fullMaterial);
  }
  return (*this);
}

Acts::BinnedSurfaceMaterial*
Acts::BinnedSurfaceMaterial::clone() const
{
  return new BinnedSurfaceMaterial(*this);
}

void
Acts::BinnedSurfaceMaterial::clearMaterial()
{
  // clear the full material
  for (auto& matMatrixIter : m_fullMaterial) {
    for (auto& matIter : matMatrixIter) 
        delete matIter;
  }
  m_fullMaterial.clear();
}

void
Acts::BinnedSurfaceMaterial::fillMaterial(
    const Acts::MaterialPropertiesMatrix& matMatrix)
{
  m_fullMaterial.reserve(m_binUtility.max(1) + 1);
  for (auto& matMatrixIter : matMatrix) {
    // the vector to be copied
    Acts::MaterialPropertiesVector matVector;
    matVector.reserve(m_binUtility.max(0) + 1);
    // reassign
    for (auto& matIter : matMatrixIter)
      matVector.push_back(matIter ? matIter->clone() : nullptr);
    m_fullMaterial.push_back(matVector);
  }
}

Acts::BinnedSurfaceMaterial&
Acts::BinnedSurfaceMaterial::operator*=(double scale)
{
  // scale the full material
  unsigned int imat1 = 0;
  for (auto& matMatrixIter : m_fullMaterial) {
    unsigned int imat0 = 0;
    // the vector iterator
    for (auto& matIter : matMatrixIter) {
      if (matIter) {
        Acts::MaterialProperties* mprop = matIter->clone();
        (*mprop) *= scale;
        delete matIter;
        m_fullMaterial[imat1][imat0] = mprop;
      }
      ++imat0;
    }
    ++imat1;
  }
  return (*this);
}

const Acts::MaterialProperties*
Acts::BinnedSurfaceMaterial::material(const Vector2D& lp) const
{
  if (!m_fullMaterial.size()) return nullptr;
  // the first bin
  size_t ibin0 = m_binUtility.bin(lp, 0);
  size_t ibin1 = m_binUtility.max(1) ? m_binUtility.bin(lp, 1) : 0;
  return m_fullMaterial[ibin1][ibin0];
}

const Acts::MaterialProperties*
Acts::BinnedSurfaceMaterial::material(const Acts::Vector3D& gp) const
{
  if (!m_fullMaterial.size()) return nullptr;
  // the first bin
  size_t ibin0 = m_binUtility.bin(gp, 0);
  size_t ibin1 = m_binUtility.max(1) ? m_binUtility.bin(gp, 1) : 0;
  return m_fullMaterial[ibin1][ibin0];
}

std::ostream&
Acts::BinnedSurfaceMaterial::dump(std::ostream& sl) const
{
  sl << "Acts::BinnedSurfaceMaterial : " << std::endl;
  sl << "   - Number of Material bins [0,1] : " << m_binUtility.max(0) + 1
     << " / " << m_binUtility.max(1) + 1 << std::endl;
  sl << "   - Parse full update material    : " << std::endl;  //
  // output  the full material
  unsigned int imat1 = 0;
  for (auto& matMatrixIter : m_fullMaterial) {
    unsigned int imat0 = 0;
    // the vector iterator
    for (auto& matIter : matMatrixIter) {
      if (matIter)
        sl << " Bin [" << imat1 << "][" << imat0 << "] - " << (*matIter)
           << std::endl;
      else
        sl << " Bin [" << imat1 << "][" << imat0 << "] -  empty " << std::endl;
      ++imat0;
    }
    ++imat1;
  }
  return sl;
}
