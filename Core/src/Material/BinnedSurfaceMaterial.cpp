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

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial()
  : Acts::SurfaceMaterial(), m_binUtility(nullptr)
{
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(Acts::BinUtility& binutility)
  : Acts::SurfaceMaterial(), m_binUtility(binutility.clone())
{
  // reserve
  m_fullMaterial.reserve(binutility.max(1) + 1);
  for (unsigned int ibin2 = 0; ibin2 < (unsigned int)binutility.max(1) + 1;
       ++ibin2) {
    // create the vector for the push_back
    Acts::MaterialPropertiesVector matVec;
    matVec.reserve(binutility.max(0) + 1);
    for (unsigned int ibin = 0; ibin < (unsigned int)binutility.max(0) + 1;
         ++ibin)
      matVec.push_back(nullptr);
    m_fullMaterial.push_back(matVec);
  }
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const Acts::BinUtility&         binutility,
    const MaterialPropertiesVector& fullProperties,
    double                          splitFactor)
  : Acts::SurfaceMaterial(splitFactor), m_binUtility(binutility.clone())
{
  // constructor from a single vector
  clearMaterial();
  // fill the material with deep copy
  m_fullMaterial.push_back(fullProperties);
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const Acts::BinUtility&         binutility,
    const MaterialPropertiesMatrix& fullProperties,
    double                          splitFactor)
  : Acts::SurfaceMaterial(splitFactor)
  , m_binUtility(binutility.clone())
  , m_fullMaterial(fullProperties)
{
}

Acts::BinnedSurfaceMaterial::~BinnedSurfaceMaterial()
{
  delete m_binUtility;
  clearMaterial();
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const Acts::BinnedSurfaceMaterial& lmp)
  : Acts::SurfaceMaterial(lmp), m_binUtility(lmp.m_binUtility->clone())
{
  // clear the material
  clearMaterial();
  // and fill the material
  fillMaterial(lmp.m_fullMaterial);
}

Acts::BinnedSurfaceMaterial&
Acts::BinnedSurfaceMaterial::operator=(const Acts::BinnedSurfaceMaterial& lmp)
{
  if (this != &lmp) {
    Acts::SurfaceMaterial::operator=(lmp);
    // first delete everything
    delete m_binUtility;
    // reassign
    m_binUtility = lmp.binUtility()->clone();
    clearMaterial();
    // reassign teh material
    fillMaterial(lmp.m_fullMaterial);
  }
  return (*this);
}

Acts::BinnedSurfaceMaterial*
Acts::BinnedSurfaceMaterial::clone() const
{
  return new Acts::BinnedSurfaceMaterial(*this);
}

void
Acts::BinnedSurfaceMaterial::clearMaterial()
{
  // clear the full material
  for (auto& matMatrixIter : m_fullMaterial) {
    for (auto& matIter : matMatrixIter) delete matIter;
  }

  m_fullMaterial.clear();
}

void
Acts::BinnedSurfaceMaterial::fillMaterial(
    const Acts::MaterialPropertiesMatrix& matMatrix)
{
  m_fullMaterial.reserve(m_binUtility->max(1) + 1);
  for (auto& matMatrixIter : matMatrix) {
    // the vector to be copied
    Acts::MaterialPropertiesVector matVector;
    matVector.reserve(m_binUtility->max(0) + 1);
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
  unsigned int imat2 = 0;
  for (auto& matMatrixIter : m_fullMaterial) {
    unsigned int imat1 = 0;
    // the vector iterator
    for (auto& matIter : matMatrixIter) {
      if (matIter) {
        Acts::MaterialProperties* mprop = matIter->clone();
        (*mprop) *= scale;
        delete matIter;
        m_fullMaterial[imat2][imat1] = mprop;
      }
      ++imat1;
    }
    ++imat2;
  }
  return (*this);
}

const Acts::MaterialProperties*
Acts::BinnedSurfaceMaterial::material(const Acts::Vector2D& lp) const
{
  if (!m_fullMaterial.size() || !m_binUtility) return nullptr;
  // the first bin
  size_t ibin1 = m_binUtility->bin(lp, 0);
  size_t ibin2 = m_binUtility->max(1) ? m_binUtility->bin(lp, 1) : 0;
  return m_fullMaterial[ibin2][ibin1];
}

const Acts::MaterialProperties*
Acts::BinnedSurfaceMaterial::material(const Acts::Vector3D& gp) const
{
  if (!m_fullMaterial.size() || !m_binUtility) return nullptr;
  // the first bin
  size_t ibin1 = m_binUtility->bin(gp, 0);
  size_t ibin2 = m_binUtility->max(1) ? m_binUtility->bin(gp, 1) : 0;
  return m_fullMaterial[ibin2][ibin1];
}

std::ostream&
Acts::BinnedSurfaceMaterial::dump(std::ostream& sl) const
{
  sl << "Acts::BinnedSurfaceMaterial : " << std::endl;
  sl << "   - Number of Material bins (1/2) : " << m_binUtility->max(0) + 1
     << " / " << m_binUtility->max(1) + 1 << std::endl;
  sl << "   - Parse full update material    : " << std::endl;  //
  // output  the full material
  unsigned int imat2 = 0;
  for (auto& matMatrixIter : m_fullMaterial) {
    unsigned int imat1 = 0;
    // the vector iterator
    for (auto& matIter : matMatrixIter) {
      if (matIter)
        sl << " Bin [" << imat2 << "][" << imat1 << "] - " << (*matIter)
           << std::endl;
      else
        sl << " Bin [" << imat2 << "][" << imat1 << "] -  empty " << std::endl;
      ++imat1;
    }
    ++imat2;
  }
  return sl;
}
