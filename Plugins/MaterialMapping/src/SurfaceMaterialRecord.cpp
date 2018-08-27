// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialRecord.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/SurfaceMaterialRecord.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::SurfaceMaterialRecord::SurfaceMaterialRecord(const Surface&    surface,
                                                   const BinUtility& binUtility)
  : m_surface(&surface)
  , m_binUtility(std::make_unique<const BinUtility>(binUtility))
  , m_mappedMaterial()
{
  // reserve
  m_mappedMaterial.reserve(m_binUtility->bins(1));
  for (size_t ibin1 = 0; ibin1 < m_binUtility->bins(1); ++ibin1) {
    // create the vector for the push_back
    RecordVector mappedVector;
    mappedVector.reserve(m_binUtility->bins(0));
    for (size_t ibin0 = 0; ibin0 < m_binUtility->bins(0); ++ibin0) {
      mappedVector.push_back(RecordBin(MaterialProperties(), 0));
    }
    m_mappedMaterial.push_back(mappedVector);
  }
}

Acts::SurfaceMaterialRecord::SurfaceMaterialRecord(
    const SurfaceMaterialRecord& lmrecord)
  : m_surface(lmrecord.m_surface)
  , m_binUtility(std::make_unique<const BinUtility>(*lmrecord.m_binUtility))
  , m_mappedMaterial(lmrecord.m_mappedMaterial)
{
}

Acts::SurfaceMaterialRecord&
Acts::SurfaceMaterialRecord::operator=(const SurfaceMaterialRecord& lmrecord)
{
  if (this != &lmrecord) {
    m_surface    = lmrecord.m_surface;
    m_binUtility = std::make_unique<const BinUtility>(*lmrecord.m_binUtility);
    m_mappedMaterial = lmrecord.m_mappedMaterial;
  }
  return (*this);
}

void
Acts::SurfaceMaterialRecord::assignEmptyStep(const Vector3D& mPosition)
{
  // get the bins corresponding to the position
  size_t bin0 = m_binUtility->bin(mPosition, 0);
  size_t bin1 = m_binUtility->bin(mPosition, 1);

  // increase the number of entries for this material bin
  // always do that, also for empty assignments
  // as they will take part in the average as well
  m_mappedMaterial[bin1][bin0].second++;
}

void
Acts::SurfaceMaterialRecord::assignMaterialStep(const MaterialStep& mStep,
                                                double pathCorrection)
{
  // get the position
  auto mPosition = mStep.position();

  // get the bins corresponding to the position
  size_t bin0 = m_binUtility->bin(mPosition, 0);
  size_t bin1 = m_binUtility->bin(mPosition, 1);

  // increase the number of entries for this material bin
  // always do that, also for empty assignments
  // as they will take part in the average as well
  m_mappedMaterial[bin1][bin0].second++;

  // it's time to correct for the path correction
  float thickness = mStep.materialProperties().thickness() / pathCorrection;
  // now simple add the material component
  float x0  = mStep.materialProperties().material().X0();
  float l0  = mStep.materialProperties().material().L0();
  float A   = mStep.materialProperties().averageA();
  float Z   = mStep.materialProperties().averageZ();
  float rho = mStep.materialProperties().averageRho();

  // access the material collected so far already
  // and simply add it
  auto materialBin = m_mappedMaterial[bin1][bin0];
  if (materialBin.first) {
    thickness += materialBin.first.thickness();
    x0 += materialBin.first.material().X0();
    l0 += materialBin.first.material().L0();
    A += materialBin.first.averageA();
    Z += materialBin.first.averageZ();
    rho += materialBin.first.averageRho();
  }

  // set the new cummulated material to the bin
  // in a final step, all parameters need to be divided
  // by the number of bin hits
  Material cummulatedMaterial(x0, l0, A, Z, rho);
  m_mappedMaterial[bin1][bin0].first
      = MaterialProperties(cummulatedMaterial, thickness);
}
