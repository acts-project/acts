// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DigitizationModule.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Digitization/DigitizationModule.hpp"

#include <cmath>
#include <cstddef>
#include <utility>

Acts::DigitizationModule::DigitizationModule(
    std::shared_ptr<const Segmentation> moduleSegmentation,
    double halfThickness, int readoutDirection, double lorentzAngle,
    double energyThreshold, bool analogue)
    :

      m_halfThickness(halfThickness),
      m_readoutDirection(readoutDirection),
      m_lorentzAngle(lorentzAngle),
      m_tanLorentzAngle(tan(lorentzAngle)),
      m_energyThreshold(energyThreshold),
      m_analogue(analogue),
      m_segmentation(std::move(moduleSegmentation)),
      m_boundarySurfaces(),
      m_segmentationSurfacesX(),
      m_segmentationSurfacesY() {
  m_segmentation->createSegmentationSurfaces(
      m_boundarySurfaces, m_segmentationSurfacesX, m_segmentationSurfacesY,
      halfThickness, readoutDirection, lorentzAngle);
}

const Acts::SurfacePtrVector Acts::DigitizationModule::segmentationSurfaces(
    const Acts::DigitizationCell& entryCids,
    const Acts::DigitizationCell& exitCids) const {
  SurfacePtrVector sSurfaces;

  auto startbinX = entryCids.channel0;
  auto endbinX = exitCids.channel0;
  // swap if needed
  if (startbinX > endbinX) {
    std::swap(startbinX, endbinX);
  }
  // now cash in the rest
  for (; startbinX <= endbinX; ++startbinX) {
    sSurfaces.push_back(m_segmentationSurfacesX[startbinX]);
  }

  // start bin, end bin
  auto startbinY = entryCids.channel1;
  auto endbinY = exitCids.channel1;
  // swap if needed
  if (startbinY > endbinY) {
    std::swap(startbinY, endbinY);
  }
  // now cash in the rest
  for (; startbinY <= endbinY; ++startbinY) {
    sSurfaces.push_back(m_segmentationSurfacesY[startbinY]);
  }

  // return what you have
  return sSurfaces;
}

Acts::SurfacePtrVector Acts::DigitizationModule::stepSurfaces(
    const Vector3& start, const Vector3& end) const {
  // prepare the return vector
  SurfacePtrVector stepSurfaces;

  const DigitizationCell startCell = m_segmentation->cell(start);
  const DigitizationCell endCell = m_segmentation->cell(end);

  // go along x - first with the naive binning (i.e. w.o lorentz angle)
  size_t sCellX = startCell.channel0;
  size_t eCellX = endCell.channel0;
  if (sCellX > eCellX) {
    std::swap(sCellX, eCellX);
  }
  // now take the boundaries as well
  if (sCellX > 0) {
    --sCellX;
  }
  ++eCellX;  // @TODO check : safe because we can assume to have eCell+1
  // the surfaces along Y are easy, just the bin surfaces
  size_t sCellY = startCell.channel1;
  size_t eCellY = endCell.channel1;
  if (sCellY > eCellY) {
    std::swap(sCellY, eCellY);
  }
  // reserve - be safe
  stepSurfaces.reserve((eCellY - sCellY) + (eCellX - sCellX) + 2);
  // now fill the x surfaces
  for (; sCellX <= eCellX && sCellX < m_segmentationSurfacesX.size();
       ++sCellX) {
    stepSurfaces.push_back(m_segmentationSurfacesX[sCellX]);
  }
  // end fill the y surfaces
  for (; sCellY <= eCellY && sCellY < m_segmentationSurfacesY.size();
       ++sCellY) {
    stepSurfaces.push_back(m_segmentationSurfacesY[sCellY]);
  }
  // return the lot
  return stepSurfaces;
}
