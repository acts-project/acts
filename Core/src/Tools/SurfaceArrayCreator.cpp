// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Utilities/MsgMacros.hpp"
#include "ACTS/Utilities/Definitions.hpp"

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Acts::Surface*>& surfaces,
    double                                   R,
    double                                   minPhi,
    double                                   maxPhi,
    double                                   halfZ,
    size_t                                   binsPhi,
    size_t                                   binsZ,
    std::shared_ptr<Acts::Transform3D>       transform) const
{
  ACTS_DEBUG("Creating a SurfaceArray on a cylinder with grid in phi x z = "
            << binsPhi
            << " x "
            << binsZ);

  // create the (plain) binUtility - with the transform
  auto arrayUtility = std::make_unique<BinUtility>(binsPhi, minPhi, maxPhi, closed, binPhi, transform);
  (*arrayUtility) += BinUtility(binsZ, -halfZ, halfZ, open, binZ);

  // the z step and the phi step
  double zStep   = (2 * halfZ / binsZ);
  double phiStep = (maxPhi - minPhi) / binsPhi;

  // prepare the surface system in phi x z
  std::vector<std::vector<std::pair<SurfacePosition, Vector3D>>> phizSystem;
  phizSystem.reserve(binsZ);
  for (size_t iZ = 0; iZ < binsZ; ++iZ) {
    // the current z value
    double currentZ = -halfZ + (iZ + 0.5) * zStep;
    // the current phi row
    std::vector<std::pair<SurfacePosition, Vector3D>> phiSystem;
    phiSystem.reserve(binsPhi);
    for (size_t iPhi = 0; iPhi < binsPhi; ++iPhi) {
      // the current phi value
      double currentPhi = minPhi + (iPhi + 0.5) * phiStep;
      // the bin position
      Vector3D binPosition = transform
          ? ((*transform)
             * Vector3D(R * cos(currentPhi), R * sin(currentPhi), currentZ))
          : Vector3D(R * cos(currentPhi), R * sin(currentPhi), currentZ);
      // the bin direction
      Vector3D binDirection = transform
          ? ((transform->linear())
             * Vector3D(cos(currentPhi), sin(currentPhi), 0.))
          : Vector3D(cos(currentPhi), sin(currentPhi), 0.);
      // push it in
      phiSystem.push_back(std::pair<SurfacePosition, Vector3D>(
          SurfacePosition(nullptr, binPosition), binDirection));
    }
    phizSystem.push_back(phiSystem);
  }

  // create and complete
  std::vector<SurfacePosition> sVector;
  // complete binning (even if bin-multipliers are in place, for x-check)
  completeBinning(surfaces, *arrayUtility, sVector, phizSystem);
  // create the surfaceArray
  auto sArray
    = std::make_unique<BinnedArrayXD<const Surface*>>(sVector, std::move(arrayUtility));
  // register the neighbours
  // @TODO -> registerNeighboursGrid(sArray->arrayObjectsOrdered(), false, true);
  // return the surface array
  return std::move(sArray);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Acts::Surface*>& surfaces,
    double                                   minR,
    double                                   maxR,
    double                                   minPhi,
    double                                   maxPhi,
    size_t                                   binsR,
    size_t                                   binsPhi,
    std::shared_ptr<Acts::Transform3D>       transform) const
{
  ACTS_DEBUG("Creating a SurfaceArray on a disc with grid in r x phi = "
            << binsR
            << " x "
            << binsPhi);

  // the z step and the phi step
  double                       phiStep = (maxPhi - minPhi) / (binsPhi - 1);
  std::vector<SurfacePosition> sVector;
  // reserve the right amount of binning
  sVector.reserve(binsPhi*binsR);
  // the array bin untility
  std::unique_ptr<BinUtility> arrayUtility = nullptr;
  // the very simple binning - nu bin multipliers are working
  if (binsR == 1 && binsPhi == surfaces.size()) {
    ACTS_DEBUG("Only one single ring in R is present - 1D surface array to be "
              "created.");
    // solve this at once
    arrayUtility = std::make_unique<BinUtility>(binsPhi, minPhi, maxPhi, closed, binPhi, transform);
    // fill the surfaces
    for (auto& surface : surfaces) {
      // fill it into the surface,position for further registration
      sVector.push_back(
          SurfacePosition(surface, surface->binningPosition(binPhi)));
    }
    // create the surface array
    auto sArray = std::make_unique< BinnedArrayXD<const Surface*> >
      (sVector, std::move(arrayUtility));
    // register the neighbours
    // copy to get const away,
    // - but take them from the surface array (and not the input vector) because
    // like this they are bin ordered
    // std::vector<const Surface*> arraySurfaces;
    // arraySurfaces.insert(arraySurfaces.begin(),
    //                      sArray->arrayObjects().begin(),
    //                      sArray->arrayObjects().end());
    //std::vector<std::vector<const Surface*>> arraySystem = {arraySurfaces};
    // prepared to run the neighbour registration now
    // @TODO -> registerNeighboursGrid(arraySystem, false, true);
    // now return
    return std::move(sArray);
  }
  // more complicated binning 2D
  double rStep = ((maxR - minR) / (binsR));
  // 1D or 2D binning - (depending on binsR)
  if (binsR == 1){
    // only phi-binning necessary
      arrayUtility = std::make_unique<BinUtility>(binsPhi, minPhi, maxPhi, closed, binPhi);
  } else {
    // r & phi binning necessary
    arrayUtility = std::make_unique<BinUtility>(binsR, minR, maxR, open, binR, transform);
    (*arrayUtility) += BinUtility(binsPhi, minPhi, maxPhi, closed, binPhi);
  }
  // prepare the surface system in r x phi
  std::vector<std::vector<std::pair<SurfacePosition, Vector3D>>> rphiSystem;
  // the bin direction 
  // @TODO should actually be +1, -1, gets important for intersection sorting method
  Vector3D binDirection(0., 0., 1);
  rphiSystem.reserve(binsPhi);
  // loop of iR and iPhi bins and fill the order positions
  for (size_t iPhi = 0; iPhi < binsPhi; ++iPhi) {
    // the current phi value
    double currentPhi = minPhi + (iPhi + 0.5) * phiStep;
    // the current phi row
    std::vector< std::pair<SurfacePosition, Vector3D> > rSystem;
    rSystem.reserve(binsR);
    for (size_t iR = 0; iR < binsR; ++iR) {
      // the current R value
      double currentR = minR + (iR + 0.5) * rStep;
      // the bin position
      Vector3D binPosition = transform
          ? ((*transform) * Vector3D(currentR * cos(currentPhi),
                                     currentR * sin(currentPhi),
                                     0.))
          : Vector3D(
                currentR * cos(currentPhi), currentR * sin(currentPhi), 0.);
      // push it in
      rSystem.push_back(std::pair<SurfacePosition, Vector3D>(
          SurfacePosition(nullptr, binPosition), binDirection));
    }
    rphiSystem.push_back(rSystem);
  }
  // create and complete
  completeBinning(surfaces, *arrayUtility, sVector, rphiSystem);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(sVector, std::move(arrayUtility));
  // register the neighbours
  // @TODO -> put it back
  // registerNeighboursGrid(sArray->arrayObjectsOrdered(), false, true);
  // return the surface array
  return std::move(sArray);
}

/** SurfaceArrayCreator interface method - create an array on a plane */
std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnPlane(
    const std::vector<const Acts::Surface*>& /*surfaces*/,
    double /*halflengthX*/,
    double /*halflengthY*/,
    size_t /*binsX*/,
    size_t /*binsY*/,
    std::shared_ptr<Acts::Transform3D> /*transform*/) const
{
  //!< @TODO implement - take from ATLAS complex TRT builder
  return nullptr;
}

void
Acts::SurfaceArrayCreator::completeBinning(
    const std::vector<const Surface*>&                  surfaces,
    const BinUtility&                                   binUtility,
    std::vector<SurfacePosition>&                       sVector,
    std::vector<std::vector<SurfacePositionDirection>>& binSystem) const
{
  ACTS_DEBUG("Complete binning by filling closest neighbour surfaces in "
             "potentially empty bins.");

  // get the number of bins
  size_t bins0 = binSystem.at(0).size();
  size_t bins1 = binSystem.size();

  ACTS_DEBUG("Prefilling a bin system with [ " << bins0 << " x " << bins1
              << " ] with " << surfaces.size() << " objects.");

  // prefill the easy ones
  for (auto& surface : surfaces) {
    // calculate the bin 
    // - we need to use binningPosition() for e.g.  
    size_t bin0 = 0; // binUtility.bin(surface->binningPosition(), 0);
    size_t bin1 = binUtility.dimensions() > 1 ?  binUtility.bin(surface->center(), 1) : 0;
    ACTS_VERBOSE("- estimated bin [ " << bin0 << " x " << bin1
                                     << " ] for surface "
                                     << surface);
    // fill it - flip if bins0 == 1, because this means phi-binning only
    if (bins0 == 1)
        binSystem.at(bin0).at(bin1).first.first = surface;
    else
        binSystem.at(bin1).at(bin0).first.first = surface;
  }

  size_t completedBins = 0;

  // now complete the system - register the neighbors if necessary
  for (size_t ibin1 = 0; ibin1 < bins1; ++ibin1) {
    for (size_t ibin0 = 0; ibin0 < bins0; ++ibin0) {
      // get the current surface
      const Surface* binSurface = binSystem.at(ibin1).at(ibin0).first.first;
      // we are done when we have a surface registerd
      if (binSurface) continue;
      // binPosition
      Vector3D binPosition = binSystem.at(ibin1).at(ibin0).first.second;
      double surfaceDist = 10e10;
      // counter
      ++completedBins;
      // brute force method
      // - find the closest surface
      // @TODO try to add intersection test
      for (auto& surface : surfaces) {
        // skip if not here
        if (!surface) continue;
        // recalculate distance
        double testDist = (binPosition - surface->center()).mag();
        if (testDist < surfaceDist) {
          binSystem.at(ibin1).at(ibin0).first.first = surface;
          surfaceDist                               = testDist;
        }
      }
    }
  }
  ACTS_DEBUG("Number of empty bins that were filled with neighbours: "
              << completedBins);

  // stream out the system for Binned Array creation & check for nullptr
  size_t emptyBins = 0;
  sVector.reserve(bins0 * bins1);
  for (auto& outer : binSystem)
    for (auto& inner : outer) {
      ACTS_VERBOSE("- bin [ " << binUtility.bin(inner.first.second, 0) << " x "
                             << (binUtility.dimensions() > 1 ?
                                 binUtility.bin(inner.first.second, 1) : 0 )
                             << " ] holds surface "
                             << inner.first.first);
      if (!inner.first.first) ++emptyBins;
      sVector.push_back(SurfacePosition(inner.first.first, inner.first.second));
    }
  // warning message in case you have empty bins
  if (emptyBins) {
    ACTS_WARNING(emptyBins << " emptyBins detected after bin completion!");
  } else
    ACTS_DEBUG("No emptyBins detected after bin completion!");
    
}
