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

  // create the (plain) binUtility - with the transform if given
  auto arrayUtility = std::make_unique<BinUtility>(binsPhi, minPhi, maxPhi, closed, binPhi, transform);
  (*arrayUtility) += BinUtility(binsZ, -halfZ, halfZ, open, binZ);

  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(binsZ, SurfaceVector( binsPhi, nullptr) ) );
  V3Matrix    v3Matrix(binsZ, V3Vector( binsPhi, Vector3D(0.,0.,0.)));

  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility->binningData();
  // create the position matrix first
  for (size_t iz = 0; iz < binsZ;  ++iz){
      // generate the z value
      double z = bdataSet[1].centerValue(iz);
      for (size_t iphi = 0; iphi < binsPhi; ++iphi){
          // generate the phi value
          double phi = bdataSet[0].centerValue(iphi);
          // fill the position
          v3Matrix[iz][iphi] = Vector3D(R*cos(phi),R*sin(phi),z);
      }
    }

  /// prefill the surfaces we have 
  for (auto& sf : surfaces){
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // get the bins and fill
    std::array<size_t,3> bTriple =
    arrayUtility->binTriple(bPosition);
    // and fill into the grid
    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }
  // complete the Binning
  completeBinning(*arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(sGrid, std::move(arrayUtility));
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

  std::unique_ptr<BinUtility> arrayUtility = nullptr;
  // 1D or 2D binning - (depending on binsR)
  if (binsR == 1){
    // only phi-binning necessary : make a 0D binning data
    BinningData r0Data(binR, minR, maxR);
    arrayUtility = std::make_unique<BinUtility>(r0Data, transform);
  } else {
    // r & phi binning necessary
    arrayUtility = std::make_unique<BinUtility>(binsR, minR, maxR, open, binR, transform);
  }
  /// add the phi binning 
  /// @TODO steer open/closed by phi values
  (*arrayUtility) += BinUtility(binsPhi, minPhi, maxPhi, closed, binPhi);

  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(binsPhi, SurfaceVector( binsR, nullptr) ) );
  V3Matrix    v3Matrix(binsPhi, V3Vector( binsR, Vector3D(0.,0.,0.)));
 
  // get the average z
  double z = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces){
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the z position
    z += bPosition.z();
    // get the bins and fill
    std::array<size_t,3> bTriple =
    arrayUtility->binTriple(bPosition);
    // and fill into the grid
    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }
  // average the z position
  z /= surfaces.size();
  
  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility->binningData();
  // create the position matrix first
  for (size_t iphi = 0; iphi < binsPhi; ++iphi){
      // generate the z value
      double phi = bdataSet[1].centerValue(iphi);
      for (size_t ir = 0; ir < binsR; ++ir){
        // generate the phi value
        double R = bdataSet[0].centerValue(ir);
        // fill the position
        v3Matrix[iphi][ir] = Vector3D(R*cos(phi),R*sin(phi),z);
    }
  }
  // complete the Binning
  completeBinning(*arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(sGrid, std::move(arrayUtility));
  // return the surface array
  return std::move(sArray);
}

/// SurfaceArrayCreator interface method - create an array on a plane 
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

/// @TODO implement nearest neighbour search 
/// - takes too long otherwise in initialization
void
Acts::SurfaceArrayCreator::completeBinning(const BinUtility&    binUtility, 
                                           const V3Matrix&      v3Matrix,
                                           const SurfaceVector& sVector,
                                           SurfaceGrid& sGrid) const
{
  ACTS_DEBUG("Complete binning by filling closest neighbour surfaces in "
             "potentially empty bins.");
  // make a copy of the surface grid
  size_t nSurfaces   = sVector.size();
  size_t nGridPoints = v3Matrix.size()*v3Matrix[0].size();
  //
  ACTS_DEBUG("- Object count : " << nSurfaces   << " number of surfaces");
  ACTS_DEBUG("- Surface grid : " << nGridPoints << " number of bins");
  ACTS_DEBUG("       to fill : " << nGridPoints-nSurfaces);
  // 
  size_t binCompleted = 0;
  //
  for (size_t io1 = 0; io1 < v3Matrix.size(); ++io1){
    for (size_t io0 = 0; io0 < v3Matrix[0].size(); ++io0){
      // screen output
      const Surface* sentry = sGrid[0][io1][io0];
        /// intersect 
        Vector3D sposition = v3Matrix[io1][io0];
        double minPath = 10e10;
        for (auto& sf : sVector){
          double testPath = (sposition-sf->binningPosition(binR)).mag();
          if (testPath < minPath){
              sGrid[0][io1][io0] = sf;
              sentry             = sf;
              minPath            = testPath;
          }
        }
        ACTS_VERBOSE("   -  bin [" << io0 << ", " << io1 << "] - surface is " << sentry);
        // increase the bin completion
        ++binCompleted;
      }
    }
    ACTS_DEBUG("       filled  : " << binCompleted);
    
}
