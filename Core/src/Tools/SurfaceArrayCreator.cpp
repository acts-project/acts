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
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//#include "ACTS/Utilities/Units.hpp"

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
  ACTS_DEBUG("Creating a SurfaceArray on a cylinder.");
  // create the 2D bin utility
  // create the (plain) binUtility - with the transform if given
  auto arrayUtility = std::make_unique<Acts::BinUtility>(createBinUtility(
      surfaces, binPhi, equidistant, binsPhi, minPhi, maxPhi, transform));
  (*arrayUtility)
      += createBinUtility(surfaces, binZ, equidistant, binsZ, -halfZ, halfZ);
  // prepare the surface matrix
  size_t      bins1 = arrayUtility->bins(1);
  size_t      bins0 = arrayUtility->bins(0);
  SurfaceGrid sGrid(1, SurfaceMatrix(bins1, SurfaceVector(bins0, nullptr)));
  V3Matrix    v3Matrix(bins1, V3Vector(bins0, Vector3D(0., 0., 0.)));
  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility->binningData();
  // create the position matrix first
  for (size_t iz = 0; iz < bins1; ++iz) {
    // generate the z value
    double z = bdataSet[1].centerValue(iz);
    for (size_t iphi = 0; iphi < bins0; ++iphi) {
      // generate the phi value
      double phi = bdataSet[0].centerValue(iphi);
      // fill the position
      v3Matrix[iz][iphi] = Vector3D(R * cos(phi), R * sin(phi), z);
    }
  }
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility->binTriple(bPosition);
    // and fill into the grid
    sGrid[bTriple.at(2)][bTriple.at(1)][bTriple.at(0)] = sf;
  }
  // complete the Binning @TODO switch on when we have a faster method for this
  completeBinning(*arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(
      sGrid, std::move(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
  // return the surface array
  return std::move(sArray);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Acts::Surface*>& surfaces,
    Acts::BinningType                        bTypePhi,
    Acts::BinningType                        bTypeZ,
    std::shared_ptr<Acts::Transform3D>       transform) const
{
  ACTS_DEBUG("Creating a SurfaceArray on a cylinder.");
  // create the 2D bin utility
  // create the (plain) binUtility - with the transform if given
  Acts::BinUtility arrayUtility;
  if (bTypePhi == equidistant)
    arrayUtility = createEquidistantBinUtility(surfaces, binPhi, transform);
  else
    arrayUtility = createArbitraryBinUtility(surfaces, binPhi, transform);
  if (bTypeZ == equidistant)
    arrayUtility += createEquidistantBinUtility(surfaces, binZ);
  else
    arrayUtility += createArbitraryBinUtility(surfaces, binZ);

  // get the number of bins
  size_t bins1 = arrayUtility.bins(1);
  size_t bins0 = arrayUtility.bins(0);
  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(bins1, SurfaceVector(bins0, nullptr)));
  V3Matrix    v3Matrix(bins1, V3Vector(bins0, Vector3D(0., 0., 0.)));
  // get the average r
  double R = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the r position
    R += bPosition.perp();
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility.binTriple(bPosition);
    // and fill into the grid
    sGrid[bTriple.at(2)][bTriple.at(1)][bTriple.at(0)] = sf;
  }
  // average the R position
  R /= surfaces.size();
  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility.binningData();
  // create the position matrix first
  for (size_t iz = 0; iz < bins1; ++iz) {
    // generate the z value
    double z = bdataSet[1].centerValue(iz);
    for (size_t iphi = 0; iphi < bins0; ++iphi) {
      // generate the phi value
      double phi = bdataSet[0].centerValue(iphi);
      // fill the position
      v3Matrix[iz][iphi] = Vector3D(R * cos(phi), R * sin(phi), z);
    }
  }
  // complete the Binning @TODO switch on when we have a faster method for this
  completeBinning(arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(
      sGrid, std::make_unique<Acts::BinUtility>(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
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
  ACTS_DEBUG("Creating a SurfaceArray on a disc.");

  auto arrayUtility = std::make_unique<Acts::BinUtility>(createBinUtility(
      surfaces, binR, equidistant, binsR, minR, maxR, transform));
  (*arrayUtility) += createBinUtility(
      surfaces, binPhi, equidistant, binsPhi, minPhi, maxPhi);

  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(binsPhi, SurfaceVector(binsR, nullptr)));
  V3Matrix    v3Matrix(binsPhi, V3Vector(binsR, Vector3D(0., 0., 0.)));

  // get the average z
  double z = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the z position
    z += bPosition.z();
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility->binTriple(bPosition);
    // and fill into the grid
    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }
  // average the z position
  z /= surfaces.size();

  ACTS_DEBUG("- z-position of disk estimated as " << z);

  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility->binningData();
  // create the position matrix first
  for (size_t iphi = 0; iphi < binsPhi; ++iphi) {
    // generate the z value
    double phi = bdataSet[1].centerValue(iphi);
    for (size_t ir = 0; ir < binsR; ++ir) {
      // generate the phi value
      double R = bdataSet[0].centerValue(ir);
      // fill the position
      v3Matrix[iphi][ir] = Vector3D(R * cos(phi), R * sin(phi), z);
    }
  }
  // complete the Binning
  completeBinning(*arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(
      sGrid, std::move(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
  // return the surface array
  return std::move(sArray);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Acts::Surface*>& surfaces,
    Acts::BinningType                        bTypeR,
    Acts::BinningType                        bTypePhi,
    std::shared_ptr<Acts::Transform3D>       transform) const
{
  ACTS_DEBUG("Creating a SurfaceArray on a disc.");
  Acts::BinUtility arrayUtility;
  if (bTypeR == equidistant)
    arrayUtility = createEquidistantBinUtility(surfaces, binR);
  else
    arrayUtility = createArbitraryBinUtility(surfaces, binR);
  if (bTypePhi == equidistant)
    arrayUtility += createEquidistantBinUtility(surfaces, binPhi, transform);
  else
    arrayUtility += createArbitraryBinUtility(surfaces, binPhi, transform);
  // get the number of bins
  size_t bins1 = arrayUtility.bins(1);
  size_t bins0 = arrayUtility.bins(0);
  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(bins1, SurfaceVector(bins0, nullptr)));
  V3Matrix    v3Matrix(bins1, V3Vector(bins0, Vector3D(0., 0., 0.)));
  // get the average z
  double z = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the z position
    z += bPosition.z();
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility.binTriple(bPosition);
    // and fill into the grid
    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }
  // average the z position
  z /= surfaces.size();

  ACTS_DEBUG("- z-position of disk estimated as " << z);

  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility.binningData();
  // create the position matrix first
  for (size_t iphi = 0; iphi < bins1; ++iphi) {
    // generate the z value
    double phi = bdataSet[1].centerValue(iphi);
    for (size_t ir = 0; ir < bins0; ++ir) {
      // generate the phi value
      double R = bdataSet[0].centerValue(ir);
      // fill the position
      v3Matrix[iphi][ir] = Vector3D(R * cos(phi), R * sin(phi), z);
    }
  }
  // complete the Binning
  completeBinning(arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  auto sArray = std::make_unique<BinnedArrayXD<const Surface*>>(
      sGrid, std::make_unique<Acts::BinUtility>(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
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

Acts::BinUtility
Acts::SurfaceArrayCreator::createArbitraryBinUtility(
    const std::vector<const Acts::Surface*>& surfaces,
    Acts::BinningValue                       bValue,
    std::shared_ptr<Acts::Transform3D>       transform) const
{
  // BinningOption is open for z and r, in case of phi binning reset later
  Acts::BinningOption bOption = Acts::open;
  // arbitrary binning - find out the binning positions and boundaries
  // of the surfaces in the specific direction
  std::vector<std::pair<float, float>> boundaries;
  // now loop through the surfaces and find out the needed information
  std::vector<float> bValues;
  if (bValue == Acts::binPhi) {
    // set the BinningOption
    bOption = closed;
    // copy the surface vector to a non constant vector
    std::vector<const Acts::Surface*> surf(surfaces);
    // the key surfaces - placed in different bins in the given binning
    // direction
    std::vector<const Acts::Surface*> keys;
    // sort first in phi
    std::stable_sort(surf.begin(),
                     surf.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().phi() < b->center().phi());
                     });
    // fill the key surfaces at the different phi positions
    std::unique_copy(begin(surf),
                     end(surf),
                     back_inserter(keys),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (fabs(a->center().phi() - b->center().phi())
                               < 10e-12);
                     });
    // loop through the key surfaces and access the internal boundaries
    for (auto& surface : keys) {
      // get the bounds
      const Acts::PlanarBounds* planarBounds
          = dynamic_cast<const Acts::PlanarBounds*>(&(surface->bounds()));
      if (!planarBounds)
        ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                   "other bounds yet! ");
      // get the vertices
      std::vector<Acts::Vector2D> vertices = planarBounds->vertices();
      if (vertices.empty()) ACTS_ERROR("Vertices of current surface empty!");

      // Phi binning
      Acts::Vector3D globPos1(0., 0., 0.);
      // @TODO think of which position to take (or mixture) for a more
      // general case (e.g. what if modules are twisted?)
      surface->localToGlobal(
          vertices.at(2), Acts::Vector3D(0., 0., 0.), globPos1);
      bValues.push_back(globPos1.phi());
      if (surface == *(keys.end() - 1)) {
        Acts::Vector3D globPos2(0., 0., 0.);
        surface->localToGlobal(
            vertices.at(1), Acts::Vector3D(0., 0., 0.), globPos1);
        bValues.push_back(globPos2.phi());
      }
    }
  } else if (bValue == Acts::binZ) {
    // copy the surface vector to a non constant vector
    std::vector<const Acts::Surface*> surf(surfaces);
    // the key surfaces - placed in different bins in the given binning
    // direction
    std::vector<const Acts::Surface*> keys;
    // sort first in z
    std::stable_sort(surf.begin(),
                     surf.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().z() < b->center().z());
                     });
    // fill the key surfaces at the different z positions
    std::unique_copy(begin(surf),
                     end(surf),
                     back_inserter(keys),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().z() == b->center().z());
                     });
    // go through key surfaces
    for (auto& surface : keys) {
      // get the bounds
      const Acts::PlanarBounds* planarBounds
          = dynamic_cast<const Acts::PlanarBounds*>(&(surface->bounds()));
      if (!planarBounds)
        ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                   "other bounds yet! ");
      // get the vertices
      std::vector<Acts::Vector2D> vertices = planarBounds->vertices();
      if (vertices.empty()) ACTS_ERROR("Vertices of current surface empty!");
      // Z binning
      // accessing any entry to access the second coordinate is sufficient
      // make local coordinate global
      // @TODO implement general approach - what if the modules are twisted
      // and it is not the second coordinate?
      double length   = fabs(vertices.front().y());
      double position = 0;
      // surfaces can be binned in z or r
      position = surface->center().z();
      bValues.push_back(position - length);
      if (surface == *(keys.end() - 1)) {
        bValues.push_back(position + length);
      }
    }
  } else {
    // copy the surface vector to a non constant vector
    std::vector<const Acts::Surface*> surf(surfaces);
    // the key surfaces - placed in different bins in the given binning
    // direction
    std::vector<const Acts::Surface*> keys;
    // sort first in r
    std::stable_sort(surf.begin(),
                     surf.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().perp() < b->center().perp());
                     });
    // fill the key surfaces at the different r positions
    std::unique_copy(begin(surf),
                     end(surf),
                     back_inserter(keys),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (fabs(a->center().perp() - b->center().perp())
                               < 10e-6);
                     });
    // go through key surfaces
    for (auto& surface : keys) {
      // get the bounds
      const Acts::PlanarBounds* planarBounds
          = dynamic_cast<const Acts::PlanarBounds*>(&(surface->bounds()));
      if (!planarBounds)
        ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                   "other bounds yet! ");
      // get the vertices
      std::vector<Acts::Vector2D> vertices = planarBounds->vertices();
      if (vertices.empty()) ACTS_ERROR("Vertices of current surface empty!");
      // R binning
      // accessing any entry to access the second coordinate is
      // sufficient
      // make local coordinate global
      // @TODO implement general approach - what if the modules are twisted
      // and it is not the second coordinate?
      double length   = fabs(vertices.front().y());
      double position = 0;
      // surfaces can be binned in z or r
      position = surface->center().perp();
      bValues.push_back(position - length);
      // eliminate double values and sort values
      if (surface == *(keys.end() - 1)) {
        bValues.push_back(position + length);
      }
    }
  }
  std::sort(bValues.begin(), bValues.end());
  ACTS_DEBUG("Create BinUtility for BinnedSurfaceArray with arbitrary "
             "BinningType");
  ACTS_DEBUG("	BinningValue: " << bValue);
  ACTS_DEBUG("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
             "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_DEBUG("	Number of bins: " << bValues.size());
  // create the BinUtility
  return (Acts::BinUtility(bValues, bOption, bValue));
}

Acts::BinUtility
Acts::SurfaceArrayCreator::createEquidistantBinUtility(
    const std::vector<const Acts::Surface*>& surfaces,
    Acts::BinningValue                       bValue,
    std::shared_ptr<Acts::Transform3D>       transform) const
{  // check the binning type first
  double minimum = 0.;
  double maximum = 0.;
  // binning option is open for z and r, in case of phi binning reset later
  Acts::BinningOption bOption = Acts::open;
  // copy the surface vector to a non constant vector
  std::vector<const Acts::Surface*> surf(surfaces);
  // the key surfaces - placed in different bins in the given binning
  // direction
  std::vector<const Acts::Surface*> keys;
  // now check the binning value
  if (bValue == Acts::binPhi) {
    // Phi binning
    // set the binning option for phi
    bOption = Acts::closed;
    // sort first in phi
    std::stable_sort(surf.begin(),
                     surf.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().phi() < b->center().phi());
                     });
    // fill the key surfaces at the different phi positions
    std::unique_copy(begin(surf),
                     end(surf),
                     back_inserter(keys),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (fabs(a->center().phi() - b->center().phi())
                               < 10e-12);
                     });
    // set the minimum and maximum
    Acts::Vector3D globPos1(0., 0., 0.);
    Acts::Vector3D globPos2(0., 0., 0.);
    // get the first and the last surface in phi
    const Acts::Surface* frontSurface = keys.front();
    const Acts::Surface* backSurface  = keys.back();
    // access the bounds of the first surface
    const Acts::PlanarBounds* frontBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(frontSurface->bounds()));
    if (!frontBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the boundaries
    std::vector<Acts::Vector2D> frontVertices = frontBounds->vertices();
    // access boundaries of last surface
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the boundaries
    std::vector<Acts::Vector2D> backVertices = backBounds->vertices();
    // @TODO think of which position to take (or mixture) for a more
    // general
    // case (e.g. what if modules are twisted?)
    // make point global
    frontSurface->localToGlobal(
        frontVertices.at(3), Acts::Vector3D(0., 0., 0.), globPos1);
    backSurface->localToGlobal(
        backVertices.at(0), Acts::Vector3D(0., 0., 0.), globPos2);
    // get minimum and maximum
    minimum = globPos1.phi();
    maximum = globPos2.phi();
    // check
    if (frontSurface->center().phi() > maximum) maximum = M_PI;
    if (backSurface->center().phi() < minimum) minimum  = -M_PI;

  } else if (bValue == Acts::binZ) {
    // Z binning
    // sort first in z
    std::stable_sort(surf.begin(),
                     surf.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().z() < b->center().z());
                     });
    // fill the key surfaces at the different z positions
    std::unique_copy(begin(surf),
                     end(surf),
                     back_inserter(keys),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().z() == b->center().z());
                     });
    // set the minimum and maximum
    // get the first and the last surface in z
    const Acts::Surface* frontSurface = keys.front();
    const Acts::Surface* backSurface  = keys.back();
    // get the boundaries of the first surface
    const Acts::PlanarBounds* frontBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(frontSurface->bounds()));
    if (!frontBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the vertices
    std::vector<Acts::Vector2D> frontVertices = frontBounds->vertices();
    // get boundaries of last surface
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the vertices
    std::vector<Acts::Vector2D> backVertices = backBounds->vertices();
    // set the minimum and maximum
    minimum = frontSurface->center().z() - fabs(frontVertices.front().y());
    maximum = backSurface->center().z() + fabs(backVertices.front().y());
  } else {
    // R binning
    // sort first in r
    std::stable_sort(surf.begin(),
                     surf.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->center().perp() < b->center().perp());
                     });
    // fill the key surfaces at the different r positions
    std::unique_copy(begin(surf),
                     end(surf),
                     back_inserter(keys),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (fabs(a->center().perp() - b->center().perp())
                               < 10e-6);
                     });
    // set the minimum and maximum
    // get the first and the last surface in phi
    const Acts::Surface* frontSurface = keys.front();
    const Acts::Surface* backSurface  = keys.back();
    // get the boundaries of the first surface
    const Acts::PlanarBounds* frontBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(keys.front()->bounds()));
    if (!frontBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the vertices
    std::vector<Acts::Vector2D> frontVertices = frontBounds->vertices();
    // get the boundaries of the last surface
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(keys.back()->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get vertices
    std::vector<Acts::Vector2D> backVertices = backBounds->vertices();
    // calculate minimum and maximum
    minimum = frontSurface->center().perp() - fabs(frontVertices.front().y());
    maximum = backSurface->center().perp() + fabs(backVertices.front().y());
  }
  // assign the bin size
  double binNumber = keys.size();
  ACTS_DEBUG("Create BinUtility for BinnedSurfaceArray with equidistant1 "
             "BinningType");
  ACTS_DEBUG("	BinningValue: " << bValue);
  ACTS_DEBUG("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
             "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_DEBUG("	Number of bins: " << binNumber);
  ACTS_DEBUG("	(Min/Max) = (" << minimum << "/" << maximum << ")");
  return (Acts::BinUtility(binNumber, minimum, maximum, bOption, bValue));
}

Acts::BinUtility
Acts::SurfaceArrayCreator::createBinUtility(
    const std::vector<const Acts::Surface*>& surfaces,
    Acts::BinningValue                       bValue,
    Acts::BinningType                        bType,
    size_t                                   bins,
    double                                   min,
    double                                   max,
    std::shared_ptr<Acts::Transform3D>       transform) const
{
  // check first
  if (surfaces.empty())
    ACTS_ERROR("No surfaces given - can not create BinUtility for "
               "BinnedSurfaceArray!");
  // introduce BinUtility to hand back
  std::unique_ptr<Acts::BinUtility> binUtility = nullptr;
  Acts::BinningOption               bOption    = open;
  // all the information already given
  if (bValue == Acts::binPhi) bOption = closed;
  ACTS_DEBUG("Create BinUtility for BinnedSurfaceArray with equidistant "
             "BinningType");
  ACTS_DEBUG("	BinningValue: " << bValue);
  ACTS_DEBUG("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
             "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_DEBUG("	Number of bins: " << bins);
  ACTS_DEBUG("	(Min/Max) = (" << min << "/" << max << ")");
  // create the BinUtility
  return (Acts::BinUtility(bins, min, max, bOption, bValue, transform));
}

/// Register the neigbourhood
void
Acts::SurfaceArrayCreator::registerNeighbourHood(
    const SurfaceArray& sArray) const
{
  ACTS_DEBUG("Register neighbours to the elements.");
  // get the grid first
  auto objectGrid = sArray.objectGrid();
  // statistics
  size_t neighboursSet = 0;
  // then go through, will respect a non-regular matrix
  size_t io2 = 0;
  for (auto& v210 : objectGrid) {
    size_t io1 = 0;
    for (auto& v10 : v210) {
      size_t io0 = 0;
      for (auto& bSurface : v10) {
        // get the member of this bin
        if (bSurface && bSurface->associatedDetectorElement()) {
          // get the bin detector element (for readability)
          auto bElement = bSurface->associatedDetectorElement();
          // get all the surfaces clustering around
          auto objectCluster = sArray.objectCluster({{io0, io1, io2}});
          // now loop and fill
          for (auto& nSurface : objectCluster) {
            // create the detector element vector with nullptr protection
            std::vector<const DetectorElementBase*> neighbourElements;
            if (nSurface && nSurface != bSurface
                && nSurface->associatedDetectorElement()) {
              // register it to the vector
              neighbourElements.push_back(
                  nSurface->associatedDetectorElement());
              // increase the counter
              ++neighboursSet;
            }
            // now register the neighbours
            bElement->registerNeighbours(neighbourElements);
          }
        }
        ++io0;
      }
      ++io1;
    }
    ++io2;
  }
  ACTS_DEBUG("Neighbours set for this layer: " << neighboursSet);
}

/// @TODO implement nearest neighbour search - this is brute force attack
/// - takes too long otherwise in initialization
void
Acts::SurfaceArrayCreator::completeBinning(const BinUtility&    binUtility,
                                           const V3Matrix&      v3Matrix,
                                           const SurfaceVector& sVector,
                                           SurfaceGrid&         sGrid) const
{
  ACTS_DEBUG("Complete binning by filling closest neighbour surfaces into "
             "empty bins.");
  // make a copy of the surface grid
  size_t nSurfaces   = sVector.size();
  size_t nGridPoints = v3Matrix.size() * v3Matrix[0].size();
  // bail out as there is nothing to do
  if (nGridPoints == nSurfaces) {
    ACTS_VERBOSE(" - Nothing to do, no empty bins present.");
    return;
  }
  // VERBOSE screen output
  ACTS_VERBOSE("- Object count : " << nSurfaces << " number of surfaces");
  ACTS_VERBOSE("- Surface grid : " << nGridPoints << " number of bins");
  ACTS_VERBOSE("       to fill : " << nGridPoints - nSurfaces);

  size_t binCompleted = 0;
  //
  for (size_t io1 = 0; io1 < v3Matrix.size(); ++io1) {
    for (size_t io0 = 0; io0 < v3Matrix[0].size(); ++io0) {
      /// intersect
      Vector3D sposition = v3Matrix[io1][io0];
      double   minPath   = 10e10;
      for (auto& sf : sVector) {
        double testPath = (sposition - sf->binningPosition(binR)).mag();
        if (testPath < minPath) {
          sGrid[0][io1][io0] = sf;
          minPath            = testPath;
        }
      }
      // increase the bin completion
      ++binCompleted;
    }
  }

  ACTS_DEBUG("       filled  : " << binCompleted);
}

std::vector<float>
Acts::SurfaceArrayCreator::createBinValues(
    std::vector<std::pair<float, float>> old) const
{
  sort(old.begin(),
       old.end(),
       [](std::pair<float, float> ap, std::pair<float, float> bp) {
         float a = (ap.second + ap.first) * 0.5;
         float b = (bp.second + bp.first) * 0.5;
         return a < b;
       });
  std::vector<float> newValues;
  std::pair<float, float> current;
  std::pair<float, float> next;
  // eliminate doubles
  std::vector<std::pair<float, float>> oldKeys;
  std::unique_copy(begin(old),
                   end(old),
                   back_inserter(oldKeys),
                   [](std::pair<float, float>& a, std::pair<float, float>& b) {
                     return a == b;
                   });
  for (std::vector<std::pair<float, float>>::iterator it = oldKeys.begin();
       it != (oldKeys.end() - 1);
       ++it) {
    current = *it;
    next    = *(it + 1);
    if (it == oldKeys.begin()) newValues.push_back(current.first);
    if (next.first <= current.second)
      newValues.push_back((current.second + next.first) * 0.5);
    else
      newValues.push_back(current.second);
    if (it == (oldKeys.end() - 2)) newValues.push_back(next.second);
  }
  return (newValues);
}
