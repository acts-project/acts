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
#include <algorithm>
#include <cmath>
#include <functional>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Surface*>& surfaces,
    size_t                             binsPhi,
    size_t                             binsZ,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << binsPhi << " x " << binsZ << " = "
                                      << binsPhi * binsZ
                                      << " bins.");

  // is definitely equidistant
  Acts::BinUtility arrayUtility = createEquidistantBinUtility(
      surfaces, binPhi, protoLayer, transform, binsPhi);
  arrayUtility += createEquidistantBinUtility(
      surfaces, binZ, protoLayer, nullptr, binsZ);

  double R = protoLayer.maxR - protoLayer.minR;

  // prepare the surface matrix
  size_t      bins1 = arrayUtility.bins(1);
  size_t      bins0 = arrayUtility.bins(0);
  SurfaceGrid sGrid(1, SurfaceMatrix(bins1, SurfaceVector(bins0, nullptr)));
  V3Matrix    v3Matrix(bins1, V3Vector(bins0, Vector3D(0., 0., 0.)));
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

  size_t nOverwrite = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility.binTriple(bPosition);
    // and fill into the grid

    if (sGrid[bTriple[2]][bTriple[1]][bTriple[0]] != nullptr) {
      nOverwrite++;
    }

    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }

  if (nOverwrite > 0) {
    ACTS_WARNING("- " << nOverwrite
                      << " bins were overwritten during bin filling");
  }

  // complete the Binning @todo switch on when we have a faster method for this
  completeBinning(arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  std::unique_ptr<Acts::SurfaceArray> sArray
      = std::make_unique<BinnedArrayXD<const Surface*>>(
          sGrid, std::make_unique<const Acts::BinUtility>(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
  // return the surface array
  return sArray;
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Surface*>& surfaces,
    BinningType                        bTypePhi,
    BinningType                        bTypeZ,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  // create the 2D bin utility
  // create the (plain) binUtility - with the transform if given
  Acts::BinUtility arrayUtility;
  if (bTypePhi == equidistant)
    arrayUtility
        = createEquidistantBinUtility(surfaces, binPhi, protoLayer, transform);
  else
    arrayUtility = createArbitraryBinUtility(surfaces, binPhi, transform);
  if (bTypeZ == equidistant)
    arrayUtility
        += createEquidistantBinUtility(surfaces, binZ, protoLayer, nullptr);
  else
    arrayUtility += createArbitraryBinUtility(surfaces, binZ, nullptr);

  // get the number of bins
  size_t bins1 = arrayUtility.bins(1);
  size_t bins0 = arrayUtility.bins(0);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");

  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(bins1, SurfaceVector(bins0, nullptr)));
  V3Matrix    v3Matrix(bins1, V3Vector(bins0, Vector3D(0., 0., 0.)));
  // get the average r
  double R          = 0;
  size_t nOverwrite = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the r position
    R += bPosition.perp();
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility.binTriple(bPosition);
    // and fill into the grid
    if (sGrid[bTriple[2]][bTriple[1]][bTriple[0]] != nullptr) {
      nOverwrite++;
    }
    sGrid[bTriple.at(2)][bTriple.at(1)][bTriple.at(0)] = sf;
  }

  if (nOverwrite > 0) {
    ACTS_WARNING("- " << nOverwrite
                      << " bins were overwritten during bin filling");
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
  // complete the Binning
  // @TODO switch on when we have a faster method for this
  completeBinning(arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  std::unique_ptr<Acts::SurfaceArray> sArray
      = std::make_unique<BinnedArrayXD<const Surface*>>(
          sGrid, std::make_unique<const Acts::BinUtility>(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
  // return the surface array
  return sArray;
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Surface*>& surfaces,
    size_t                             binsR,
    size_t                             binsPhi,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  // is definitely equidistant
  Acts::BinUtility arrayUtility
      = createEquidistantBinUtility(surfaces, binR, protoLayer, nullptr, binsR);
  arrayUtility += createEquidistantBinUtility(
      surfaces, binPhi, protoLayer, transform, binsPhi);

  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(binsPhi, SurfaceVector(binsR, nullptr)));
  V3Matrix    v3Matrix(binsPhi, V3Vector(binsR, Vector3D(0., 0., 0.)));

  // get the average z
  double z          = 0;
  size_t nOverwrite = 0;

  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the z position
    z += bPosition.z();
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility.binTriple(bPosition);
    // and fill into the grid
    if (sGrid[bTriple[2]][bTriple[1]][bTriple[0]] != nullptr) {
      nOverwrite++;
    }

    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }
  // average the z position
  z /= surfaces.size();

  if (nOverwrite > 0) {
    ACTS_WARNING("- " << nOverwrite
                      << " bins were overwritten during bin filling");
  }

  ACTS_VERBOSE("- z-position of disk estimated as " << z);

  // get access to the binning data
  const std::vector<BinningData>& bdataSet = arrayUtility.binningData();
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
  completeBinning(arrayUtility, v3Matrix, surfaces, sGrid);
  // create the surfaceArray
  std::unique_ptr<Acts::SurfaceArray> sArray
      = std::make_unique<BinnedArrayXD<const Surface*>>(
          sGrid, std::make_unique<const BinUtility>(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
  // return the surface array
  return sArray;
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Surface*>& surfaces,
    BinningType                        bTypeR,
    BinningType                        bTypePhi,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");
  Acts::BinUtility arrayUtility;
  if (bTypeR == equidistant)
    arrayUtility
        = createEquidistantBinUtility(surfaces, binR, protoLayer, nullptr);
  else
    arrayUtility = createArbitraryBinUtility(surfaces, binR, nullptr);
  if (bTypePhi == equidistant)
    arrayUtility
        += createEquidistantBinUtility(surfaces, binPhi, protoLayer, transform);
  else
    arrayUtility += createArbitraryBinUtility(surfaces, binPhi, transform);
  // get the number of bins
  size_t bins1 = arrayUtility.bins(1);
  size_t bins0 = arrayUtility.bins(0);
  // prepare the surface matrix
  SurfaceGrid sGrid(1, SurfaceMatrix(bins1, SurfaceVector(bins0, nullptr)));
  V3Matrix    v3Matrix(bins1, V3Vector(bins0, Vector3D(0., 0., 0.)));
  // get the average z
  double z          = 0;
  size_t nOverwrite = 0;
  /// prefill the surfaces we have
  for (auto& sf : surfaces) {
    /// get the binning position
    Vector3D bPosition = sf->binningPosition(binR);
    // record the z position
    z += bPosition.z();
    // get the bins and fill
    std::array<size_t, 3> bTriple = arrayUtility.binTriple(bPosition);
    // and fill into the grid
    if (sGrid[bTriple[2]][bTriple[1]][bTriple[0]] != nullptr) {
      nOverwrite++;
    }

    sGrid[bTriple[2]][bTriple[1]][bTriple[0]] = sf;
  }
  // average the z position
  z /= surfaces.size();

  if (nOverwrite > 0) {
    ACTS_WARNING("- " << nOverwrite
                      << " bins were overwritten during bin filling");
  }

  ACTS_VERBOSE("- z-position of disk estimated as " << z);

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
  std::unique_ptr<Acts::SurfaceArray> sArray
      = std::make_unique<BinnedArrayXD<const Surface*>>(
          sGrid, std::make_unique<const Acts::BinUtility>(arrayUtility));
  // define neigbourhood
  registerNeighbourHood(*sArray);
  // return the surface array
  return sArray;
}

/// SurfaceArrayCreator interface method - create an array on a plane
std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnPlane(
    const std::vector<const Surface*>& /*surfaces*/,
    double /*halflengthX*/,
    double /*halflengthY*/,
    size_t /*binsX*/,
    size_t /*binsY*/,
    std::shared_ptr<const Transform3D> /*transform*/) const
{
  //!< @todo implement - take from ATLAS complex TRT builder
  return nullptr;
}

size_t
Acts::SurfaceArrayCreator::determineBinCount(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue) const
{

  std::function<bool(const Acts::Surface*, const Acts::Surface*)> sorter;

  // bind matcher with binning type
  auto matcher = std::bind(m_cfg.surfaceMatcher,
                           bValue,
                           std::placeholders::_1,
                           std::placeholders::_2);

  if (bValue == Acts::binPhi) {
    sorter = [](const Acts::Surface* a, const Acts::Surface* b) {
      return (a->center().phi() < b->center().phi());
    };
  } else if (bValue == Acts::binZ) {
    // Z binning
    // sort first in z
    sorter = [](const Acts::Surface* a, const Acts::Surface* b) {
      return (a->center().z() < b->center().z());
    };
  } else {
    // R binning
    // sort first in r
    sorter = [](const Acts::Surface* a, const Acts::Surface* b) {
      return (a->center().perp() < b->center().perp());
    };
  }

  std::vector<const Acts::Surface*> surf(surfaces);
  std::stable_sort(surf.begin(), surf.end(), sorter);
  // fill the key surfaces
  std::vector<const Acts::Surface*> keys;
  std::unique_copy(begin(surf), end(surf), back_inserter(keys), matcher);
  return keys.size();
}

Acts::BinUtility
Acts::SurfaceArrayCreator::createArbitraryBinUtility(
    const std::vector<const Surface*>& surfaces,
    Acts::BinningValue                 bValue,
    std::shared_ptr<const Transform3D> transform) const
{
  if (!surfaces.size())
    throw std::logic_error(
        "No surfaces handed over for creating arbitrary bin utility!");
  // BinningOption is open for z and r, in case of phi binning reset later
  Acts::BinningOption bOption = Acts::open;
  // the vector with the binning Values (boundaries for each bin)

  // bind matcher with binning type
  auto matcher = std::bind(m_cfg.surfaceMatcher,
                           bValue,
                           std::placeholders::_1,
                           std::placeholders::_2);

  std::vector<float> bValues;
  if (bValue == Acts::binPhi) {
    // set the BinningOption closed for binning in phi
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
    std::unique_copy(begin(surf), end(surf), back_inserter(keys), matcher);

    // the phi-center position of the previous surface
    bool phiCorrected = false;

    // get minimum and maximum
    // the first boundary (minimum) and the last boundary (maximum) need to
    // be calculated separately with the vertices of the first and last
    // surface in the binning direction

    // get the bounds of the first surfaces
    const Acts::Surface*      frontSurface = keys.front();
    const Acts::PlanarBounds* frontBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(frontSurface->bounds()));
    if (!frontBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> frontVertices
        = makeGlobalVertices(*frontSurface, frontBounds->vertices());
    double minBValue = std::min_element(frontVertices.begin(),
                                        frontVertices.end(),
                                        [](const Acts::Vector3D& a,
                                           const Acts::Vector3D& b) {
                                          return a.phi() < b.phi();
                                        })
                           ->phi();
    // phi correction
    if (frontSurface->center().phi() < minBValue) {
      bValues.push_back(-M_PI);
      bValues.push_back(M_PI);
      phiCorrected = true;
    }
    bValues.push_back(minBValue);

    // get the bounds of the last surfaces
    const Acts::Surface*      backSurface = keys.back();
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> backVertices
        = makeGlobalVertices(*backSurface, backBounds->vertices());
    double maxBValue = std::max_element(backVertices.begin(),
                                        backVertices.end(),
                                        [](const Acts::Vector3D& a,
                                           const Acts::Vector3D& b) {
                                          return a.phi() < b.phi();
                                        })
                           ->phi();
    // phi correction
    if (backSurface->center().phi() > maxBValue && !phiCorrected) {
      bValues.push_back(M_PI);
      bValues.push_back(-M_PI);
    }
    bValues.push_back(maxBValue);

    double previous = frontSurface->center().phi();
    // go through key surfaces
    for (auto surface = keys.begin(); surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(0.5 * (previous + (*surface)->center().phi()));
      previous = (*surface)->center().phi();
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
    std::unique_copy(begin(surf), end(surf), back_inserter(keys), matcher);

    // get minimum and maximum
    // the first boundary (minimum) and the last boundary (maximum) need to
    // be calculated separately with the vertices of the first and last
    // surface in the binning direction

    // get the bounds of the first surfaces
    const Acts::Surface*      frontSurface = keys.front();
    const Acts::PlanarBounds* frontBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(frontSurface->bounds()));
    if (!frontBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> frontVertices
        = makeGlobalVertices(*frontSurface, frontBounds->vertices());
    double minBValue = std::min_element(frontVertices.begin(),
                                        frontVertices.end(),
                                        [](const Acts::Vector3D& a,
                                           const Acts::Vector3D& b) {
                                          return a.z() < b.z();
                                        })
                           ->z();

    bValues.push_back(minBValue);

    // get the bounds of the last surfaces
    const Acts::Surface*      backSurface = keys.back();
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> backVertices
        = makeGlobalVertices(*backSurface, backBounds->vertices());
    double maxBValue = std::max_element(backVertices.begin(),
                                        backVertices.end(),
                                        [](const Acts::Vector3D& a,
                                           const Acts::Vector3D& b) {
                                          return a.z() < b.z();
                                        })
                           ->z();

    bValues.push_back(maxBValue);

    // the z-center position of the previous surface
    double previous = frontSurface->center().z();
    // go through key surfaces
    for (auto surface = keys.begin(); surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(0.5 * (previous + (*surface)->center().z()));
      previous = (*surface)->center().z();
    }
  } else {  // binR
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
    std::unique_copy(begin(surf), end(surf), back_inserter(keys), matcher);

    // get minimum and maximum
    // the first boundary (minimum) and the last boundary (maximum) need to
    // be calculated separately with the vertices of the first and last
    // surface in the binning direction

    // get the bounds of the first surfaces
    const Acts::Surface*      frontSurface = keys.front();
    const Acts::PlanarBounds* frontBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(frontSurface->bounds()));
    if (!frontBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> frontVertices
        = makeGlobalVertices(*frontSurface, frontBounds->vertices());
    double minBValue = std::min_element(frontVertices.begin(),
                                        frontVertices.end(),
                                        [](const Acts::Vector3D& a,
                                           const Acts::Vector3D& b) {
                                          return (a.perp() < b.perp());
                                        })
                           ->perp();

    bValues.push_back(minBValue);

    // get the bounds of the last surfaces
    const Acts::Surface*      backSurface = keys.back();
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> backVertices
        = makeGlobalVertices(*backSurface, backBounds->vertices());
    double maxBValue = std::max_element(backVertices.begin(),
                                        backVertices.end(),
                                        [](const Acts::Vector3D& a,
                                           const Acts::Vector3D& b) {
                                          return (a.perp() < b.perp());
                                        })
                           ->perp();

    bValues.push_back(maxBValue);

    // the r-center position of the previous surface
    double previous = frontSurface->center().perp();

    // go through key surfaces
    for (auto surface = keys.begin(); surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(0.5 * (previous + (*surface)->center().perp()));
      previous = (*surface)->center().perp();
    }
  }
  std::sort(bValues.begin(), bValues.end());
  ACTS_VERBOSE("Create BinUtility for BinnedSurfaceArray with arbitrary "
               "BinningType");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
               "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << bValues.size());
  // create the BinUtility
  return (Acts::BinUtility(bValues, bOption, bValue));
}

Acts::BinUtility
Acts::SurfaceArrayCreator::createEquidistantBinUtility(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue,
    ProtoLayer                         protoLayer,
    std::shared_ptr<const Transform3D> transform,
    size_t                             nBins) const
{
  if (!surfaces.size())
    throw std::logic_error(
        "No surfaces handed over for creating equidistant bin utility!");
  // check the binning type first
  double minimum = 0.;
  double maximum = 0.;
  // binning option is open for z and r, in case of phi binning reset later
  Acts::BinningOption bOption = Acts::open;
  // copy the surface vector to a non constant vector
  std::vector<const Acts::Surface*> surf(surfaces);
  // the key surfaces - placed in different bins in the given binning
  // direction
  std::vector<const Acts::Surface*> keys;

  size_t binNumber;
  if (nBins == 0) {
    // determine bin count
    binNumber = determineBinCount(surfaces, bValue);
  } else {
    // use bin count
    binNumber = nBins;
  }

  // bind matcher with binning type
  auto matcher = std::bind(m_cfg.surfaceMatcher,
                           bValue,
                           std::placeholders::_1,
                           std::placeholders::_2);

  // now check the binning value
  if (bValue == Acts::binPhi) {
    // Phi binning
    // set the binning option for phi
    // sort first in phi
    const Acts::Surface* maxElem
        = *std::max_element(surfaces.begin(),
                            surfaces.end(),
                            [](const Acts::Surface* a, const Acts::Surface* b) {
                              return a->center().phi() < b->center().phi();
                            });

    // fill the key surfaces at the different phi positions
    std::unique_copy(begin(surf), end(surf), back_inserter(keys), matcher);

    // multiple surfaces, we bin from -pi to pi closed
    if (keys.size() > 1) {
      bOption = Acts::closed;

      minimum = -M_PI;
      maximum = M_PI;

      // double step = 2 * M_PI / keys.size();
      double step = 2 * M_PI / binNumber;
      // rotate to max phi module plus one half step
      // this should make sure that phi wrapping at +- pi
      // never falls on a module center
      double max   = maxElem->center().phi();
      double angle = M_PI - (max + 0.5 * step);

      if (transform == nullptr) {
        transform = std::make_shared<const Transform3D>(
            AngleAxis3D(angle, Eigen::Vector3d::UnitZ()));
      } else {
        transform = std::make_shared<const Transform3D>(
            (*transform) * AngleAxis3D(angle, Eigen::Vector3d::UnitZ()));
      }

    } else {
      // only one surface

      // binning is open if only one element in phi
      bOption = Acts::open;

      minimum = protoLayer.minPhi;
      maximum = protoLayer.maxPhi;

      // we do not need a transform in this case
    }

  } else if (bValue == Acts::binZ) {
    // Z binning

    // just use maximum and minimum of all surfaces
    // we do not need key surfaces here
    maximum = protoLayer.maxZ;
    minimum = protoLayer.minZ;

  } else {
    // R binning

    // just use maximum and minimum of all surfaces
    // we do not need key surfaces here
    maximum = protoLayer.maxR;
    minimum = protoLayer.minR;
  }
  // assign the bin size
  ACTS_VERBOSE("Create BinUtility for BinnedSurfaceArray with equidistant "
               "BinningType");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
               "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << binNumber);
  ACTS_VERBOSE("	(Min/Max) = (" << minimum << "/" << maximum << ")");

  return (Acts::BinUtility(
      binNumber, minimum, maximum, bOption, bValue, transform));
}

Acts::BinUtility
Acts::SurfaceArrayCreator::createBinUtility(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue,
    BinningType                        bType,
    size_t                             bins,
    double                             min,
    double                             max,
    std::shared_ptr<const Transform3D> transform) const
{
  // check first
  if (surfaces.empty())
    ACTS_ERROR("No surfaces given - can not create BinUtility for "
               "BinnedSurfaceArray!");
  // introduce BinUtility to hand back
  Acts::BinningOption bOption = open;
  // all the information already given
  if (bValue == Acts::binPhi) bOption = closed;
  ACTS_VERBOSE("Create BinUtility for BinnedSurfaceArray with equidistant "
               "BinningType");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
               "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << bins);
  ACTS_VERBOSE("	(Min/Max) = (" << min << "/" << max << ")");
  // create the BinUtility
  return (Acts::BinUtility(bins, min, max, bOption, bValue, transform));
}

/// Register the neigbourhood
void
Acts::SurfaceArrayCreator::registerNeighbourHood(SurfaceArray& sArray) const
{
  ACTS_VERBOSE("Register neighbours to the elements.");
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
            // now, mutate bElement to register the neighbours
            auto mutableBElement = const_cast<DetectorElementBase*>(bElement);
            mutableBElement->registerNeighbours(neighbourElements);
          }
        }
        ++io0;
      }
      ++io1;
    }
    ++io2;
  }
  ACTS_VERBOSE("Neighbours set for this layer: " << neighboursSet);
}

/// @todo implement nearest neighbour search - this is brute force attack
/// - takes too long otherwise in initialization
void
Acts::SurfaceArrayCreator::completeBinning(const BinUtility&    binUtility,
                                           const V3Matrix&      v3Matrix,
                                           const SurfaceVector& sVector,
                                           SurfaceGrid&         sGrid) const
{
  ACTS_VERBOSE("Complete binning by filling closest neighbour surfaces into "
               "empty bins.");
  // make a copy of the surface grid
  size_t nSurfaces   = sVector.size();
  size_t nGridPoints = v3Matrix.size() * v3Matrix[0].size();

  // VERBOSE screen output
  ACTS_VERBOSE("- Object count : " << nSurfaces << " number of surfaces");
  ACTS_VERBOSE("- Surface grid : " << nGridPoints << " number of bins");

  size_t binCompleted = 0;
  //
  for (size_t io1 = 0; io1 < v3Matrix.size(); ++io1) {
    for (size_t io0 = 0; io0 < v3Matrix[0].size(); ++io0) {
      // only loop if the bin is empty, else skip
      if (sGrid[0][io1][io0] != nullptr) continue;

      Vector3D sposition = v3Matrix[io1][io0];
      double   minPath   = 10e10;

      for (auto& sf : sVector) {
        double testPath = (sposition - sf->binningPosition(binR)).mag();
        if (testPath < minPath) {
          sGrid[0][io1][io0] = sf;
          minPath            = testPath;
        }
      }

      // inc number of newly filled bins
      ++binCompleted;
    }
  }

  ACTS_VERBOSE("       filled  : " << binCompleted);
}

std::vector<Acts::Vector3D>
Acts::SurfaceArrayCreator::makeGlobalVertices(
    const Acts::Surface&               surface,
    const std::vector<Acts::Vector2D>& locVertices) const
{
  std::vector<Acts::Vector3D> globVertices;
  for (auto& vertex : locVertices) {
    Acts::Vector3D globVertex(0., 0., 0.);
    surface.localToGlobal(vertex, Acts::Vector3D(), globVertex);
    globVertices.push_back(globVertex);
  }
  return globVertices;
}
