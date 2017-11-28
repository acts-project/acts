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
#include "ACTS/Surfaces/SurfaceArray.hpp"
#include "ACTS/Surfaces/concept/AnySurfaceGridLookup.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ACTS/Utilities/detail/Axis.hpp"

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Surface*>& surfaces,
    size_t                             binsPhi,
    size_t                             binsZ,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> _transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << binsPhi << " x " << binsZ << " = "
                                      << binsPhi * binsZ
                                      << " bins.");

  Transform3D transform
      = _transform != nullptr ? *_transform : Transform3D::Identity();

  ProtoAxis pAxisPhi
      = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, binsPhi);
  ProtoAxis pAxisZ
      = createEquidistantAxis(surfaces, binZ, protoLayer, transform, binsZ);

  double R = protoLayer.maxR - protoLayer.minR;

  Transform3D itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(loc.phi(), loc.z());
  };
  auto localToGlobal = [itransform, R](const Vector2D& loc) {
    return itransform
        * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  SurfaceArray::SurfaceGridLookup2D sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Closed,
                                detail::AxisWrapping::Open>(
          surfaces, globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

  sl.fill(surfaces);
  completeBinning(sl, surfaces);

  return std::make_unique<SurfaceArray>(sl, surfaces);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Surface*>& surfaces,
    BinningType                        bTypePhi,
    BinningType                        bTypeZ,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> _transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  double      R = protoLayer.maxR - protoLayer.minR;
  Transform3D transform
      = _transform != nullptr ? *_transform : Transform3D::Identity();

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisZ;

  if (bTypePhi == equidistant)
    pAxisPhi
        = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, 0);
  else
    pAxisPhi = createVariableAxis(surfaces, binPhi, transform);

  if (bTypeZ == equidistant)
    pAxisZ = createEquidistantAxis(surfaces, binZ, protoLayer, transform);
  else
    pAxisZ = createVariableAxis(surfaces, binZ, transform);

  Transform3D itransform = transform.inverse();
  auto globalToLocal     = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(loc.phi(), loc.z());
  };
  auto localToGlobal = [itransform, R](const Vector2D& loc) {
    return itransform
        * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  SurfaceArray::SurfaceGridLookup2D sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Closed,
                                detail::AxisWrapping::Open>(
          surfaces, globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

  sl.fill(surfaces);
  completeBinning(sl, surfaces);

  // get the number of bins
  auto   axes  = sl.getAxes();
  size_t bins0 = axes.at(0).getNBins();
  size_t bins1 = axes.at(1).getNBins();

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");

  return std::make_unique<SurfaceArray>(sl, surfaces);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Surface*>& surfaces,
    size_t                             binsR,
    size_t                             binsPhi,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> _transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3D transform
      = _transform != nullptr ? *_transform : Transform3D::Identity();

  ProtoAxis pAxisR
      = createEquidistantAxis(surfaces, binR, protoLayer, transform, binsR);
  ProtoAxis pAxisPhi
      = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, binsPhi);

  double Z = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  ACTS_VERBOSE("- z-position of disk estimated as " << Z);

  Transform3D itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(loc.phi(), loc.perp());
  };
  auto localToGlobal = [itransform, Z](const Vector2D& loc) {
    return itransform
        * Vector3D(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  SurfaceArray::SurfaceGridLookup2D sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Open,
                                detail::AxisWrapping::Closed>(
          surfaces, globalToLocal, localToGlobal, pAxisR, pAxisPhi);

  // get the number of bins
  auto   axes  = sl.getAxes();
  size_t bins0 = axes.at(0).getNBins();
  size_t bins1 = axes.at(1).getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");
  sl.fill(surfaces);
  completeBinning(sl, surfaces);

  return std::make_unique<SurfaceArray>(sl, surfaces);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Surface*>& surfaces,
    BinningType                        bTypeR,
    BinningType                        bTypePhi,
    boost::optional<ProtoLayer>        _protoLayer,
    std::shared_ptr<const Transform3D> _transform) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3D transform
      = _transform != nullptr ? *_transform : Transform3D::Identity();

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisR;

  if (bTypePhi == equidistant)
    pAxisPhi
        = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, 0);
  else
    pAxisPhi = createVariableAxis(surfaces, binPhi, transform);

  if (bTypeR == equidistant)
    pAxisR = createEquidistantAxis(surfaces, binR, protoLayer, transform);
  else
    pAxisR = createVariableAxis(surfaces, binR, transform);

  double Z = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  ACTS_VERBOSE("- z-position of disk estimated as " << Z);

  Transform3D itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(loc.phi(), loc.perp());
  };
  auto localToGlobal = [itransform, Z](const Vector2D& loc) {
    return itransform
        * Vector3D(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  SurfaceArray::SurfaceGridLookup2D sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Open,
                                detail::AxisWrapping::Closed>(
          surfaces, globalToLocal, localToGlobal, pAxisR, pAxisPhi);

  // get the number of bins
  auto   axes  = sl.getAxes();
  size_t bins0 = axes.at(0).getNBins();
  size_t bins1 = axes.at(1).getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");

  sl.fill(surfaces);
  completeBinning(sl, surfaces);

  return std::make_unique<SurfaceArray>(sl, surfaces);
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

Acts::SurfaceArrayCreator::ProtoAxis
Acts::SurfaceArrayCreator::createVariableAxis(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue,
    Transform3D&                       transform) const
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

  std::vector<double> bValues;
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

  ProtoAxis pAxis;
  pAxis.bType    = arbitrary;
  pAxis.bValue   = bValue;
  pAxis.binEdges = bValues;

  return pAxis;
  // return (Acts::BinUtility(bValues, bOption, bValue));
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

Acts::SurfaceArrayCreator::ProtoAxis
Acts::SurfaceArrayCreator::createEquidistantAxis(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue,
    ProtoLayer                         protoLayer,
    Transform3D&                       transform,
    size_t                             nBins) const
{
  if (!surfaces.size())
    throw std::logic_error(
        "No surfaces handed over for creating equidistant bin utility!");
  // check the binning type first

  double minimum = 0.;
  double maximum = 0.;

  // binning option is open for z and r, in case of phi binning reset later
  // Acts::BinningOption bOption = Acts::open;

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
      // bOption = Acts::closed;

      minimum = -M_PI;
      maximum = M_PI;

      // double step = 2 * M_PI / keys.size();
      double step = 2 * M_PI / binNumber;
      // rotate to max phi module plus one half step
      // this should make sure that phi wrapping at +- pi
      // never falls on a module center
      double max   = maxElem->center().phi();
      double angle = M_PI - (max + 0.5 * step);

      // replace given transform ref
      transform = (transform)*AngleAxis3D(angle, Vector3D::UnitZ());

    } else {
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
  ACTS_VERBOSE("Create equidistant binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
               "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << binNumber);
  ACTS_VERBOSE("	(Min/Max) = (" << minimum << "/" << maximum << ")");

  ProtoAxis pAxis;
  pAxis.max    = maximum;
  pAxis.min    = minimum;
  pAxis.bType  = equidistant;
  pAxis.bValue = bValue;
  pAxis.nBins  = binNumber;

  return pAxis;

  // return detail::Axis<detail::AxisType::Equidistant, wrap>(
  // minimum, maximum, binNumber);
  // return (Acts::BinUtility(
  // binNumber, minimum, maximum, bOption, bValue, transform));
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
  /*auto objectGrid = sArray.objectGrid();
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
  }*/
  // ACTS_VERBOSE("Neighbours set for this layer: " << neighboursSet);
  ACTS_WARNING("registerNeighbourHood does nothing right now");
}

/// @todo implement nearest neighbour search - this is brute force attack
/// - takes too long otherwise in initialization
void
Acts::SurfaceArrayCreator::completeBinning(const BinUtility&    binUtility,
                                           const V3Matrix&      v3Matrix,
                                           const SurfaceVector& sVector,
                                           SurfaceGrid_old&     sGrid) const
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
