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

  SurfaceArray::AnySurfaceGridLookup_t sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Closed,
                                detail::AxisWrapping::Open>(
          globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

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

  SurfaceArray::AnySurfaceGridLookup_t sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Closed,
                                detail::AxisWrapping::Open>(
          globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

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
    return Vector2D(loc.perp(), loc.phi());
  };
  auto localToGlobal = [itransform, Z](const Vector2D& loc) {
    return itransform
        * Vector3D(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  SurfaceArray::AnySurfaceGridLookup_t sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Open,
                                detail::AxisWrapping::Closed>(
          globalToLocal, localToGlobal, pAxisR, pAxisPhi);

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
    return Vector2D(loc.perp(), loc.phi());
  };
  auto localToGlobal = [itransform, Z](const Vector2D& loc) {
    return itransform
        * Vector3D(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  SurfaceArray::AnySurfaceGridLookup_t sl
      = makeSurfaceGridLookup2D<detail::AxisWrapping::Open,
                                detail::AxisWrapping::Closed>(
          globalToLocal, localToGlobal, pAxisR, pAxisPhi);

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
  // the vector with the binning Values (boundaries for each bin)

  // bind matcher with binning type
  auto matcher = std::bind(m_cfg.surfaceMatcher,
                           bValue,
                           std::placeholders::_1,
                           std::placeholders::_2);

  std::vector<double> bValues;
  if (bValue == Acts::binPhi) {
    // set the BinningOption closed for binning in phi
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

    // get the bounds of the last surfaces
    const Acts::Surface*      backSurface = keys.back();
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (!backBounds)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");

    double maxPhi
        = 0.5 * (keys.at(0)->center().phi() + keys.at(1)->center().phi());

    // create rotation, so that maxPhi is +pi
    double angle = -(M_PI + maxPhi);
    transform    = (transform)*AngleAxis3D(angle, Vector3D::UnitZ());

    // iterate over all key surfaces, and use their mean position as bValues,
    // but
    // rotate using transform from before
    double previous = keys.at(0)->center().phi();
    // go through key surfaces
    for (size_t i = 1; i < keys.size(); i++) {
      const Surface* surface = keys.at(i);
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      double edge = 0.5 * (previous + surface->center().phi()) + angle;
      bValues.push_back(edge);
      previous = surface->center().phi();
    }

    bValues.push_back(M_PI);
    ;

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
  ACTS_VERBOSE("Create variable binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
               "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << (bValues.size() - 1));

  ProtoAxis pAxis;
  pAxis.bType    = arbitrary;
  pAxis.bValue   = bValue;
  pAxis.binEdges = bValues;

  return pAxis;
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
        "No surfaces handed over for creating equidistant axis!");
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
