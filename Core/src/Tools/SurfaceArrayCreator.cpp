// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include <algorithm>
#include <cmath>
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/detail/Axis.hpp"

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Surface*>&        surfaces,
    size_t                                    binsPhi,
    size_t                                    binsZ,
    boost::optional<ProtoLayer>               protoLayerOpt,
    const std::shared_ptr<const Transform3D>& transformOpt) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = protoLayerOpt ? *protoLayerOpt : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << binsPhi << " x " << binsZ << " = "
                                      << binsPhi * binsZ
                                      << " bins.");

  Transform3D transform
      = transformOpt != nullptr ? *transformOpt : Transform3D::Identity();

  ProtoAxis pAxisPhi
      = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, binsPhi);
  ProtoAxis pAxisZ
      = createEquidistantAxis(surfaces, binZ, protoLayer, transform, binsZ);

  double R = protoLayer.maxR - protoLayer.minR;

  Transform3D itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(LA::phi(loc), loc.z());
  };
  auto localToGlobal = [itransform, R](const Vector2D& loc) {
    return itransform
        * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl
      = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Closed,
                                detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

  sl->fill(surfaces);
  completeBinning(*sl, surfaces);

  return std::make_unique<SurfaceArray>(
      std::move(sl), surfaces, std::make_shared<const Transform3D>(transform));
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const std::vector<const Surface*>&        surfaces,
    BinningType                               bTypePhi,
    BinningType                               bTypeZ,
    boost::optional<ProtoLayer>               protoLayerOpt,
    const std::shared_ptr<const Transform3D>& transformOpt) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = protoLayerOpt ? *protoLayerOpt : ProtoLayer(surfaces);

  double      R = 0.5 * (protoLayer.maxR - protoLayer.minR);
  Transform3D transform
      = transformOpt != nullptr ? *transformOpt : Transform3D::Identity();

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisZ;

  if (bTypePhi == equidistant) {
    pAxisPhi
        = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, 0);
  } else {
    pAxisPhi = createVariableAxis(surfaces, binPhi, protoLayer, transform);
  }

  if (bTypeZ == equidistant) {
    pAxisZ = createEquidistantAxis(surfaces, binZ, protoLayer, transform);
  } else {
    pAxisZ = createVariableAxis(surfaces, binZ, protoLayer, transform);
  }

  Transform3D itransform = transform.inverse();
  auto globalToLocal     = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(LA::phi(loc), loc.z());
  };
  auto localToGlobal = [itransform, R](const Vector2D& loc) {
    return itransform
        * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl
      = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Closed,
                                detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

  sl->fill(surfaces);
  completeBinning(*sl, surfaces);

  // get the number of bins
  auto   axes  = sl->getAxes();
  size_t bins0 = axes.at(0)->getNBins();
  size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");

  return std::make_unique<SurfaceArray>(
      std::move(sl), surfaces, std::make_shared<const Transform3D>(transform));
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Surface*>&        surfaces,
    size_t                                    binsR,
    size_t                                    binsPhi,
    boost::optional<ProtoLayer>               protoLayerOpt,
    const std::shared_ptr<const Transform3D>& transformOpt) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = protoLayerOpt ? *protoLayerOpt : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3D transform
      = transformOpt != nullptr ? *transformOpt : Transform3D::Identity();

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
    return Vector2D(LA::perp(loc), LA::phi(loc));
  };
  auto localToGlobal = [itransform, Z](const Vector2D& loc) {
    return itransform
        * Vector3D(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl
      = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                                detail::AxisBoundaryType::Closed>(
          globalToLocal, localToGlobal, pAxisR, pAxisPhi);

  // get the number of bins
  auto   axes  = sl->getAxes();
  size_t bins0 = axes.at(0)->getNBins();
  size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");
  sl->fill(surfaces);
  completeBinning(*sl, surfaces);

  return std::make_unique<SurfaceArray>(
      std::move(sl), surfaces, std::make_shared<const Transform3D>(transform));
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const std::vector<const Surface*>&        surfaces,
    BinningType                               bTypeR,
    BinningType                               bTypePhi,
    boost::optional<ProtoLayer>               protoLayerOpt,
    const std::shared_ptr<const Transform3D>& transformOpt) const
{
  // check if we have proto layer, else build it
  ProtoLayer protoLayer = protoLayerOpt ? *protoLayerOpt : ProtoLayer(surfaces);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3D transform
      = transformOpt != nullptr ? *transformOpt : Transform3D::Identity();

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisR;

  if (bTypeR == equidistant) {
    pAxisR = createEquidistantAxis(surfaces, binR, protoLayer, transform);
  } else {
    pAxisR = createVariableAxis(surfaces, binR, protoLayer, transform);
  }

  // if we have more than one R ring, we need to figure out
  // the number of phi bins.
  if (pAxisR.nBins > 1) {
    // more than one R-Ring, we need to adjust
    // this FORCES equidistant binning
    std::vector<std::vector<const Surface*>> phiModules(pAxisR.nBins);
    for (const auto& srf : surfaces) {
      Vector3D bpos = srf->binningPosition(binR);
      size_t   bin  = pAxisR.getBin(LA::perp(bpos));
      phiModules.at(bin).push_back(srf);
    }

    std::vector<size_t> nPhiModules;
    auto                matcher = m_cfg.surfaceMatcher;
    auto equal = [&matcher](const Surface* a, const Surface* b) {
      return matcher(binPhi, a, b);
    };

    std::transform(
        phiModules.begin(),
        phiModules.end(),
        std::back_inserter(nPhiModules),
        [&equal, this](std::vector<const Surface*> surfaces_) -> size_t {
          return this->findKeySurfaces(surfaces_, equal).size();
        });

    // @FIXME: Problem: phi binning runs rotation to optimize
    // for bin edges. This FAILS after this modification, since
    // the bin count is the one from the lowest module-count bin,
    // but the rotation is done considering all bins.
    // This might be resolved through bin completion, but not sure.
    // @TODO: check in extrapolation
    size_t nBinsPhi
        = (*std::min_element(nPhiModules.begin(), nPhiModules.end()));
    pAxisPhi = createEquidistantAxis(
        surfaces, binPhi, protoLayer, transform, nBinsPhi);

  } else {
    // use regular determination
    if (bTypePhi == equidistant) {
      pAxisPhi
          = createEquidistantAxis(surfaces, binPhi, protoLayer, transform, 0);
    } else {
      pAxisPhi = createVariableAxis(surfaces, binPhi, protoLayer, transform);
    }
  }

  double Z = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  ACTS_VERBOSE("- z-position of disk estimated as " << Z);

  Transform3D itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [transform](const Vector3D& pos) {
    Vector3D loc = transform * pos;
    return Vector2D(LA::perp(loc), LA::phi(loc));
  };
  auto localToGlobal = [itransform, Z](const Vector2D& loc) {
    return itransform
        * Vector3D(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl
      = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                                detail::AxisBoundaryType::Closed>(
          globalToLocal, localToGlobal, pAxisR, pAxisPhi);

  // get the number of bins
  auto   axes  = sl->getAxes();
  size_t bins0 = axes.at(0)->getNBins();
  size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1
                                      << " bins.");

  sl->fill(surfaces);
  completeBinning(*sl, surfaces);

  return std::make_unique<SurfaceArray>(
      std::move(sl), surfaces, std::make_shared<const Transform3D>(transform));
}

/// SurfaceArrayCreator interface method - create an array on a plane
std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnPlane(
    const std::vector<const Surface*>& /*surfaces*/,
    double /*halflengthX*/,
    double /*halflengthY*/,
    size_t /*binsX*/,
    size_t /*binsY*/,
    const std::shared_ptr<const Transform3D>& /*transform*/) const
{
  //!< @todo implement - take from ATLAS complex TRT builder
  return nullptr;
}

std::vector<const Acts::Surface*>
Acts::SurfaceArrayCreator::findKeySurfaces(
    const std::vector<const Surface*>& surfaces,
    const std::function<bool(const Surface*, const Surface*)>& equal) const
{
  std::vector<const Surface*> keys;
  for (const auto& srfA : surfaces) {
    bool exists = false;
    for (const auto& srfB : keys) {
      if (equal(srfA, srfB)) {
        exists = true;
        break;
      }
    }
    if (!exists) {
      keys.push_back(srfA);
    }
  }

  return keys;
}

size_t
Acts::SurfaceArrayCreator::determineBinCount(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue) const
{

  auto matcher = m_cfg.surfaceMatcher;
  auto equal   = [&bValue, &matcher](const Surface* a, const Surface* b) {
    return matcher(bValue, a, b);
  };
  std::vector<const Surface*> keys = findKeySurfaces(surfaces, equal);

  return keys.size();
}

Acts::SurfaceArrayCreator::ProtoAxis
Acts::SurfaceArrayCreator::createVariableAxis(
    const std::vector<const Surface*>& surfaces,
    BinningValue                       bValue,
    ProtoLayer                         protoLayer,
    Transform3D&                       transform) const
{
  if (surfaces.empty()) {
    throw std::logic_error(
        "No surfaces handed over for creating arbitrary bin utility!");
  }
  // BinningOption is open for z and r, in case of phi binning reset later
  // the vector with the binning Values (boundaries for each bin)

  // bind matcher with binning type
  auto matcher = m_cfg.surfaceMatcher;
  // find the key surfaces
  auto equal = [&bValue, &matcher](const Surface* a, const Surface* b) {
    return matcher(bValue, a, b);
  };
  std::vector<const Acts::Surface*> keys = findKeySurfaces(surfaces, equal);

  std::vector<double> bValues;
  if (bValue == Acts::binPhi) {
    std::stable_sort(keys.begin(),
                     keys.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (LA::phi(a->binningPosition(binPhi))
                               < LA::phi(b->binningPosition(binPhi)));
                     });

    double maxPhi = 0.5 * (LA::phi(keys.at(0)->binningPosition(binPhi))
                           + LA::phi(keys.at(1)->binningPosition(binPhi)));

    // create rotation, so that maxPhi is +pi
    double angle = -(M_PI + maxPhi);
    transform    = (transform)*AngleAxis3D(angle, Vector3D::UnitZ());

    // iterate over all key surfaces, and use their mean position as bValues,
    // but
    // rotate using transform from before
    double previous = LA::phi(keys.at(0)->binningPosition(binPhi));
    // go through key surfaces
    for (size_t i = 1; i < keys.size(); i++) {
      const Surface* surface = keys.at(i);
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      double edge = 0.5 * (previous + LA::phi(surface->binningPosition(binPhi)))
          + angle;
      bValues.push_back(edge);
      previous = LA::phi(surface->binningPosition(binPhi));
    }

    // get the bounds of the last surfaces
    const Acts::Surface*      backSurface = keys.back();
    const Acts::PlanarBounds* backBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (backBounds == nullptr)
      ACTS_ERROR("Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the global vertices
    std::vector<Acts::Vector3D> backVertices
        = makeGlobalVertices(*backSurface, backBounds->vertices());
    double maxBValue = LA::phi(
        *std::max_element(backVertices.begin(),
                          backVertices.end(),
                          [](const Acts::Vector3D& a, const Acts::Vector3D& b) {
                            return LA::phi(a) < LA::phi(b);
                          }));

    bValues.push_back(maxBValue);

    bValues.push_back(M_PI);

  } else if (bValue == Acts::binZ) {
    std::stable_sort(keys.begin(),
                     keys.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->binningPosition(binZ).z()
                               < b->binningPosition(binZ).z());
                     });

    bValues.push_back(protoLayer.minZ);
    bValues.push_back(protoLayer.maxZ);

    // the z-center position of the previous surface
    double previous = keys.front()->binningPosition(binZ).z();
    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(0.5
                        * (previous + (*surface)->binningPosition(binZ).z()));
      previous = (*surface)->binningPosition(binZ).z();
    }
  } else {  // binR
    std::stable_sort(keys.begin(),
                     keys.end(),
                     [](const Acts::Surface* a, const Acts::Surface* b) {
                       return (LA::perp(a->binningPosition(binR))
                               < LA::perp(b->binningPosition(binR)));
                     });

    bValues.push_back(protoLayer.minR);
    bValues.push_back(protoLayer.maxR);

    // the r-center position of the previous surface
    double previous = LA::perp(keys.front()->binningPosition(binR));

    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(
          0.5 * (previous + LA::perp((*surface)->binningPosition(binR))));
      previous = LA::perp((*surface)->binningPosition(binR));
    }
  }
  std::sort(bValues.begin(), bValues.end());
  ACTS_VERBOSE("Create variable binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE("	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
               "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << (bValues.size() - 1));
  ACTS_VERBOSE("	(Min/Max) = (" << bValues.front() << "/" << bValues.back()
                                 << ")");

  ProtoAxis pAxis;
  pAxis.bType    = arbitrary;
  pAxis.bValue   = bValue;
  pAxis.binEdges = bValues;
  pAxis.nBins    = bValues.size() - 1;

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
  if (surfaces.empty()) {
    throw std::logic_error(
        "No surfaces handed over for creating equidistant axis!");
  }
  // check the binning type first

  double minimum = 0.;
  double maximum = 0.;

  // binning option is open for z and r, in case of phi binning reset later
  // Acts::BinningOption bOption = Acts::open;

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
  auto matcher = m_cfg.surfaceMatcher;

  // now check the binning value
  if (bValue == Acts::binPhi) {

    if (m_cfg.doPhiBinningOptimization) {
      // Phi binning
      // set the binning option for phi
      // sort first in phi
      const Acts::Surface* maxElem = *std::max_element(
          surfaces.begin(),
          surfaces.end(),
          [](const Acts::Surface* a, const Acts::Surface* b) {
            return LA::phi(a->binningPosition(binR))
                < LA::phi(b->binningPosition(binR));
          });

      // get the key surfaces at the different phi positions
      auto equal = [&bValue, &matcher](const Surface* a, const Surface* b) {
        return matcher(bValue, a, b);
      };
      keys = findKeySurfaces(surfaces, equal);

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
        double max   = LA::phi(maxElem->binningPosition(binR));
        double angle = M_PI - (max + 0.5 * step);

        // replace given transform ref
        transform = (transform)*AngleAxis3D(angle, Vector3D::UnitZ());

      } else {
        minimum = protoLayer.minPhi;
        maximum = protoLayer.maxPhi;

        // we do not need a transform in this case
      }
    } else {
      minimum = -M_PI;
      maximum = M_PI;
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
