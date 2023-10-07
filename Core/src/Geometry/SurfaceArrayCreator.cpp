// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/SurfaceArrayCreator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, size_t binsPhi,
    size_t binsZ, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpack_shared_vector(surfaces);
  // Check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << binsPhi << " x " << binsZ << " = "
                                      << binsPhi * binsZ << " bins.");

  Transform3 ftransform = transform;
  ProtoAxis pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, binPhi,
                                             protoLayer, ftransform, binsPhi);
  ProtoAxis pAxisZ = createEquidistantAxis(gctx, surfacesRaw, binZ, protoLayer,
                                           ftransform, binsZ);

  double R = protoLayer.medium(binR, true);

  Transform3 itransform = ftransform.inverse();
  // Transform lambda captures the transform matrix
  auto globalToLocal = [ftransform](const Vector3& pos) {
    Vector3 loc = ftransform * pos;
    return Vector2(phi(loc), loc.z());
  };
  auto localToGlobal = [itransform, R](const Vector2& loc) {
    return itransform *
           Vector3(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<detail::AxisBoundaryType::Closed,
                              detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

  sl->fill(gctx, surfacesRaw);
  completeBinning(gctx, *sl, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        ftransform);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypePhi,
    BinningType bTypeZ, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpack_shared_vector(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  double R = protoLayer.medium(binR, true);

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisZ;

  Transform3 ftransform = transform;

  if (bTypePhi == equidistant) {
    pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, binPhi, protoLayer,
                                     ftransform, 0);
  } else {
    pAxisPhi =
        createVariableAxis(gctx, surfacesRaw, binPhi, protoLayer, ftransform);
  }

  if (bTypeZ == equidistant) {
    pAxisZ =
        createEquidistantAxis(gctx, surfacesRaw, binZ, protoLayer, ftransform);
  } else {
    pAxisZ =
        createVariableAxis(gctx, surfacesRaw, binZ, protoLayer, ftransform);
  }

  Transform3 itransform = ftransform.inverse();
  auto globalToLocal = [ftransform](const Vector3& pos) {
    Vector3 loc = ftransform * pos;
    return Vector2(phi(loc), loc.z());
  };
  auto localToGlobal = [itransform, R](const Vector2& loc) {
    return itransform *
           Vector3(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<detail::AxisBoundaryType::Closed,
                              detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxisPhi, pAxisZ);

  sl->fill(gctx, surfacesRaw);
  completeBinning(gctx, *sl, surfacesRaw);

  // get the number of bins
  auto axes = sl->getAxes();
  size_t bins0 = axes.at(0)->getNBins();
  size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with phi x z  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        ftransform);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, size_t binsR,
    size_t binsPhi, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpack_shared_vector(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3 ftransform = transform;
  ProtoAxis pAxisR = createEquidistantAxis(gctx, surfacesRaw, binR, protoLayer,
                                           ftransform, binsR);
  ProtoAxis pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, binPhi,
                                             protoLayer, ftransform, binsPhi);

  double Z = protoLayer.medium(binZ, true);
  ACTS_VERBOSE("- z-position of disk estimated as " << Z);

  Transform3 itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [ftransform](const Vector3& pos) {
    Vector3 loc = ftransform * pos;
    return Vector2(perp(loc), phi(loc));
  };
  auto localToGlobal = [itransform, Z](const Vector2& loc) {
    return itransform *
           Vector3(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                              detail::AxisBoundaryType::Closed>(
          globalToLocal, localToGlobal, pAxisR, pAxisPhi);

  // get the number of bins
  auto axes = sl->getAxes();
  size_t bins0 = axes.at(0)->getNBins();
  size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");
  sl->fill(gctx, surfacesRaw);
  completeBinning(gctx, *sl, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        ftransform);
}

std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnDisc(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypeR,
    BinningType bTypePhi, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpack_shared_vector(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisR;

  Transform3 ftransform = transform;

  if (bTypeR == equidistant) {
    pAxisR =
        createEquidistantAxis(gctx, surfacesRaw, binR, protoLayer, ftransform);
  } else {
    pAxisR =
        createVariableAxis(gctx, surfacesRaw, binR, protoLayer, ftransform);
  }

  // if we have more than one R ring, we need to figure out
  // the number of phi bins.
  if (pAxisR.nBins > 1) {
    // more than one R-Ring, we need to adjust
    // this FORCES equidistant binning
    std::vector<std::vector<const Surface*>> phiModules(pAxisR.nBins);
    for (const auto& srf : surfacesRaw) {
      Vector3 bpos = srf->binningPosition(gctx, binR);
      size_t bin = pAxisR.getBin(perp(bpos));
      phiModules.at(bin).push_back(srf);
    }

    std::vector<size_t> nPhiModules;
    auto matcher = m_cfg.surfaceMatcher;
    auto equal = [&gctx, &matcher](const Surface* a, const Surface* b) {
      return matcher(gctx, binPhi, a, b);
    };

    std::transform(
        phiModules.begin(), phiModules.end(), std::back_inserter(nPhiModules),
        [&equal, this](const std::vector<const Surface*>& surfaces_) -> size_t {
          return this->findKeySurfaces(surfaces_, equal).size();
        });

    // @FIXME: Problem: phi binning runs rotation to optimize
    // for bin edges. This FAILS after this modification, since
    // the bin count is the one from the lowest module-count bin,
    // but the rotation is done considering all bins.
    // This might be resolved through bin completion, but not sure.
    // @TODO: check in extrapolation
    size_t nBinsPhi =
        (*std::min_element(nPhiModules.begin(), nPhiModules.end()));
    pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, binPhi, protoLayer,
                                     ftransform, nBinsPhi);

  } else {
    // use regular determination
    if (bTypePhi == equidistant) {
      pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, binPhi, protoLayer,
                                       ftransform, 0);
    } else {
      pAxisPhi =
          createVariableAxis(gctx, surfacesRaw, binPhi, protoLayer, ftransform);
    }
  }

  double Z = protoLayer.medium(binZ, true);
  ACTS_VERBOSE("- z-position of disk estimated as " << Z);

  Transform3 itransform = ftransform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [ftransform](const Vector3& pos) {
    Vector3 loc = ftransform * pos;
    return Vector2(perp(loc), phi(loc));
  };
  auto localToGlobal = [itransform, Z](const Vector2& loc) {
    return itransform *
           Vector3(loc[0] * std::cos(loc[1]), loc[0] * std::sin(loc[1]), Z);
  };

  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                              detail::AxisBoundaryType::Closed>(
          globalToLocal, localToGlobal, pAxisR, pAxisPhi);

  // get the number of bins
  auto axes = sl->getAxes();
  size_t bins0 = axes.at(0)->getNBins();
  size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  sl->fill(gctx, surfacesRaw);
  completeBinning(gctx, *sl, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        ftransform);
}

/// SurfaceArrayCreator interface method - create an array on a plane
std::unique_ptr<Acts::SurfaceArray>
Acts::SurfaceArrayCreator::surfaceArrayOnPlane(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, size_t bins1,
    size_t bins2, BinningValue bValue, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpack_shared_vector(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a plance");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.")
  ACTS_VERBOSE(" -- with " << bins1 << " x " << bins2 << " = " << bins1 * bins2
                           << " bins.");
  Transform3 ftransform = transform;
  Transform3 itransform = transform.inverse();
  // transform lambda captures the transform matrix
  auto globalToLocal = [ftransform](const Vector3& pos) {
    Vector3 loc = ftransform * pos;
    return Vector2(loc.x(), loc.y());
  };
  auto localToGlobal = [itransform](const Vector2& loc) {
    return itransform * Vector3(loc.x(), loc.y(), 0.);
  };
  // Build the grid
  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl;

  // Axis along the binning
  switch (bValue) {
    case BinningValue::binX: {
      ProtoAxis pAxis1 = createEquidistantAxis(gctx, surfacesRaw, binY,
                                               protoLayer, ftransform, bins1);
      ProtoAxis pAxis2 = createEquidistantAxis(gctx, surfacesRaw, binZ,
                                               protoLayer, ftransform, bins2);
      sl = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                                   detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxis1, pAxis2);
      break;
    }
    case BinningValue::binY: {
      ProtoAxis pAxis1 = createEquidistantAxis(gctx, surfacesRaw, binX,
                                               protoLayer, ftransform, bins1);
      ProtoAxis pAxis2 = createEquidistantAxis(gctx, surfacesRaw, binZ,
                                               protoLayer, ftransform, bins2);
      sl = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                                   detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxis1, pAxis2);
      break;
    }
    case BinningValue::binZ: {
      ProtoAxis pAxis1 = createEquidistantAxis(gctx, surfacesRaw, binX,
                                               protoLayer, ftransform, bins1);
      ProtoAxis pAxis2 = createEquidistantAxis(gctx, surfacesRaw, binY,
                                               protoLayer, ftransform, bins2);
      sl = makeSurfaceGridLookup2D<detail::AxisBoundaryType::Bound,
                                   detail::AxisBoundaryType::Bound>(
          globalToLocal, localToGlobal, pAxis1, pAxis2);
      break;
    }
    default: {
      throw std::invalid_argument(
          "Acts::SurfaceArrayCreator::"
          "surfaceArrayOnPlane: Invalid binning "
          "direction");
    }
  }

  sl->fill(gctx, surfacesRaw);
  completeBinning(gctx, *sl, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        ftransform);
  //!< @todo implement - take from ATLAS complex TRT builder
}

std::vector<const Acts::Surface*> Acts::SurfaceArrayCreator::findKeySurfaces(
    const std::vector<const Surface*>& surfaces,
    const std::function<bool(const Surface*, const Surface*)>& equal) const {
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

size_t Acts::SurfaceArrayCreator::determineBinCount(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    BinningValue bValue) const {
  auto matcher = m_cfg.surfaceMatcher;
  auto equal = [&gctx, &bValue, &matcher](const Surface* a, const Surface* b) {
    return matcher(gctx, bValue, a, b);
  };
  std::vector<const Surface*> keys = findKeySurfaces(surfaces, equal);

  return keys.size();
}

Acts::SurfaceArrayCreator::ProtoAxis
Acts::SurfaceArrayCreator::createVariableAxis(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    BinningValue bValue, const ProtoLayer& protoLayer,
    Transform3& transform) const {
  if (surfaces.empty()) {
    throw std::logic_error(
        "No surfaces handed over for creating arbitrary bin utility!");
  }
  // BinningOption is open for z and r, in case of phi binning reset later
  // the vector with the binning Values (boundaries for each bin)

  // bind matcher with binning type
  auto matcher = m_cfg.surfaceMatcher;
  // find the key surfaces
  auto equal = [&gctx, &bValue, &matcher](const Surface* a, const Surface* b) {
    return matcher(gctx, bValue, a, b);
  };
  std::vector<const Acts::Surface*> keys = findKeySurfaces(surfaces, equal);

  std::vector<AxisScalar> bValues;
  if (bValue == Acts::binPhi) {
    std::stable_sort(keys.begin(), keys.end(),
                     [&gctx](const Acts::Surface* a, const Acts::Surface* b) {
                       return (phi(a->binningPosition(gctx, binPhi)) <
                               phi(b->binningPosition(gctx, binPhi)));
                     });

    AxisScalar maxPhi = 0.5 * (phi(keys.at(0)->binningPosition(gctx, binPhi)) +
                               phi(keys.at(1)->binningPosition(gctx, binPhi)));

    // create rotation, so that maxPhi is +pi
    AxisScalar angle = -(M_PI + maxPhi);
    transform = (transform)*AngleAxis3(angle, Vector3::UnitZ());

    // iterate over all key surfaces, and use their mean position as bValues,
    // but
    // rotate using transform from before
    AxisScalar previous = phi(keys.at(0)->binningPosition(gctx, binPhi));
    // go through key surfaces
    for (size_t i = 1; i < keys.size(); i++) {
      const Surface* surface = keys.at(i);
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      AxisScalar edge =
          0.5 * (previous + phi(surface->binningPosition(gctx, binPhi))) +
          angle;
      bValues.push_back(edge);
      previous = phi(surface->binningPosition(gctx, binPhi));
    }

    // segments
    unsigned int segments = 72;

    // get the bounds of the last surfaces
    const Acts::Surface* backSurface = keys.back();
    const Acts::PlanarBounds* backBounds =
        dynamic_cast<const Acts::PlanarBounds*>(&(backSurface->bounds()));
    if (backBounds == nullptr) {
      ACTS_ERROR(
          "Given SurfaceBounds are not planar - not implemented for "
          "other bounds yet! ");
    }
    // get the global vertices
    std::vector<Acts::Vector3> backVertices =
        makeGlobalVertices(gctx, *backSurface, backBounds->vertices(segments));
    AxisScalar maxBValue = phi(
        *std::max_element(backVertices.begin(), backVertices.end(),
                          [](const Acts::Vector3& a, const Acts::Vector3& b) {
                            return phi(a) < phi(b);
                          }));

    bValues.push_back(maxBValue);

    bValues.push_back(M_PI);

  } else if (bValue == Acts::binZ) {
    std::stable_sort(keys.begin(), keys.end(),
                     [&gctx](const Acts::Surface* a, const Acts::Surface* b) {
                       return (a->binningPosition(gctx, binZ).z() <
                               b->binningPosition(gctx, binZ).z());
                     });

    bValues.push_back(protoLayer.min(binZ));
    bValues.push_back(protoLayer.max(binZ));

    // the z-center position of the previous surface
    AxisScalar previous = keys.front()->binningPosition(gctx, binZ).z();
    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(
          0.5 * (previous + (*surface)->binningPosition(gctx, binZ).z()));
      previous = (*surface)->binningPosition(gctx, binZ).z();
    }
  } else {  // binR
    std::stable_sort(keys.begin(), keys.end(),
                     [&gctx](const Acts::Surface* a, const Acts::Surface* b) {
                       return (perp(a->binningPosition(gctx, binR)) <
                               perp(b->binningPosition(gctx, binR)));
                     });

    bValues.push_back(protoLayer.min(binR));
    bValues.push_back(protoLayer.max(binR));

    // the r-center position of the previous surface
    AxisScalar previous = perp(keys.front()->binningPosition(gctx, binR));

    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      bValues.push_back(
          0.5 * (previous + perp((*surface)->binningPosition(gctx, binR))));
      previous = perp((*surface)->binningPosition(gctx, binR));
    }
  }
  std::sort(bValues.begin(), bValues.end());
  ACTS_VERBOSE("Create variable binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE(
      "	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
      "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << (bValues.size() - 1));
  ACTS_VERBOSE("	(Min/Max) = (" << bValues.front() << "/"
                                       << bValues.back() << ")");

  ProtoAxis pAxis;
  pAxis.bType = arbitrary;
  pAxis.bValue = bValue;
  pAxis.binEdges = bValues;
  pAxis.nBins = bValues.size() - 1;

  return pAxis;
}

Acts::SurfaceArrayCreator::ProtoAxis
Acts::SurfaceArrayCreator::createEquidistantAxis(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    BinningValue bValue, const ProtoLayer& protoLayer, Transform3& transform,
    size_t nBins) const {
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

  size_t binNumber = 0;
  if (nBins == 0) {
    // determine bin count
    binNumber = determineBinCount(gctx, surfaces, bValue);
  } else {
    // use bin count
    binNumber = nBins;
  }

  // bind matcher & context with binning type
  auto matcher = m_cfg.surfaceMatcher;

  // now check the binning value
  if (bValue == binPhi) {
    if (m_cfg.doPhiBinningOptimization) {
      // Phi binning
      // set the binning option for phi
      // sort first in phi
      const Acts::Surface* maxElem = *std::max_element(
          surfaces.begin(), surfaces.end(),
          [&gctx](const Acts::Surface* a, const Acts::Surface* b) {
            return phi(a->binningPosition(gctx, binR)) <
                   phi(b->binningPosition(gctx, binR));
          });

      // get the key surfaces at the different phi positions
      auto equal = [&gctx, &bValue, &matcher](const Surface* a,
                                              const Surface* b) {
        return matcher(gctx, bValue, a, b);
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
        double max = phi(maxElem->binningPosition(gctx, binR));
        double angle = M_PI - (max + 0.5 * step);

        // replace given transform ref
        transform = (transform)*AngleAxis3(angle, Vector3::UnitZ());

      } else {
        minimum = protoLayer.min(binPhi, true);
        maximum = protoLayer.max(binPhi, true);
        // we do not need a transform in this case
      }
    } else {
      minimum = -M_PI;
      maximum = M_PI;
    }
  } else {
    maximum = protoLayer.max(bValue, false);
    minimum = protoLayer.min(bValue, false);
  }

  // assign the bin size
  ACTS_VERBOSE("Create equidistant binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	BinningValue: " << bValue);
  ACTS_VERBOSE(
      "	(binX = 0, binY = 1, binZ = 2, binR = 3, binPhi = 4, "
      "binRPhi = 5, binH = 6, binEta = 7)");
  ACTS_VERBOSE("	Number of bins: " << binNumber);
  ACTS_VERBOSE("	(Min/Max) = (" << minimum << "/" << maximum << ")");

  ProtoAxis pAxis;
  pAxis.max = maximum;
  pAxis.min = minimum;
  pAxis.bType = equidistant;
  pAxis.bValue = bValue;
  pAxis.nBins = binNumber;

  return pAxis;
}

std::vector<Acts::Vector3> Acts::SurfaceArrayCreator::makeGlobalVertices(
    const GeometryContext& gctx, const Acts::Surface& surface,
    const std::vector<Acts::Vector2>& locVertices) const {
  std::vector<Acts::Vector3> globVertices;
  for (auto& vertex : locVertices) {
    Acts::Vector3 globVertex =
        surface.localToGlobal(gctx, vertex, Acts::Vector3());
    globVertices.push_back(globVertex);
  }
  return globVertices;
}
