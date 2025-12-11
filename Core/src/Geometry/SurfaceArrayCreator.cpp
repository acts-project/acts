// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/SurfaceArrayCreator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <numbers>
#include <stdexcept>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

std::unique_ptr<SurfaceArray> SurfaceArrayCreator::surfaceArrayOnCylinder(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, std::size_t binsPhi,
    std::size_t binsZ, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // Check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with phi x z  = " << binsPhi << " x " << binsZ << " = "
                                      << binsPhi * binsZ << " bins.");

  Transform3 fullTransform = transform;
  ProtoAxis pAxisPhi =
      createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                            protoLayer, fullTransform, binsPhi);
  ProtoAxis pAxisZ =
      createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisZ, protoLayer,
                            fullTransform, binsZ);

  const double R = protoLayer.medium(AxisDirection::AxisR, true);
  const double halfZ = protoLayer.range(AxisDirection::AxisZ, true) * 0.5;
  const double layerTolerance = protoLayer.range(AxisDirection::AxisR) * 0.5;

  auto surface = Surface::makeShared<CylinderSurface>(fullTransform, R, halfZ);
  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<AxisBoundaryType::Closed,
                              AxisBoundaryType::Bound>(
          std::move(surface), layerTolerance, pAxisPhi, pAxisZ);

  sl->fill(gctx, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        fullTransform);
}

std::unique_ptr<SurfaceArray> SurfaceArrayCreator::surfaceArrayOnCylinder(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypePhi,
    BinningType bTypeZ, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  const double R = protoLayer.medium(AxisDirection::AxisR, true);
  const double halfZ = protoLayer.range(AxisDirection::AxisZ, true) * 0.5;
  const double layerTolerance = protoLayer.range(AxisDirection::AxisR) * 0.5;

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisZ;

  Transform3 fullTransform = transform;

  if (bTypePhi == equidistant) {
    pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                                     protoLayer, fullTransform, 0);
  } else {
    pAxisPhi = createVariableAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                                  protoLayer, fullTransform);
  }

  if (bTypeZ == equidistant) {
    pAxisZ = createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisZ,
                                   protoLayer, fullTransform);
  } else {
    pAxisZ = createVariableAxis(gctx, surfacesRaw, AxisDirection::AxisZ,
                                protoLayer, fullTransform);
  }

  auto surface = Surface::makeShared<CylinderSurface>(fullTransform, R, halfZ);
  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<AxisBoundaryType::Closed,
                              AxisBoundaryType::Bound>(
          std::move(surface), layerTolerance, pAxisPhi, pAxisZ);

  sl->fill(gctx, surfacesRaw);

  // get the number of bins
  auto axes = sl->getAxes();
  const std::size_t bins0 = axes.at(0)->getNBins();
  const std::size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with phi x z  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        fullTransform);
}

std::unique_ptr<SurfaceArray> SurfaceArrayCreator::surfaceArrayOnDisc(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, std::size_t binsR,
    std::size_t binsPhi, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3 fullTransform = transform;
  ProtoAxis pAxisR =
      createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisR, protoLayer,
                            fullTransform, binsR);
  ProtoAxis pAxisPhi =
      createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                            protoLayer, fullTransform, binsPhi);

  const double Z = protoLayer.medium(AxisDirection::AxisZ, true);
  const double Rmin = protoLayer.min(AxisDirection::AxisR, true);
  const double Rmax = protoLayer.max(AxisDirection::AxisR, true);
  const double layerThickness = protoLayer.range(AxisDirection::AxisZ) * 0.5;
  ACTS_VERBOSE("- z-position of disc estimated as " << Z);

  auto surface = Surface::makeShared<DiscSurface>(fullTransform, Rmin, Rmax);
  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<AxisBoundaryType::Bound,
                              AxisBoundaryType::Closed>(
          std::move(surface), layerThickness, pAxisR, pAxisPhi);

  // get the number of bins
  auto axes = sl->getAxes();
  std::size_t bins0 = axes.at(0)->getNBins();
  std::size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  sl->fill(gctx, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        fullTransform);
}

std::unique_ptr<SurfaceArray> SurfaceArrayCreator::surfaceArrayOnDisc(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypeR,
    BinningType bTypePhi, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  ProtoAxis pAxisPhi;
  ProtoAxis pAxisR;

  Transform3 fullTransform = transform;
  Transform3 inverseTransform = transform.inverse();

  if (bTypeR == equidistant) {
    pAxisR = createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisR,
                                   protoLayer, fullTransform);
  } else {
    pAxisR = createVariableAxis(gctx, surfacesRaw, AxisDirection::AxisR,
                                protoLayer, fullTransform);
  }

  // if we have more than one R ring, we need to figure out
  // the number of phi bins.
  if (pAxisR.nBins > 1) {
    // more than one R-Ring, we need to adjust
    // this FORCES equidistant binning
    std::vector<std::vector<const Surface*>> phiModules(pAxisR.nBins);
    for (const auto& srf : surfacesRaw) {
      Vector3 bpos =
          inverseTransform * srf->referencePosition(gctx, AxisDirection::AxisR);
      std::size_t bin = pAxisR.getBin(perp(bpos));
      phiModules.at(bin).push_back(srf);
    }

    std::vector<std::size_t> nPhiModules;
    auto matcher = m_cfg.surfaceMatcher;
    auto equal = [&gctx, &matcher](const Surface* a, const Surface* b) {
      return matcher(gctx, AxisDirection::AxisPhi, a, b);
    };

    std::transform(
        phiModules.begin(), phiModules.end(), std::back_inserter(nPhiModules),
        [&equal,
         this](const std::vector<const Surface*>& surfaces_) -> std::size_t {
          return this->findKeySurfaces(surfaces_, equal).size();
        });

    // @FIXME: Problem: phi binning runs rotation to optimize
    // for bin edges. This FAILS after this modification, since
    // the bin count is the one from the lowest module-count bin,
    // but the rotation is done considering all bins.
    // This might be resolved through bin completion, but not sure.
    // @TODO: check in extrapolation
    std::size_t nBinsPhi =
        (*std::min_element(nPhiModules.begin(), nPhiModules.end()));
    pAxisPhi = createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                                     protoLayer, fullTransform, nBinsPhi);

  } else {
    // use regular determination
    if (bTypePhi == equidistant) {
      pAxisPhi =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                                protoLayer, fullTransform, 0);
    } else {
      pAxisPhi = createVariableAxis(gctx, surfacesRaw, AxisDirection::AxisPhi,
                                    protoLayer, fullTransform);
    }
  }

  const double Z = protoLayer.medium(AxisDirection::AxisZ, true);
  const double Rmin = protoLayer.min(AxisDirection::AxisR, true);
  const double Rmax = protoLayer.max(AxisDirection::AxisR, true);
  const double layerThickness = protoLayer.range(AxisDirection::AxisZ) * 0.5;
  ACTS_VERBOSE("- z-position of disc estimated as " << Z);

  auto surface = Surface::makeShared<DiscSurface>(fullTransform, Rmin, Rmax);
  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl =
      makeSurfaceGridLookup2D<AxisBoundaryType::Bound,
                              AxisBoundaryType::Closed>(
          std::move(surface), layerThickness, pAxisR, pAxisPhi);

  // get the number of bins
  auto axes = sl->getAxes();
  const std::size_t bins0 = axes.at(0)->getNBins();
  const std::size_t bins1 = axes.at(1)->getNBins();

  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  sl->fill(gctx, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        fullTransform);
}

/// SurfaceArrayCreator interface method - create an array on a plane
std::unique_ptr<SurfaceArray> SurfaceArrayCreator::surfaceArrayOnPlane(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, std::size_t bins1,
    std::size_t bins2, AxisDirection aDir,
    std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform) const {
  std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a plance");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with " << bins1 << " x " << bins2 << " = " << bins1 * bins2
                           << " bins.");
  Transform3 fullTransform = transform;
  // Build the grid
  std::unique_ptr<SurfaceArray::ISurfaceGridLookup> sl;

  const double layerTolerance = protoLayer.range(aDir) * 0.5;

  // Axis along the binning
  switch (aDir) {
    case AxisDirection::AxisX: {
      ProtoAxis pAxis1 =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisY,
                                protoLayer, fullTransform, bins1);
      ProtoAxis pAxis2 =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisZ,
                                protoLayer, fullTransform, bins2);
      auto surface = Surface::makeShared<PlaneSurface>(
          fullTransform, std::make_shared<RectangleBounds>(
                             Vector2(protoLayer.min(AxisDirection::AxisY),
                                     protoLayer.min(AxisDirection::AxisZ)),
                             Vector2(protoLayer.max(AxisDirection::AxisY),
                                     protoLayer.max(AxisDirection::AxisZ))));
      sl = makeSurfaceGridLookup2D<AxisBoundaryType::Bound,
                                   AxisBoundaryType::Bound>(
          std::move(surface), layerTolerance, pAxis1, pAxis2);
      break;
    }
    case AxisDirection::AxisY: {
      ProtoAxis pAxis1 =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisX,
                                protoLayer, fullTransform, bins1);
      ProtoAxis pAxis2 =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisZ,
                                protoLayer, fullTransform, bins2);
      auto surface = Surface::makeShared<PlaneSurface>(
          fullTransform, std::make_shared<RectangleBounds>(
                             Vector2(protoLayer.min(AxisDirection::AxisX),
                                     protoLayer.min(AxisDirection::AxisY)),
                             Vector2(protoLayer.max(AxisDirection::AxisX),
                                     protoLayer.max(AxisDirection::AxisY))));
      sl = makeSurfaceGridLookup2D<AxisBoundaryType::Bound,
                                   AxisBoundaryType::Bound>(
          std::move(surface), layerTolerance, pAxis1, pAxis2);
      break;
    }
    case AxisDirection::AxisZ: {
      ProtoAxis pAxis1 =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisX,
                                protoLayer, fullTransform, bins1);
      ProtoAxis pAxis2 =
          createEquidistantAxis(gctx, surfacesRaw, AxisDirection::AxisY,
                                protoLayer, fullTransform, bins2);
      auto surface = Surface::makeShared<PlaneSurface>(
          fullTransform, std::make_shared<RectangleBounds>(
                             Vector2(protoLayer.min(AxisDirection::AxisX),
                                     protoLayer.min(AxisDirection::AxisY)),
                             Vector2(protoLayer.max(AxisDirection::AxisX),
                                     protoLayer.max(AxisDirection::AxisY))));
      sl = makeSurfaceGridLookup2D<AxisBoundaryType::Bound,
                                   AxisBoundaryType::Bound>(
          std::move(surface), layerTolerance, pAxis1, pAxis2);
      break;
    }
    default: {
      throw std::invalid_argument(
          "SurfaceArrayCreator::"
          "surfaceArrayOnPlane: Invalid binning "
          "direction");
    }
  }

  sl->fill(gctx, surfacesRaw);

  return std::make_unique<SurfaceArray>(std::move(sl), std::move(surfaces),
                                        fullTransform);
  //!< @todo implement - take from ATLAS complex TRT builder
}

std::vector<const Surface*> SurfaceArrayCreator::findKeySurfaces(
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

std::size_t SurfaceArrayCreator::determineBinCount(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    AxisDirection aDir) const {
  auto matcher = m_cfg.surfaceMatcher;
  auto equal = [&gctx, &aDir, &matcher](const Surface* a, const Surface* b) {
    return matcher(gctx, aDir, a, b);
  };
  std::vector<const Surface*> keys = findKeySurfaces(surfaces, equal);

  return keys.size();
}

SurfaceArrayCreator::ProtoAxis SurfaceArrayCreator::createVariableAxis(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    AxisDirection aDir, const ProtoLayer& protoLayer,
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
  auto equal = [&gctx, &aDir, &matcher](const Surface* a, const Surface* b) {
    return matcher(gctx, aDir, a, b);
  };
  std::vector<const Surface*> keys = findKeySurfaces(surfaces, equal);

  std::vector<AxisScalar> aDirs;
  if (aDir == AxisDirection::AxisPhi) {
    std::stable_sort(
        keys.begin(), keys.end(), [&gctx](const Surface* a, const Surface* b) {
          return (phi(a->referencePosition(gctx, AxisDirection::AxisPhi)) <
                  phi(b->referencePosition(gctx, AxisDirection::AxisPhi)));
        });

    AxisScalar maxPhi =
        0.5 *
        (phi(keys.at(0)->referencePosition(gctx, AxisDirection::AxisPhi)) +
         phi(keys.at(1)->referencePosition(gctx, AxisDirection::AxisPhi)));

    // create rotation, so that maxPhi is +pi
    AxisScalar angle = -(std::numbers::pi + maxPhi);
    transform = transform * AngleAxis3(angle, Vector3::UnitZ());

    // iterate over all key surfaces, and use their mean position as aDirs,
    // but
    // rotate using transform from before
    AxisScalar previous =
        phi(keys.at(0)->referencePosition(gctx, AxisDirection::AxisPhi));
    // go through key surfaces
    for (std::size_t i = 1; i < keys.size(); i++) {
      const Surface* surface = keys.at(i);
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      AxisScalar edge = 0.5 * (previous + phi(surface->referencePosition(
                                              gctx, AxisDirection::AxisPhi))) +
                        angle;
      aDirs.push_back(edge);
      previous = phi(surface->referencePosition(gctx, AxisDirection::AxisPhi));
    }

    // segments
    unsigned int segments = 72;

    // get the bounds of the last surfaces
    const Surface* backSurface = keys.back();
    const PlanarBounds* backBounds =
        dynamic_cast<const PlanarBounds*>(&(backSurface->bounds()));
    if (backBounds == nullptr) {
      ACTS_ERROR(
          "Given SurfaceBounds are not planar - not implemented for "
          "other bounds yet! ");
    }
    // get the global vertices
    std::vector<Vector3> backVertices =
        makeGlobalVertices(gctx, *backSurface, backBounds->vertices(segments));
    AxisScalar maxBValue = phi(*std::max_element(
        backVertices.begin(), backVertices.end(),
        [](const Vector3& a, const Vector3& b) { return phi(a) < phi(b); }));

    aDirs.push_back(maxBValue);

    aDirs.push_back(std::numbers::pi_v<AxisScalar>);

  } else if (aDir == AxisDirection::AxisZ) {
    std::stable_sort(
        keys.begin(), keys.end(), [&gctx](const Surface* a, const Surface* b) {
          return (a->referencePosition(gctx, AxisDirection::AxisZ).z() <
                  b->referencePosition(gctx, AxisDirection::AxisZ).z());
        });

    aDirs.push_back(protoLayer.min(AxisDirection::AxisZ));
    aDirs.push_back(protoLayer.max(AxisDirection::AxisZ));

    // the z-center position of the previous surface
    AxisScalar previous =
        keys.front()->referencePosition(gctx, AxisDirection::AxisZ).z();
    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      aDirs.push_back(
          0.5 *
          (previous +
           (*surface)->referencePosition(gctx, AxisDirection::AxisZ).z()));
      previous = (*surface)->referencePosition(gctx, AxisDirection::AxisZ).z();
    }
  } else {  // AxisDirection::AxisR
    std::stable_sort(
        keys.begin(), keys.end(), [&gctx](const Surface* a, const Surface* b) {
          return (perp(a->referencePosition(gctx, AxisDirection::AxisR)) <
                  perp(b->referencePosition(gctx, AxisDirection::AxisR)));
        });

    aDirs.push_back(protoLayer.min(AxisDirection::AxisR));
    aDirs.push_back(protoLayer.max(AxisDirection::AxisR));

    // the r-center position of the previous surface
    AxisScalar previous =
        perp(keys.front()->referencePosition(gctx, AxisDirection::AxisR));

    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      aDirs.push_back(0.5 * (previous + perp((*surface)->referencePosition(
                                            gctx, AxisDirection::AxisR))));
      previous =
          perp((*surface)->referencePosition(gctx, AxisDirection::AxisR));
    }
  }
  std::ranges::sort(aDirs);
  ACTS_VERBOSE("Create variable binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	AxisDirection: " << aDir);
  ACTS_VERBOSE("	Number of bins: " << (aDirs.size() - 1));
  ACTS_VERBOSE("	(Min/Max) = (" << aDirs.front() << "/" << aDirs.back()
                                       << ")");

  ProtoAxis pAxis;
  pAxis.bType = arbitrary;
  pAxis.axisDir = aDir;
  pAxis.binEdges = aDirs;
  pAxis.nBins = aDirs.size() - 1;

  return pAxis;
}

SurfaceArrayCreator::ProtoAxis SurfaceArrayCreator::createEquidistantAxis(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    AxisDirection aDir, const ProtoLayer& protoLayer, Transform3& transform,
    std::size_t nBins) const {
  if (surfaces.empty()) {
    throw std::logic_error(
        "No surfaces handed over for creating equidistant axis!");
  }
  // check the binning type first

  double minimum = protoLayer.min(aDir, false);
  double maximum = protoLayer.max(aDir, false);

  std::size_t binNumber = 0;
  if (nBins == 0) {
    // determine bin count
    binNumber = determineBinCount(gctx, surfaces, aDir);
  } else {
    // use bin count
    binNumber = nBins;
  }

  // bind matcher & context with binning type
  auto matcher = m_cfg.surfaceMatcher;

  // now check the binning value
  if (aDir == AxisDirection::AxisPhi) {
    minimum = protoLayer.min(AxisDirection::AxisPhi, true);
    maximum = protoLayer.max(AxisDirection::AxisPhi, true);

    if (m_cfg.doPhiBinningOptimization) {
      minimum = -std::numbers::pi;
      maximum = std::numbers::pi;

      // Phi binning
      // set the binning option for phi
      // sort first in phi
      const Surface* maxElem = *std::max_element(
          surfaces.begin(), surfaces.end(),
          [&gctx](const Surface* a, const Surface* b) {
            return phi(a->referencePosition(gctx, AxisDirection::AxisR)) <
                   phi(b->referencePosition(gctx, AxisDirection::AxisR));
          });

      // rotate to max phi module plus one half step
      // this should make sure that phi wrapping at +- pi
      // never falls on a module center
      double surfaceMax =
          phi(maxElem->referencePosition(gctx, AxisDirection::AxisR));
      double gridStep = 2 * std::numbers::pi / binNumber;
      double gridMax = std::numbers::pi - 0.5 * gridStep;
      double angle = gridMax - surfaceMax;

      // replace given transform ref
      transform = transform * AngleAxis3(angle, Vector3::UnitZ());
    }
  }

  // assign the bin size
  ACTS_VERBOSE("Create equidistant binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	AxisDirection: " << aDir);
  ACTS_VERBOSE("	Number of bins: " << binNumber);
  ACTS_VERBOSE("	(Min/Max) = (" << minimum << "/" << maximum << ")");

  ProtoAxis pAxis;
  pAxis.max = maximum;
  pAxis.min = minimum;
  pAxis.bType = equidistant;
  pAxis.axisDir = aDir;
  pAxis.nBins = binNumber;

  return pAxis;
}

std::vector<Vector3> SurfaceArrayCreator::makeGlobalVertices(
    const GeometryContext& gctx, const Surface& surface,
    const std::vector<Vector2>& locVertices) const {
  std::vector<Vector3> globVertices;
  for (auto& vertex : locVertices) {
    Vector3 globVertex = surface.localToGlobal(gctx, vertex, Vector3());
    globVertices.push_back(globVertex);
  }
  return globVertices;
}

}  // namespace Acts
