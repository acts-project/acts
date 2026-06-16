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

SurfaceArray SurfaceArrayCreator::surfaceArrayOnCylinder(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, std::size_t binsPhi,
    std::size_t binsZ, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform, std::uint8_t maxNeighborDistance) const {
  using enum AxisDirection;

  const std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // Check if we have proto layer, else build it
  const ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with phi x z  = " << binsPhi << " x " << binsZ << " = "
                                      << binsPhi * binsZ << " bins.");

  Transform3 fullTransform = transform;
  const auto pAxisPhi =
      createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                            AxisPhi, protoLayer, fullTransform, binsPhi);
  const auto pAxisZ =
      createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound, AxisZ,
                            protoLayer, fullTransform, binsZ);

  const double R = protoLayer.medium(AxisR, true);
  const double halfZ = protoLayer.range(AxisZ, true) * 0.5;
  const double layerTolerance = protoLayer.range(AxisR) * 0.5;

  auto surface = Surface::makeShared<CylinderSurface>(fullTransform, R, halfZ);
  ACTS_VERBOSE("- projection surface is: " << surface->toString(gctx));

  return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                      layerTolerance, {*pAxisPhi, *pAxisZ},
                      maxNeighborDistance);
}

SurfaceArray SurfaceArrayCreator::surfaceArrayOnCylinder(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypePhi,
    BinningType bTypeZ, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform, std::uint8_t maxNeighborDistance) const {
  using enum AxisDirection;

  const std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  const ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  const double R = protoLayer.medium(AxisR, true);
  const double halfZ = protoLayer.range(AxisZ, true) * 0.5;
  const double layerTolerance = protoLayer.range(AxisR) * 0.5;

  std::unique_ptr<const IAxis> pAxisPhi;
  std::unique_ptr<const IAxis> pAxisZ;

  Transform3 fullTransform = transform;

  if (bTypePhi == equidistant) {
    pAxisPhi =
        createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                              AxisPhi, protoLayer, fullTransform, 0);
  } else {
    pAxisPhi = createVariableAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                                  AxisPhi, protoLayer, fullTransform);
  }

  if (bTypeZ == equidistant) {
    pAxisZ = createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                   AxisZ, protoLayer, fullTransform);
  } else {
    pAxisZ = createVariableAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisZ, protoLayer, fullTransform);
  }

  auto surface = Surface::makeShared<CylinderSurface>(fullTransform, R, halfZ);

  const std::size_t bins0 = pAxisPhi->getNBins();
  const std::size_t bins1 = pAxisZ->getNBins();
  ACTS_VERBOSE("Creating a SurfaceArray on a cylinder");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with phi x z  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                      layerTolerance, {*pAxisPhi, *pAxisZ},
                      maxNeighborDistance);
}

SurfaceArray SurfaceArrayCreator::surfaceArrayOnDisc(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, std::size_t binsR,
    std::size_t binsPhi, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform, std::uint8_t maxNeighborDistance) const {
  using enum AxisDirection;

  const std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  const ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  Transform3 fullTransform = transform;
  const auto pAxisR =
      createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound, AxisR,
                            protoLayer, fullTransform, binsR);
  const auto pAxisPhi =
      createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                            AxisPhi, protoLayer, fullTransform, binsPhi);

  const double Z = protoLayer.medium(AxisZ, true);
  const double Rmin = protoLayer.min(AxisR, true);
  const double Rmax = protoLayer.max(AxisR, true);
  const double layerThickness = protoLayer.range(AxisZ) * 0.5;
  ACTS_VERBOSE("- z-position of disc estimated as " << Z);
  ACTS_VERBOSE("- full transform is \n" << fullTransform.matrix());

  if (fullTransform.translation().norm() < s_transformEquivalentTolerance) {
    ACTS_VERBOSE(
        "input transform does not have translation: putting projection surface "
        "at center of gravity in z");
    fullTransform.translate(Vector3::UnitZ() * Z);
  }

  auto surface = Surface::makeShared<DiscSurface>(fullTransform, Rmin, Rmax);
  ACTS_VERBOSE("- projection surface is: " << surface->toString(gctx));

  const std::size_t bins0 = pAxisR->getNBins();
  const std::size_t bins1 = pAxisPhi->getNBins();
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                      layerThickness, {*pAxisR, *pAxisPhi},
                      maxNeighborDistance);
}

SurfaceArray SurfaceArrayCreator::surfaceArrayOnDisc(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypeR,
    BinningType bTypePhi, std::optional<ProtoLayer> protoLayerOpt,
    const Transform3& transform, std::uint8_t maxNeighborDistance) const {
  using enum AxisDirection;

  const std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  const ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a disc");

  std::unique_ptr<const IAxis> pAxisPhi;
  std::unique_ptr<const IAxis> pAxisR;

  Transform3 fullTransform = transform;
  Transform3 inverseTransform = transform.inverse();

  if (bTypeR == equidistant) {
    pAxisR = createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                   AxisR, protoLayer, fullTransform);
  } else {
    pAxisR = createVariableAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisR, protoLayer, fullTransform);
  }

  // if we have more than one R ring, we need to figure out
  // the number of phi bins.
  if (pAxisR->getNBins() > 1) {
    // more than one R-Ring, we need to adjust
    // this FORCES equidistant binning
    std::vector<std::vector<const Surface*>> phiModules(pAxisR->getNBins());
    for (const auto& srf : surfacesRaw) {
      const Vector3 bpos =
          inverseTransform * srf->referencePosition(gctx, AxisR);
      const std::size_t bin =
          pAxisR->getBin(perp(bpos)) - 1;  // subtract underflow bin
      phiModules.at(bin).push_back(srf);
    }

    std::vector<std::size_t> nPhiModules;
    const auto& matcher = m_cfg.surfaceMatcher;
    const auto equal = [&gctx, &matcher](const Surface* a, const Surface* b) {
      return matcher(gctx, AxisPhi, a, b);
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
    const std::size_t nBinsPhi =
        *std::min_element(nPhiModules.begin(), nPhiModules.end());
    pAxisPhi =
        createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                              AxisPhi, protoLayer, fullTransform, nBinsPhi);

  } else {
    // use regular determination
    if (bTypePhi == equidistant) {
      pAxisPhi =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                                AxisPhi, protoLayer, fullTransform, 0);
    } else {
      pAxisPhi = createVariableAxis(gctx, surfacesRaw, AxisBoundaryType::Closed,
                                    AxisPhi, protoLayer, fullTransform);
    }
  }

  const double Z = protoLayer.medium(AxisZ, true);
  const double Rmin = protoLayer.min(AxisR, true);
  const double Rmax = protoLayer.max(AxisR, true);
  const double layerThickness = protoLayer.range(AxisZ) * 0.5;
  ACTS_VERBOSE("- z-position of disc estimated as " << Z);

  if (fullTransform.translation().norm() < s_transformEquivalentTolerance) {
    ACTS_VERBOSE(
        "input transform does not have translation: putting projection surface "
        "at center of gravity in z");
    fullTransform.translate(Vector3::UnitZ() * Z);
  }

  auto surface = Surface::makeShared<DiscSurface>(fullTransform, Rmin, Rmax);

  // get the number of bins
  const std::size_t bins0 = pAxisR->getNBins();
  const std::size_t bins1 = pAxisPhi->getNBins();
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with r x phi  = " << bins0 << " x " << bins1 << " = "
                                      << bins0 * bins1 << " bins.");

  return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                      layerThickness, {*pAxisR, *pAxisPhi},
                      maxNeighborDistance);
}

/// SurfaceArrayCreator interface method - create an array on a plane
SurfaceArray SurfaceArrayCreator::surfaceArrayOnPlane(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, std::size_t bins1,
    std::size_t bins2, AxisDirection aDir,
    std::optional<ProtoLayer> protoLayerOpt, const Transform3& transform,
    std::uint8_t maxNeighborDistance) const {
  using enum AxisDirection;

  const std::vector<const Surface*> surfacesRaw = unpackSmartPointers(surfaces);
  // check if we have proto layer, else build it
  const ProtoLayer protoLayer =
      protoLayerOpt ? *protoLayerOpt : ProtoLayer(gctx, surfacesRaw);

  ACTS_VERBOSE("Creating a SurfaceArray on a plance");
  ACTS_VERBOSE(" -- with " << surfaces.size() << " surfaces.");
  ACTS_VERBOSE(" -- with " << bins1 << " x " << bins2 << " = " << bins1 * bins2
                           << " bins.");
  Transform3 fullTransform = transform;

  const double layerTolerance = protoLayer.range(aDir) * 0.5;

  // Axis along the binning
  switch (aDir) {
    case AxisX: {
      const auto pAxis1 =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisY, protoLayer, fullTransform, bins1);
      const auto pAxis2 =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisZ, protoLayer, fullTransform, bins2);
      auto surface = Surface::makeShared<PlaneSurface>(
          fullTransform,
          std::make_shared<RectangleBounds>(
              Vector2(protoLayer.min(AxisY), protoLayer.min(AxisZ)),
              Vector2(protoLayer.max(AxisY), protoLayer.max(AxisZ))));
      return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                          layerTolerance, {*pAxis1, *pAxis2},
                          maxNeighborDistance);
    }
    case AxisY: {
      const auto pAxis1 =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisX, protoLayer, fullTransform, bins1);
      const auto pAxis2 =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisZ, protoLayer, fullTransform, bins2);
      auto surface = Surface::makeShared<PlaneSurface>(
          fullTransform,
          std::make_shared<RectangleBounds>(
              Vector2(protoLayer.min(AxisX), protoLayer.min(AxisY)),
              Vector2(protoLayer.max(AxisX), protoLayer.max(AxisY))));
      return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                          layerTolerance, {*pAxis1, *pAxis2},
                          maxNeighborDistance);
    }
    case AxisZ: {
      const auto pAxis1 =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisX, protoLayer, fullTransform, bins1);
      const auto pAxis2 =
          createEquidistantAxis(gctx, surfacesRaw, AxisBoundaryType::Bound,
                                AxisY, protoLayer, fullTransform, bins2);
      auto surface = Surface::makeShared<PlaneSurface>(
          fullTransform,
          std::make_shared<RectangleBounds>(
              Vector2(protoLayer.min(AxisX), protoLayer.min(AxisY)),
              Vector2(protoLayer.max(AxisX), protoLayer.max(AxisY))));
      return SurfaceArray(gctx, std::move(surfaces), std::move(surface),
                          layerTolerance, {*pAxis1, *pAxis2},
                          maxNeighborDistance);
    }
    default:
      break;
  }

  throw std::invalid_argument(
      "SurfaceArrayCreator::surfaceArrayOnPlane: Invalid binning direction");
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
  const auto& matcher = m_cfg.surfaceMatcher;
  const auto equal = [&gctx, &aDir, &matcher](const Surface* a,
                                              const Surface* b) {
    return matcher(gctx, aDir, a, b);
  };
  const std::vector<const Surface*> keys = findKeySurfaces(surfaces, equal);

  return keys.size();
}

std::unique_ptr<const IAxis> SurfaceArrayCreator::createVariableAxis(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    AxisBoundaryType aBoundaryType, AxisDirection aDir,
    const ProtoLayer& protoLayer, Transform3& transform) const {
  using enum AxisDirection;

  if (surfaces.empty()) {
    throw std::logic_error(
        "No surfaces handed over for creating arbitrary bin utility!");
  }
  // BinningOption is open for z and r, in case of phi binning reset later
  // the vector with the binning Values (boundaries for each bin)

  // bind matcher with binning type
  const auto& matcher = m_cfg.surfaceMatcher;
  // find the key surfaces
  const auto equal = [&gctx, &aDir, &matcher](const Surface* a,
                                              const Surface* b) {
    return matcher(gctx, aDir, a, b);
  };
  std::vector<const Surface*> keys = findKeySurfaces(surfaces, equal);

  std::vector<double> binEdges;
  if (aDir == AxisPhi) {
    std::stable_sort(keys.begin(), keys.end(),
                     [&gctx](const Surface* a, const Surface* b) {
                       return (phi(a->referencePosition(gctx, AxisPhi)) <
                               phi(b->referencePosition(gctx, AxisPhi)));
                     });

    const double maxPhi =
        0.5 * (phi(keys.at(0)->referencePosition(gctx, AxisPhi)) +
               phi(keys.at(1)->referencePosition(gctx, AxisPhi)));

    // create rotation, so that maxPhi is +pi
    const double angle = -(std::numbers::pi + maxPhi);
    transform = transform * AngleAxis3(angle, Vector3::UnitZ());

    // iterate over all key surfaces, and use their mean position as aDirs,
    // but
    // rotate using transform from before
    double previous = phi(keys.at(0)->referencePosition(gctx, AxisPhi));
    // go through key surfaces
    for (std::size_t i = 1; i < keys.size(); i++) {
      const Surface* surface = keys.at(i);
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      const double edge =
          0.5 * (previous + phi(surface->referencePosition(gctx, AxisPhi))) +
          angle;
      binEdges.push_back(edge);
      previous = phi(surface->referencePosition(gctx, AxisPhi));
    }

    // segments
    constexpr unsigned int segments = 72;

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
    const std::vector<Vector3> backVertices =
        makeGlobalVertices(gctx, *backSurface, backBounds->vertices(segments));
    const double maxBValue = phi(*std::ranges::max_element(
        backVertices, {}, [](const Vector3& v) { return phi(v); }));

    binEdges.push_back(maxBValue);

    binEdges.push_back(std::numbers::pi);

  } else if (aDir == AxisZ) {
    std::ranges::stable_sort(keys, {}, [&gctx](const Surface* s) {
      return s->referencePosition(gctx, AxisZ).z();
    });

    binEdges.push_back(protoLayer.min(AxisZ));
    binEdges.push_back(protoLayer.max(AxisZ));

    // the z-center position of the previous surface
    double previous = keys.front()->referencePosition(gctx, AxisZ).z();
    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      binEdges.push_back(
          0.5 * (previous + (*surface)->referencePosition(gctx, AxisZ).z()));
      previous = (*surface)->referencePosition(gctx, AxisZ).z();
    }
  } else {  // AxisR
    std::ranges::stable_sort(keys, {}, [&gctx](const Surface* s) {
      return perp(s->referencePosition(gctx, AxisR));
    });

    binEdges.push_back(protoLayer.min(AxisR));
    binEdges.push_back(protoLayer.max(AxisR));

    // the r-center position of the previous surface
    double previous = perp(keys.front()->referencePosition(gctx, AxisR));

    // go through key surfaces
    for (auto surface = keys.begin() + 1; surface != keys.end(); surface++) {
      // create central binning values which is the mean of the center
      // positions in the binning direction of the current and previous
      // surface
      binEdges.push_back(
          0.5 * (previous + perp((*surface)->referencePosition(gctx, AxisR))));
      previous = perp((*surface)->referencePosition(gctx, AxisR));
    }
  }
  std::ranges::sort(binEdges);
  ACTS_VERBOSE("Create variable binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	AxisDirection: " << aDir);
  ACTS_VERBOSE("	Number of bins: " << (binEdges.size() - 1));
  ACTS_VERBOSE("	(Min/Max) = (" << binEdges.front() << "/"
                                       << binEdges.back() << ")");

  return IAxis::createVariable(aBoundaryType, binEdges, aDir);
}

std::unique_ptr<const IAxis> SurfaceArrayCreator::createEquidistantAxis(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    AxisBoundaryType aBoundaryType, AxisDirection aDir,
    const ProtoLayer& protoLayer, Transform3& transform,
    std::size_t nBins) const {
  using enum AxisDirection;

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
  if (aDir == AxisPhi) {
    minimum = protoLayer.min(AxisPhi, true);
    maximum = protoLayer.max(AxisPhi, true);

    if (m_cfg.doPhiBinningOptimization) {
      minimum = -std::numbers::pi;
      maximum = std::numbers::pi;

      // Phi binning
      // set the binning option for phi
      // sort first in phi
      const Surface* maxElem =
          *std::ranges::max_element(surfaces, {}, [&gctx](const Surface* s) {
            return phi(s->referencePosition(gctx, AxisR));
          });

      // rotate to max phi module plus one half step
      // this should make sure that phi wrapping at +- pi
      // never falls on a module center
      const double surfaceMax = phi(maxElem->referencePosition(gctx, AxisR));
      const double gridStep =
          2 * std::numbers::pi / static_cast<double>(binNumber);
      const double gridMax = std::numbers::pi - 0.5 * gridStep;
      const double angle = gridMax - surfaceMax;

      // replace given transform ref
      transform = transform * AngleAxis3(angle, Vector3::UnitZ());
    }
  }

  // assign the bin size
  ACTS_VERBOSE("Create equidistant binning Axis for binned SurfaceArray");
  ACTS_VERBOSE("	AxisDirection: " << aDir);
  ACTS_VERBOSE("	Number of bins: " << binNumber);
  ACTS_VERBOSE("	(Min/Max) = (" << minimum << "/" << maximum << ")");

  return IAxis::createEquidistant(aBoundaryType, minimum, maximum, binNumber,
                                  aDir);
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
