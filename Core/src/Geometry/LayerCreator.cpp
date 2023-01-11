// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/LayerCreator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <iterator>
#include <set>
#include <utility>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::LayerCreator::LayerCreator(const Acts::LayerCreator::Config& lcConfig,
                                 std::unique_ptr<const Logger> logger)
    : m_cfg(lcConfig), m_logger(std::move(logger)) {}

void Acts::LayerCreator::setConfiguration(
    const Acts::LayerCreator::Config& lcConfig) {
  // @todo check consistency
  // copy the configuration
  m_cfg = lcConfig;
}

void Acts::LayerCreator::setLogger(std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

Acts::MutableLayerPtr Acts::LayerCreator::cylinderLayer(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, size_t binsPhi,
    size_t binsZ, std::optional<ProtoLayer> _protoLayer,
    const Transform3& transform, std::unique_ptr<ApproachDescriptor> ad) const {
  ProtoLayer protoLayer =
      _protoLayer ? *_protoLayer : ProtoLayer(gctx, surfaces);

  // Remaining layer parameters - they include the envelopes
  double layerR = protoLayer.medium(binR);
  double layerZ = protoLayer.medium(binZ);
  double layerHalfZ = 0.5 * protoLayer.range(binZ);
  double layerThickness = protoLayer.range(binR);

  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R     = " << layerR);
  ACTS_VERBOSE(" - from R min/max   = " << protoLayer.min(binR, false) << " / "
                                        << protoLayer.max(binR, false));
  ACTS_VERBOSE(" - with R thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envelope[binR][0u]
                                        << " / "
                                        << protoLayer.envelope[binR][1u]);

  ACTS_VERBOSE(" - with z min/max   = "
               << protoLayer.min(binZ, false) << " (-"
               << protoLayer.envelope[binZ][0u] << ") / "
               << protoLayer.max(binZ, false) << " (+"
               << protoLayer.envelope[binZ][1u] << ")");

  ACTS_VERBOSE(" - z center         = " << layerZ);
  ACTS_VERBOSE(" - halflength z     = " << layerHalfZ);

  // create the layer transforms if not given
  // we need to transform in case layerZ != 0, so that the layer will be
  // correctly defined using the halflength
  Translation3 addTranslation(0., 0., 0.);
  if (transform.isApprox(Transform3::Identity())) {
    // double shift = -(layerZ + envZShift);
    addTranslation = Translation3(0., 0., layerZ);
    ACTS_VERBOSE(" - layer z shift  = " << -layerZ);
  }

  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.min(binPhi, false)
                                        << " / "
                                        << protoLayer.max(binPhi, false));
  ACTS_VERBOSE(" - # of modules     = " << surfaces.size() << " ordered in ( "
                                        << binsPhi << " x " << binsZ << ")");
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
        gctx, std::move(surfaces), binsPhi, binsZ, protoLayer);

    checkBinning(gctx, *sArray);
  }

  // create the layer and push it back
  std::shared_ptr<const CylinderBounds> cBounds(
      new CylinderBounds(layerR, layerHalfZ));

  // create the layer
  MutableLayerPtr cLayer = CylinderLayer::create(
      addTranslation * transform, cBounds, std::move(sArray), layerThickness,
      std::move(ad), active);

  if (!cLayer) {
    ACTS_ERROR("Creation of cylinder layer did not succeed!");
  }
  associateSurfacesToLayer(*cLayer);

  // now return
  return cLayer;
}

Acts::MutableLayerPtr Acts::LayerCreator::cylinderLayer(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypePhi,
    BinningType bTypeZ, std::optional<ProtoLayer> _protoLayer,
    const Transform3& transform, std::unique_ptr<ApproachDescriptor> ad) const {
  ProtoLayer protoLayer =
      _protoLayer ? *_protoLayer : ProtoLayer(gctx, surfaces);

  // remaining layer parameters
  double layerR = protoLayer.medium(binR);
  double layerZ = protoLayer.medium(binZ);
  double layerHalfZ = 0.5 * protoLayer.range(binZ);
  double layerThickness = protoLayer.range(binR);

  // adjust the layer radius
  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R     = " << layerR);
  ACTS_VERBOSE(" - from R min/max   = " << protoLayer.min(binR, false) << " / "
                                        << protoLayer.max(binR, false));
  ACTS_VERBOSE(" - with R thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envelope[binR][0u]
                                        << " / "
                                        << protoLayer.envelope[binR][1u]);
  ACTS_VERBOSE(" - with z min/max   = "
               << protoLayer.min(binZ, false) << " (-"
               << protoLayer.envelope[binZ][0u] << ") / "
               << protoLayer.max(binZ, false) << " (+"
               << protoLayer.envelope[binZ][1u] << ")");
  ACTS_VERBOSE(" - z center         = " << layerZ);
  ACTS_VERBOSE(" - halflength z     = " << layerHalfZ);

  // create the layer transforms if not given
  // we need to transform in case layerZ != 0, so that the layer will be
  // correctly defined using the halflength
  // create the layer transforms if not given
  Translation3 addTranslation(0., 0., 0.);
  if (transform.isApprox(Transform3::Identity()) && bTypeZ == equidistant) {
    addTranslation = Translation3(0., 0., layerZ);
    ACTS_VERBOSE(" - layer z shift    = " << -layerZ);
  }

  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.min(binPhi, false)
                                        << " / "
                                        << protoLayer.max(binPhi, false));
  ACTS_VERBOSE(" - # of modules     = " << surfaces.size() << "");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
        gctx, std::move(surfaces), bTypePhi, bTypeZ, protoLayer);

    checkBinning(gctx, *sArray);
  }

  // create the layer and push it back
  std::shared_ptr<const CylinderBounds> cBounds(
      new CylinderBounds(layerR, layerHalfZ));

  // create the layer
  MutableLayerPtr cLayer = CylinderLayer::create(
      addTranslation * transform, cBounds, std::move(sArray), layerThickness,
      std::move(ad), active);

  if (!cLayer) {
    ACTS_ERROR("Creation of cylinder layer did not succeed!");
  }
  associateSurfacesToLayer(*cLayer);

  // now return
  return cLayer;
}

Acts::MutableLayerPtr Acts::LayerCreator::discLayer(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, size_t binsR,
    size_t binsPhi, std::optional<ProtoLayer> _protoLayer,
    const Transform3& transform, std::unique_ptr<ApproachDescriptor> ad) const {
  ProtoLayer protoLayer =
      _protoLayer ? *_protoLayer : ProtoLayer(gctx, surfaces);

  double layerZ = protoLayer.medium(binZ);
  double layerThickness = protoLayer.range(binZ);

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position    = " << layerZ);
  ACTS_VERBOSE(" - from Z min/max   = " << protoLayer.min(binZ, false) << " / "
                                        << protoLayer.max(binZ, false));
  ACTS_VERBOSE(" - with Z thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envelope[binZ][0u]
                                        << " / "
                                        << protoLayer.envelope[binZ][1u]);
  ACTS_VERBOSE(" - with R min/max   = "
               << protoLayer.min(binR, false) << " (-"
               << protoLayer.envelope[binR][0u] << ") / "
               << protoLayer.max(binR, false) << " (+"
               << protoLayer.envelope[binR][1u] << ")");
  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.min(binPhi, false)
                                        << " / "
                                        << protoLayer.max(binPhi, false));
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << " ordered in ( "
                                       << binsR << " x " << binsPhi << ")");

  // create the layer transforms if not given
  Translation3 addTranslation(0., 0., 0.);
  if (transform.isApprox(Transform3::Identity())) {
    addTranslation = Translation3(0., 0., layerZ);
  }
  // create the surface array
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
        gctx, std::move(surfaces), binsR, binsPhi, protoLayer, transform);

    checkBinning(gctx, *sArray);
  }

  // create the share disc bounds
  auto dBounds = std::make_shared<const RadialBounds>(protoLayer.min(binR),
                                                      protoLayer.max(binR));

  // create the layers
  // we use the same transform here as for the layer itself
  // for disk this is fine since we don't bin in Z, so does not matter
  MutableLayerPtr dLayer =
      DiscLayer::create(addTranslation * transform, dBounds, std::move(sArray),
                        layerThickness, std::move(ad), active);

  if (!dLayer) {
    ACTS_ERROR("Creation of disc layer did not succeed!");
  }
  associateSurfacesToLayer(*dLayer);
  // return the layer
  return dLayer;
}

Acts::MutableLayerPtr Acts::LayerCreator::discLayer(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, BinningType bTypeR,
    BinningType bTypePhi, std::optional<ProtoLayer> _protoLayer,
    const Transform3& transform, std::unique_ptr<ApproachDescriptor> ad) const {
  ProtoLayer protoLayer =
      _protoLayer ? *_protoLayer : ProtoLayer(gctx, surfaces);

  double layerZ = protoLayer.medium(binZ);
  double layerThickness = protoLayer.range(binZ);

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position    = " << layerZ);
  ACTS_VERBOSE(" - from Z min/max   = " << protoLayer.min(binZ, false) << " / "
                                        << protoLayer.max(binZ, false));
  ACTS_VERBOSE(" - with Z thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envelope[binZ][0u]
                                        << " / "
                                        << protoLayer.envelope[binZ][1u]);
  ACTS_VERBOSE(" - with R min/max   = "
               << protoLayer.min(binR, false) << " (-"
               << protoLayer.envelope[binR][0u] << ") / "
               << protoLayer.max(binR, false) << " (+"
               << protoLayer.envelope[binR][1u] << ")");
  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.min(binPhi, false)
                                        << " / "
                                        << protoLayer.max(binPhi, false));
  ACTS_VERBOSE(" - # of modules     = " << surfaces.size());

  // create the layer transforms if not given
  Translation3 addTranslation(0., 0., 0.);
  if (transform.isApprox(Transform3::Identity())) {
    addTranslation = Translation3(0., 0., layerZ);
  }

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
        gctx, std::move(surfaces), bTypeR, bTypePhi, protoLayer, transform);

    checkBinning(gctx, *sArray);
  }

  // create the shared disc bounds
  auto dBounds = std::make_shared<const RadialBounds>(protoLayer.min(binR),
                                                      protoLayer.max(binR));

  // create the layers
  MutableLayerPtr dLayer =
      DiscLayer::create(addTranslation * transform, dBounds, std::move(sArray),
                        layerThickness, std::move(ad), active);
  if (!dLayer) {
    ACTS_ERROR("Creation of disc layer did not succeed!");
  }
  associateSurfacesToLayer(*dLayer);
  // return the layer
  return dLayer;
}

Acts::MutableLayerPtr Acts::LayerCreator::planeLayer(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Surface>> surfaces, size_t bins1,
    size_t bins2, BinningValue bValue, std::optional<ProtoLayer> _protoLayer,
    const Transform3& transform, std::unique_ptr<ApproachDescriptor> ad) const {
  ProtoLayer protoLayer =
      _protoLayer ? *_protoLayer : ProtoLayer(gctx, surfaces);

  // remaining layer parameters
  double layerHalf1 = 0, layerHalf2 = 0, layerThickness = 0;
  switch (bValue) {
    case BinningValue::binX: {
      layerHalf1 = 0.5 * (protoLayer.max(binY) - protoLayer.min(binY));
      layerHalf2 = 0.5 * (protoLayer.max(binZ) - protoLayer.min(binZ));
      layerThickness = (protoLayer.max(binX) - protoLayer.min(binX));
      break;
    }
    case BinningValue::binY: {
      layerHalf1 = 0.5 * (protoLayer.max(binX) - protoLayer.min(binX));
      layerHalf2 = 0.5 * (protoLayer.max(binZ) - protoLayer.min(binZ));
      layerThickness = (protoLayer.max(binY) - protoLayer.min(binY));
      break;
    }
    default: {
      layerHalf1 = 0.5 * (protoLayer.max(binX) - protoLayer.min(binX));
      layerHalf2 = 0.5 * (protoLayer.max(binY) - protoLayer.min(binY));
      layerThickness = (protoLayer.max(binZ) - protoLayer.min(binZ));
    }
  }

  double centerX = 0.5 * (protoLayer.max(binX) + protoLayer.min(binX));
  double centerY = 0.5 * (protoLayer.max(binY) + protoLayer.min(binY));
  double centerZ = 0.5 * (protoLayer.max(binZ) + protoLayer.min(binZ));

  ACTS_VERBOSE("Creating a plane Layer:");
  ACTS_VERBOSE(" - with layer center     = "
               << "(" << centerX << ", " << centerY << ", " << centerZ << ")");
  ACTS_VERBOSE(" - from X min/max   = " << protoLayer.min(binX) << " / "
                                        << protoLayer.max(binX));
  ACTS_VERBOSE(" - from Y min/max   = " << protoLayer.min(binY) << " / "
                                        << protoLayer.max(binY));
  ACTS_VERBOSE(" - with Z thickness = " << layerThickness);

  // create the layer transforms if not given
  // we need to transform in case centerX/centerY/centerZ != 0, so that the
  // layer will be correctly defined
  Translation3 addTranslation(0., 0., 0.);
  if (transform.isApprox(Transform3::Identity())) {
    // double shift = (layerZ + envZShift);
    addTranslation = Translation3(centerX, centerY, centerZ);
    ACTS_VERBOSE(" - layer shift  = "
                 << "(" << centerX << ", " << centerY << ", " << centerZ
                 << ")");
  }

  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnPlane(
        gctx, std::move(surfaces), bins1, bins2, bValue, protoLayer, transform);

    checkBinning(gctx, *sArray);
  }

  // create the layer and push it back
  std::shared_ptr<const PlanarBounds> pBounds(
      new RectangleBounds(layerHalf1, layerHalf2));

  // create the layer
  MutableLayerPtr pLayer =
      PlaneLayer::create(addTranslation * transform, pBounds, std::move(sArray),
                         layerThickness, std::move(ad), active);

  if (!pLayer) {
    ACTS_ERROR("Creation of plane layer did not succeed!");
  }
  associateSurfacesToLayer(*pLayer);

  // now return
  return pLayer;
}

void Acts::LayerCreator::associateSurfacesToLayer(Layer& layer) const {
  if (layer.surfaceArray() != nullptr) {
    auto surfaces = layer.surfaceArray()->surfaces();

    for (auto& surface : surfaces) {
      auto mutableSurface = const_cast<Surface*>(surface);
      mutableSurface->associateLayer(layer);
    }
  }
}

bool Acts::LayerCreator::checkBinning(const GeometryContext& gctx,
                                      const SurfaceArray& sArray) const {
  // do consistency check: can we access all sensitive surfaces
  // through the binning? If not, surfaces get lost and the binning does not
  // work

  ACTS_VERBOSE("Performing consistency check")

  std::vector<const Surface*> surfaces = sArray.surfaces();
  std::set<const Surface*> sensitiveSurfaces(surfaces.begin(), surfaces.end());
  std::set<const Surface*> accessibleSurfaces;
  size_t nEmptyBins = 0;
  size_t nBinsChecked = 0;

  // iterate over all bins
  size_t size = sArray.size();
  for (size_t b = 0; b < size; ++b) {
    std::vector<const Surface*> binContent = sArray.at(b);
    // we don't check under/overflow bins
    if (!sArray.isValidBin(b)) {
      continue;
    }
    for (const auto& srf : binContent) {
      accessibleSurfaces.insert(srf);
    }
    if (binContent.empty()) {
      nEmptyBins++;
    }
    nBinsChecked++;
  }

  std::vector<const Acts::Surface*> diff;
  std::set_difference(sensitiveSurfaces.begin(), sensitiveSurfaces.end(),
                      accessibleSurfaces.begin(), accessibleSurfaces.end(),
                      std::inserter(diff, diff.begin()));

  ACTS_VERBOSE(" - Checked " << nBinsChecked << " valid bins");

  if (nEmptyBins > 0) {
    ACTS_ERROR(" -- Not all bins point to surface. " << nEmptyBins << " empty");
  } else {
    ACTS_VERBOSE(" -- All bins point to a surface");
  }

  if (!diff.empty()) {
    ACTS_ERROR(
        " -- Not all sensitive surfaces are accessible through binning. "
        "sensitive: "
        << sensitiveSurfaces.size()
        << "    accessible: " << accessibleSurfaces.size());

    // print all inaccessibles
    ACTS_ERROR(" -- Inaccessible surfaces: ");
    for (const auto& srf : diff) {
      // have to choose BinningValue here
      Vector3 ctr = srf->binningPosition(gctx, binR);
      ACTS_ERROR(" Surface(x=" << ctr.x() << ", y=" << ctr.y()
                               << ", z=" << ctr.z() << ", r=" << perp(ctr)
                               << ", phi=" << phi(ctr) << ")");
    }

  } else {
    ACTS_VERBOSE(" -- All sensitive surfaces are accessible through binning.");
  }

  return nEmptyBins == 0 && diff.empty();
}
