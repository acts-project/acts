// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

using namespace Acts;

namespace ActsTests {

std::vector<Surface*> CylindricalTrackingGeometry::surfacesRing(
    DetectorStore& detStore, double moduleHalfXminY, double moduleHalfXmaxY,
    double moduleHalfY, double moduleThickness, double moduleTilt,
    double ringRadius, double ringZ, double zStagger, int nPhi) {
  std::vector<Surface*> layerSurfaces;

  // Module material from input
  MaterialSlab moduleMaterial(makeSilicon(), moduleThickness);

  // Create a new surface material
  std::shared_ptr<const ISurfaceMaterial> moduleMaterialPtr =
      std::shared_ptr<const ISurfaceMaterial>(
          new HomogeneousSurfaceMaterial(moduleMaterial));

  // The rectangle/trapezoid bounds for all modules
  std::shared_ptr<PlanarBounds> mBounds = nullptr;
  if (moduleHalfXminY == moduleHalfXmaxY) {
    mBounds = std::make_shared<RectangleBounds>(moduleHalfXminY, moduleHalfY);
  } else {
    mBounds = std::make_shared<TrapezoidBounds>(moduleHalfXminY,
                                                moduleHalfXmaxY, moduleHalfY);
  }

  double phiStep = 2 * std::numbers::pi / nPhi;

  for (int im = 0; im < nPhi; ++im) {
    // Get the moduleTransform
    double phi = -std::numbers::pi + im * phiStep;
    auto mModuleTransform = Transform3(
        Translation3(ringRadius * std::cos(phi), ringRadius * std::sin(phi),
                     ringZ + (im % 2) * zStagger) *
        AngleAxis3(phi - std::numbers::pi / 2., Vector3::UnitZ()) *
        AngleAxis3(moduleTilt, Vector3::UnitY()));

    // Create the detector element
    auto detElement = std::make_unique<DetectorElementStub>(
        mModuleTransform, mBounds, moduleThickness, moduleMaterialPtr);

    layerSurfaces.push_back(&detElement->surface());
    detStore.push_back(std::move(detElement));
  }

  return layerSurfaces;
}

std::vector<Surface*> CylindricalTrackingGeometry::surfacesCylinder(
    DetectorStore& detStore, double moduleHalfX, double moduleHalfY,
    double moduleThickness, double moduleTiltPhi, double layerRadius,
    double radialStagger, double longitudinalOverlap,
    const std::pair<int, int>& binningSchema) {
  std::vector<Surface*> layerSurfaces;

  // Module material from input
  MaterialSlab moduleMaterial(makeSilicon(), moduleThickness);

  // Create a new surface material
  std::shared_ptr<const ISurfaceMaterial> moduleMaterialPtr =
      std::shared_ptr<const ISurfaceMaterial>(
          new HomogeneousSurfaceMaterial(moduleMaterial));

  // The rectangle bounds for all modules
  auto mBounds = std::make_shared<RectangleBounds>(moduleHalfX, moduleHalfY);

  // Create the module centers
  auto moduleCenters =
      modulePositionsCylinder(layerRadius, radialStagger, moduleHalfY,
                              longitudinalOverlap, binningSchema);

  for (auto& mCenter : moduleCenters) {
    // The association transform
    double modulePhi = VectorHelpers::phi(mCenter);
    // Local z axis is the normal vector
    Vector3 moduleLocalZ(std::cos(modulePhi + moduleTiltPhi),
                         std::sin(modulePhi + moduleTiltPhi), 0.);
    // Local y axis is the global z axis
    Vector3 moduleLocalY(0., 0., 1);
    // Local x axis the normal to local y,z
    Vector3 moduleLocalX(-std::sin(modulePhi + moduleTiltPhi),
                         std::cos(modulePhi + moduleTiltPhi), 0.);
    // Create the RotationMatrix
    RotationMatrix3 moduleRotation;
    moduleRotation.col(0) = moduleLocalX;
    moduleRotation.col(1) = moduleLocalY;
    moduleRotation.col(2) = moduleLocalZ;
    // Get the moduleTransform
    auto mModuleTransform = Transform3(Translation3(mCenter) * moduleRotation);
    // Create the detector element
    auto detElement = std::make_unique<DetectorElementStub>(
        mModuleTransform, mBounds, moduleThickness, moduleMaterialPtr);

    layerSurfaces.push_back(&detElement->surface());
    detStore.push_back(std::move(detElement));
  }
  return layerSurfaces;
}

/// Helper method for cylinder layer
/// create the positions for module surfaces on a cylinder
std::vector<Vector3> CylindricalTrackingGeometry::modulePositionsCylinder(
    double radius, double zStagger, double moduleHalfLength, double lOverlap,
    const std::pair<int, int>& binningSchema) {
  int nPhiBins = binningSchema.first;
  int nZbins = binningSchema.second;
  // prepare the return value
  std::vector<Vector3> mPositions;
  mPositions.reserve(nPhiBins * nZbins);
  // prep work
  double phiStep = 2 * std::numbers::pi / (nPhiBins);
  double minPhi = -std::numbers::pi + phiStep / 2.;
  double zStart = -0.5 * (nZbins - 1) * (2 * moduleHalfLength - lOverlap);
  double zStep = 2 * std::abs(zStart) / (nZbins - 1);
  // loop over the bins
  for (std::size_t zBin = 0; zBin < static_cast<std::size_t>(nZbins); ++zBin) {
    // prepare z and r
    double moduleZ = zStart + zBin * zStep;
    double moduleR =
        (zBin % 2) != 0u ? radius - 0.5 * zStagger : radius + 0.5 * zStagger;
    for (std::size_t phiBin = 0; phiBin < static_cast<std::size_t>(nPhiBins);
         ++phiBin) {
      // calculate the current phi value
      double modulePhi = minPhi + phiBin * phiStep;
      mPositions.push_back(Vector3(moduleR * std::cos(modulePhi),
                                   moduleR * std::sin(modulePhi), moduleZ));
    }
  }
  return mPositions;
}

std::shared_ptr<TrackingGeometry> CylindricalTrackingGeometry::buildGen1(
    const Logger& /*logger*/) {
  using namespace UnitLiterals;

  Logging::Level surfaceLLevel = Logging::INFO;
  Logging::Level layerLLevel = Logging::INFO;
  Logging::Level volumeLLevel = Logging::INFO;

  // configure surface array creator
  auto surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>(
      getDefaultLogger("SurfaceArrayCreator", surfaceLLevel));
  // configure the layer creator that uses the surface array creator
  LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const LayerCreator>(
      lcConfig, getDefaultLogger("LayerCreator", layerLLevel));
  // configure the layer array creator
  LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
      lacConfig, getDefaultLogger("LayerArrayCreator", layerLLevel));

  // tracking volume array creator
  TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
      tvacConfig, getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));

  // ----------------- build a beam pipe -----------------------------------
  MaterialSlab beamPipeMaterial(makeBeryllium(), 0.8_mm);
  PassiveLayerBuilder::Config bplConfig;
  bplConfig.layerIdentification = "BeamPipe";
  bplConfig.centralLayerRadii = std::vector<double>(1, kBeamPipeRadius);
  bplConfig.centralLayerHalflengthZ =
      std::vector<double>(1, kBeamPipeHalfLengthZ);
  bplConfig.centralLayerThickness = std::vector<double>(1, kBeamPipeThickness);
  bplConfig.centralLayerMaterial = {
      std::make_shared<const HomogeneousSurfaceMaterial>(beamPipeMaterial)};
  auto beamPipeBuilder = std::make_shared<const PassiveLayerBuilder>(
      bplConfig, getDefaultLogger("BeamPipeLayerBuilder", layerLLevel));
  // create the volume for the beam pipe
  CylinderVolumeBuilder::Config bpvConfig;
  bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
  bpvConfig.volumeName = "BeamPipe";
  bpvConfig.layerBuilder = beamPipeBuilder;
  bpvConfig.layerEnvelopeR = {1_mm, 1_mm};
  bpvConfig.buildToRadiusZero = true;
  auto beamPipeVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
      bpvConfig, getDefaultLogger("BeamPipeVolumeBuilder", volumeLLevel));

  // create the bounds and the volume
  auto beamPipeBounds =
      std::make_shared<const CylinderVolumeBounds>(0., 25., 1100.);
  auto beamPipeVolume = beamPipeVolumeBuilder->trackingVolume(
      geoContext, nullptr, beamPipeBounds);

  //-------------------------------------------------------------------------------------
  // some prep work for the material
  // Layer material properties - thickness, X0, L0, A, Z, Rho
  MaterialSlab lProperties(makeSilicon(), 1.5_mm);

  auto layerMaterialPtr =
      std::make_shared<HomogeneousSurfaceMaterial>(lProperties);

  std::vector<LayerPtr> pLayers;

  for (std::size_t ilp = 0; ilp < kLayerRadii.size(); ++ilp) {
    std::vector<Surface*> layerSurfacesMutable =
        surfacesCylinder(detectorStore, kModuleHalfX[ilp], kModuleHalfY[ilp],
                         kModuleThickness[ilp], kModuleTiltPhi[ilp],
                         kLayerRadii[ilp], 2_mm, 5_mm, kLayerBinning[ilp]);

    std::vector<const Surface*> layerSurfaces;
    std::ranges::transform(layerSurfacesMutable,
                           std::back_inserter(layerSurfaces),
                           [](Surface* s) { return s; });

    // Make a shared version out of it
    std::vector<std::shared_ptr<const Surface>> layerSurfacePtrs;
    layerSurfacePtrs.reserve(layerSurfaces.size());
    for (auto& sf : layerSurfaces) {
      layerSurfacePtrs.push_back(sf->getSharedPtr());
    }

    // create the layer and store it
    ProtoLayer protoLayer(geoContext, layerSurfaces);
    protoLayer.envelope[AxisDirection::AxisR] = {0.5, 0.5};
    auto pLayer = layerCreator->cylinderLayer(
        geoContext, std::move(layerSurfacePtrs), kLayerBinning[ilp].first,
        kLayerBinning[ilp].second, protoLayer);
    auto approachSurfaces = pLayer->approachDescriptor()->containedSurfaces();
    auto mutableOuterSurface = const_cast<Surface*>(approachSurfaces.at(1));
    mutableOuterSurface->assignSurfaceMaterial(layerMaterialPtr);
    /// now push back the layer
    pLayers.push_back(pLayer);

  }  // loop over layers

  // layer array
  auto pLayerArray = layerArrayCreator->layerArray(
      geoContext, pLayers, 25., 300., arbitrary, AxisDirection::AxisR);
  auto pVolumeBounds = std::make_shared<CylinderVolumeBounds>(25., 300., 1100.);
  // create the Tracking volume
  auto pVolume = std::make_shared<TrackingVolume>(
      Transform3::Identity(), pVolumeBounds, nullptr, std::move(pLayerArray),
      nullptr, MutableTrackingVolumeVector{}, "Pixel::Barrel");

  // The combined volume
  auto detectorVolume = cylinderVolumeHelper->createContainerTrackingVolume(
      geoContext, {beamPipeVolume, pVolume});

  // create and return the geometry
  return std::make_shared<TrackingGeometry>(detectorVolume);
}

std::shared_ptr<TrackingGeometry> CylindricalTrackingGeometry::buildGen3(
    const Logger& logger) {
  using namespace Experimental;
  using namespace UnitLiterals;
  using enum CylinderVolumeBounds::Face;
  using enum AxisDirection;
  using LayerType = LayerBlueprintNode::LayerType;

  const MaterialSlab lProperties(makeSilicon(), 1.5_mm);

  // Create a binned material in 2 bins - irregularly in z, 2 bins in phi
  std::vector<float> binEdges = {// empirical bin edges. these are not checked!
                                 -476.5, 0, 476.5};

  BinUtility binUtility(2u, -std::numbers::pi, std::numbers::pi,
                        BinningOption::closed, AxisPhi);

  binUtility += Acts::BinUtility(binEdges, BinningOption::open, AxisZ);

  std::vector<MaterialSlab> materialSlabs0 = {lProperties, lProperties};
  std::vector<MaterialSlab> materialSlabs1 = {lProperties, lProperties};

  auto binnedMaterial = std::make_shared<BinnedSurfaceMaterial>(
      binUtility, std::vector{materialSlabs0, materialSlabs1}, 0.,
      MappingType::Default);

  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .z = {20_mm, 20_mm},
      .r = {0_mm, 20_mm},
  }};
  Blueprint root{cfg};

  root.addCylinderContainer("Detector", AxisR, [&](auto& detector) {
    auto beampipeBounds =
        std::make_unique<CylinderVolumeBounds>(0_mm, kBeamPipeRadius, 100_mm);
    auto beampipe = std::make_unique<TrackingVolume>(
        Transform3::Identity(), std::move(beampipeBounds), "Beampipe");

    detector.addMaterial("BeampipeMaterial", [&](auto& bpMat) {
      MaterialSlab beamPipeMaterial(makeBeryllium(), kBeamPipeThickness);
      bpMat.configureFace(OuterCylinder, binnedMaterial);
      bpMat.addStaticVolume(std::move(beampipe));
    });

    auto layerMaterialPtr =
        std::make_shared<HomogeneousSurfaceMaterial>(lProperties);

    for (std::size_t ilp = 0; ilp < kLayerRadii.size(); ++ilp) {
      std::vector<Surface*> layerSurfaces =
          surfacesCylinder(detectorStore, kModuleHalfX[ilp], kModuleHalfY[ilp],
                           kModuleThickness[ilp], kModuleTiltPhi[ilp],
                           kLayerRadii[ilp], 2_mm, 5_mm, kLayerBinning[ilp]);

      // Make a shared version out of it
      std::vector<std::shared_ptr<Surface>> layerSurfacePtrs;
      layerSurfacePtrs.reserve(layerSurfaces.size());
      for (auto& sf : layerSurfaces) {
        layerSurfacePtrs.push_back(sf->getSharedPtr());
      }

      // create the layer and store it
      MutableProtoLayer protoLayer(geoContext, layerSurfaces);
      protoLayer.envelope[AxisR] = {0.5, 0.5};

      std::string layerName = "L" + std::to_string(ilp);
      detector.addMaterial(layerName + "_Mat", [&](auto& mat) {
        mat.configureFace(OuterCylinder, layerMaterialPtr);
        mat.addLayer(layerName, [&](auto& layer) {
          layer.setProtoLayer(protoLayer);
          layer.setLayerType(LayerType::Cylinder);
          layer.setEnvelope(
              ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {5_mm, 5_mm}}});

          using SrfArrayNavPol = SurfaceArrayNavigationPolicy;

          layer.setNavigationPolicyFactory(
              NavigationPolicyFactory{}
                  .add<TryAllNavigationPolicy>(
                      TryAllNavigationPolicy::Config{.sensitives = false})
                  .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                      .layerType = SrfArrayNavPol::LayerType::Cylinder,
                      .bins = kLayerBinning[ilp]})
                  .asUniquePtr());
        });
      });

    }  // loop over layers
  });

  BlueprintOptions opts;
  return root.construct(opts, geoContext, logger);
}

std::shared_ptr<TrackingGeometry> CylindricalTrackingGeometry::operator()(
    const Logger& logger) {
  if (gen3) {
    return buildGen3(logger);
  } else {
    return buildGen1(logger);
  }
}
}  // namespace ActsTests
