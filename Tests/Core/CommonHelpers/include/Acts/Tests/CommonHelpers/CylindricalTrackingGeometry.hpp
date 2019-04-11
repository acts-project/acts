// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <vector>

#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace Test {

struct CylindricalTrackingGeometry {
  std::reference_wrapper<const GeometryContext> geoContext;

  /// Only allowed constructor with reference wrapper
  CylindricalTrackingGeometry(
      std::reference_wrapper<const GeometryContext> gctx)
      : geoContext(gctx) {}

  /// The detector store for memory management
  std::vector<std::unique_ptr<const DetectorElementStub>> detectorStore = {};

  /// helper method for cylinder layer
  /// create the positions for module surfaces on a cylinder
  std::vector<Vector3D> modulePositionsCylinder(
      double radius, double zStagger, double moduleHalfLength, double lOverlap,
      const std::pair<int, int>& binningSchema) {
    int nPhiBins = binningSchema.first;
    int nZbins = binningSchema.second;
    // prepare the return value
    std::vector<Vector3D> mPositions;
    mPositions.reserve(nPhiBins * nZbins);
    // prep work
    double phiStep = 2 * M_PI / (nPhiBins);
    double minPhi = -M_PI + 0.5 * phiStep;
    double zStart = -0.5 * (nZbins - 1) * (2 * moduleHalfLength - lOverlap);
    double zStep = 2 * std::abs(zStart) / (nZbins - 1);
    // loop over the bins
    for (size_t zBin = 0; zBin < size_t(nZbins); ++zBin) {
      // prepare z and r
      double moduleZ = zStart + zBin * zStep;
      double moduleR =
          (zBin % 2) != 0u ? radius - 0.5 * zStagger : radius + 0.5 * zStagger;
      for (size_t phiBin = 0; phiBin < size_t(nPhiBins); ++phiBin) {
        // calculate the current phi value
        double modulePhi = minPhi + phiBin * phiStep;
        mPositions.push_back(Vector3D(moduleR * cos(modulePhi),
                                      moduleR * sin(modulePhi), moduleZ));
      }
    }
    return mPositions;
  }

  // @brief Call operator for the creation method of the tracking geometry
  std::shared_ptr<const TrackingGeometry> operator()() {
    using namespace Acts::UnitLiterals;

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
    auto tVolumeArrayCreator =
        std::make_shared<const TrackingVolumeArrayCreator>(
            tvacConfig,
            getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
    // configure the cylinder volume helper
    CylinderVolumeHelper::Config cvhConfig;
    cvhConfig.layerArrayCreator = layerArrayCreator;
    cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
    auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
        cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));

    // ----------------- build a beam pipe -----------------------------------
    MaterialProperties beamPipeMaterial{352.8, 407., 9.012, 4., 1.848e-3, 0.8};
    PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = "BeamPipe";
    bplConfig.centralLayerRadii = std::vector<double>(1, 19.);
    bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 1000.);
    bplConfig.centralLayerThickness = std::vector<double>(1, 0.8);
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
    bpvConfig.volumeSignature = 0;
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
    MaterialProperties lProperties(95.7, 465.2, 28.03, 14., 2.32e-3, 1.5_mm);

    std::shared_ptr<const ISurfaceMaterial> layerMaterialPtr =
        std::shared_ptr<const ISurfaceMaterial>(
            new Acts::HomogeneousSurfaceMaterial(lProperties));

    // Module material - X0, L0, A, Z, Rho
    Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

    std::vector<double> pLayerRadii = {32., 72., 116., 172.};
    std::vector<std::pair<int, int>> pLayerBinning = {
        {16, 14}, {32, 14}, {52, 14}, {78, 14}};
    std::vector<double> pModuleTiltPhi = {0.145, 0.145, 0.145, 0.145};
    std::vector<double> pModuleHalfX = {8.4, 8.4, 8.4, 8.4};
    std::vector<double> pModuleHalfY = {36., 36., 36., 36.};
    std::vector<double> pModuleThickness = {0.15, 0.15, 0.15, 0.15};

    std::vector<LayerPtr> pLayers;

    for (size_t ilp = 0; ilp < pLayerRadii.size(); ++ilp) {
      std::vector<std::shared_ptr<const Surface>> layerModules;

      // Module material from input
      MaterialProperties moduleMaterialProperties(pcMaterial,
                                                  pModuleThickness[ilp]);
      // Create a new surface material
      std::shared_ptr<const ISurfaceMaterial> moduleMaterialPtr =
          std::shared_ptr<const ISurfaceMaterial>(
              new Acts::HomogeneousSurfaceMaterial(moduleMaterialProperties));

      // The rectangle bounds for all modules
      auto mBounds = std::make_shared<RectangleBounds>(pModuleHalfX[ilp],
                                                       pModuleHalfY[ilp]);
      // Create the module centers
      auto moduleCenters = modulePositionsCylinder(
          pLayerRadii[ilp], 2_mm, pModuleHalfY[ilp], 5_mm, pLayerBinning[ilp]);

      for (auto& mCenter : moduleCenters) {
        // The association transform
        double modulePhi = VectorHelpers::phi(mCenter);
        // Local z axis is the normal vector
        Vector3D moduleLocalZ(cos(modulePhi + pModuleTiltPhi[ilp]),
                              sin(modulePhi + pModuleTiltPhi[ilp]), 0.);
        // Local y axis is the global z axis
        Vector3D moduleLocalY(0., 0., 1);
        // Local x axis the normal to local y,z
        Vector3D moduleLocalX(-sin(modulePhi + pModuleTiltPhi[ilp]),
                              cos(modulePhi + pModuleTiltPhi[ilp]), 0.);
        // Create the RotationMatrix
        RotationMatrix3D moduleRotation;
        moduleRotation.col(0) = moduleLocalX;
        moduleRotation.col(1) = moduleLocalY;
        moduleRotation.col(2) = moduleLocalZ;
        // Get the moduleTransform
        std::shared_ptr<Transform3D> mModuleTransform =
            std::make_shared<Transform3D>(Translation3D(mCenter) *
                                          moduleRotation);
        // Create the detector element
        auto detElement = std::make_unique<const DetectorElementStub>(
            mModuleTransform, mBounds, pModuleThickness[ilp],
            moduleMaterialPtr);

        layerModules.push_back(detElement->surface().getSharedPtr());
        detectorStore.push_back(std::move(detElement));
      }
      // create the layer and store it
      ProtoLayer protoLayer(geoContext, layerModules);
      protoLayer.envR = {0.5, 0.5};
      auto pLayer = layerCreator->cylinderLayer(
          geoContext, std::move(layerModules), pLayerBinning[ilp].first,
          pLayerBinning[ilp].second, protoLayer);
      auto approachSurfaces = pLayer->approachDescriptor()->containedSurfaces();
      auto mutableOuterSurface =
          const_cast<Acts::Surface*>(approachSurfaces.at(1));
      mutableOuterSurface->assignSurfaceMaterial(layerMaterialPtr);
      /// now push back the layer
      pLayers.push_back(pLayer);

    }  // loop over layers

    // layer array
    auto pLayerArray = layerArrayCreator->layerArray(geoContext, pLayers, 25.,
                                                     300., arbitrary, binR);
    auto pVolumeBounds =
        std::make_shared<const CylinderVolumeBounds>(25., 300., 1100.);
    // create the Tracking volume
    auto pVolume = TrackingVolume::create(nullptr, pVolumeBounds, nullptr,
                                          std::move(pLayerArray), nullptr,
                                          {}, "Pixel::Barrel");

    // The combined volume
    auto detectorVolume = cylinderVolumeHelper->createContainerTrackingVolume(
        geoContext, {beamPipeVolume, pVolume});

    // create and return the geometry
    return std::make_shared<const TrackingGeometry>(detectorVolume);
  }
};

}  // namespace Test
}  // namespace Acts
