// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTRAPOLATION_TEST_GEOMETRY_HPP
#define ACTS_EXTRAPOLATION_TEST_GEOMETRY_HPP

#include <vector>
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/ProtoLayer.hpp"
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/LayerArrayCreator.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/PassiveLayerBuilder.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Tools/TrackingVolumeArrayCreator.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

namespace Acts {

namespace Test {

  // a non-free module surface to test (fully) connected geometry
  // i.e. it avoids surface copying
  class ModuleSurface : public PlaneSurface
  {
  public:
    // constructor to forward to PlaneSurface
    ModuleSurface(std::shared_ptr<const Transform3D>  htrans,
                  std::shared_ptr<const PlanarBounds> pbounds)
      : PlaneSurface(htrans, pbounds)
    {
    }

    // causes pointer assignment
    virtual bool
    isFree() const override final
    {
      return false;
    }
  };

  /// helper method for cylinder layer
  /// create the positions for module surfaces on a cylinder
  static std::vector<Vector3D>
  modulePositionsCylinder(double radius,
                          double zStagger,
                          double moduleHalfLength,
                          double lOverlap,
                          const std::pair<int, int>& binningSchema)
  {
    int nPhiBins = binningSchema.first;
    int nZbins   = binningSchema.second;
    // prepare the return value
    std::vector<Vector3D> mPositions;
    mPositions.reserve(nPhiBins * nZbins);
    // prep work
    double phiStep = 2 * M_PI / (nPhiBins);
    double minPhi  = -M_PI + 0.5 * phiStep;
    double zStart  = -0.5 * (nZbins - 1) * (2 * moduleHalfLength - lOverlap);
    double zStep   = 2 * std::abs(zStart) / (nZbins - 1);
    // loop over the bins
    for (size_t zBin = 0; zBin < size_t(nZbins); ++zBin) {
      // prepare z and r
      double moduleZ = zStart + zBin * zStep;
      double moduleR
          = (zBin % 2) ? radius - 0.5 * zStagger : radius + 0.5 * zStagger;
      for (size_t phiBin = 0; phiBin < size_t(nPhiBins); ++phiBin) {
        // calculate the current phi value
        double modulePhi = minPhi + phiBin * phiStep;
        mPositions.push_back(Vector3D(
            moduleR * cos(modulePhi), moduleR * sin(modulePhi), moduleZ));
      }
    }
    return mPositions;
  }

  // the creation method for the simple tracking geometry
  //
  // @param surface_cache is the surface cache to be filled
  //        for memory management
  template <typename Plane_type>
  std::shared_ptr<const TrackingGeometry>
  testGeometry(std::vector<std::unique_ptr<const Surface>>& surface_cache)
  {

    Logging::Level surfaceLLevel = Logging::INFO;
    Logging::Level layerLLevel   = Logging::INFO;
    Logging::Level volumeLLevel  = Logging::INFO;

    // configure surface array creator
    auto surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>(
        getDefaultLogger("SurfaceArrayCreator", surfaceLLevel));
    // configure the layer creator that uses the surface array creator
    LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator            = std::make_shared<const LayerCreator>(
        lcConfig, getDefaultLogger("LayerCreator", layerLLevel));
    // configure the layer array creator
    auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
        getDefaultLogger("LayerArrayCreator", layerLLevel));

    // tracking volume array creator
    auto tVolumeArrayCreator
        = std::make_shared<const TrackingVolumeArrayCreator>(
            getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
    // configure the cylinder volume helper
    CylinderVolumeHelper::Config cvhConfig;
    cvhConfig.layerArrayCreator          = layerArrayCreator;
    cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
    auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
        cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));

    // ----------------- build a beam pipe -----------------------------------
    PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification     = "BeamPipe";
    bplConfig.centralLayerRadii       = std::vector<double>(1, 19.);
    bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 1000.);
    bplConfig.centralLayerThickness   = std::vector<double>(1, 0.8);
    bplConfig.centralLayerMaterial
        = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
    auto beamPipeBuilder = std::make_shared<const PassiveLayerBuilder>(
        bplConfig, getDefaultLogger("BeamPipeLayerBuilder", layerLLevel));
    // create the volume for the beam pipe
    CylinderVolumeBuilder::Config bpvConfig;
    bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
    bpvConfig.volumeName           = "BeamPipe";
    bpvConfig.layerBuilder         = beamPipeBuilder;
    bpvConfig.layerEnvelopeR       = {1. * units::_mm, 1. * units::_mm};
    bpvConfig.buildToRadiusZero    = true;
    bpvConfig.volumeSignature      = 0;
    auto beamPipeVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
        bpvConfig, getDefaultLogger("BeamPipeVolumeBuilder", volumeLLevel));

    // create the bounds and the volume
    auto beamPipeBounds
        = std::make_shared<const CylinderVolumeBounds>(0., 25., 1100.);
    auto beamPipeVolume
        = beamPipeVolumeBuilder->trackingVolume(nullptr, beamPipeBounds);

    //-------------------------------------------------------------------------------------
    // some prep work for the material
    // Layer material properties - thickness, X0, L0, A, Z, Rho
    MaterialProperties lProperties(
        95.7, 465.2, 28.03, 14., 2.32e-3, 1.5 * units::_mm);

    std::shared_ptr<const SurfaceMaterial> layerMaterialPtr
        = std::shared_ptr<const SurfaceMaterial>(
            new Acts::HomogeneousSurfaceMaterial(lProperties));

    // Module material - X0, L0, A, Z, Rho
    Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

    std::vector<double> pLayerRadii = {32., 72., 116., 172.};
    std::vector<std::pair<int, int>> pLayerBinning
        = {{16, 14}, {32, 14}, {52, 14}, {78, 14}};
    std::vector<double> pModuleTiltPhi   = {0.145, 0.145, 0.145, 0.145};
    std::vector<double> pModuleHalfX     = {8.4, 8.4, 8.4, 8.4};
    std::vector<double> pModuleHalfY     = {36., 36., 36., 36.};
    std::vector<double> pModuleThickness = {0.15, 0.15, 0.15, 0.15};

    std::vector<LayerPtr> pLayers;

    for (size_t ilp = 0; ilp < pLayerRadii.size(); ++ilp) {

      std::vector<const Surface*> layerModules;

      // Module material from input
      MaterialProperties moduleMaterialProperties(pcMaterial,
                                                  pModuleThickness[ilp]);
      // create a new surface material
      std::shared_ptr<const SurfaceMaterial> moduleMaterialPtr
          = std::shared_ptr<const SurfaceMaterial>(
              new Acts::HomogeneousSurfaceMaterial(moduleMaterialProperties));

      // The rectangle bounds for all modules
      auto mBounds = std::make_shared<RectangleBounds>(pModuleHalfX[ilp],
                                                       pModuleHalfY[ilp]);
      // create the module centers
      auto moduleCenters = modulePositionsCylinder(pLayerRadii[ilp],
                                                   2. * units::_mm,
                                                   pModuleHalfY[ilp],
                                                   5. * units::_mm,
                                                   pLayerBinning[ilp]);

      for (auto& mCenter : moduleCenters) {
        // create the association transform
        double modulePhi = mCenter.phi();
        // the local z axis is the normal vector
        Vector3D moduleLocalZ(cos(modulePhi + pModuleTiltPhi[ilp]),
                              sin(modulePhi + pModuleTiltPhi[ilp]),
                              0.);
        // the local y axis is the global z axis
        Vector3D moduleLocalY(0., 0., 1);
        // the local x axis the normal to local y,z
        Vector3D moduleLocalX(-sin(modulePhi + pModuleTiltPhi[ilp]),
                              cos(modulePhi + pModuleTiltPhi[ilp]),
                              0.);
        // create the RotationMatrix
        RotationMatrix3D moduleRotation;
        moduleRotation.col(0) = moduleLocalX;
        moduleRotation.col(1) = moduleLocalY;
        moduleRotation.col(2) = moduleLocalZ;
        // get the moduleTransform
        std::shared_ptr<Transform3D> mModuleTransform(new Transform3D(
            getTransformFromRotTransl(moduleRotation, mCenter)));

        Plane_type* mSurface = new Plane_type(mModuleTransform, mBounds);
        // let's assign the material
        mSurface->setAssociatedMaterial(moduleMaterialPtr);

        layerModules.push_back(mSurface);
        surface_cache.push_back(std::unique_ptr<const Surface>(mSurface));
      }
      // create the layer and store it
      ProtoLayer protoLayer(layerModules);
      protoLayer.envR = {0.5, 0.5};
      auto pLayer     = layerCreator->cylinderLayer(layerModules,
                                                pLayerBinning[ilp].first,
                                                pLayerBinning[ilp].second,
                                                protoLayer);
      auto approachSurfaces = pLayer->approachDescriptor()->containedSurfaces();
      auto mutableOuterSurface
          = const_cast<Acts::Surface*>(approachSurfaces.at(1));
      mutableOuterSurface->setAssociatedMaterial(layerMaterialPtr);
      /// now push back the layer
      pLayers.push_back(pLayer);

    }  // loop over layers

    // layer array
    auto pLayerArray
        = layerArrayCreator->layerArray(pLayers, 25., 300., arbitrary, binR);
    auto pVolumeBounds
        = std::make_shared<const CylinderVolumeBounds>(25., 300., 1100.);
    // create the Tracking volume
    auto pVolume = TrackingVolume::create(nullptr,
                                          pVolumeBounds,
                                          nullptr,
                                          std::move(pLayerArray),
                                          {},
                                          {},
                                          {},
                                          "Pixel::Barrel");

    // the combined volume
    auto detectorVolume = cylinderVolumeHelper->createContainerTrackingVolume(
        {beamPipeVolume, pVolume});

    // create and return the geometry
    return std::make_shared<const TrackingGeometry>(detectorVolume);
  }

}  // end of namespace Test

}  // end of namespace Acts

#endif