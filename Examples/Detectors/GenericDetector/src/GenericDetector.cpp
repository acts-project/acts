// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"
#include "ActsExamples/GenericDetector/LayerBuilder.hpp"

#include <fstream>

#include "./GenericDetectorBuilder.hpp"

namespace ActsExamples {

namespace Generic {
namespace {

class Gen1GenericDetectorBuilder : public GenericDetectorBuilder {
 public:
  using GenericDetectorBuilder::GenericDetectorBuilder;

  std::unique_ptr<const Acts::TrackingGeometry> buildTrackingGeometry(
      const Acts::GeometryContext& gctx,
      std::shared_ptr<const Acts::IMaterialDecorator> matDecorator,
      Acts::Logging::Level surfaceLLevel, Acts::Logging::Level volumeLLevel);
};

std::unique_ptr<const Acts::TrackingGeometry>
Gen1GenericDetectorBuilder::buildTrackingGeometry(
    const Acts::GeometryContext& gctx,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator,
    Acts::Logging::Level surfaceLLevel, Acts::Logging::Level volumeLLevel) {
  ACTS_INFO("Building tracking geometry for Generic Detector in Gen1 mode");
  using namespace Acts::UnitLiterals;

  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig, logger().clone("SurfaceArrayCreator", surfaceLLevel));
  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, logger().clone("LayerCreator", m_cfg.layerLogLevel));
  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, logger().clone("LayerArrayCreator", m_cfg.layerLogLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          logger().clone("TrackingVolumeArrayCreator", volumeLLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig, logger().clone("CylinderVolumeHelper", volumeLLevel));
  //-------------------------------------------------------------------------------------
  // vector of the volume builders
  std::vector<std::shared_ptr<const Acts::ITrackingVolumeBuilder>>
      volumeBuilders;

  ACTS_DEBUG("Building BeamPipe");
  //-------------------------------------------------------------------------------------
  // Beam Pipe
  //-------------------------------------------------------------------------------------

  // configure the beam pipe layer builder
  Acts::PassiveLayerBuilder::Config bplConfig;
  bplConfig.layerIdentification = "BeamPipe";
  bplConfig.centralLayerRadii = std::vector<double>(1, kBeamPipeRadius);
  bplConfig.centralLayerHalflengthZ =
      std::vector<double>(1, kBeamPipeHalfLengthZ);
  bplConfig.centralLayerThickness = std::vector<double>(1, kBeamPipeThickness);
  bplConfig.centralLayerMaterial = {m_beamPipeMaterial};
  auto beamPipeBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
      bplConfig, logger().clone("BeamPipeLayerBuilder", m_cfg.layerLogLevel));
  // create the volume for the beam pipe
  Acts::CylinderVolumeBuilder::Config bpvConfig;
  bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
  bpvConfig.volumeName = "BeamPipe";
  bpvConfig.layerBuilder = beamPipeBuilder;
  bpvConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                              1. * Acts::UnitConstants::mm};
  bpvConfig.buildToRadiusZero = true;
  auto beamPipeVolumeBuilder =
      std::make_shared<const Acts::CylinderVolumeBuilder>(
          bpvConfig, logger().clone("BeamPipeVolumeBuilder", volumeLLevel));
  // add to the list of builders
  volumeBuilders.push_back(beamPipeVolumeBuilder);

  ACTS_DEBUG("Building Pixel");
  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // Pixel detector
  //-------------------------------------------------------------------------------------
  // some prep work

  ProtoLayerCreator pplCreator = createPixelProtoLayerCreator();

  // configure pixel layer builder
  LayerBuilder::Config plbConfig;
  plbConfig.layerCreator = layerCreator;
  plbConfig.layerIdentification = "Pixel";
  // material concentration alsways outside the modules
  plbConfig.centralProtoLayers = pplCreator.centralProtoLayers(gctx);
  plbConfig.centralLayerMaterialConcentration = {1, 1, 1, 1};
  plbConfig.centralLayerMaterial = {
      m_pixelCentralMaterial, m_pixelCentralMaterial, m_pixelCentralMaterial,
      m_pixelCentralMaterial};
  if (m_cfg.buildLevel > 0) {
    // material concentration is always behind the layer in the pixels
    plbConfig.posnegLayerMaterialConcentration = std::vector<int>(7, 0);
    // layer structure surface has pixel material properties
    plbConfig.posnegLayerMaterial = {
        m_pixelEndcapMaterial, m_pixelEndcapMaterial, m_pixelEndcapMaterial,
        m_pixelEndcapMaterial, m_pixelEndcapMaterial, m_pixelEndcapMaterial,
        m_pixelEndcapMaterial};
    // negative proto layers
    plbConfig.negativeProtoLayers = pplCreator.negativeProtoLayers(gctx);
    plbConfig.positiveProtoLayers = pplCreator.positiveProtoLayers(gctx);
  }
  // define the builder
  auto pixelLayerBuilder = std::make_shared<const LayerBuilder>(
      plbConfig, logger().clone("PixelLayerBuilder", m_cfg.layerLogLevel));
  //-------------------------------------------------------------------------------------
  // build the pixel volume
  Acts::CylinderVolumeBuilder::Config pvbConfig;
  pvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  pvbConfig.volumeName = "Pixel";
  pvbConfig.buildToRadiusZero = false;
  pvbConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                              5. * Acts::UnitConstants::mm};
  pvbConfig.layerBuilder = pixelLayerBuilder;
  auto pixelVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
      pvbConfig, logger().clone("PixelVolumeBuilder", volumeLLevel));
  // add to the list of builders
  volumeBuilders.push_back(pixelVolumeBuilder);

  if (m_cfg.buildLevel > 1) {
    ACTS_DEBUG("Building PST");
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------
    // Pixel Support Tybe (PST)
    //-------------------------------------------------------------------------------------

    // Configuration
    Acts::PassiveLayerBuilder::Config pstConfig;
    pstConfig.layerIdentification = "PST";
    pstConfig.centralLayerRadii = std::vector<double>(1, kPstRadius);
    pstConfig.centralLayerHalflengthZ = std::vector<double>(1, kPstHalfLengthZ);
    pstConfig.centralLayerThickness = std::vector<double>(1, kPstThickness);
    pstConfig.centralLayerMaterial = {m_pstMaterial};
    auto pstBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        pstConfig, logger().clone("PSTLayerBuilder", m_cfg.layerLogLevel));
    // create the volume for the beam pipe
    Acts::CylinderVolumeBuilder::Config pstvolConfig;
    pstvolConfig.trackingVolumeHelper = cylinderVolumeHelper;
    pstvolConfig.volumeName = "PST";
    pstvolConfig.buildToRadiusZero = false;
    pstvolConfig.layerBuilder = pstBuilder;
    auto pstVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        pstvolConfig, logger().clone("PSTVolumeBuilder", volumeLLevel));
    // add to the detector builds
    volumeBuilders.push_back(pstVolumeBuilder);

    ACTS_DEBUG("Building SStrip");
    //-------------------------------------------------------------------------------------
    // SHORT strip detector
    //-------------------------------------------------------------------------------------
    // first add a Pixel Support Tube
    // STRIPS
    //
    // fill necessary vectors for configuration
    //-------------------------------------------------------------------------------------
    // some prep work

    ProtoLayerCreator ssplCreator = createShortStripProtoLayerCreator();

    std::size_t nposnegs = ssplCreator.config().posnegLayerPositionsZ.size();

    // configure short strip layer builder
    LayerBuilder::Config sslbConfig;
    sslbConfig.layerCreator = layerCreator;
    sslbConfig.layerIdentification = "SStrip";

    sslbConfig.centralProtoLayers = ssplCreator.centralProtoLayers(gctx);
    sslbConfig.centralLayerMaterialConcentration = {-1, -1, -1, -1};
    sslbConfig.centralLayerMaterial = {
        m_shortStripCentralMaterial, m_shortStripCentralMaterial,
        m_shortStripCentralMaterial, m_shortStripCentralMaterial};

    if (m_cfg.buildLevel > 2) {
      sslbConfig.negativeProtoLayers = ssplCreator.negativeProtoLayers(gctx);
      sslbConfig.positiveProtoLayers = ssplCreator.positiveProtoLayers(gctx);

      sslbConfig.posnegLayerMaterialConcentration =
          std::vector<int>(nposnegs, 0);
      sslbConfig.posnegLayerMaterial =
          std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>(
              nposnegs, m_shortStripEndcapMaterial);
    }

    // define the builder
    auto sstripLayerBuilder = std::make_shared<const LayerBuilder>(
        sslbConfig, logger().clone("SStripLayerBuilder", m_cfg.layerLogLevel));
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config ssvbConfig;
    ssvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    ssvbConfig.volumeName = "SStrip";
    ssvbConfig.buildToRadiusZero = false;
    ssvbConfig.layerBuilder = sstripLayerBuilder;
    auto sstripVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            ssvbConfig, logger().clone("SStripVolumeBuilder", volumeLLevel));

    //-------------------------------------------------------------------------------------
    // add to the list of builders
    volumeBuilders.push_back(sstripVolumeBuilder);
    //-------------------------------------------------------------------------------------
    ACTS_DEBUG("Building LStrip");
    //-------------------------------------------------------------------------------------
    // LONG strip detector
    //-------------------------------------------------------------------------------------
    // fill necessary vectors for configuration
    //-------------------------------------------------------------------------------------

    // some prep work

    ProtoLayerCreator lsplCreator = createLongStripProtoLayerCreator();

    // configure short strip layer builder
    Generic::LayerBuilder::Config lslbConfig;
    lslbConfig.layerCreator = layerCreator;
    lslbConfig.layerIdentification = "LStrip";
    lslbConfig.centralLayerMaterialConcentration = {-1, -1};
    lslbConfig.centralLayerMaterial = {m_longStripCentralMaterial,
                                       m_longStripCentralMaterial};
    lslbConfig.centralProtoLayers = lsplCreator.centralProtoLayers(gctx);

    if (m_cfg.buildLevel > 2) {
      lslbConfig.posnegLayerMaterialConcentration =
          std::vector<int>(nposnegs, 0);
      lslbConfig.posnegLayerMaterial =
          std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>(
              nposnegs, m_longStripEndcapMaterial);
      lslbConfig.negativeProtoLayers = lsplCreator.negativeProtoLayers(gctx);
      lslbConfig.positiveProtoLayers = lsplCreator.positiveProtoLayers(gctx);
    }

    // define the builder
    auto lstripLayerBuilder = std::make_shared<const LayerBuilder>(
        lslbConfig, logger().clone("LStripLayerBuilder", m_cfg.layerLogLevel));
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config lsvbConfig;
    lsvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    lsvbConfig.volumeName = "LStrip";
    lsvbConfig.buildToRadiusZero = false;
    lsvbConfig.layerBuilder = lstripLayerBuilder;
    auto lstripVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            lsvbConfig, logger().clone("LStripVolumeBuilder", volumeLLevel));
    // add to the list of builders
    volumeBuilders.push_back(lstripVolumeBuilder);
  }

  ACTS_DEBUG("Building TrackingGeometry");
  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  Acts::TrackingGeometryBuilder::Config tgConfig;
  // Add the build call functions
  for (auto& vb : volumeBuilders) {
    tgConfig.trackingVolumeBuilders.emplace_back(
        [=](const auto& context, const auto& inner, const auto&) {
          return vb->trackingVolume(context, inner);
        });
  }
  tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  tgConfig.materialDecorator = std::move(matDecorator);

  auto cylinderGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(
          tgConfig,
          logger().clone("CylinderGeometryBuilder", Acts::Logging::INFO));
  // get the geometry
  auto trackingGeometry = cylinderGeometryBuilder->trackingGeometry(gctx);
  // return the tracking geometry
  return trackingGeometry;
}

class Gen3GenericDetectorBuilder : public GenericDetectorBuilder {
 public:
  using GenericDetectorBuilder::GenericDetectorBuilder;

  std::unique_ptr<const Acts::TrackingGeometry> buildTrackingGeometry(
      const Acts::GeometryContext& gctx,
      std::optional<std::filesystem::path> graphvizFile);

  void buildPixel(Acts::Experimental::BlueprintNode& parent,
                  const Acts::GeometryContext& gctx);

  void buildShortStrip(Acts::Experimental::BlueprintNode& parent,
                       const Acts::GeometryContext& gctx);

  void buildLongStrip(Acts::Experimental::BlueprintNode& parent,
                      const Acts::GeometryContext& gctx);

  static std::shared_ptr<const Acts::HomogeneousSurfaceMaterial> asHomogeneous(
      std::string_view debugLabel,
      const std::shared_ptr<const Acts::ISurfaceMaterial>& material) {
    auto hm = std::dynamic_pointer_cast<const Acts::HomogeneousSurfaceMaterial>(
        material);
    if (hm == nullptr) {
      throw std::runtime_error(std::string{debugLabel} +
                               " material is not homogeneous");
    }
    return hm;
  }

  using LayerType = Acts::SurfaceArrayNavigationPolicy::LayerType;
  std::unique_ptr<Acts::NavigationPolicyFactory> createNavigationPolicy(
      LayerType layerType,
      std::pair<std::size_t, std::size_t> bins = {0, 0}) const {
    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    return Acts::NavigationPolicyFactory::make()
        .add<Acts::TryAllNavigationPolicy>(
            Acts::TryAllNavigationPolicy::Config{.sensitives = false})
        .add<SrfArrayNavPol>(
            SrfArrayNavPol::Config{.layerType = layerType, .bins = bins})
        .asUniquePtr();
  }
};

std::unique_ptr<const Acts::TrackingGeometry>
Gen3GenericDetectorBuilder::buildTrackingGeometry(
    const Acts::GeometryContext& gctx,
    std::optional<std::filesystem::path> graphvizFile) {
  using enum Acts::AxisDirection;
  using namespace Acts::UnitLiterals;
  using namespace Acts::Experimental;
  using enum Acts::AxisBoundaryType;
  using enum Acts::CylinderVolumeBounds::Face;
  ACTS_INFO("GenericDetector construction in  Gen3 mode");

  Blueprint::Config cfg;
  cfg.envelope = Acts::ExtentEnvelope{{
      .z = {20_mm, 20_mm},
      .r = {0_mm, 20_mm},
  }};
  Blueprint root{cfg};

  root.addCylinderContainer("Detector", AxisR, [&](auto& detector) {
    auto beampipeBounds = std::make_unique<Acts::CylinderVolumeBounds>(
        0_mm, kBeamPipeRadius, kBeamPipeHalfLengthZ);
    auto beampipe = std::make_unique<Acts::TrackingVolume>(
        Acts::Transform3::Identity(), std::move(beampipeBounds), "Beampipe");

    detector.addMaterial("BeampipeMaterial", [&](auto& bpMat) {
      bpMat.configureFace(OuterCylinder,
                          asHomogeneous("Beam pipe", m_beamPipeMaterial));
      bpMat.addStaticVolume(std::move(beampipe));
    });

    buildPixel(detector, gctx);

    if (m_cfg.buildLevel > 1) {
      buildShortStrip(detector, gctx);
      buildLongStrip(detector, gctx);
    }
  });

  if (graphvizFile) {
    std::ofstream file(*graphvizFile);
    root.graphviz(file);
  }

  BlueprintOptions options;
  auto trackingGeometry =
      root.construct(options, gctx, *logger().clone("Blprnt"));
  return trackingGeometry;
}

void Gen3GenericDetectorBuilder::buildPixel(
    Acts::Experimental::BlueprintNode& parent,
    const Acts::GeometryContext& gctx) {
  using enum Acts::AxisDirection;
  using namespace Acts::Experimental;
  using namespace Acts::UnitLiterals;
  using enum Acts::CylinderVolumeBounds::Face;
  using AttachmentStrategy = Acts::VolumeAttachmentStrategy;
  using ResizeStrategy = Acts::VolumeResizeStrategy;

  ACTS_DEBUG("Building Pixel");

  ProtoLayerCreator pplCreator = createPixelProtoLayerCreator();

  auto& pixel = parent.addCylinderContainer("Pixel", AxisR);

  if (m_cfg.buildLevel > 1) {
    ACTS_DEBUG("Building PST");
    auto pstVolume = std::make_unique<Acts::TrackingVolume>(
        Acts::Transform3::Identity(),
        std::make_unique<Acts::CylinderVolumeBounds>(
            kPstRadius - kPstThickness, kPstRadius, kPstHalfLengthZ),
        "PST");

    pixel.addMaterial("PSTMaterial", [&](auto& pstMat) {
      pstMat.configureFace(OuterCylinder, asHomogeneous("PST", m_pstMaterial));
      pstMat.addStaticVolume(std::move(pstVolume));
    });
  }

  auto& sensitive = pixel.addCylinderContainer("Pixel", AxisZ);
  sensitive.addCylinderContainer("Pixel_Barrel", AxisR, [&](auto& barrel) {
    barrel.setAttachmentStrategy(AttachmentStrategy::Gap);
    barrel.setResizeStrategies(ResizeStrategy::Expand, ResizeStrategy::Gap);

    auto centralProtoLayerSurfaces = pplCreator.centralProtoLayers(gctx);
    ACTS_DEBUG("Adding " << centralProtoLayerSurfaces.size()
                         << " central proto layers to "
                         << "Pixel_Barrel");
    for (const auto& [idx, temp] :

         Acts::enumerate(centralProtoLayerSurfaces)) {
      auto& pl = temp.protoLayer;
      std::string layerName = "Pixel_Barrel_L" + std::to_string(idx);
      barrel.addMaterial(layerName + "_Mat", [&](auto& mat) {
        mat.configureFace(OuterCylinder,
                          asHomogeneous("Pixel", m_pixelCentralMaterial));
        mat.addLayer(layerName, [&](auto& layer) {
          ACTS_VERBOSE("Adding layer " << layer.name());
          layer.setProtoLayer(pl);
          layer.setEnvelope(
              Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {5_mm, 5_mm}}});
          layer.setNavigationPolicyFactory(
              createNavigationPolicy(LayerType::Cylinder));

          ACTS_VERBOSE("-> Number of surfaces: " << layer.surfaces().size());
        });
      });
    }
  });

  sensitive.addCylinderContainer("Pixel_nEC", AxisZ, [&](auto& endcap) {
    endcap.setAttachmentStrategy(AttachmentStrategy::Gap);
    endcap.setResizeStrategy(ResizeStrategy::Gap);

    auto protoLayerSurfaces = pplCreator.negativeProtoLayers(gctx);
    ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                         << " negative proto layers to "
                         << "Pixel_nEC");
    for (const auto& [idx, temp] :

         Acts::enumerate(protoLayerSurfaces)) {
      auto& pl = temp.protoLayer;
      std::string layerName = "Pixel_nEC_L" + std::to_string(idx);
      endcap.addMaterial(layerName + "_Mat", [&](auto& mat) {
        mat.configureFace(NegativeDisc,
                          asHomogeneous("Pixel", m_pixelEndcapMaterial));
        mat.addLayer(layerName, [&](auto& layer) {
          layer.setNavigationPolicyFactory(
              createNavigationPolicy(LayerType::Disc));
          ACTS_VERBOSE("Adding layer " << layer.name());
          layer.setProtoLayer(pl);
          layer.setEnvelope(
              Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {5_mm, 5_mm}}});
        });
      });
    }
  });

  sensitive.addCylinderContainer("Pixel_pEC", AxisZ, [&](auto& endcap) {
    endcap.setAttachmentStrategy(AttachmentStrategy::Gap);
    endcap.setResizeStrategy(ResizeStrategy::Gap);
    auto protoLayerSurfaces = pplCreator.positiveProtoLayers(gctx);
    ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                         << " positive proto layers to "
                         << "Pixel_pEC");
    for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
      auto& pl = temp.protoLayer;
      std::string layerName = "Pixel_pEC_L" + std::to_string(idx);
      endcap.addMaterial(layerName + "_Mat", [&](auto& mat) {
        mat.configureFace(PositiveDisc,
                          asHomogeneous("Pixel", m_pixelEndcapMaterial));
        mat.addLayer(layerName, [&](auto& layer) {
          layer.setNavigationPolicyFactory(
              createNavigationPolicy(LayerType::Disc));
          ACTS_VERBOSE("Adding layer " << layer.name());
          layer.setProtoLayer(pl);
          layer.setEnvelope(
              Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {5_mm, 5_mm}}});
        });
      });
    }
  });
}

void Gen3GenericDetectorBuilder::buildShortStrip(
    Acts::Experimental::BlueprintNode& parent,
    const Acts::GeometryContext& gctx) {
  using enum Acts::AxisDirection;
  using namespace Acts::Experimental;
  using namespace Acts::UnitLiterals;
  using enum Acts::CylinderVolumeBounds::Face;
  using AttachmentStrategy = Acts::VolumeAttachmentStrategy;
  using ResizeStrategy = Acts::VolumeResizeStrategy;

  ACTS_DEBUG("Building Short Strip");

  ProtoLayerCreator ssplCreator = createShortStripProtoLayerCreator();

  auto& sstrip = parent.addCylinderContainer("ShortStrip", AxisZ);
  sstrip.setAttachmentStrategy(AttachmentStrategy::Gap);
  sstrip.setResizeStrategy(ResizeStrategy::Gap);

  sstrip.addCylinderContainer("ShortStrip_Barrel", AxisR, [&](auto& barrel) {
    barrel.setAttachmentStrategy(AttachmentStrategy::Gap);
    barrel.setResizeStrategy(ResizeStrategy::Gap);

    auto protoLayerSurfaces = ssplCreator.centralProtoLayers(gctx);
    ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                         << " central proto layers to "
                         << "ShortStrip");
    for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
      auto& pl = temp.protoLayer;
      std::string layerName = "ShortStrip_Barrel_L" + std::to_string(idx);
      barrel.addMaterial(layerName + "_Mat", [&](auto& mat) {
        mat.configureFace(
            InnerCylinder,
            asHomogeneous("ShortStrip", m_shortStripCentralMaterial));
        mat.addLayer(layerName, [&](auto& layer) {
          layer.setNavigationPolicyFactory(
              createNavigationPolicy(LayerType::Cylinder));
          ACTS_VERBOSE("Adding layer " << layer.name());
          layer.setProtoLayer(pl);
          layer.setEnvelope(
              Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {5_mm, 5_mm}}});
        });
      });
    }
  });

  if (m_cfg.buildLevel > 2) {
    sstrip.addCylinderContainer("ShortStrip_nEC", AxisZ, [&](auto& endcap) {
      endcap.setAttachmentStrategy(AttachmentStrategy::Gap);
      endcap.setResizeStrategy(ResizeStrategy::Gap);

      auto protoLayerSurfaces = ssplCreator.negativeProtoLayers(gctx);

      ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                           << " negative proto layers to "
                           << "ShortStrip_nEC");
      for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
        auto& pl = temp.protoLayer;
        std::string layerName = "ShortStrip_nEC_L" + std::to_string(idx);
        endcap.addMaterial(layerName + "_Mat", [&](auto& mat) {
          mat.configureFace(
              NegativeDisc,
              asHomogeneous("ShortStrip", m_shortStripEndcapMaterial));
          mat.addLayer(layerName, [&](auto& layer) {
            layer.setNavigationPolicyFactory(
                createNavigationPolicy(LayerType::Disc));
            ACTS_VERBOSE("Adding layer " << layer.name());
            layer.setProtoLayer(pl);
            layer.setEnvelope(
                Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {0_mm, 0_mm}}});
          });
        });
      }
    });

    sstrip.addCylinderContainer("ShortStrip_pEC", AxisZ, [&](auto& endcap) {
      endcap.setAttachmentStrategy(AttachmentStrategy::Gap);
      endcap.setResizeStrategy(ResizeStrategy::Gap);

      auto protoLayerSurfaces = ssplCreator.positiveProtoLayers(gctx);

      ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                           << " positive proto layers to "
                           << "ShortStrip_pEC");
      for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
        auto& pl = temp.protoLayer;
        std::string layerName = "ShortStrip_pEC_L" + std::to_string(idx);
        endcap.addMaterial(layerName + "_Mat", [&](auto& mat) {
          mat.configureFace(
              PositiveDisc,
              asHomogeneous("ShortStrip", m_shortStripEndcapMaterial));
          mat.addLayer(layerName, [&](auto& layer) {
            layer.setNavigationPolicyFactory(
                createNavigationPolicy(LayerType::Disc));
            ACTS_VERBOSE("Adding layer " << layer.name());
            layer.setProtoLayer(pl);
            layer.setEnvelope(
                Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {0_mm, 0_mm}}});
          });
        });
      }
    });
  }
}

void Gen3GenericDetectorBuilder::buildLongStrip(
    Acts::Experimental::BlueprintNode& parent,
    const Acts::GeometryContext& gctx) {
  using enum Acts::AxisDirection;
  using namespace Acts::Experimental;
  using namespace Acts::UnitLiterals;
  using enum Acts::CylinderVolumeBounds::Face;
  using AttachmentStrategy = Acts::VolumeAttachmentStrategy;
  using ResizeStrategy = Acts::VolumeResizeStrategy;

  ACTS_DEBUG("Building Long Strip");

  ProtoLayerCreator lsplCreator = createLongStripProtoLayerCreator();

  auto& lstrip = parent.addCylinderContainer("LongStrip", AxisZ);
  lstrip.setAttachmentStrategy(AttachmentStrategy::Gap);
  lstrip.setResizeStrategy(ResizeStrategy::Gap);

  lstrip.addCylinderContainer("LongStrip_Barrel", AxisR, [&](auto& barrel) {
    barrel.setAttachmentStrategy(AttachmentStrategy::Gap);
    barrel.setResizeStrategy(ResizeStrategy::Gap);

    auto protoLayerSurfaces = lsplCreator.centralProtoLayers(gctx);
    ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                         << " central proto layers to "
                         << "LongStrip_Barrel");
    for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
      auto& pl = temp.protoLayer;
      std::string layerName = "LongStrip_Barrel_L" + std::to_string(idx);
      barrel.addMaterial(layerName + "_Mat", [&](auto& mat) {
        mat.configureFace(
            InnerCylinder,
            asHomogeneous("LongStrip", m_longStripCentralMaterial));
        mat.addLayer(layerName, [&](auto& layer) {
          layer.setNavigationPolicyFactory(
              createNavigationPolicy(LayerType::Cylinder));
          ACTS_VERBOSE("Adding layer " << layer.name());
          layer.setProtoLayer(pl);
          layer.setEnvelope(
              Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {5_mm, 5_mm}}});
        });
      });
    }
  });

  if (m_cfg.buildLevel > 2) {
    lstrip.addCylinderContainer("LongStrip_nEC", AxisZ, [&](auto& endcap) {
      endcap.setAttachmentStrategy(AttachmentStrategy::Gap);
      endcap.setResizeStrategy(ResizeStrategy::Gap);

      auto protoLayerSurfaces = lsplCreator.negativeProtoLayers(gctx);

      ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                           << " negative proto layers to "
                           << "LongStrip_nEC");
      for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
        auto& pl = temp.protoLayer;
        std::string layerName = "LongStrip_nEC_L" + std::to_string(idx);
        endcap.addMaterial(layerName + "_Mat", [&](auto& mat) {
          mat.configureFace(
              NegativeDisc,
              asHomogeneous("LongStrip", m_longStripEndcapMaterial));
          mat.addLayer(layerName, [&](auto& layer) {
            layer.setNavigationPolicyFactory(
                createNavigationPolicy(LayerType::Disc));
            ACTS_VERBOSE("Adding layer " << layer.name());
            layer.setProtoLayer(pl);
            layer.setEnvelope(
                Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {0_mm, 0_mm}}});
          });
        });
      }
    });

    lstrip.addCylinderContainer("LongStrip_pEC", AxisZ, [&](auto& endcap) {
      endcap.setAttachmentStrategy(AttachmentStrategy::Gap);
      endcap.setResizeStrategy(ResizeStrategy::Gap);

      auto protoLayerSurfaces = lsplCreator.positiveProtoLayers(gctx);

      ACTS_DEBUG("Adding " << protoLayerSurfaces.size()
                           << " positive proto layers to "
                           << "LongStrip_pEC");
      for (const auto& [idx, temp] : Acts::enumerate(protoLayerSurfaces)) {
        auto& pl = temp.protoLayer;
        std::string layerName = "LongStrip_pEC_L" + std::to_string(idx);
        endcap.addMaterial(layerName + "_Mat", [&](auto& mat) {
          mat.configureFace(
              PositiveDisc,
              asHomogeneous("LongStrip", m_longStripEndcapMaterial));
          mat.addLayer(layerName, [&](auto& layer) {
            layer.setNavigationPolicyFactory(
                createNavigationPolicy(LayerType::Disc));
            ACTS_VERBOSE("Adding layer " << layer.name());
            layer.setProtoLayer(pl);
            layer.setEnvelope(
                Acts::ExtentEnvelope{{.z = {5_mm, 5_mm}, .r = {0_mm, 0_mm}}});
          });
        });
      }
    });
  }
}
}  // namespace
}  // namespace Generic

GenericDetector::GenericDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("GenericDetector", cfg.logLevel)),
      m_cfg(cfg) {
  m_nominalGeometryContext = Acts::GeometryContext();
  auto detectorElementFactory =
      [this](std::shared_ptr<const Acts::Transform3> transform,
             std::shared_ptr<const Acts::PlanarBounds> bounds, double thickness,
             std::shared_ptr<const Acts::ISurfaceMaterial> material)
      -> std::shared_ptr<GenericDetectorElement> {
    auto id =
        static_cast<GenericDetectorElement::Identifier>(m_detectorStore.size());
    auto detElem = std::make_shared<GenericDetectorElement>(
        id, std::move(transform), std::move(bounds), thickness,
        std::move(material));
    m_detectorStore.push_back(detElem);
    return detElem;
  };
  buildTrackingGeometry(detectorElementFactory);
}

GenericDetector::GenericDetector(const Config& cfg, NoBuildTag /*tag*/)
    : Detector(Acts::getDefaultLogger("GenericDetector", cfg.logLevel)),
      m_cfg(cfg) {}

void GenericDetector::buildTrackingGeometry(
    const Generic::ProtoLayerCreator::DetectorElementFactory&
        detectorElementFactory) {
  ACTS_INFO("Building tracking geometry");
  if (m_trackingGeometry != nullptr) {
    throw std::runtime_error("Tracking geometry already built");
  }

  Generic::GenericDetectorBuilder::Config cfg;
  cfg.detectorElementFactory = detectorElementFactory;
  cfg.protoMaterial = m_cfg.buildProto;
  cfg.layerLogLevel = m_cfg.layerLogLevel;
  cfg.buildLevel = m_cfg.buildLevel;

  if (m_cfg.gen3) {
    m_trackingGeometry =
        Generic::Gen3GenericDetectorBuilder(cfg, logger().clone())
            .buildTrackingGeometry(m_nominalGeometryContext,
                                   m_cfg.graphvizFile);
  } else {
    m_trackingGeometry =
        Generic::Gen1GenericDetectorBuilder(cfg, logger().clone())
            .buildTrackingGeometry(m_nominalGeometryContext,
                                   m_cfg.materialDecorator,
                                   m_cfg.surfaceLogLevel, m_cfg.volumeLogLevel);
  }

  if (m_trackingGeometry == nullptr) {
    throw std::runtime_error("Error building tracking geometry");
  }

  ACTS_INFO("Tracking geometry built");
}

}  // namespace ActsExamples
