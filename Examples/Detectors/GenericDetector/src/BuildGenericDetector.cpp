// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <numbers>

#include "./GenericDetectorBuilder.hpp"

using namespace Acts::UnitLiterals;

namespace ActsExamples {

namespace Generic {

class Gen1GenericDetectorBuilder : public GenericDetectorBuilder {
 public:
  using GenericDetectorBuilder::GenericDetectorBuilder;

  std::unique_ptr<const Acts::TrackingGeometry> buildTrackingGeometry(
      const Acts::GeometryContext& gctx, int level,
      std::shared_ptr<const Acts::IMaterialDecorator> matDecorator,
      Acts::Logging::Level surfaceLLevel, Acts::Logging::Level layerLLevel,
      Acts::Logging::Level volumeLLevel);
};

std::unique_ptr<const Acts::TrackingGeometry>
Gen1GenericDetectorBuilder::buildTrackingGeometry(
    const Acts::GeometryContext& gctx, int level,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator,
    Acts::Logging::Level surfaceLLevel, Acts::Logging::Level layerLLevel,
    Acts::Logging::Level volumeLLevel) {
  using namespace Acts::UnitLiterals;

  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig, logger().clone("SurfaceArrayCreator", surfaceLLevel));
  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, logger().clone("LayerCreator", layerLLevel));
  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, logger().clone("LayerArrayCreator", layerLLevel));
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

  //-------------------------------------------------------------------------------------
  // Beam Pipe
  //-------------------------------------------------------------------------------------

  // configure the beam pipe layer builder
  Acts::PassiveLayerBuilder::Config bplConfig;
  bplConfig.layerIdentification = "BeamPipe";
  bplConfig.centralLayerRadii = std::vector<double>(1, 19.);
  bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 3000.);
  bplConfig.centralLayerThickness = std::vector<double>(1, 0.8);
  bplConfig.centralLayerMaterial = {m_beamPipeMaterial};
  auto beamPipeBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
      bplConfig, logger().clone("BeamPipeLayerBuilder", layerLLevel));
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
  if (level > 0) {
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
      plbConfig, logger().clone("PixelLayerBuilder", layerLLevel));
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

  if (level > 1) {
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------
    // Pixel Support Tybe (PST)
    //-------------------------------------------------------------------------------------

    // Configuration
    Acts::PassiveLayerBuilder::Config pstConfig;
    pstConfig.layerIdentification = "PST";
    pstConfig.centralLayerRadii = std::vector<double>(1, 200.);
    pstConfig.centralLayerHalflengthZ = std::vector<double>(1, 2800.);
    pstConfig.centralLayerThickness = std::vector<double>(1, 1.8);
    pstConfig.centralLayerMaterial = {m_pstMaterial};
    auto pstBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        pstConfig, logger().clone("PSTLayerBuilder", layerLLevel));
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

    if (level > 2) {
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
        sslbConfig, logger().clone("SStripLayerBuilder", layerLLevel));
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
    //-------------------------------------------------------------------------------------
    // LONG strip detector
    //-------------------------------------------------------------------------------------
    // fill necessary vectors for configuration
    //-------------------------------------------------------------------------------------

    // some prep work

    ProtoLayerCreator lsplCreator = createLongStripProtoLayerCreator();

    // configure short strip layer builder
    LayerBuilder::Config lslbConfig;
    lslbConfig.layerCreator = layerCreator;
    lslbConfig.layerIdentification = "LStrip";
    lslbConfig.centralLayerMaterialConcentration = {-1, -1};
    lslbConfig.centralLayerMaterial = {m_longStripCentralMaterial,
                                       m_longStripCentralMaterial};
    lslbConfig.centralProtoLayers = lsplCreator.centralProtoLayers(gctx);

    if (level > 2) {
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
        lslbConfig, logger().clone("LStripLayerBuilder", layerLLevel));
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
  Gen3GenericDetectorBuilder(const Config& cfg,
                             std::unique_ptr<const Acts::Logger> logger =
                                 Acts::getDefaultLogger("Gen3GenDetBldr",
                                                        Acts::Logging::INFO))
      : GenericDetectorBuilder(cfg), m_logger(std::move(logger)) {}

  std::unique_ptr<const Acts::TrackingGeometry> buildTrackingGeometry(
      const Acts::GeometryContext& gctx);

  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }
};

std::unique_ptr<const Acts::TrackingGeometry>
Gen3GenericDetectorBuilder::buildTrackingGeometry(
    const Acts::GeometryContext& /*gctx*/) {
  ACTS_INFO("GenericDetector construction in  Gen3 mode");

  using namespace Acts::Experimental;
  Blueprint::Config cfg;
  Blueprint root{cfg};

  // @TODO: Add beampipe passive layer
  // @TODO: Add Pixel

  ProtoLayerCreator pplCreator = createPixelProtoLayerCreator();

  // @TODO: Add Pixel Support Tube
  // @TODO: Add Short Strip
  // @TODO: Add Long Strip

  return nullptr;
}

}  // namespace Generic

std::unique_ptr<const Acts::TrackingGeometry> Generic::buildDetector(
    const Acts::GeometryContext& gctxIn,
    const ProtoLayerCreator::DetectorElementFactory& detectorElementFactory,
    std::size_t level,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator,
    bool protoMaterial, Acts::Logging::Level surfaceLLevel,
    Acts::Logging::Level layerLLevel, Acts::Logging::Level volumeLLevel,
    bool gen3) {
  Gen1GenericDetectorBuilder::Config cfg;
  cfg.detectorElementFactory = detectorElementFactory;
  cfg.protoMaterial = protoMaterial;
  cfg.layerLogLevel = layerLLevel;

  if (gen3) {
    return Gen3GenericDetectorBuilder(cfg).buildTrackingGeometry(gctxIn);
  } else {
    return Gen1GenericDetectorBuilder(cfg).buildTrackingGeometry(
        gctxIn, level, std::move(matDecorator), surfaceLLevel, layerLLevel,
        volumeLLevel);
  }
}

}  // namespace ActsExamples
