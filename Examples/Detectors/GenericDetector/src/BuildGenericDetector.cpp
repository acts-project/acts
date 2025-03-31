// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"

#include <cmath>
#include <numbers>

namespace ActsExamples::Generic {

std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const Acts::GeometryContext& gctxIn,
    const ProtoLayerCreator::DetectorElementFactory& detectorElementFactory,
    std::size_t level,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator,
    bool protoMaterial, Acts::Logging::Level surfaceLLevel,
    Acts::Logging::Level layerLLevel, Acts::Logging::Level volumeLLevel) {
  using namespace Acts::UnitLiterals;

  //   auto gctx = Acts::GeometryContext::make(gctxIn);
  Acts::GeometryContext gctx{gctxIn};

  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig, Acts::getDefaultLogger("SurfaceArrayCreator", surfaceLLevel));
  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, Acts::getDefaultLogger("LayerCreator", layerLLevel));
  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, Acts::getDefaultLogger("LayerArrayCreator", layerLLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          Acts::getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig,
          Acts::getDefaultLogger("CylinderVolumeHelper", volumeLLevel));
  //-------------------------------------------------------------------------------------
  // vector of the volume builders
  std::vector<std::shared_ptr<const Acts::ITrackingVolumeBuilder>>
      volumeBuilders;

  // Prepare the proto material - in case it's designed to do so
  // - cylindrical
  Acts::BinUtility pCylinderUtility(10, -1, 1, Acts::closed,
                                    Acts::AxisDirection::AxisPhi);
  pCylinderUtility +=
      Acts::BinUtility(10, -1, 1, Acts::open, Acts::AxisDirection::AxisZ);
  auto pCylinderMaterial =
      std::make_shared<const Acts::ProtoSurfaceMaterial>(pCylinderUtility);
  // - disc
  Acts::BinUtility pDiscUtility(10, 0, 1, Acts::open,
                                Acts::AxisDirection::AxisR);
  pDiscUtility +=
      Acts::BinUtility(10, -1, 1, Acts::closed, Acts::AxisDirection::AxisPhi);
  auto pDiscMaterial =
      std::make_shared<const Acts::ProtoSurfaceMaterial>(pDiscUtility);
  // - plane
  Acts::BinUtility pPlaneUtility(1, -1, 1, Acts::open,
                                 Acts::AxisDirection::AxisX);
  auto pPlaneMaterial =
      std::make_shared<const Acts::ProtoSurfaceMaterial>(pPlaneUtility);

  //-------------------------------------------------------------------------------------
  // Beam Pipe
  //-------------------------------------------------------------------------------------
  // BeamPipe material
  const auto beryllium = Acts::Material::fromMassDensity(
      352.8_mm, 407_mm, 9.012, 4.0, 1.848_g / 1_cm3);
  std::shared_ptr<const Acts::ISurfaceMaterial> beamPipeMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          Acts::MaterialSlab(beryllium, 0.8_mm));
  if (protoMaterial) {
    beamPipeMaterial = pCylinderMaterial;
  }

  // configure the beam pipe layer builder
  Acts::PassiveLayerBuilder::Config bplConfig;
  bplConfig.layerIdentification = "BeamPipe";
  bplConfig.centralLayerRadii = std::vector<double>(1, 19.);
  bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 3000.);
  bplConfig.centralLayerThickness = std::vector<double>(1, 0.8);
  bplConfig.centralLayerMaterial = {beamPipeMaterial};
  auto beamPipeBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
      bplConfig, Acts::getDefaultLogger("BeamPipeLayerBuilder", layerLLevel));
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
          bpvConfig,
          Acts::getDefaultLogger("BeamPipeVolumeBuilder", volumeLLevel));
  // add to the list of builders
  volumeBuilders.push_back(beamPipeVolumeBuilder);

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // Pixel detector
  //-------------------------------------------------------------------------------------
  // some prep work
  // envelope for layers
  std::pair<double, double> pcEnvelope(2., 2.);

  double pCentralModuleT = 0.15;
  double pEndcapModuleT = 0.15;

  // Module material properties - X0, L0, A, Z, Rho
  // Acts::Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
  const auto silicon = Acts::Material::fromMassDensity(95.7_mm, 465.2_mm, 28.03,
                                                       14., 2.32_g / 1_cm3);
  Acts::MaterialSlab pcModuleMaterial(silicon, pCentralModuleT);
  Acts::MaterialSlab peModuleMaterial(silicon, pEndcapModuleT);
  // Layer material properties - thickness, X0, L0, A, Z, Rho
  Acts::MaterialSlab pcmbProperties(silicon, 1.5_mm);
  Acts::MaterialSlab pcmecProperties(silicon, 1.5_mm);

  // Module, central and disc material
  std::shared_ptr<const Acts::ISurfaceMaterial> pCentralMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(pcmbProperties);
  std::shared_ptr<const Acts::ISurfaceMaterial> pEndcapMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(pcmecProperties);
  std::shared_ptr<const Acts::ISurfaceMaterial> pCentralModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          pcModuleMaterial);
  std::shared_ptr<const Acts::ISurfaceMaterial> pEndcapModuleMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          peModuleMaterial);
  if (protoMaterial) {
    pCentralMaterial = pCylinderMaterial;
    pCentralModuleMaterial = pPlaneMaterial;
    pEndcapMaterial = pDiscMaterial;
    pEndcapModuleMaterial = pPlaneMaterial;
  }

  // configure the pixel proto layer builder
  typename ProtoLayerCreator::Config pplConfig;
  pplConfig.detectorElementFactory = detectorElementFactory;

  // standard, an approach envelope
  pplConfig.approachSurfaceEnvelope = 1.;
  // BARREL :
  // 4 pixel layers
  // configure the central barrel
  pplConfig.centralLayerBinMultipliers = {1, 1};
  pplConfig.centralLayerRadii = {32., 72., 116., 172.};
  pplConfig.centralLayerEnvelopes = {pcEnvelope, pcEnvelope, pcEnvelope,
                                     pcEnvelope};
  pplConfig.centralModuleBinningSchema = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};
  pplConfig.centralModuleTiltPhi = {0.14, 0.14, 0.14, 0.14};
  pplConfig.centralModuleHalfX = {8.4, 8.4, 8.4, 8.4};
  pplConfig.centralModuleHalfY = {36., 36., 36., 36.};
  pplConfig.centralModuleThickness = {pCentralModuleT, pCentralModuleT,
                                      pCentralModuleT, pCentralModuleT};
  pplConfig.centralModuleMaterial = {
      pCentralModuleMaterial, pCentralModuleMaterial, pCentralModuleMaterial,
      pCentralModuleMaterial};

  // no frontside/backside
  pplConfig.centralModuleFrontsideStereo = {};
  pplConfig.centralModuleBacksideStereo = {};
  pplConfig.centralModuleBacksideGap = {};
  // mPositions
  std::vector<std::vector<Acts::Vector3>> pplCentralModulePositions;
  for (std::size_t plb = 0; plb < pplConfig.centralLayerRadii.size(); ++plb) {
    // call the helper function
    pplCentralModulePositions.push_back(
        modulePositionsCylinder(pplConfig.centralLayerRadii[plb],
                                0.5,  // 1 mm stagger
                                pplConfig.centralModuleHalfY[plb],
                                2.,  // 4 mm module overlap in z
                                pplConfig.centralModuleBinningSchema[plb]));
  }
  pplConfig.centralModulePositions = pplCentralModulePositions;
  // ENDCAP :
  // 7 pixel layers each side
  // configure the endcaps
  pplConfig.posnegLayerBinMultipliers = {1, 1};

  pplConfig.posnegLayerPositionsZ = {
      600. * Acts::UnitConstants::mm, 700. * Acts::UnitConstants::mm,
      820. * Acts::UnitConstants::mm, 960. * Acts::UnitConstants::mm,
      1100 * Acts::UnitConstants::mm, 1300 * Acts::UnitConstants::mm,
      1500 * Acts::UnitConstants::mm};

  pplConfig.posnegLayerEnvelopeR = {
      1. * Acts::UnitConstants::mm, 1. * Acts::UnitConstants::mm,
      1. * Acts::UnitConstants::mm, 1. * Acts::UnitConstants::mm,
      1. * Acts::UnitConstants::mm, 1. * Acts::UnitConstants::mm,
      1. * Acts::UnitConstants::mm};
  std::vector<double> perHX = {8.4, 8.4};     // half length x
  std::vector<double> perHY = {36., 36.};     // half length y
  std::vector<std::size_t> perBP = {40, 68};  // bins in phi
  std::vector<double> perT = {pEndcapModuleT,
                              pEndcapModuleT};    // module thickness
  std::vector<std::size_t> perBX = {336, 336};    // bins in x
  std::vector<std::size_t> perBY = {1280, 1280};  // bins in y
  std::vector<int> perRS = {-1, -1};              // readout side
  std::vector<double> perLA = {0., 0.};           // lorentz angle
  std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>> perM = {
      pEndcapModuleMaterial, pEndcapModuleMaterial};  // material

  pplConfig.posnegModuleMinHalfX = std::vector<std::vector<double>>(7, perHX);
  pplConfig.posnegModuleMaxHalfX = {};
  pplConfig.posnegModuleHalfY = std::vector<std::vector<double>>(7, perHY);
  pplConfig.posnegModulePhiBins =
      std::vector<std::vector<std::size_t>>(7, perBP);
  pplConfig.posnegModuleThickness = std::vector<std::vector<double>>(7, perT);
  pplConfig.posnegModuleMaterial =
      std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>(
          7, perM);

  // no frontside/backside
  pplConfig.posnegModuleFrontsideStereo = {};
  pplConfig.posnegModuleBacksideStereo = {};
  pplConfig.posnegModuleBacksideGap = {};
  // mPositions
  std::vector<std::vector<std::vector<Acts::Vector3>>> pplPosnegModulePositions;
  for (std::size_t id = 0; id < pplConfig.posnegLayerPositionsZ.size(); ++id) {
    pplPosnegModulePositions.push_back(modulePositionsDisc(
        pplConfig.posnegLayerPositionsZ[id], 0.0, {4.0, 4.0}, {0.5, 0.}, 30.,
        176., pplConfig.posnegModulePhiBins[id],
        pplConfig.posnegModuleHalfY[id]));
  }
  pplConfig.posnegModulePositions = pplPosnegModulePositions;

  /// The ProtoLayer creator
  ProtoLayerCreator pplCreator(
      pplConfig, Acts::getDefaultLogger("PixelProtoLayerCreator", layerLLevel));

  // configure pixel layer builder
  typename LayerBuilder::Config plbConfig;
  plbConfig.layerCreator = layerCreator;
  plbConfig.layerIdentification = "Pixel";
  // material concentration alsways outside the modules
  plbConfig.centralProtoLayers = pplCreator.centralProtoLayers(gctx);
  plbConfig.centralLayerMaterialConcentration = {1, 1, 1, 1};
  plbConfig.centralLayerMaterial = {pCentralMaterial, pCentralMaterial,
                                    pCentralMaterial, pCentralMaterial};
  if (level > 0) {
    // material concentration is always behind the layer in the pixels
    plbConfig.posnegLayerMaterialConcentration = std::vector<int>(7, 0);
    // layer structure surface has pixel material properties
    plbConfig.posnegLayerMaterial = {
        pEndcapMaterial, pEndcapMaterial, pEndcapMaterial, pEndcapMaterial,
        pEndcapMaterial, pEndcapMaterial, pEndcapMaterial};
    // negative proto layers
    plbConfig.negativeProtoLayers = pplCreator.negativeProtoLayers(gctx);
    plbConfig.positiveProtoLayers = pplCreator.positiveProtoLayers(gctx);
  }
  // define the builder
  auto pixelLayerBuilder = std::make_shared<const LayerBuilder>(
      plbConfig, Acts::getDefaultLogger("PixelLayerBuilder", layerLLevel));
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
      pvbConfig, Acts::getDefaultLogger("PixelVolumeBuilder", volumeLLevel));
  // add to the list of builders
  volumeBuilders.push_back(pixelVolumeBuilder);

  if (level > 1) {
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------
    // Pixel Support Tybe (PST)
    //-------------------------------------------------------------------------------------
    // Material
    std::shared_ptr<const Acts::ISurfaceMaterial> pstMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            Acts::MaterialSlab(beryllium, 1.8_mm));
    if (protoMaterial) {
      pstMaterial = pCylinderMaterial;
    }

    // Configuration
    Acts::PassiveLayerBuilder::Config pstConfig;
    pstConfig.layerIdentification = "PST";
    pstConfig.centralLayerRadii = std::vector<double>(1, 200.);
    pstConfig.centralLayerHalflengthZ = std::vector<double>(1, 2800.);
    pstConfig.centralLayerThickness = std::vector<double>(1, 1.8);
    pstConfig.centralLayerMaterial = {pstMaterial};
    auto pstBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        pstConfig, Acts::getDefaultLogger("PSTBuilder", layerLLevel));
    // create the volume for the beam pipe
    Acts::CylinderVolumeBuilder::Config pstvolConfig;
    pstvolConfig.trackingVolumeHelper = cylinderVolumeHelper;
    pstvolConfig.volumeName = "PST";
    pstvolConfig.buildToRadiusZero = false;
    pstvolConfig.layerBuilder = pstBuilder;
    auto pstVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        pstvolConfig, Acts::getDefaultLogger("PSTVolumeBuilder", volumeLLevel));
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

    double ssCentralModuleT = 0.25;
    double ssEndcapModuleT = 0.25;
    // envelope double
    std::pair<double, double> ssEnvelope(2., 2.);

    // Module material properties - X0, L0, A, Z, Rho
    // Acts::Material sscMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
    Acts::MaterialSlab sscModuleMaterial(silicon, ssCentralModuleT);
    Acts::MaterialSlab sseModuleMaterial(silicon, ssEndcapModuleT);

    // Layer material properties - thickness, X0, L0, A, Z, Rho
    Acts::MaterialSlab ssbmProperties(silicon, 2_mm);
    Acts::MaterialSlab ssecmProperties(silicon, 2.5_mm);

    // Module, central and disc material
    std::shared_ptr<const Acts::ISurfaceMaterial> ssCentralMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            ssbmProperties);
    std::shared_ptr<const Acts::ISurfaceMaterial> ssEndcapMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            ssecmProperties);
    std::shared_ptr<const Acts::ISurfaceMaterial> ssCentralModuleMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            sscModuleMaterial);
    std::shared_ptr<const Acts::ISurfaceMaterial> ssEndcapModuleMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            sseModuleMaterial);
    if (protoMaterial) {
      ssCentralMaterial = pCylinderMaterial;
      ssCentralModuleMaterial = pPlaneMaterial;
      ssEndcapMaterial = pDiscMaterial;
      ssEndcapModuleMaterial = pPlaneMaterial;
    }

    // ----------------------------------------------------------------------------
    // Configure the short strip proto layer builder
    typename ProtoLayerCreator::Config ssplConfig;
    ssplConfig.detectorElementFactory = detectorElementFactory;

    // configure the central barrel
    ssplConfig.centralLayerBinMultipliers = {1, 1};
    ssplConfig.centralLayerRadii = {260., 360., 500., 660.};
    ssplConfig.centralLayerEnvelopes = {ssEnvelope, ssEnvelope, ssEnvelope,
                                        ssEnvelope};

    ssplConfig.centralModuleBinningSchema = {
        {40, 21}, {56, 21}, {78, 21}, {102, 21}};
    ssplConfig.centralModuleTiltPhi = {-0.15, -0.15, -0.15, -0.15};
    ssplConfig.centralModuleHalfX = {24., 24., 24., 24.};
    ssplConfig.centralModuleHalfY = {54., 54., 54., 54.};
    ssplConfig.centralModuleThickness = {ssCentralModuleT, ssCentralModuleT,
                                         ssCentralModuleT, ssCentralModuleT};

    ssplConfig.centralModuleMaterial = {
        ssCentralModuleMaterial, ssCentralModuleMaterial,
        ssCentralModuleMaterial, ssCentralModuleMaterial};
    ssplConfig.centralModuleFrontsideStereo = {};
    ssplConfig.centralModuleBacksideStereo = {};
    ssplConfig.centralModuleBacksideGap = {};
    // mPositions
    std::vector<std::vector<Acts::Vector3>> ssplCentralModulePositions;
    for (std::size_t sslb = 0; sslb < ssplConfig.centralLayerRadii.size();
         ++sslb) {
      // call the helper function
      ssplCentralModulePositions.push_back(
          modulePositionsCylinder(ssplConfig.centralLayerRadii[sslb],
                                  3.,  // 3 mm stagger
                                  ssplConfig.centralModuleHalfY[sslb],
                                  5.,  // 5 mm module overlap
                                  ssplConfig.centralModuleBinningSchema[sslb]));
    }
    ssplConfig.centralModulePositions = ssplCentralModulePositions;

    // configure the endcaps
    std::vector<double> mrMinHx = {16.4, 24.2, 32.2};
    std::vector<double> mrMaxHx = {24.2, 32.2, 40.0};
    std::vector<double> mrHy = {78., 78., 78.};

    std::vector<std::size_t> mPhiBins = {54, 56, 60};
    std::vector<double> mThickness = {ssEndcapModuleT, ssEndcapModuleT,
                                      ssEndcapModuleT};
    std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>> mMaterial = {
        ssEndcapModuleMaterial, ssEndcapModuleMaterial, ssEndcapModuleMaterial};

    ssplConfig.posnegLayerBinMultipliers = {1, 2};

    ssplConfig.posnegLayerPositionsZ = {1220., 1500., 1800.,
                                        2150., 2550., 2950.};
    std::size_t nposnegs = ssplConfig.posnegLayerPositionsZ.size();
    ssplConfig.posnegLayerEnvelopeR = std::vector<double>(nposnegs, 5.);

    ssplConfig.posnegModuleMinHalfX =
        std::vector<std::vector<double>>(nposnegs, mrMinHx);
    ssplConfig.posnegModuleMaxHalfX =
        std::vector<std::vector<double>>(nposnegs, mrMaxHx);
    ssplConfig.posnegModuleHalfY =
        std::vector<std::vector<double>>(nposnegs, mrHy);
    ssplConfig.posnegModulePhiBins =
        std::vector<std::vector<std::size_t>>(nposnegs, mPhiBins);
    ssplConfig.posnegModuleThickness =
        std::vector<std::vector<double>>(nposnegs, mThickness);

    ssplConfig.posnegModuleMaterial =
        std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>(
            nposnegs, mMaterial);

    ssplConfig.posnegModuleFrontsideStereo = {};
    ssplConfig.posnegModuleBacksideStereo = {};
    ssplConfig.posnegModuleBacksideGap = {};

    // mPositions
    std::vector<std::vector<std::vector<Acts::Vector3>>>
        ssplPosnegModulePositions;
    for (std::size_t id = 0; id < ssplConfig.posnegLayerPositionsZ.size();
         ++id) {
      ssplPosnegModulePositions.push_back(modulePositionsDisc(
          ssplConfig.posnegLayerPositionsZ[id], 6.0, {3., 3., 3.}, {0., 0., 0.},
          240., 700., ssplConfig.posnegModulePhiBins[id],
          ssplConfig.posnegModuleHalfY[id]));
    }
    ssplConfig.posnegModulePositions = ssplPosnegModulePositions;

    // The ProtoLayer creator
    ProtoLayerCreator ssplCreator(
        ssplConfig,
        Acts::getDefaultLogger("SStripProtoLayerCreator", layerLLevel));

    // configure short strip layer builder
    typename LayerBuilder::Config sslbConfig;
    sslbConfig.layerCreator = layerCreator;
    sslbConfig.layerIdentification = "SStrip";

    sslbConfig.centralProtoLayers = ssplCreator.centralProtoLayers(gctx);
    sslbConfig.centralLayerMaterialConcentration = {-1, -1, -1, -1};
    sslbConfig.centralLayerMaterial = {ssCentralMaterial, ssCentralMaterial,
                                       ssCentralMaterial, ssCentralMaterial};

    if (level > 2) {
      sslbConfig.negativeProtoLayers = ssplCreator.negativeProtoLayers(gctx);
      sslbConfig.positiveProtoLayers = ssplCreator.positiveProtoLayers(gctx);

      sslbConfig.posnegLayerMaterialConcentration =
          std::vector<int>(nposnegs, 0);
      sslbConfig.posnegLayerMaterial =
          std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>(
              nposnegs, ssEndcapMaterial);
    }

    // define the builder
    auto sstripLayerBuilder = std::make_shared<const LayerBuilder>(
        sslbConfig, Acts::getDefaultLogger("SStripLayerBuilder", layerLLevel));
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config ssvbConfig;
    ssvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    ssvbConfig.volumeName = "SStrip";
    ssvbConfig.buildToRadiusZero = false;
    ssvbConfig.layerBuilder = sstripLayerBuilder;
    auto sstripVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            ssvbConfig,
            Acts::getDefaultLogger("SStripVolumeBuilder", volumeLLevel));

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
    // envelope double
    std::pair<double, double> lsEnvelope(2., 2.);

    double lsCentralModuleT = 0.35;
    double lsEndcapModuleT = 0.35;

    // Module material properties - X0, L0, A, Z, Rho
    // Acts::Material lsMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
    Acts::MaterialSlab lscModuleMaterial(silicon, lsCentralModuleT);
    Acts::MaterialSlab lseModuleMaterial(silicon, lsEndcapModuleT);

    // Layer material properties - thickness, X0, L0, A, Z, Rho - barrel
    Acts::MaterialSlab lsbmProperties(silicon, 2.5_mm);
    Acts::MaterialSlab lsecmProperties(silicon, 3.5_mm);

    // Module, central and disc material
    std::shared_ptr<const Acts::ISurfaceMaterial> lsCentralMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            lsbmProperties);
    std::shared_ptr<const Acts::ISurfaceMaterial> lsEndcapMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            lsecmProperties);
    std::shared_ptr<const Acts::ISurfaceMaterial> lsCentralModuleMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            lscModuleMaterial);
    std::shared_ptr<const Acts::ISurfaceMaterial> lsEndcapModuleMaterial =
        std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
            lseModuleMaterial);
    if (protoMaterial) {
      lsCentralMaterial = pCylinderMaterial;
      lsCentralModuleMaterial = pPlaneMaterial;
      lsEndcapMaterial = pDiscMaterial;
      lsEndcapModuleMaterial = pPlaneMaterial;
    }

    // The proto layer creator
    typename ProtoLayerCreator::Config lsplConfig;
    lsplConfig.detectorElementFactory = detectorElementFactory;

    // configure the central barrel
    lsplConfig.centralLayerBinMultipliers = {1, 1};
    lsplConfig.centralLayerRadii = {820., 1020.};
    lsplConfig.centralLayerEnvelopes = {lsEnvelope, lsEnvelope};

    lsplConfig.centralModuleBinningSchema = {{120, 21}, {152, 21}};
    lsplConfig.centralModuleTiltPhi = {-0.15, -0.15};
    lsplConfig.centralModuleHalfX = {24., 24.};
    lsplConfig.centralModuleHalfY = {54., 54.};
    lsplConfig.centralModuleThickness = {lsCentralModuleT, lsCentralModuleT};
    lsplConfig.centralModuleMaterial = {lsCentralModuleMaterial,
                                        lsCentralModuleMaterial};

    lsplConfig.centralModuleFrontsideStereo = {};
    lsplConfig.centralModuleBacksideStereo = {};
    lsplConfig.centralModuleBacksideGap = {};
    // mPositions
    std::vector<std::vector<Acts::Vector3>> lslbCentralModulePositions;
    for (std::size_t lslb = 0; lslb < lsplConfig.centralLayerRadii.size();
         ++lslb) {
      // call the helper function
      lslbCentralModulePositions.push_back(
          modulePositionsCylinder(lsplConfig.centralLayerRadii[lslb],
                                  3.,  // 3 mm stagger
                                  lsplConfig.centralModuleHalfY[lslb],
                                  5.,  // 5 mm module overlap
                                  lsplConfig.centralModuleBinningSchema[lslb]));
    }

    lsplConfig.centralModulePositions = lslbCentralModulePositions;
    // configure the endcaps
    mrMinHx = {54., 66.};
    mrMaxHx = {64.2, 72.};
    mrHy = {78., 78.};
    mPhiBins = {48, 50};
    mThickness = {lsEndcapModuleT, lsEndcapModuleT};
    mMaterial = {lsEndcapModuleMaterial, lsEndcapModuleMaterial};

    // endcap
    lsplConfig.posnegLayerBinMultipliers = {1, 2};
    lsplConfig.posnegLayerPositionsZ = {1220., 1500., 1800.,
                                        2150., 2550., 2950.};
    nposnegs = lsplConfig.posnegLayerPositionsZ.size();
    lsplConfig.posnegLayerEnvelopeR = std::vector<double>(nposnegs, 5.);

    lsplConfig.posnegModuleMinHalfX =
        std::vector<std::vector<double>>(nposnegs, mrMinHx);
    lsplConfig.posnegModuleMaxHalfX =
        std::vector<std::vector<double>>(nposnegs, mrMaxHx);
    lsplConfig.posnegModuleHalfY =
        std::vector<std::vector<double>>(nposnegs, mrHy);
    lsplConfig.posnegModulePhiBins =
        std::vector<std::vector<std::size_t>>(nposnegs, mPhiBins);
    lsplConfig.posnegModuleThickness =
        std::vector<std::vector<double>>(nposnegs, mThickness);

    lsplConfig.posnegModuleMaterial =
        std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>(
            nposnegs, mMaterial);
    lsplConfig.posnegModuleFrontsideStereo = {};
    lsplConfig.posnegModuleBacksideStereo = {};
    lsplConfig.posnegModuleBacksideGap = {};

    // mPositions
    std::vector<std::vector<std::vector<Acts::Vector3>>>
        lssbPosnegModulePositions;
    for (std::size_t id = 0; id < lsplConfig.posnegLayerPositionsZ.size();
         ++id) {
      lssbPosnegModulePositions.push_back(modulePositionsDisc(
          lsplConfig.posnegLayerPositionsZ[id],
          8.0,  // staggering of rings, we put the disk structure in between
          {3., 3.}, {0., 0.}, 750., 1020., lsplConfig.posnegModulePhiBins[id],
          lsplConfig.posnegModuleHalfY[id]));
    }
    lsplConfig.posnegModulePositions = lssbPosnegModulePositions;

    // The ProtoLayer creator
    ProtoLayerCreator lsplCreator(
        lsplConfig,
        Acts::getDefaultLogger("LStripProtoLayerCreator", layerLLevel));

    // configure short strip layer builder
    typename LayerBuilder::Config lslbConfig;
    lslbConfig.layerCreator = layerCreator;
    lslbConfig.layerIdentification = "LStrip";
    lslbConfig.centralLayerMaterialConcentration = {-1, -1};
    lslbConfig.centralLayerMaterial = {lsCentralMaterial, lsCentralMaterial};
    lslbConfig.centralProtoLayers = lsplCreator.centralProtoLayers(gctx);

    if (level > 2) {
      lslbConfig.posnegLayerMaterialConcentration =
          std::vector<int>(nposnegs, 0);
      lslbConfig.posnegLayerMaterial =
          std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>(
              nposnegs, lsEndcapMaterial);
      lslbConfig.negativeProtoLayers = lsplCreator.negativeProtoLayers(gctx);
      lslbConfig.positiveProtoLayers = lsplCreator.positiveProtoLayers(gctx);
    }

    // define the builder
    auto lstripLayerBuilder = std::make_shared<const LayerBuilder>(
        lslbConfig, Acts::getDefaultLogger("LStripLayerBuilder", layerLLevel));
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config lsvbConfig;
    lsvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    lsvbConfig.volumeName = "LStrip";
    lsvbConfig.buildToRadiusZero = false;
    lsvbConfig.layerBuilder = lstripLayerBuilder;
    auto lstripVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            lsvbConfig,
            Acts::getDefaultLogger("LStripVolumeBuilder", volumeLLevel));
    // add to the list of builders
    volumeBuilders.push_back(lstripVolumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  Acts::TrackingGeometryBuilder::Config tgConfig;
  // Add the build call functions
  for (auto& vb : volumeBuilders) {
    tgConfig.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto&) {
          return vb->trackingVolume(context, inner);
        });
  }
  tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  tgConfig.materialDecorator = std::move(matDecorator);

  auto cylinderGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(
          tgConfig,
          Acts::getDefaultLogger("TrackerGeometryBuilder", volumeLLevel));
  // get the geometry
  auto trackingGeometry = cylinderGeometryBuilder->trackingGeometry(gctx);
  // return the tracking geometry
  return trackingGeometry;
}

/// helper method for cylinder
std::vector<Acts::Vector3> modulePositionsCylinder(
    double radius, double zStagger, double moduleHalfLength, double lOverlap,
    const std::pair<int, int>& binningSchema) {
  int nPhiBins = binningSchema.first;
  int nZbins = binningSchema.second;
  // prepare the return value
  std::vector<Acts::Vector3> mPositions;
  mPositions.reserve(nPhiBins * nZbins);
  // prep work
  double phiStep = 2 * std::numbers::pi / nPhiBins;
  double minPhi = -std::numbers::pi + 0.5 * phiStep;
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
      mPositions.push_back(Acts::Vector3(moduleR * cos(modulePhi),
                                         moduleR * sin(modulePhi), moduleZ));
    }
  }
  return mPositions;
}

/// helper method for disc
std::vector<std::vector<Acts::Vector3>> modulePositionsDisc(
    double z, double ringStagger, std::vector<double> phiStagger,
    std::vector<double> phiSubStagger, double innerRadius, double outerRadius,
    const std::vector<std::size_t>& discBinning,
    const std::vector<double>& moduleHalfLength) {
  // calculate the radii
  std::vector<double> radii;
  // the radial span of the disc
  double deltaR = outerRadius - innerRadius;
  // quick exits
  if (discBinning.size() == 1) {
    radii.push_back(0.5 * (innerRadius + outerRadius));
  } else {
    double totalLength = 0;
    // sum up the total length
    for (auto& mhlength : moduleHalfLength) {
      totalLength += 2 * mhlength;
    }
    // now calculate the overlap (equal pay)
    double rOverlap = (totalLength - deltaR) / (moduleHalfLength.size() - 1);
    // and now fill the radii and gaps
    double lastR = innerRadius;
    double lastHl = 0.;
    double lastOl = 0.;
    // now calculate
    for (auto& mhlength : moduleHalfLength) {
      // calculate the radius
      radii.push_back(lastR + lastHl - lastOl + mhlength);
      lastR = radii[radii.size() - 1];
      lastOl = rOverlap;
      lastHl = mhlength;
    }
  }
  // now prepare the return method
  std::vector<std::vector<Acts::Vector3>> mPositions;
  for (std::size_t ir = 0; ir < radii.size(); ++ir) {
    // generate the z value
    // convention inner ring is closer to origin : makes sense
    double rz = radii.size() == 1 ? z
                                  : ((ir % 2) != 0u ? z + 0.5 * ringStagger
                                                    : z - 0.5 * ringStagger);
    // fill the ring positions
    double psStagger = phiSubStagger.empty() ? 0. : phiSubStagger[ir];
    mPositions.push_back(modulePositionsRing(rz, radii[ir], phiStagger[ir],
                                             psStagger, discBinning[ir]));
  }
  return mPositions;
}

/// Helper method for positioning
std::vector<Acts::Vector3> modulePositionsRing(double z, double radius,
                                               double phiStagger,
                                               double phiSubStagger,
                                               int nPhiBins) {
  // create and fill the positions
  std::vector<Acts::Vector3> rPositions;
  rPositions.reserve(nPhiBins);
  // prep work
  double phiStep = 2 * std::numbers::pi / nPhiBins;
  double minPhi = -std::numbers::pi + 0.5 * phiStep;
  // phi loop
  for (std::size_t iphi = 0; iphi < static_cast<std::size_t>(nPhiBins);
       ++iphi) {
    // if we have a phi sub stagger presents
    double rzs = 0.;
    // phi stagger affects 0 vs 1, 2 vs 3 ... etc
    // -> only works if it is a %4
    // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
    if (phiSubStagger != 0. && ((nPhiBins % 4) == 0)) {
      // switch sides
      if ((iphi % 4) == 0u) {
        rzs = phiSubStagger;
      } else if (((iphi + 1) % 4) == 0u) {
        rzs = -phiSubStagger;
      }
    }
    // the module phi
    double phi = minPhi + iphi * phiStep;
    // main z position depending on phi bin
    double rz = (iphi % 2) != 0u ? z - 0.5 * phiStagger : z + 0.5 * phiStagger;
    rPositions.push_back(
        Acts::Vector3(radius * cos(phi), radius * sin(phi), rz + rzs));
  }
  return rPositions;
}

}  // namespace ActsExamples::Generic
