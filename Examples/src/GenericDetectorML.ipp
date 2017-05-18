// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//-------------------------------------------------------------------------------------
// Beam Pipe
//-------------------------------------------------------------------------------------
// configure the beam pipe layer builder
PassiveLayerBuilder::Config bplConfig;
bplConfig.layerIdentification     = "BeamPipe";
bplConfig.centralLayerRadii       = std::vector<double>(1, 19.);
bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 3000.);
bplConfig.centralLayerThickness   = std::vector<double>(1, 0.8);
bplConfig.centralLayerMaterial = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
auto beamPipeBuilder              = std::make_shared<const PassiveLayerBuilder>(
    bplConfig,
    getDefaultLogger("BeamPipeLayerBuilder", layerLLevel));
// create the volume for the beam pipe
CylinderVolumeBuilder::Config bpvConfig;
bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
bpvConfig.volumeName           = "BeamPipe";
bpvConfig.layerBuilder         = beamPipeBuilder;
bpvConfig.layerEnvelopeR       = {1. * Acts::units::_mm, 1. * Acts::units::_mm};
bpvConfig.buildToRadiusZero    = true;
bpvConfig.volumeSignature      = 0;
auto beamPipeVolumeBuilder     = std::make_shared<const CylinderVolumeBuilder>(
    bpvConfig,
    getDefaultLogger("BeamPipeVolumeBuilder", volumeLLevel));
// add to the list of builders
volumeBuilders.push_back(beamPipeVolumeBuilder);
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
// Pixel detector
//-------------------------------------------------------------------------------------
// configure pixel layer builder
GenericLayerBuilder::Config plbConfig;
plbConfig.layerCreator        = layerCreator;
plbConfig.layerIdentification = "Pixel";
// fill necessary vectors for configuration
//-------------------------------------------------------------------------------------
// some prep work
// envelope for layers
std::pair<double, double> pcEnvelope(2., 2.);
// Layer material properties - thickness, X0, L0, A, Z, Rho
MaterialProperties pcmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
// Module material - X0, L0, A, Z, Rho
Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

// standard, an approach envelope
plbConfig.approachSurfaceEnvelope = 1.;

// BARREL :
// 4 pixel layers
// configure the central barrel
plbConfig.centralLayerBinMultipliers = {1, 1};
plbConfig.centralLayerRadii          = {32., 72., 116., 172.};
plbConfig.centralLayerEnvelopes
    = {pcEnvelope, pcEnvelope, pcEnvelope, pcEnvelope};
// material concentration alsways outside the modules
plbConfig.centralLayerMaterialConcentration = {1, 1, 1, 1};
plbConfig.centralLayerMaterialProperties
    = {pcmProperties, pcmProperties, pcmProperties, pcmProperties};
plbConfig.centralModuleBinningSchema = {{16, 14}, {32, 14}, {52, 14}, {78, 14}};
plbConfig.centralModuleTiltPhi       = {0.14, 0.14, 0.14, 0.14};
plbConfig.centralModuleHalfX         = {8.4, 8.4, 8.4, 8.4};
plbConfig.centralModuleHalfY         = {36., 36., 36., 36.};
plbConfig.centralModuleThickness     = {0.15, 0.15, 0.15, 0.15};
plbConfig.centralModuleMaterial
    = {pcMaterial, pcMaterial, pcMaterial, pcMaterial};
// pitch definitions
plbConfig.centralModuleReadoutBinsX = {336, 336, 336, 336};
plbConfig.centralModuleReadoutBinsY = {1280, 1280, 1280, 1280};
plbConfig.centralModuleReadoutSide  = {-1, -1, -1, -1};
plbConfig.centralModuleLorentzAngle = {0.12, 0.12, 0.12, 0.12};

// no frontside/backside
plbConfig.centralModuleFrontsideStereo = {};
plbConfig.centralModuleBacksideStereo  = {};
plbConfig.centralModuleBacksideGap     = {};
// mPositions
std::vector<std::vector<Vector3D>> plbCentralModulePositions;
for (size_t plb = 0; plb < plbConfig.centralLayerRadii.size(); ++plb) {
  // call the helper function
  plbCentralModulePositions.push_back(
      modulePositionsCylinder(plbConfig.centralLayerRadii[plb],
                              0.5,  // 1 mm stagger
                              plbConfig.centralModuleHalfY[plb],
                              2.,  // 4 mm module overlap in z
                              plbConfig.centralModuleBinningSchema[plb]));
}
plbConfig.centralModulePositions = plbCentralModulePositions;

// ENDCAP :
// 4 pixel layers
// configure the endcaps
plbConfig.posnegLayerBinMultipliers = {1, 1};

plbConfig.posnegLayerPositionsZ = {600. * Acts::units::_mm,
                                   700. * Acts::units::_mm,
                                   820. * Acts::units::_mm,
                                   960. * Acts::units::_mm,
                                   1100 * Acts::units::_mm,
                                   1500 * Acts::units::_mm};

plbConfig.posnegLayerEnvelopeR = {1. * Acts::units::_mm,
                                  1. * Acts::units::_mm,
                                  1. * Acts::units::_mm,
                                  1. * Acts::units::_mm,
                                  1. * Acts::units::_mm,
                                  1. * Acts::units::_mm};
// material concentration is always behing
plbConfig.posnegLayerMaterialConcentration = std::vector<int>(6, 0);
// layer structure surface has pixel material properties
plbConfig.posnegLayerMaterialProperties = {pcmProperties,
                                           pcmProperties,
                                           pcmProperties,
                                           pcmProperties,
                                           pcmProperties,
                                           pcmProperties};
std::vector<double>         perHX = {8.4, 8.4};    // half length x
std::vector<double>         perHY = {36., 36.};    // half length y
std::vector<size_t>         perBP = {40, 68};      // bins in phi
std::vector<double>         perT  = {0.15, 0.15};  // module thickness
std::vector<size_t>         perBX = {336, 336};    // bins in x
std::vector<size_t>         perBY = {1280, 1280};  // bins in y
std::vector<int>            perRS = {-1, -1};      // readout side
std::vector<double>         perLA = {0., 0.};      // lorentz angle
std::vector<Acts::Material> perM  = {pcMaterial, pcMaterial};  // material

plbConfig.posnegModuleMinHalfX     = std::vector<std::vector<double>>(6, perHX);
plbConfig.posnegModuleMaxHalfX     = {};
plbConfig.posnegModuleHalfY        = std::vector<std::vector<double>>(6, perHY);
plbConfig.posnegModulePhiBins      = std::vector<std::vector<size_t>>(6, perBP);
plbConfig.posnegModuleThickness    = std::vector<std::vector<double>>(6, perT);
plbConfig.posnegModuleReadoutBinsX = std::vector<std::vector<size_t>>(6, perBX);
plbConfig.posnegModuleReadoutBinsY = std::vector<std::vector<size_t>>(6, perBY);
plbConfig.posnegModuleReadoutSide  = std::vector<std::vector<int>>(6, perRS);
plbConfig.posnegModuleLorentzAngle = std::vector<std::vector<double>>(6, perLA);
plbConfig.posnegModuleMaterial
    = std::vector<std::vector<Acts::Material>>(6, perM);

// no frontside/backside
plbConfig.posnegModuleFrontsideStereo = {};
plbConfig.posnegModuleBacksideStereo  = {};
plbConfig.posnegModuleBacksideGap     = {};
// mPositions
std::vector<std::vector<std::vector<Vector3D>>> plbPosnegModulePositions;
for (size_t id = 0; id < plbConfig.posnegLayerPositionsZ.size(); ++id) {
  plbPosnegModulePositions.push_back(
      modulePositionsDisc(plbConfig.posnegLayerPositionsZ[id],
                          0.0,
                          {4.0, 4.0},
                          {0.5, 0.},
                          30.,
                          176.,
                          plbConfig.posnegModulePhiBins[id],
                          plbConfig.posnegModuleHalfY[id]));
}
plbConfig.posnegModulePositions = plbPosnegModulePositions;

// define the builder
auto pixelLayerBuilder = std::make_shared<const GenericLayerBuilder>(
    plbConfig,
    getDefaultLogger("PixelLayerBuilder", layerLLevel));
//-------------------------------------------------------------------------------------
// build the pixel volume
CylinderVolumeBuilder::Config pvbConfig;
pvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
pvbConfig.volumeName           = "Pixel";
pvbConfig.buildToRadiusZero    = false;
pvbConfig.layerEnvelopeR       = {1. * Acts::units::_mm, 5. * Acts::units::_mm};
pvbConfig.layerBuilder         = pixelLayerBuilder;
pvbConfig.volumeSignature      = 0;
auto pixelVolumeBuilder        = std::make_shared<const CylinderVolumeBuilder>(
    pvbConfig,
    getDefaultLogger("PixelVolumeBuilder", volumeLLevel));
// add to the list of builders
volumeBuilders.push_back(pixelVolumeBuilder);

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
// Pixel Support Tybe (PST)
//-------------------------------------------------------------------------------------
PassiveLayerBuilder::Config pstConfig;
pstConfig.layerIdentification     = "PST";
pstConfig.centralLayerRadii       = std::vector<double>(1, 200.);
pstConfig.centralLayerHalflengthZ = std::vector<double>(1, 2800.);
pstConfig.centralLayerThickness   = std::vector<double>(1, 1.8);
pstConfig.centralLayerMaterial = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
auto pstBuilder                   = std::make_shared<const PassiveLayerBuilder>(
    pstConfig,
    getDefaultLogger("PSTBuilder", layerLLevel));
// create the volume for the beam pipe
CylinderVolumeBuilder::Config pstvolConfig;
pstvolConfig.trackingVolumeHelper = cylinderVolumeHelper;
pstvolConfig.volumeName           = "PST";
pstvolConfig.buildToRadiusZero    = false;
pstvolConfig.layerBuilder         = pstBuilder;
pstvolConfig.volumeSignature      = 0;
auto pstVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
    pstvolConfig,
    getDefaultLogger("PSTVolumeBuilder", volumeLLevel));
// add to the detector builds
volumeBuilders.push_back(pstVolumeBuilder);

//-------------------------------------------------------------------------------------
// SHORT strip detector
//-------------------------------------------------------------------------------------
// first add a Pixel Support Tube
// STRIPS
// ----------------------------------------------------------------------------
// configure short strip layer builder
GenericLayerBuilder::Config sslbConfig;
sslbConfig.layerCreator        = layerCreator;
sslbConfig.layerIdentification = "SStrip";
// fill necessary vectors for configuration
//-------------------------------------------------------------------------------------
// some prep work
// envelope double
std::pair<double, double> ssEnvelope(2., 2.);
// Layer material properties - thickness, X0, L0, A, Z, Rho
MaterialProperties ssmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
// Module material - X0, L0, A, Z, Rho
Material ssMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

// configure the central barrel
sslbConfig.centralLayerBinMultipliers = {1, 1};
sslbConfig.centralLayerRadii          = {260., 360., 500., 660.};
sslbConfig.centralLayerEnvelopes
    = {ssEnvelope, ssEnvelope, ssEnvelope, ssEnvelope};
sslbConfig.centralLayerMaterialConcentration = {-1, -1, -1, -1};
sslbConfig.centralLayerMaterialProperties
    = {ssmProperties, ssmProperties, ssmProperties, ssmProperties};
sslbConfig.centralModuleBinningSchema
    = {{40, 21}, {56, 21}, {78, 21}, {102, 21}};
sslbConfig.centralModuleTiltPhi   = {-0.15, -0.15, -0.15, -0.15};
sslbConfig.centralModuleHalfX     = {24., 24., 24., 24.};
sslbConfig.centralModuleHalfY     = {54., 54., 54., 54.};
sslbConfig.centralModuleThickness = {0.25, 0.25, 0.25, 0.25};

sslbConfig.centralModuleReadoutBinsX = {800, 800, 800, 800};  // 60 um pitch
sslbConfig.centralModuleReadoutBinsY = {9, 9, 9, 9};          // 1.2 cm strixels
sslbConfig.centralModuleReadoutSide  = {1, 1, 1, 1};
sslbConfig.centralModuleLorentzAngle = {0.12, 0.12, 0.12, 0.12};

sslbConfig.centralModuleMaterial
    = {ssMaterial, ssMaterial, ssMaterial, ssMaterial};
sslbConfig.centralModuleFrontsideStereo = {-0.02, -0.02, -0.02, -0.02};
sslbConfig.centralModuleBacksideStereo  = {0.02, 0.02, 0.02, 0.02};
sslbConfig.centralModuleBacksideGap     = {2., 2., 2., 2.};
// mPositions
std::vector<std::vector<Vector3D>> sslbCentralModulePositions;
for (size_t sslb = 0; sslb < sslbConfig.centralLayerRadii.size(); ++sslb) {
  // call the helper function
  sslbCentralModulePositions.push_back(
      modulePositionsCylinder(sslbConfig.centralLayerRadii[sslb],
                              3.,  // 3 mm stagger
                              sslbConfig.centralModuleHalfY[sslb],
                              5.,  // 5 mm module overlap
                              sslbConfig.centralModuleBinningSchema[sslb]));
}
sslbConfig.centralModulePositions = sslbCentralModulePositions;

// configure the endcaps
std::vector<double> mrMinHx = {16.4, 24.2, 32.2};
std::vector<double> mrMaxHx = {24.2, 32.2, 40.0};
std::vector<double> mrHy    = {78., 78., 78.};

// simplified strixels readout
std::vector<size_t> mrReadoutBinsX = {968, 1288, 1600};  // 50 um pitch
std::vector<size_t> mrReadoutBinsY = {60, 60, 60};       // 1.6 mm strixels
std::vector<int>    mrReadoutSide  = {1, 1, 1};
std::vector<double> mrLorentzAngle = {0., 0., 0.};

std::vector<size_t>   mPhiBins   = {52, 56, 60};
std::vector<double>   mThickness = {0.2, 0.2, 0.2};
std::vector<Material> mMaterial  = {ssMaterial, ssMaterial, ssMaterial};
std::vector<double>   mfStereo   = {-0.02, -0.02, -0.02};
std::vector<double>   mbStereo   = {0.02, 0.02, 0.02};
std::vector<double>   mfbGap     = {2., 2., 2.};

sslbConfig.posnegLayerBinMultipliers = {1, 2};

sslbConfig.posnegLayerPositionsZ = {1220., 1500., 1800., 2150., 2550., 2950.};
size_t nposnegs                  = sslbConfig.posnegLayerPositionsZ.size();
sslbConfig.posnegLayerEnvelopeR  = std::vector<double>(nposnegs, 5.);
sslbConfig.posnegLayerMaterialConcentration = std::vector<int>(nposnegs, 0);
sslbConfig.posnegLayerMaterialProperties
    = std::vector<MaterialProperties>(nposnegs, ssmProperties);
sslbConfig.posnegModuleMinHalfX
    = std::vector<std::vector<double>>(nposnegs, mrMinHx);
sslbConfig.posnegModuleMaxHalfX
    = std::vector<std::vector<double>>(nposnegs, mrMaxHx);
sslbConfig.posnegModuleHalfY = std::vector<std::vector<double>>(nposnegs, mrHy);
sslbConfig.posnegModulePhiBins
    = std::vector<std::vector<size_t>>(nposnegs, mPhiBins);
sslbConfig.posnegModuleThickness
    = std::vector<std::vector<double>>(nposnegs, mThickness);

sslbConfig.posnegModuleReadoutBinsX
    = std::vector<std::vector<size_t>>(nposnegs, mrReadoutBinsX);
sslbConfig.posnegModuleReadoutBinsY
    = std::vector<std::vector<size_t>>(nposnegs, mrReadoutBinsY);
sslbConfig.posnegModuleReadoutSide
    = std::vector<std::vector<int>>(nposnegs, mrReadoutSide);
sslbConfig.posnegModuleLorentzAngle
    = std::vector<std::vector<double>>(nposnegs, mrLorentzAngle);

sslbConfig.posnegModuleMaterial
    = std::vector<std::vector<Material>>(nposnegs, mMaterial);
sslbConfig.posnegModuleFrontsideStereo
    = std::vector<std::vector<double>>(nposnegs, mfStereo);
sslbConfig.posnegModuleBacksideStereo
    = std::vector<std::vector<double>>(nposnegs, mbStereo);
sslbConfig.posnegModuleBacksideGap
    = std::vector<std::vector<double>>(nposnegs, mfbGap);
// mPositions
std::vector<std::vector<std::vector<Vector3D>>> sslbPosnegModulePositions;
for (size_t id = 0; id < sslbConfig.posnegLayerPositionsZ.size(); ++id) {
  sslbPosnegModulePositions.push_back(
      modulePositionsDisc(sslbConfig.posnegLayerPositionsZ[id],
                          6.0,
                          {3., 3., 3.},
                          {0., 0., 0.},
                          240.,
                          700.,
                          sslbConfig.posnegModulePhiBins[id],
                          sslbConfig.posnegModuleHalfY[id]));
}
sslbConfig.posnegModulePositions = sslbPosnegModulePositions;

// define the builder
auto sstripLayerBuilder = std::make_shared<const GenericLayerBuilder>(
    sslbConfig,
    getDefaultLogger("SStripLayerBuilder", layerLLevel));
//-------------------------------------------------------------------------------------
// build the pixel volume
CylinderVolumeBuilder::Config ssvbConfig;
ssvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
ssvbConfig.volumeName           = "SStrip";
ssvbConfig.buildToRadiusZero    = false;
ssvbConfig.layerBuilder         = sstripLayerBuilder;
ssvbConfig.volumeSignature      = 0;
auto sstripVolumeBuilder        = std::make_shared<const CylinderVolumeBuilder>(
    ssvbConfig,
    getDefaultLogger("SStripVolumeBuilder", volumeLLevel));

//-------------------------------------------------------------------------------------
// add to the list of builders
volumeBuilders.push_back(sstripVolumeBuilder);

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
// LONG strip detector
//-------------------------------------------------------------------------------------

// configure short strip layer builder
GenericLayerBuilder::Config lslbConfig;
lslbConfig.layerCreator        = layerCreator;
lslbConfig.layerIdentification = "LStrip";
// fill necessary vectors for configuration
//-------------------------------------------------------------------------------------
// some prep work
// envelope double
std::pair<double, double> lsEnvelope(2., 2.);
// Layer material properties - thickness, X0, L0, A, Z, Rho
MaterialProperties lsmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
// Module material - X0, L0, A, Z, Rho
Material lsMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

// configure the central barrel
lslbConfig.centralLayerBinMultipliers        = {1, 1};
lslbConfig.centralLayerRadii                 = {820., 1020.};
lslbConfig.centralLayerEnvelopes             = {lsEnvelope, lsEnvelope};
lslbConfig.centralLayerMaterialConcentration = {-1, -1};
lslbConfig.centralLayerMaterialProperties    = {lsmProperties, lsmProperties};
lslbConfig.centralModuleBinningSchema        = {{120, 21}, {152, 21}};
lslbConfig.centralModuleTiltPhi              = {-0.15, -0.15};
lslbConfig.centralModuleHalfX                = {24., 24.};
lslbConfig.centralModuleHalfY                = {54., 54.};
lslbConfig.centralModuleThickness            = {0.25, 0.25};
lslbConfig.centralModuleMaterial             = {lsMaterial, lsMaterial};

lslbConfig.centralModuleReadoutBinsX = {840, 840};  // 100 um pitch
lslbConfig.centralModuleReadoutBinsY = {1, 1};      // 1 strip
lslbConfig.centralModuleReadoutSide  = {1, 1};
lslbConfig.centralModuleLorentzAngle = {0.08, 0.08};

lslbConfig.centralModuleFrontsideStereo = {-0.04, -0.04};
lslbConfig.centralModuleBacksideStereo  = {0.04, 0.04};
lslbConfig.centralModuleBacksideGap     = {2., 2.};
// mPositions
std::vector<std::vector<Vector3D>> lslbCentralModulePositions;
for (size_t lslb = 0; lslb < lslbConfig.centralLayerRadii.size(); ++lslb) {
  // call the helper function
  lslbCentralModulePositions.push_back(
      modulePositionsCylinder(lslbConfig.centralLayerRadii[lslb],
                              3.,  // 3 mm stagger
                              lslbConfig.centralModuleHalfY[lslb],
                              5.,  // 5 mm module overlap
                              lslbConfig.centralModuleBinningSchema[lslb]));
}

lslbConfig.centralModulePositions = lslbCentralModulePositions;

// configure the endcaps
mrMinHx    = {40., 78.};
mrMaxHx    = {64.};
mrHy       = {78., 78.};
mPhiBins   = {58, 46};
mThickness = {0.25, 0.25};
mMaterial  = {lsMaterial, lsMaterial};
mfStereo   = {-0.04, -0.04};
mbStereo   = {0.04, 0.04, 0.04};
mfbGap     = {2., 2.};

mrReadoutBinsX = {1120, 1120};  // 100 um pitch
mrReadoutBinsY = {1, 1};        // 2 strips
mrReadoutSide  = {1, 1};
mrLorentzAngle = {0., 0.};

// endcap
lslbConfig.posnegLayerBinMultipliers = {1, 2};
lslbConfig.posnegLayerPositionsZ = {1220., 1500., 1800., 2150., 2550., 2950.};
nposnegs                         = lslbConfig.posnegLayerPositionsZ.size();
lslbConfig.posnegLayerEnvelopeR  = std::vector<double>(nposnegs, 5.);
lslbConfig.posnegLayerMaterialConcentration = std::vector<int>(nposnegs, 0);
lslbConfig.posnegLayerMaterialProperties
    = std::vector<MaterialProperties>(nposnegs, lsmProperties);
lslbConfig.posnegModuleMinHalfX
    = std::vector<std::vector<double>>(nposnegs, mrMinHx);
lslbConfig.posnegModuleMaxHalfX
    = std::vector<std::vector<double>>(nposnegs, mrMaxHx);
lslbConfig.posnegModuleHalfY = std::vector<std::vector<double>>(nposnegs, mrHy);
lslbConfig.posnegModulePhiBins
    = std::vector<std::vector<size_t>>(nposnegs, mPhiBins);
lslbConfig.posnegModuleThickness
    = std::vector<std::vector<double>>(nposnegs, mThickness);

lslbConfig.posnegModuleReadoutBinsX
    = std::vector<std::vector<size_t>>(nposnegs, mrReadoutBinsX);
lslbConfig.posnegModuleReadoutBinsY
    = std::vector<std::vector<size_t>>(nposnegs, mrReadoutBinsY);
lslbConfig.posnegModuleReadoutSide
    = std::vector<std::vector<int>>(nposnegs, mrReadoutSide);
lslbConfig.posnegModuleLorentzAngle
    = std::vector<std::vector<double>>(nposnegs, mrLorentzAngle);

lslbConfig.posnegModuleMaterial
    = std::vector<std::vector<Material>>(nposnegs, mMaterial);
lslbConfig.posnegModuleFrontsideStereo
    = std::vector<std::vector<double>>(nposnegs, mfStereo);
lslbConfig.posnegModuleBacksideStereo
    = std::vector<std::vector<double>>(nposnegs, mbStereo);
lslbConfig.posnegModuleBacksideGap
    = std::vector<std::vector<double>>(nposnegs, mfbGap);
// mPositions
std::vector<std::vector<std::vector<Vector3D>>> lssbPosnegModulePositions;
for (size_t id = 0; id < lslbConfig.posnegLayerPositionsZ.size(); ++id) {
  lssbPosnegModulePositions.push_back(modulePositionsDisc(
      lslbConfig.posnegLayerPositionsZ[id],
      8.0,  // staggering of rings, we put the disk strucutre in between
      {3., 3.},
      {0., 0.},
      750.,
      1020.,
      lslbConfig.posnegModulePhiBins[id],
      lslbConfig.posnegModuleHalfY[id]));
}
lslbConfig.posnegModulePositions = lssbPosnegModulePositions;

// define the builder
auto lstripLayerBuilder = std::make_shared<const GenericLayerBuilder>(
    lslbConfig,
    getDefaultLogger("LStripLayerBuilder", layerLLevel));
//-------------------------------------------------------------------------------------
// build the pixel volume
CylinderVolumeBuilder::Config lsvbConfig;
lsvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
lsvbConfig.volumeName           = "LStrip";
lsvbConfig.buildToRadiusZero    = false;
lsvbConfig.layerBuilder         = lstripLayerBuilder;
lsvbConfig.volumeSignature      = 0;
auto lstripVolumeBuilder        = std::make_shared<const CylinderVolumeBuilder>(
    lsvbConfig,
    getDefaultLogger("LStripVolumeBuilder", volumeLLevel));
// add to the list of builders
volumeBuilders.push_back(lstripVolumeBuilder);
