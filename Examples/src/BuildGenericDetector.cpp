// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Examples/BuildGenericDetector.hpp"
#include <array>
#include <iostream>
#include <vector>
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Examples/GenericLayerBuilder.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/CylinderGeometryBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/LayerArrayCreator.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/PassiveLayerBuilder.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Tools/TrackingVolumeArrayCreator.hpp"

namespace Acts {

std::unique_ptr<const Acts::TrackingGeometry>
trackingGeometry(Logging::Level lvl, size_t stage)
{
  // configure surface array creator
  SurfaceArrayCreator::Config sacConfig;
  sacConfig.logger         = getDefaultLogger("SurfaceArrayCreator", lvl);
  auto surfaceArrayCreator = std::make_shared<SurfaceArrayCreator>(sacConfig);
  // configure the layer creator that uses the surface array creator
  LayerCreator::Config lcConfig;
  lcConfig.logger              = getDefaultLogger("LayerCreator", lvl);
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator            = std::make_shared<LayerCreator>(lcConfig);
  // configure the layer array creator
  LayerArrayCreator::Config lacConfig;
  lacConfig.logger       = getDefaultLogger("LayerArrayCreator", lvl);
  auto layerArrayCreator = std::make_shared<LayerArrayCreator>(lacConfig);
  // tracking volume array creator
  TrackingVolumeArrayCreator::Config tvacConfig;
  tvacConfig.logger = getDefaultLogger("TrackingVolumeArrayCreator", lvl);
  auto tVolumeArrayCreator
      = std::make_shared<TrackingVolumeArrayCreator>(tvacConfig);
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator          = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper            = std::make_shared<CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", lvl));
  //-------------------------------------------------------------------------------------
  // beam pipe
  //-------------------------------------------------------------------------------------
  // configure the beam pipe layer builder
  PassiveLayerBuilder::Config bplConfig;
  bplConfig.logger              = getDefaultLogger("PassiveLayerBuilder", lvl);
  bplConfig.layerIdentification = "BeamPipe";
  bplConfig.centralLayerRadii   = std::vector<double>(1, 19.);
  bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 200.);
  bplConfig.centralLayerThickness   = std::vector<double>(1, 0.8);
  bplConfig.centralLayerMaterial = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
  auto beamPipeBuilder = std::make_shared<PassiveLayerBuilder>(bplConfig);
  // create the volume for the beam pipe
  CylinderVolumeBuilder::Config bpvConfig;
  bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
  bpvConfig.volumeName           = "BeamPipe";
  bpvConfig.layerBuilder         = beamPipeBuilder;
  bpvConfig.layerEnvelopeR       = 1.;
  bpvConfig.layerEnvelopeZ       = 1.;
  bpvConfig.volumeSignature      = 0;
  auto beamPipeVolumeBuilder     = std::make_shared<CylinderVolumeBuilder>(
      bpvConfig, getDefaultLogger("CylinderVolumeBuilder", lvl));
  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // pixel detector
  //-------------------------------------------------------------------------------------
  // configure pixel layer builder
  GenericLayerBuilder::Config plbConfig;
  plbConfig.layerCreator        = layerCreator;
  plbConfig.layerIdentification = "Pixel";
  // fill necessary vectors for configuration
  //-------------------------------------------------------------------------------------
  // some prep work
  // envelope double
  std::pair<double, double> pcEnvelope(2., 2.);
  // Layer material properties - thickness, X0, L0, A, Z, Rho
  MaterialProperties pcmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
  // Module material - X0, L0, A, Z, Rho
  Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
  // configure the central barrel
  plbConfig.centralLayerBinMultipliers = {2, 2};
  plbConfig.centralLayerRadii          = {29., 55., 88., 120.};
  plbConfig.centralLayerEnvelopes
      = {pcEnvelope, pcEnvelope, pcEnvelope, pcEnvelope};
  plbConfig.centralLayerMaterialConcentration = {1, 1, 1, 1};
  plbConfig.centralLayerMaterialProperties
      = {pcmProperties, pcmProperties, pcmProperties, pcmProperties};
  plbConfig.centralModuleBinningSchema
      = {{16, 13}, {24, 13}, {38, 13}, {60, 13}};
  plbConfig.centralModuleTiltPhi   = {0.18, 0.18, 0.2, 0.2};
  plbConfig.centralModuleHalfX     = {8.4, 8.4, 8.4, 8.4};
  plbConfig.centralModuleHalfY     = {32., 32., 32., 32.};
  plbConfig.centralModuleThickness = {0.15, 0.15, 0.15, 0.15};
  plbConfig.centralModuleMaterial
      = {pcMaterial, pcMaterial, pcMaterial, pcMaterial};
  plbConfig.centralModuleFrontsideStereo      = {};
  plbConfig.centralModuleBacksideStereo       = {};
  plbConfig.centralModuleBacksideGap          = {};
  // mPositions
  std::vector<std::vector<Vector3D>> centralModulePositions;
  for (size_t plb = 0; plb < plbConfig.centralLayerRadii.size(); ++plb) {
    // call the helper function
    centralModulePositions.push_back(
        modulePositionsCylinder(plbConfig.centralLayerRadii[plb],
                                0.5,  // 1 mm stagger
                                plbConfig.centralModuleHalfY[plb],
                                2.,  // 2 mm module overlap
                                plbConfig.centralModuleBinningSchema[plb]));
  }
  plbConfig.centralModulePositions            = centralModulePositions;

  // configure the endcaps
  plbConfig.posnegLayerBinMultipliers          = { 1, 1 };
  plbConfig.posnegLayerPositionsZ              = {500., 580., 680., 700.};
  plbConfig.posnegLayerEnvelopeR               = {1., 1., 1., 1.};
  plbConfig.posnegLayerMaterialConcentration   = {1, 1, 1, 1};
  plbConfig.posnegLayerMaterialProperties
      = {pcmProperties, pcmProperties, pcmProperties, pcmProperties};
  plbConfig.posnegModuleMinHalfX               = {{8.4}, {8.4}, {8.4}, {8.4}};
  plbConfig.posnegModuleMaxHalfX               = {};
  plbConfig.posnegModuleHalfY                  = {{42.}, {42.}, {42.}, {42.}};
  plbConfig.posnegModulePhiBins                = {{24}, {24}, {24}, {24}};
  plbConfig.posnegModuleThickness = {{0.15}, {0.15}, {0.15}, {0.15}};
  plbConfig.posnegModuleMaterial
      = {{pcMaterial}, {pcMaterial}, {pcMaterial}, {pcMaterial}};
  plbConfig.posnegModuleFrontsideStereo        = {};
  plbConfig.posnegModuleBacksideStereo         = {};
  plbConfig.posnegModuleBacksideGap            = {};
  // mPositions  
  std::vector< std::vector< std::vector<Vector3D> > > posnegModulePositions;
  for (size_t id = 0; id < plbConfig.posnegLayerPositionsZ.size(); ++id ){
    posnegModulePositions.push_back(
        modulePositionsDisc(plbConfig.posnegLayerPositionsZ[id],
                            2.0,
                            0.5,
                            29.,
                            120.,
                            plbConfig.posnegModulePhiBins[id],
                            plbConfig.posnegModuleHalfY[id]));
  }
  plbConfig.posnegModulePositions = posnegModulePositions;
  // define the builder
  auto pixelLayerBuilder = std::make_shared<GenericLayerBuilder>(
      plbConfig, getDefaultLogger("GenericLayerBuilder", lvl));
  //-------------------------------------------------------------------------------------
  // build the pixel volume
  CylinderVolumeBuilder::Config pvbConfig;
  pvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  pvbConfig.volumeName           = "Pixel";
  pvbConfig.volumeToBeamPipe     = false;
  pvbConfig.layerBuilder         = pixelLayerBuilder;
  pvbConfig.layerEnvelopeR       = 1.;
  pvbConfig.layerEnvelopeZ       = 10.;
  pvbConfig.volumeSignature      = 0;
  auto pixelVolumeBuilder        = std::make_shared<CylinderVolumeBuilder>(
      pvbConfig, getDefaultLogger("CylinderVolumeBuilder", lvl));

  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<ITrackingVolumeBuilder>> detectorBuilders;
  detectorBuilders.push_back(pixelVolumeBuilder);

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // strip detector
  //-------------------------------------------------------------------------------------
  if (stage > 0) {
    // configure short strip layer builder
    GenericLayerBuilder::Config sslbConfig;
    sslbConfig.logger       = getDefaultLogger("SStripLayerBuilder", lvl);
    sslbConfig.layerCreator = layerCreator;
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
    sslbConfig.centralLayerRadii          = {220., 320., 420.};
    sslbConfig.centralLayerEnvelopes = {ssEnvelope, ssEnvelope, ssEnvelope};
    sslbConfig.centralLayerMaterialConcentration = {1, 1, 1};
    sslbConfig.centralLayerMaterialProperties
        = {ssmProperties, ssmProperties, ssmProperties};
    sslbConfig.centralModuleBinningSchema = {{42, 12}, {64, 12}, {78, 12}};
    sslbConfig.centralModuleTiltPhi       = {-0.15, -0.15, -0.15};
    sslbConfig.centralModuleHalfX         = {18.2, 18.2, 18.2};
    sslbConfig.centralModuleHalfY         = {68., 68., 68.};
    sslbConfig.centralModuleThickness     = {0.25, 0.25, 0.25};
    sslbConfig.centralModuleMaterial
        = {ssMaterial, ssMaterial, ssMaterial, ssMaterial};
    sslbConfig.centralModuleFrontsideStereo = {-0.2, 0.2, -0.2};
    sslbConfig.centralModuleBacksideStereo  = {0.2, -0.2, 0.2};
    sslbConfig.centralModuleBacksideGap     = {2., 2., 2.};
    // mPositions
    std::vector<std::vector<Vector3D>> centralModulePositions;
    for (size_t sslb = 0; sslb < sslbConfig.centralLayerRadii.size(); ++sslb) {
      // call the helper function
      centralModulePositions.push_back(
          modulePositionsCylinder(sslbConfig.centralLayerRadii[sslb],
                                  0.5,  // 1 mm stagger
                                  sslbConfig.centralModuleHalfY[sslb],
                                  2.,  // 2 mm module overlap
                                  sslbConfig.centralModuleBinningSchema[sslb]));
    }
    sslbConfig.centralModulePositions = centralModulePositions;

    // configure the endcaps
    sslbConfig.posnegLayerBinMultipliers        = {1, 1};
    sslbConfig.posnegLayerPositionsZ            = {880.};
    sslbConfig.posnegLayerEnvelopeR             = {5.};
    sslbConfig.posnegLayerMaterialConcentration = {1};
    sslbConfig.posnegLayerMaterialProperties    = {ssmProperties};
    sslbConfig.posnegModuleMinHalfX             = {{16.4, 18.2}};
    sslbConfig.posnegModuleMaxHalfX             = {{24.2, 32.2}};
    sslbConfig.posnegModuleHalfY                = {{52.5, 52.5}};
    sslbConfig.posnegModulePhiBins              = {{42, 42}};
    sslbConfig.posnegModuleThickness            = {{0.25, 0.25}};
    sslbConfig.posnegModuleMaterial             = {{ssMaterial, ssMaterial}};
    sslbConfig.posnegModuleFrontsideStereo      = {};
    sslbConfig.posnegModuleBacksideStereo       = {};
    sslbConfig.posnegModuleBacksideGap          = {};
    // mPositions
    std::vector<std::vector<std::vector<Vector3D>>> posnegModulePositions;
    for (size_t id = 0; id < sslbConfig.posnegLayerPositionsZ.size(); ++id) {
      posnegModulePositions.push_back(
          modulePositionsDisc(sslbConfig.posnegLayerPositionsZ[id],
                              2.0,
                              0.5,
                              220.,
                              420.,
                              sslbConfig.posnegModulePhiBins[id],
                              sslbConfig.posnegModuleHalfY[id]));
    }
    sslbConfig.posnegModulePositions = posnegModulePositions;

    // define the builder
    auto sstripLayerBuilder = std::make_shared<GenericLayerBuilder>(sslbConfig);
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    CylinderVolumeBuilder::Config ssvbConfig;
    ssvbConfig.logger = Acts::getDefaultLogger("SStripVolumeBuilder", lvl);
    ssvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    ssvbConfig.volumeName           = "SStrip";
    ssvbConfig.volumeToBeamPipe     = false;
    ssvbConfig.layerBuilder         = sstripLayerBuilder;
    ssvbConfig.layerEnvelopeR       = 1.;
    ssvbConfig.layerEnvelopeZ       = 10.;
    ssvbConfig.volumeSignature      = 0;
    auto sstripVolumeBuilder
        = std::make_shared<CylinderVolumeBuilder>(ssvbConfig);

    //-------------------------------------------------------------------------------------
    // list the volume builders
    detectorBuilders.push_back(sstripVolumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  CylinderGeometryBuilder::Config tgConfig;
  tgConfig.beamPipeBuilder        = beamPipeVolumeBuilder;
  tgConfig.trackingVolumeBuilders = detectorBuilders;
  tgConfig.trackingVolumeHelper   = cylinderVolumeHelper;
  auto cylinderGeometryBuilder
      = std::make_shared<const CylinderGeometryBuilder>(
          tgConfig, getDefaultLogger("CylinderGeometryBuilder", lvl));
  return cylinderGeometryBuilder->trackingGeometry();
}

/// helper method for cylinder
std::vector<Acts::Vector3D>
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
  double zStep   = 2 * fabs(zStart) / (nZbins - 1);
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

/// helper method for disc
std::vector<std::vector<Acts::Vector3D>>
modulePositionsDisc(double                     z,
                    double                     ringStagger,
                    double                     phiStagger,
                    double                     innerRadius,
                    double                     outerRadius,
                    const std::vector<int>&    discBinning,
                    const std::vector<double>& moduleHalfLength)
{
  // calculate the radii
  std::vector<double> radii;
  // calculate the radial borders
  std::vector<double> radialBoarders;
  // the radial span of the disc
  double deltaR = outerRadius - innerRadius;
  // quick exits
  if (discBinning.size() == 1) {
    radii.push_back(0.5 * (innerRadius + outerRadius));
    radialBoarders = {innerRadius, outerRadius};
  } else {
    double totalLength = 0;
    // sum up the total length
    for (auto& mhlength : moduleHalfLength) totalLength += 2 * mhlength;
    // now calculate the overlap (equal pay)
    double rOverlap = (totalLength - deltaR) / (moduleHalfLength.size() - 1);
    // and now fill the radii and gaps
    double lastR  = innerRadius;
    double lastHl = 0.;
    double lastOl = 0.;
    // remember the radial boarders
    radialBoarders.push_back(innerRadius);
    // now calculate
    for (auto& mhlength : moduleHalfLength) {
      // calculate the radius
      radii.push_back(lastR + lastHl - lastOl + mhlength);
      lastR  = radii[radii.size() - 1];
      lastOl = rOverlap;
      lastHl = mhlength;
      // and register the radial boarder
      radialBoarders.push_back(lastR + 2 * lastHl - 0.5 * lastOl);
    }
  }
  // now prepare the return method
  std::vector<std::vector<Vector3D>> mPositions;
  for (size_t ir = 0; ir < radii.size(); ++ir) {
    // generate the z value
    double rz = radii.size() == 1 ? z : (ir % 2 ? z - 0.5 * ringStagger
                                                : z + 0.5 * ringStagger);
    // fill the ring positions
    mPositions.push_back(
        modulePositionsRing(rz, radii[ir], phiStagger, discBinning[ir]));
  }
  return mPositions;
}

/// Helper method for positioning
std::vector<Acts::Vector3D>
modulePositionsRing(double z, double radius, double phiStagger, int nPhiBins)
{
  // create and fill the positions
  std::vector<Vector3D> rPositions;
  rPositions.reserve(nPhiBins);
  // prep work
  double phiStep = 2 * M_PI / (nPhiBins);
  double minPhi  = -M_PI + 0.5 * phiStep;
  // phi loop
  for (size_t iphi = 0; iphi < size_t(nPhiBins); ++iphi) {
    double phi = minPhi + iphi * phiStep;
    double rz  = iphi % 2 ? z - 0.5 * phiStagger : z + 0.5 * phiStagger;
    rPositions.push_back(Vector3D(radius * cos(phi), radius * sin(phi), rz));
  }
  return rPositions;
}

}  // end of namespace Acts
