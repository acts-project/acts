// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Examples/BuildGenericDetector.hpp"
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
#include <array>
#include <iostream>
#include <vector>

namespace Acts {
  
std::unique_ptr<const Acts::TrackingGeometry>
trackingGeometry(Logging::Level lvl, size_t version)
{
  // configure surface array creator
  SurfaceArrayCreator::Config sacConfig;
  sacConfig.logger             = getDefaultLogger("SurfaceArrayCreator", lvl);
    auto surfaceArrayCreator = std::make_shared<SurfaceArrayCreator>(sacConfig);
  // configure the layer creator that uses the surface array creator
  LayerCreator::Config lcConfig;
  lcConfig.logger              = getDefaultLogger("LayerCreator", lvl);
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<LayerCreator>(lcConfig);
  // configure the layer array creator 
  LayerArrayCreator::Config lacConfig;
  lacConfig.logger            = getDefaultLogger("LayerArrayCreator", lvl);  
  auto layerArrayCreator = std::make_shared<LayerArrayCreator>(lacConfig);
  // tracking volume array creator
  TrackingVolumeArrayCreator::Config tvacConfig;
  tvacConfig.logger           = getDefaultLogger("TrackingVolumeArrayCreator", lvl);  
  auto tVolumeArrayCreator = std::make_shared<TrackingVolumeArrayCreator>(tvacConfig);
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.logger            = getDefaultLogger("CylinderVolumeHelper", lvl);
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper
      = std::make_shared<CylinderVolumeHelper>(cvhConfig);
  //-------------------------------------------------------------------------------------
  // beam pipe
  //-------------------------------------------------------------------------------------
  // configure the beam pipe layer builder
  PassiveLayerBuilder::Config bplConfig;
  bplConfig.logger = getDefaultLogger("PassiveLayerBuilder", lvl);
  bplConfig.layerIdentification     = "BeamPipe";
  bplConfig.centralLayerRadii       = std::vector<double>(1, 21.);
  bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 200.);
  bplConfig.centralLayerThickness   = std::vector<double>(1, 0.8);
    bplConfig.centralLayerMaterial    = { Material(352.8,407.,9.012, 4., 1.848e-3) };
  auto beamPipeBuilder = std::make_shared<PassiveLayerBuilder>(bplConfig);
  // create the volume for the beam pipe
  CylinderVolumeBuilder::Config bpvConfig;
  bpvConfig.logger = getDefaultLogger("CylinderVolumeBuilder", lvl);
  bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
  bpvConfig.volumeName           = "BeamPipe";
  bpvConfig.layerBuilder         = beamPipeBuilder;
  bpvConfig.layerEnvelopeR       = 1.;
  bpvConfig.layerEnvelopeZ       = 1.;
  bpvConfig.volumeSignature      = 0;
  auto beamPipeVolumeBuilder
      = std::make_shared<CylinderVolumeBuilder>(bpvConfig);
  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // pixel detector
  //-------------------------------------------------------------------------------------
  // configure pixel layer builder
  GenericLayerBuilder::Config plbConfig;
  plbConfig.logger = getDefaultLogger("GenericLayerBuilder", lvl);
  plbConfig.layerCreator        = layerCreator;
  plbConfig.layerIdentification = "Pixel";
  // fill necessary vectors for configuration
  //-------------------------------------------------------------------------------------
  // some prep work
  // envelope double
  std::pair<double, double> pcEnvelope(2.,2.); 
  // Layer material properties - thickness, X0, L0, A, Z, Rho
  MaterialProperties pcmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);  
  // Module material - X0, L0, A, Z, Rho
  Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);
  // configure the central barrel
  plbConfig.centralLayerBinMultipliers        = { 1, 1 };
  plbConfig.centralLayerRadii                 = { 29., 55., 88. };
  plbConfig.centralLayerEnvelopes             = { pcEnvelope, pcEnvelope, pcEnvelope };
  plbConfig.centralLayerMaterialConcentration = { 1, 1, 1 };
  plbConfig.centralLayerMaterialProperties    = { pcmProperties, pcmProperties, pcmProperties };
  plbConfig.centralModuleBinningSchema        = { {15,13} , {24,3}, {40,13} };
  plbConfig.centralModuleTiltPhi              = { 0.18, 0.18, 0.2 };
  plbConfig.centralModuleHalfX                = { 8.4, 8.4, 8.4 }; 
  plbConfig.centralModuleHalfY                = { 32., 32., 32.};
  plbConfig.centralModuleThickness            = { 0.15 , 0.15 , 0.15 };
  plbConfig.centralModuleMaterial             = { pcMaterial, pcMaterial, pcMaterial };
  plbConfig.centralModuleFrontsideStereo      = {};
  plbConfig.centralModuleBacksideStereo       = {};
  plbConfig.centralModuleBacksideGap          = {};
  // mPositions
  std::vector< std::vector< Vector3D> > centralModulePositions;
  for (size_t plb = 0; plb < plbConfig.centralLayerRadii.size(); ++plb){
    // call the helper function
    centralModulePositions.push_back(modulePositionsCylinder(plbConfig.centralLayerRadii[plb],
                                                             1., // 1 mm stagger
                                                             2*plbConfig.centralModuleHalfY[plb],
                                                             2., // 2 mm module overlap
                                                             plbConfig.centralModuleBinningSchema[plb]));
    
  }
            
  plbConfig.centralModulePositions            = centralModulePositions;

///
///
///
///
///
///
///  // pixel radii
///  std::vector<double>   pcLayerRadii = {29., 55., 88.};
///  // the envelope
///  std::vector<double>   pcEnvelopes(pcLayers, std::make_pair<double,double>(2.,2.));
///  // the material concentration : outside
///  std::vector<double>   pcMaterialConcentration(pcLayers, 1.);
///  // thickness,X0,L0,A,Z,Rho
///  std::vector<size_t>   pcNumPhi = {15, 24, 40};
///  std::vector<double>   pcLayerModulesTiltPhi = {};
///  
///  std::vector<std::pair<int,int>>       pcNumZ     = {13, 13, 13};
///  
///  
///  std::vector<std::vector<double>> pcLayerModulesPositionPhi;
///  for (int i = 0; i < pcLayers; i++) {
///    std::vector<double> positionsPhi;
///    double              phiStep = 2. * M_PI / pcNumPhi[i];
///    for (int j = 0; j < pcNumPhi[i]; j++) {
///      positionsPhi.push_back(-M_PI + j * phiStep);
///    }
///    pcLayerModulesPositionPhi.push_back(positionsPhi);
///  }
///  // z properties centrallayers
///  double              pcNumZ[pcLayers]     = {13, 13, 13};
///  double              pcOverlapZ[pcLayers] = {2., 2., 2.};
///  std::vector<double> pcLayerModulesHalfY  = {32., 32., 32.};
///
///  std::vector<std::vector<double>> pcLayerModulesPositionZ;
///  for (int i = 0; i < pcLayers; i++) {
///    std::vector<double> positionsZ;
///    double              firstToLast
///        = (pcNumZ[i] - 1) * (2. * pcLayerModulesHalfY.at(i) - pcOverlapZ[i]);
///    double zStep = firstToLast / (pcNumZ[i] - 1);
///    for (int j = 0; j < pcNumZ[i]; j++) {
///      positionsZ.push_back(-0.5 * firstToLast + j * zStep);
///    }
///    pcLayerModulesPositionZ.push_back(positionsZ);
///  }
///  std::vector<double> pcLayerModulesStaggerZ(pcLayers, 0.5);
///  std::vector<double> pcLayerModulesMinHalfX;
///  std::vector<double> pcLayerModulesMaxHalfX(pcLayers, 8.4);
///  std::vector<double> pcLayerModulesThickness(pcLayers, 0.15);
///  std::vector<double> pcModuleMaterialParameters
///      = {95.7, 465.2, 28.03, 14., 2.32e-3};  // X0,L0,A,Z,Rho
///  std::vector<std::vector<double>> pcLayerModulesMaterial(
///      pcLayers, pcModuleMaterialParameters);
///  //-------------------------------------------------------------------------------------
///  // posneg
///  const int           ppnLayers       = 3;
///  const int           ppnRings        = 1;
///  std::vector<double> pixelPositionsZ = {500., 580., 650.};
///
///  std::vector<double>              pixelEnvelopesR(ppnLayers, 5.);
///  std::vector<double>              ppnRadii(ppnRings, 65.);
///  std::vector<std::vector<double>> ppnLayerModulesRadii(ppnLayers, ppnRadii);
///  std::vector<double>              ppnLayerModulesStaggerR(ppnLayers, 3.);
///  // phi properties posneglayers
///  std::vector<double>              ppnInPhi = {24., 24., 24.};
///  std::vector<std::vector<double>> ppnLayerModulesInPhi(ppnLayers, ppnInPhi);
///  std::vector<std::vector<double>> ppnPositionPhi;
///  for (int i = 0; i < ppnRings; i++) {
///    std::vector<double> positionsPhi;
///    double              phiStep = 2. * M_PI / ppnInPhi.at(i);
///    for (int j = 0; j < ppnInPhi.at(i); j++) {
///      positionsPhi.push_back(-M_PI + j * phiStep);
///    }
///    ppnPositionPhi.push_back(positionsPhi);
///  }
///  std::vector<std::vector<std::vector<double>>> ppnLayerModulesPositionPhi(
///      ppnLayers, ppnPositionPhi);
///  std::vector<double>              ppnStaggerPhi(ppnRings, 0.5);
///  std::vector<std::vector<double>> ppnLayerModulesStaggerPhi(ppnLayers,
///                                                             ppnStaggerPhi);
///  // Module properties
///  std::vector<std::vector<double>> ppnLayerModulesMinHalfX(ppnLayers);
///  std::vector<double>              ppnMaxHalfX(ppnRings, 8.4);
///  std::vector<std::vector<double>> ppnLayerModulesMaxHalfX(ppnLayers,
///                                                           ppnMaxHalfX);
///  std::vector<double>              ppnMaxHalfY(ppnRings, 32.);
///  std::vector<std::vector<double>> ppnLayerModulesHalfY(ppnLayers, ppnMaxHalfY);
///  std::vector<double>              ppnThickness(ppnRings, 0.15);
///  std::vector<std::vector<double>> ppnLayerModulesThickness(pcLayers,
///                                                            ppnThickness);
///  std::vector<double> ppnModuleMaterialParameters
///      = {95.7, 465.2, 28.03, 14., 2.32e-3};  // X0,L0,A,Z,Rho
///  std::vector<std::vector<double>> ppnMaterial(ppnLayers,
///                                               ppnModuleMaterialParameters);
///  std::vector<std::vector<std::vector<double>>> ppnLayerModulesMaterial(
///      ppnLayers, ppnMaterial);
///
///  //-------------------------------------------------------------------------------------
///
///  // configure the central barrel
///  plbConfig.centralLayerRadii                 = pixelLayerRadii;
///  plbConfig.centralLayerEnvelopeZ             = pixelEnvelopeZ;
///  plbConfig.centralLayerMaterialConcentration = pixelMaterialConcentration;
///  plbConfig.centralLayerMaterialProperties    = pixelMaterialProperties;
///  plbConfig.centralModulePositionPhi          = pcLayerModulesPositionPhi;
///  plbConfig.centralModuleTiltPhi              = pcLayerModulesTiltPhi;
///  plbConfig.centralModulePositionZ            = pcLayerModulesPositionZ;
///  plbConfig.centralModuleStaggerZ             = pcLayerModulesStaggerZ;
///  plbConfig.centralModuleHalfX                = pcLayerModulesMaxHalfX;
///  plbConfig.centralModuleHalfY                = pcLayerModulesHalfY;
///  plbConfig.centralModuleThickness            = pcLayerModulesThickness;
///  plbConfig.centralModuleMaterial             = pcLayerModulesMaterial;
///  //-------------------------------------------------------------------------------------
///  // configure end caps
///  plbConfig.posnegLayerPositionsZ            = pixelPositionsZ;
///  plbConfig.posnegLayerEnvelopeR             = pixelEnvelopesR;
///  plbConfig.posnegLayerMaterialConcentration = pixelMaterialConcentration;
///  plbConfig.posnegLayerMaterialProperties    = pixelMaterialProperties;
///  plbConfig.posnegModuleRadii                = ppnLayerModulesRadii;
///  plbConfig.posnegModuleStaggerR             = ppnLayerModulesStaggerR;
///  plbConfig.posnegModuleInPhi                = ppnLayerModulesInPhi;
///  plbConfig.posnegModulePositionPhi          = ppnLayerModulesPositionPhi;
///  plbConfig.posnegModuleStaggerPhi           = ppnLayerModulesStaggerPhi;
///  plbConfig.posnegModuleMinHalfX             = ppnLayerModulesMinHalfX;
///  plbConfig.posnegModuleMaxHalfX             = ppnLayerModulesMaxHalfX;
///  plbConfig.posnegModuleHalfY                = ppnLayerModulesHalfY;
///  plbConfig.posnegModuleThickness            = ppnLayerModulesThickness;
///  plbConfig.posnegModuleMaterial             = ppnLayerModulesMaterial;
  auto pixelLayerBuilder
      = std::make_shared<GenericLayerBuilder>(plbConfig);
  //-------------------------------------------------------------------------------------
  // build the pixel volume
  CylinderVolumeBuilder::Config pvbConfig;
  pvbConfig.logger = Acts::getDefaultLogger("CylinderVolumeBuilder", lvl);
  pvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  pvbConfig.volumeName           = "Pixel";
  pvbConfig.volumeToBeamPipe     = false;
  pvbConfig.layerBuilder         = pixelLayerBuilder;
  pvbConfig.layerEnvelopeR       = 1.;
  pvbConfig.layerEnvelopeZ       = 10.;
  pvbConfig.volumeSignature      = 0;
  auto pixelVolumeBuilder
      = std::make_shared<CylinderVolumeBuilder>(pvbConfig);

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  CylinderGeometryBuilder::Config tgConfig;
  tgConfig.logger = getDefaultLogger("CylinderGeometryBuilder", lvl);
  tgConfig.beamPipeBuilder = beamPipeVolumeBuilder;
  tgConfig.trackingVolumeBuilders = {pixelVolumeBuilder};
  tgConfig.trackingVolumeHelper   = cylinderVolumeHelper;
  auto cylinderGeometryBuilder
      = std::make_shared<const CylinderGeometryBuilder>(tgConfig);
  return cylinderGeometryBuilder->trackingGeometry();
}

std::vector< Acts::Vector3D >
modulePositionsCylinder(double radius,
                        double zStagger,
                        double lModule,
                        double lOverlap,
                        const std::pair<int,int>& binningSchema)
{
  int nPhiBins = binningSchema.first;
  int nZbins = binningSchema.second;
  // prepare the return value
  std::vector< Vector3D > mPositions;
  mPositions.reserve(nPhiBins*nZbins);
  // prep work
  double phiStep = 2*M_PI/(nPhiBins+1);
  double minPhi  = -M_PI + 0.5*phiStep;
  double zStart  = -0.5*(nZbins-1)*(lModule-lOverlap); 
  double zStep   =  fabs(zStart)/(nZbins-1);
  // loop over the bins
  for (size_t zBin = 0; zBin < size_t(nZbins); ++zBin){
    // prepare z and r
    double moduleZ = zStart + zBin*zStep;
    double moduleR = (zBin%2) ? radius : radius + zStagger;
    for (size_t phiBin = 0; phiBin < size_t(nPhiBins); ++phiBin){
      // calculate the current phi value
      double modulePhi = minPhi + phiBin*phiStep;
      mPositions.push_back( Vector3D(moduleR*cos(modulePhi),moduleR*sin(modulePhi),moduleZ));
      }
   }
}
                          

}  // end of namespace Acts
