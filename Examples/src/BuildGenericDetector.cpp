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
  bplConfig.centralLayerRadii       = std::vector<double>(1, 19.);
  bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 200.);
  bplConfig.centralLayerThickness   = std::vector<double>(1, 0.8);
  bplConfig.centralLayerMaterial = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
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
  plbConfig.logger              = getDefaultLogger("GenericLayerBuilder", lvl);
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
  plbConfig.centralLayerBinMultipliers = {1, 1};
  plbConfig.centralLayerRadii          = {29., 55., 88.};
  plbConfig.centralLayerEnvelopes      = {pcEnvelope, pcEnvelope, pcEnvelope};
  plbConfig.centralLayerMaterialConcentration = {1, 1, 1};
  plbConfig.centralLayerMaterialProperties
      = {pcmProperties, pcmProperties, pcmProperties};
  plbConfig.centralModuleBinningSchema = {{15, 13}, {24, 13}, {40, 13}};
  plbConfig.centralModuleTiltPhi       = {0.18, 0.18, 0.2};
  plbConfig.centralModuleHalfX         = {8.4, 8.4, 8.4};
  plbConfig.centralModuleHalfY         = {32., 32., 32.};
  plbConfig.centralModuleThickness     = {0.15, 0.15, 0.15};
  plbConfig.centralModuleMaterial      = {pcMaterial, pcMaterial, pcMaterial};
  plbConfig.centralModuleFrontsideStereo      = {};
  plbConfig.centralModuleBacksideStereo       = {};
  plbConfig.centralModuleBacksideGap          = {};
  // mPositions
  std::vector< std::vector< Vector3D> > centralModulePositions;
  for (size_t plb = 0; plb < plbConfig.centralLayerRadii.size(); ++plb){
    // call the helper function
    centralModulePositions.push_back(modulePositionsCylinder(plbConfig.centralLayerRadii[plb],
                                                             0.5, // 1 mm stagger
                                                             plbConfig.centralModuleHalfY[plb],
                                                             2., // 2 mm module overlap
                                                             plbConfig.centralModuleBinningSchema[plb]));
    
  }
  plbConfig.centralModulePositions            = centralModulePositions;
  //
  plbConfig.posnegLayerBinMultipliers          = { 1, 1 };
  plbConfig.posnegLayerPositionsZ              = {500., 580., 680.};
  plbConfig.posnegLayerEnvelopeR               = {5., 5., 5.};
  plbConfig.posnegLayerMaterialConcentration   = {1, 1, 1};
  plbConfig.posnegLayerMaterialProperties
      = {pcmProperties, pcmProperties, pcmProperties};
  plbConfig.posnegModuleMinHalfX               = {{8.4}, {8.4}, {8.4}};
  plbConfig.posnegModuleMaxHalfX               = {};
  plbConfig.posnegModuleHalfY                  = {{32.}, {32.}, {32.}};
  plbConfig.posnegModulePhiBins                = {{24}, {24}, {24}};
  plbConfig.posnegModuleThickness              = {{0.15}, {0.15}, {0.15}};
  plbConfig.posnegModuleMaterial = {{pcMaterial}, {pcMaterial}, {pcMaterial}};
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
                            93.,
                            plbConfig.posnegModulePhiBins[id],
                            plbConfig.posnegModuleHalfY[id]));
  }
  plbConfig.posnegModulePositions = posnegModulePositions;
  // define the builder
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
  tgConfig.beamPipeBuilder        = beamPipeVolumeBuilder;
  tgConfig.trackingVolumeBuilders = {pixelVolumeBuilder};
  tgConfig.trackingVolumeHelper   = cylinderVolumeHelper;
  auto cylinderGeometryBuilder
      = std::make_shared<const CylinderGeometryBuilder>(tgConfig);
  return cylinderGeometryBuilder->trackingGeometry();
}

    
/// helper method for cylinder
std::vector< Acts::Vector3D >
modulePositionsCylinder(double radius,
                        double zStagger,
                        double moduleHalfLength,
                        double lOverlap,
                        const std::pair<int,int>& binningSchema)
{
  int nPhiBins = binningSchema.first;
  int nZbins = binningSchema.second;
  // prepare the return value
  std::vector< Vector3D > mPositions;
  mPositions.reserve(nPhiBins*nZbins);
  // prep work
  double phiStep = 2*M_PI/(nPhiBins);
  double minPhi  = -M_PI + 0.5*phiStep;
  double zStart  = -0.5*(nZbins-1)*(2*moduleHalfLength-lOverlap); 
  double zStep   =  2*fabs(zStart)/(nZbins-1);
  // loop over the bins
  for (size_t zBin = 0; zBin < size_t(nZbins); ++zBin){
    // prepare z and r
    double moduleZ = zStart + zBin*zStep;
    double moduleR = (zBin%2) ? radius - 0.5*zStagger : radius + 0.5*zStagger;
    for (size_t phiBin = 0; phiBin < size_t(nPhiBins); ++phiBin){
      // calculate the current phi value
      double modulePhi = minPhi + phiBin*phiStep;
      mPositions.push_back( Vector3D(moduleR*cos(modulePhi),moduleR*sin(modulePhi),moduleZ));
     }
   }
    return mPositions;
}

/// helper method for disc 
std::vector<  std::vector< Acts::Vector3D > >
modulePositionsDisc(double z, 
                    double ringStagger, double phiStagger,
                    double innerRadius,
                    double outerRadius,
                    const std::vector< int >& discBinning,
                    const std::vector< double >& moduleHalfLength)
{
   
  // calculate the radii
  std::vector<double> radii;  
  // calculate the radial borders
  std::vector<double> radialBoarders;
  // the radial span of the disc
  double deltaR = outerRadius-innerRadius;
  // quick exits
  if (discBinning.size()==1){
    radii.push_back(0.5*(innerRadius+outerRadius));
    radialBoarders = { innerRadius, outerRadius };
  } else {
    double totalLength = 0;
    // sum up the total length
    for (auto& mhlength : moduleHalfLength)
      totalLength += 2*mhlength;
    // now calculate the overlap (equal pay)
    double rOverlap = (totalLength-deltaR)/(moduleHalfLength.size()-1);
    // and now fill the radii and gaps
    double lastR   = innerRadius;
    double lastHl  = 0.;
    double lastOl  = 0.;
    // remember the radial boarders   
    radialBoarders.push_back(innerRadius);
    // now calculate 
    for (auto& mhlength : moduleHalfLength){
      // calculate the radius
      radii.push_back(lastR+lastHl-lastOl+mhlength);
      lastR   = radii[radii.size()-1];
      lastOl  = rOverlap;
      lastHl  = mhlength;
      // and register the radial boarder
      radialBoarders.push_back(lastR+2*lastHl-0.5*lastOl);
    }
  }
  // now prepare the return method
  std::vector< std::vector< Vector3D > > mPositions;
  for (size_t ir = 0; ir < radii.size(); ++ir){
    // generate the z value
    double rz = radii.size() == 1 ? z : 
          ( ir%2 ? z-0.5*ringStagger : z+0.5*ringStagger);
    // fill the ring positions
    mPositions.push_back(modulePositionsRing(rz,radii[ir],phiStagger,discBinning[ir]));
  }
  return mPositions;
}  

/// Helper method for positioning
std::vector<  Acts::Vector3D >
modulePositionsRing(double z,
                    double radius,
                    double phiStagger,
                    int nPhiBins)
{
  // create and fill the positions
  std::vector < Vector3D > rPositions;
  rPositions.reserve(nPhiBins);
  // prep work
  double phiStep = 2*M_PI/(nPhiBins);
  double minPhi  = -M_PI + 0.5*phiStep;
  // phi loop
  for (size_t iphi = 0; iphi < size_t(nPhiBins); ++iphi){
    double phi = minPhi + iphi*phiStep;
    double rz  = iphi%2 ? z - 0.5*phiStagger : z + 0.5*phiStagger;
    rPositions.push_back(Vector3D(radius*cos(phi),radius*sin(phi),rz));
  }
  return rPositions;
}    

                          
}  // end of namespace Acts
