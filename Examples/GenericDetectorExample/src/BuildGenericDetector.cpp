// ACTS include(s)
#include "ACTS/Examples/GenericDetectorExample/BuildGenericDetector.hpp"
#include "ACTS/Examples/GenericDetectorExample/GenericLayerBuilder.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/PassiveLayerBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/CylinderGeometryBuilder.hpp"
#include "ACTS/Tools/TrackingVolumeArrayCreator.hpp"
#include "ACTS/Tools/LayerArrayCreator.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"

#include <iostream>

namespace Acts
{
  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry()
  {  
      // configure the layer creator with the surface array creator
      Acts::LayerCreator::Config lcConfig;
      lcConfig.surfaceArrayCreator = std::make_shared<Acts::SurfaceArrayCreator>();
      auto LayerCreator = std::make_shared<Acts::LayerCreator>(lcConfig);
      // configure the cylinder volume helper
      Acts::CylinderVolumeHelper::Config cvhConfig;
      cvhConfig.layerArrayCreator             = std::make_shared<Acts::LayerArrayCreator>();
      cvhConfig.trackingVolumeArrayCreator    = std::make_shared<Acts::TrackingVolumeArrayCreator>();
      auto CylinderVolumeHelper = std::make_shared<Acts::CylinderVolumeHelper>(cvhConfig);
      //-------------------------------------------------------------------------------------
      //beam pipe
      //-------------------------------------------------------------------------------------
      // configure the beam pipe layer builder
      Acts::PassiveLayerBuilder::Config bplConfig;
      bplConfig.layerIdentification           = "BeamPipe";
      bplConfig.centralLayerRadii             = std::vector<double>(1,21.);
      bplConfig.centralLayerHalflengthZ       = std::vector<double>(1,200.);
      bplConfig.centralLayerThickness         = std::vector<double>(1,0.8);
      bplConfig.centralLayerMaterialX0        = std::vector<double>(1,352.8);
      bplConfig.centralLayerMaterialL0        = std::vector<double>(1,407.);
      bplConfig.centralLayerMaterialA         = std::vector<double>(1,9.012);
      bplConfig.centralLayerMaterialZ         = std::vector<double>(1,4.);
      bplConfig.centralLayerMaterialRho       = std::vector<double>(1,1.848e-3);
      auto BeamPipeBuilder = std::make_shared<Acts::PassiveLayerBuilder>(bplConfig);
      // create the volume for the beam pipe
      Acts::CylinderVolumeBuilder::Config bpvConfig;
      bpvConfig.trackingVolumeHelper          = CylinderVolumeHelper;
      bpvConfig.volumeName                    = "BeamPipe";
      bpvConfig.layerEnvelopeR                = 1.;
      bpvConfig.layerEnvelopeZ                = 1.;
      auto BeamPipeVolumeBuilder = std::make_shared<Acts::CylinderVolumeBuilder>(bpvConfig);
      //-------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------
      //pixel detector
      //-------------------------------------------------------------------------------------
      // configure pixel layer builder
      Acts::GenericLayerBuilder::Config plbConfig;
      plbConfig.layerCreator                      = LayerCreator;
      plbConfig.layerIdentification               = "Pixel";
      //fill necessary vectors for configuration
      //-------------------------------------------------------------------------------------
      //central
      int pcLayers = 3;
      std::vector<double> pixelLayerRadii {29.,55.,88.};
      
      std::vector<double> pixelEnvelopeZ(pcLayers,2.);
      
      std::vector<double> pixelMaterialConcentration(pcLayers,1.);
      
      std::vector<double> pixelMaterialParamters = {1.,95.7,465.2,28.03,14.,2.32e-3}; //thickness,X0,L0,A,Z,Rho
      std::vector<std::vector<double>> pixelMaterialProperties(pcLayers,pixelMaterialParamters);
      //phi properties centrallayers
      double pcNumPhi[pcLayers] = {15,24,40};
      std::vector<double> pcLayerModulesTiltPhi = {0.18,0.18,0.2};
      std::vector<std::vector<double>> pcLayerModulesPositionPhi;
      for (int i=0; i<pcLayers;i++) {
          std::vector<double> positionsPhi;
          double phiStep = 2.*M_PI/pcNumPhi[i];
          for (int j=0; j<pcNumPhi[i]; j++) {
              positionsPhi.push_back(-M_PI+j*phiStep);
          }
          pcLayerModulesPositionPhi.push_back(positionsPhi);
      }
      //z properties centrallayers
      double pcNumZ[pcLayers]                   = {13,13,13};
      double pcOverlapZ[pcLayers]               = {2.,2.,2.};
      std::vector<double> pcLayerModulesHalfY   = {32.,32.,32.};
      
      std::vector<std::vector<double>> pcLayerModulesPositionZ;
      for (int i=0; i<pcLayers; i++) {
          std::vector<double> positionsZ;
          double firstToLast = (pcNumZ[i]-1)*(2.*pcLayerModulesHalfY.at(i)-pcOverlapZ[i]);
          double zStep = firstToLast/(pcNumZ[i]-1);
          for (int j=0; j<pcNumZ[i]; j++) {
              positionsZ.push_back(-0.5*firstToLast + j*zStep);
          }
          pcLayerModulesPositionZ.push_back(positionsZ);
      }
      std::vector<double> pcLayerModulesStaggerZ(pcLayers,0.5);
      std::vector<double> pcLayerModulesMinHalfX();
      std::vector<double> pcLayerModulesMaxHalfX(pcLayers,8.4);
      std::vector<double> pcLayerModulesThickness(pcLayers,0.15);
      std::vector<double> pcModuleMaterialParameters = {95.7,465.2,28.03,14.,2.32e-3}; //X0,L0,A,Z,Rho
      std::vector<std::vector<double>> pcLayerModulesMaterial(pcLayers,pcModuleMaterialParameters);
      //-------------------------------------------------------------------------------------
      //posneg
      int ppnLayers=3;
      int ppnRings =1;
      std::vector<double> pixelPositionsZ = {500.,580.,650.};
      
      std::vector<double> pixelEnvelopesR(ppnLayers,5.);
      std::vector<double> ppnRadii(ppnRings,65.);
      std::vector<std::vector<double>> ppnLayerModulesRadii(ppnLayers,ppnRadii);
      std::vector<double> ppnLayerModulesStaggerR(ppnLayers,3.);
      //phi properties posneglayers
      std::vector<double> ppnInPhi = {24.,24.,24.};
      std::vector<std::vector<double>> ppnLayerModulesInPhi(ppnLayers,ppnInPhi);
      std::vector<std::vector<double>> ppnPositionPhi;
      for (int i=0; i<ppnRings;i++) {
          std::vector<double> positionsPhi;
          double phiStep = 2.*M_PI/ppnInPhi.at(i);
          for (int j=0; j<ppnInPhi.at(i); j++) {
              positionsPhi.push_back(-M_PI+j*phiStep);
          }
          ppnPositionPhi.push_back(positionsPhi);
      }
      std::vector<std::vector<std::vector<double>>> ppnLayerModulesPositionPhi(ppnLayers,ppnPositionPhi);
      std::vector<double> ppnStaggerPhi(ppnRings,0.5);
      std::vector<std::vector<double>> ppnLayerModulesStaggerPhi(ppnLayers,ppnStaggerPhi);
      //Module properties
      std::vector<std::vector<double>> ppnLayerModulesMinHalfX(ppnLayers);
      std::vector<double> ppnMaxHalfX(ppnRings,8.4);
      std::vector<std::vector<double>> ppnLayerModulesMaxHalfX(ppnLayers,ppnMaxHalfX);
      std::vector<double> ppnMaxHalfY(ppnRings,32.);
      std::vector<std::vector<double>> ppnLayerModulesHalfY(ppnLayers,ppnMaxHalfY);
      std::vector<double> ppnThickness(ppnRings,0.15);
      std::vector<std::vector<double>> ppnLayerModulesThickness(pcLayers,ppnThickness);
      std::vector<double> ppnModuleMaterialParameters = {95.7,465.2,28.03,14.,2.32e-3}; //X0,L0,A,Z,Rho
      std::vector<std::vector<double>> ppnMaterial(ppnLayers,ppnModuleMaterialParameters);
      std::vector<std::vector<std::vector<double>>> ppnLayerModulesMaterial(ppnLayers,ppnMaterial);
      
      //-------------------------------------------------------------------------------------
      
      //configure the central barrel
      plbConfig.centralLayerRadii                 = pixelLayerRadii;
      plbConfig.centralLayerEnvelopeZ             = pixelEnvelopeZ;
      plbConfig.centralLayerMaterialConcentration = pixelMaterialConcentration;
      plbConfig.centralLayerMaterialProperties    = pixelMaterialProperties;
      plbConfig.centralModulePositionPhi          = pcLayerModulesPositionPhi;
      plbConfig.centralModuleTiltPhi              = pcLayerModulesTiltPhi;
      plbConfig.centralModulePositionZ            = pcLayerModulesPositionZ;
      plbConfig.centralModuleStaggerZ             = pcLayerModulesStaggerZ;
      plbConfig.centralModuleHalfX                = pcLayerModulesMaxHalfX;
      plbConfig.centralModuleHalfY                = pcLayerModulesHalfY;
      plbConfig.centralModuleThickness            = pcLayerModulesThickness;
      plbConfig.centralModuleMaterial             = pcLayerModulesMaterial;
      //-------------------------------------------------------------------------------------
      //configure end caps
      plbConfig.posnegLayerPositionsZ             = pixelPositionsZ;
      plbConfig.posnegLayerEnvelopeR              = pixelEnvelopesR;
      plbConfig.posnegLayerMaterialConcentration  = pixelMaterialConcentration;
      plbConfig.posnegLayerMaterialProperties     = pixelMaterialProperties;
      plbConfig.posnegModuleRadii                 = ppnLayerModulesRadii;
      plbConfig.posnegModuleStaggerR              = ppnLayerModulesStaggerR;
      plbConfig.posnegModuleInPhi                 = ppnLayerModulesInPhi;
      plbConfig.posnegModulePositionPhi           = ppnLayerModulesPositionPhi;
      plbConfig.posnegModuleStaggerPhi            = ppnLayerModulesStaggerPhi;
      plbConfig.posnegModuleMinHalfX              = ppnLayerModulesMinHalfX;
      plbConfig.posnegModuleMaxHalfX              = ppnLayerModulesMaxHalfX;
      plbConfig.posnegModuleHalfY                 = ppnLayerModulesHalfY;
      plbConfig.posnegModuleThickness             = ppnLayerModulesThickness;
      plbConfig.posnegModuleMaterial              = ppnLayerModulesMaterial;
      auto PixelLayerBuilder = std::make_shared<Acts::GenericLayerBuilder>(plbConfig);
      //-------------------------------------------------------------------------------------
      //build the pixel volume
      Acts::CylinderVolumeBuilder::Config pvbConfig;
      pvbConfig.trackingVolumeHelper              = CylinderVolumeHelper;
      pvbConfig.volumeName                        = "Pixel";
      pvbConfig.volumeToBeamPipe                  = false;
      pvbConfig.layerBuilder                      = PixelLayerBuilder;
      pvbConfig.layerEnvelopeR                    = 1.;
      pvbConfig.layerEnvelopeZ                    = 10.;
      auto PixelVolumeBuilder = std::make_shared<Acts::CylinderVolumeBuilder>(pvbConfig);
      //-------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------
      //create the tracking geometry
      Acts::CylinderGeometryBuilder::Config tgConfig;
      tgConfig.beamPipeBuilder                    = BeamPipeVolumeBuilder;
      std::list<std::shared_ptr<Acts::ITrackingVolumeBuilder> > trackingVolumeBuilders;
      trackingVolumeBuilders.push_back(PixelVolumeBuilder);
      tgConfig.trackingVolumeBuilders             = trackingVolumeBuilders;
      tgConfig.trackingVolumeHelper               = CylinderVolumeHelper;
      auto CylinderGeometryBuilder = std::make_shared<const Acts::CylinderGeometryBuilder>(tgConfig);
      return CylinderGeometryBuilder->trackingGeometry();
  }
} // end of namespace Acts

