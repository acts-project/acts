// ACTS include(s)
#include "ACTS/Examples/BuildGenericDetector.h"
#include "ACTS/Detector/TrackingGeometry.h"
#include "ACTS/Material/Material.h"
#include "ACTS/Tools/LayerCreator.h"
#include "ACTS/Tools/PassiveLayerBuilder.h"
#include "ACTS/Tools/CylinderVolumeBuilder.h"
#include "ACTS/Tools/CylinderVolumeHelper.h"
#include "ACTS/Tools/CylinderGeometryBuilder.h"
#include "ACTS/Tools/TrackingVolumeArrayCreator.h"
#include "ACTS/Tools/LayerArrayCreator.h"
#include "ACTS/Tools/SurfaceArrayCreator.h"
#include "ACTS/Tools/CylinderVolumeHelper.h"

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
    pLayerRadii[pcLayers]    = {29.,55.,88.};
    std::vector<double> pixelLayerRadii(pLayerRadii, pLayerRadii+sizeof(pLayerRadii)/sizeof(double));
    
    std::vector<double> pixelEnvelopeZ(pcLayers,2.);
    
    std::vector<double> pixelMaterialConcentration(pcLayers,1.);
    
    std::vector<std::vector<double>> pixelMaterialProperties;
    std::vector<double> pixelMaterialParamters;
    pixelMaterialParamters.push_back(1.);       //thickness
    pixelMaterialParamters.push_back(95.7);     //X0
    pixelMaterialParamters.push_back(465.2);    //L0
    pixelMaterialParamters.push_back(28.03);    //A
    pixelMaterialParamters.push_back(14.);      //Z
    pixelMaterialParamters.push_back(2.32e-3);  //Rho
    pixelMaterialProperties.push_back(pixelMaterialParamters);
    
    //phi properties centrallayers
    double[pcLayers] pcNumPhi = {15,24,40};
    double[pcLayers] pcTiltPhi = {0.18,0.18,0.2};
    std::vector<double> pcLayerModulesTiltPhi(pcTiltPhi, pcTiltPhi+sizeof(tiltPhi)/sizeof(double));
    std::vector<double> pcLayerModulesPositionPhi;
    for (int i=0; i<pcLayers;i++) {
      std::vector<double> positionsPhi;
      double phiStep = 2.*M_PI/pcNumPhi[i];
      for (int j=0; j<pcNumPhi[i]; j++) {
	positionsPhi.push_back(-M_PI+j*phiStep);
      }
      pcLayerModulesPositionPhi.push_back(positionsPhi);
    }
    
    //z properties centrallayers
    double[pcLayers] pcNumZ                   = {13,13,13};
    double[pcLayers] pcOverlapZ               = {2.,2.,2.};
    double[pcLayers] pcModuleHalfY            = {32.,32.,32.};
    
    std::vector<double> pcLayerModulesPositionZ;
    for (int i=0; i<pcLayers; i++) {
      std::vector<double> positionsZ;
      double firstToLast = (pcNumZ[i]-1)*(2.*pcModuleHalfY[i]-pcOverlapZ[i]);
      double zStep = firstToLast/(pcNumZ[i]-1);
      for (int j=0; j<pcNumZ[i]; j++) {
	positionsZ.push_back(-0.5*firstToLast + j*zStep);
      }
      pcLayerModulesPositionZ.push_back(positionsZ);
    }
    
    std::vector<double> pcLayerModulesStaggerZ(pcLayers,0.5);
    std::vector<double> ppnLayerModulesMinHalfX();
    std::vector<double> pcLayerModulesMaxHalfX(pcLayers,8.4);
    std::vector<double> pcLayerModulesHalfY(pcModuleHalfY, pcModuleHalfY+sizeof(pcModuleHalfY)/sizeof(double));
    std::vector<double> pcLayerModulesThickness(pcLayers,0.15);
    std::vector<double> pcLayerModulesMaterial;
    std::vector<double> pcModuleMaterialParamters;
    pcMaterialParamters.push_back(95.7);     //X0
    pcMaterialParamters.push_back(465.2);    //L0
    pcMaterialParamters.push_back(28.03);    //A
    pcMaterialParamters.push_back(14.);      //Z
    pcMaterialParamters.push_back(2.32e-3);  //Rho
    pcLayerModulesMaterial.push_back(pcModuleMaterialParamters);
    
    //-------------------------------------------------------------------------------------
    //posneg
    int ppnLayers=3;
    ppnlayerPosZ[ppnLayers] = {500.,580.,650.};
    std::vector<double> layerPositionsZ(ppnlayerPosZ,ppnlayerPosZ+sizeof(ppnlayerPosZ)/sizeof(double));
    
    std::vector<double> pixelEnvelopesR(ppnLayers,5.);
    std::vector<double> ppnLayerModulesRadii(ppnLayers,65.);
    std::vector<double> ppnLayerModulesStaggerR(ppnLayers,3.);
    //phi properties posneglayers
    double[ppnLayers] ppnNumPhi = {24.,24.,24.};
    std::vector<double> ppnLayerModulesInPhi(ppnNumPhi,ppnNumPhi + sizeof(ppnNumPhi)/sizeof(double));
    std::vector<double> ppnLayerModulesPositionPhi;
    for (int i=0; i<ppnLayers;i++) {
      std::vector<double> positionsPhi;
      double phiStep = 2.*M_PI/ppnNumPhi[i];
      for (int j=0; j<ppnNumPhi[i]; j++) {
	positionsPhi.push_back(-M_PI+j*phiStep);
      }
      ppnLayerModulesPositionPhi.push_back(positionsPhi);
    }
    std::vector<double> ppnLayerModulesStaggerPhi(ppnLayers,0.5);
    
    //Module properties
    std::vector<double> ppnLayerModulesMinHalfX();
    std::vector<double> ppnLayerModulesMaxHalfX(ppnLayers,8.4);
    std::vector<double> ppnLayerModulesHalfY(ppnLayers,32.);
    std::vector<double> ppnLayerModulesThickness(pcLayers,0.15);
    std::vector<double> ppnLayerModulesMaterial;
    std::vector<double> ppnModuleMaterialParamters;
    ppnMaterialParamters.push_back(95.7);     //X0
    ppnMaterialParamters.push_back(465.2);    //L0
    ppnMaterialParamters.push_back(28.03);    //A
    ppnMaterialParamters.push_back(14.);      //Z
    ppnMaterialParamters.push_back(2.32e-3);  //Rho
    ppnLayerModulesMaterial.push_back(ppnModuleMaterialParamters);

    
    //-------------------------------------------------------------------------------------
    
    //configure the central barrel
    plbConfig.centralLayerRadii                 = pixelLayerRadii;
    plbConfig.centralLayerEnvelopeZ             = pixelEnvelopeZ;
    plbConfig.centralLayerMaterialConcentration = pixelMaterialConcentration;
    plbConfig.centralLayerMaterialProperties    = pixelMaterialProperties;
    plbConfig.centralModulePositionPhi          = pcLayerModulesPositionPhi;
    plbConfig.centralMoudleTiltPhi              = pcLayerModulesTiltPhi;
    plbConfig.centralModulePositionZ            = pcLayerModulesPositionZ;
    plbConfig.centralModuleStaggerZ             = pcLayerModulesStaggerZ:
    plbConfig.centralModuleHalfX                = pcLayerModulesHalfX:
    plbConfig.centralModuleHalfY                = pcLayerModulesHalfY:
    plbConfig.centralModuleThickness            = pcLayerModulesThickness:
    plbConfig.centralModuleMaterial             = pcLayerModulesMaterial:
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
    vbConfig.volumeToBeamPipe                   = false;
    pvbConfig.layerBuilder                      = PixelLayerBuilder;
    pvbConfig.layerEnvelopeR                    = 1.;
    pvbConfig.layerEnvelopeZ                    = 10.;
    auto PixelVolumeBuilder = std::make_shared<Acts::CylinderVolumeBuilder>(pvbConfig);
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------
    
    //create the tracking geometry
    Acts::CylinderGeometryBuilder::Config tgConfig;
    tgConfig.beamPipeBuilder                    = BeamPipeBuilder;
    std::list<std::shared_ptr<Acts::CylinderGeometryBuilder>> trackingVolumeBuilders;
    trackingVolumeBuilders.push_back(PixelVolumeBuilder);
    tgConfig.trackingVolumeBuilders             = trackingVolumeBuilders;
    tgConfig.trackingVolumeHelper               = CylinderVolumeHelper;
    auto CylinderGeometryBuilder = std::make_shared<Acts::CylinderGeometryBuilder>(tgConfig);
    
    return CylinderGeometryBuilder->trackingGeometry();

  }
} // end of namespace Acts

