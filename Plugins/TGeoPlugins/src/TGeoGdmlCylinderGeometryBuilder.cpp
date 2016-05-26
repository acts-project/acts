//
//  TGeoGdmlCylinderGeometryBuilder.cpp
//  ACTS-Development
//
//  Created by Andreas Salzburger on 25/05/16.
//
//

#include <stdio.h>

#include "ACTS/Plugins/TGeoPlugins/TGeoGdmlCylinderGeometryBuilder.hpp"
#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
// Geometry module
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Utilities/MsgMacros.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"


#include "TGeoManager.h"

Acts::TGeoGdmlCylinderGeometryBuilder::TGeoGdmlCylinderGeometryBuilder(const Acts::TGeoGdmlCylinderGeometryBuilder::Config& dgbConfig) :
Acts::ITrackingGeometryBuilder(),
m_tgManager(new TGeoManager("TGeoManagerGDML","TGeo GDML Manager"))
{
    setConfiguration(dgbConfig);
}

Acts::TGeoGdmlCylinderGeometryBuilder::~TGeoGdmlCylinderGeometryBuilder()
{}

// configuration
void Acts::TGeoGdmlCylinderGeometryBuilder::setConfiguration(const Acts::TGeoGdmlCylinderGeometryBuilder::Config& dgbConfig)
{
    // @TODO check consistency
    // copy the configuration
    m_config = dgbConfig;
}


std::unique_ptr<Acts::TrackingGeometry> Acts::TGeoGdmlCylinderGeometryBuilder::trackingGeometry() const
{
    
    // get the TGeoManager
    TGeoManager::Import(m_config.gdmlFile.c_str());
    // get the volume
    
    std::vector<Acts::TGeoDetectorElement*> tgDetectorElement;

    //  now get the volumes chosen for translating
    for (auto volumeCfg : m_config.volumeConfigs){
        // translating layer volume ID
        MSG_INFO("Translating Volume with name " << volumeCfg.volumeName );
        for (auto layerCfg : volumeCfg.layerConfigs ){
            MSG_INFO("- layer configuration found for layer " << layerCfg.layerName << " with sensor " << layerCfg.sensorName );
            TGeoVolume* volume = gGeoManager->GetVolume(layerCfg.layerName.c_str());
            if (volume){
                // prepare the vector for the sensitive nodes
                std::vector<TGeoNode*> sensitiveNodes;
                // run recursive collection
                collectSensitive(volume,nullptr,layerCfg.sensorName,sensitiveNodes);
                MSG_INFO("- layer found to have " << sensitiveNodes.size() << " sensitive sensors ");
                // create teh detector surface vector
                std::vector<const Acts::Surface*> detSurfaces;
                detSurfaces.reserve(sensitiveNodes.size());
                // loop and fill
                for (auto& sNode : sensitiveNodes){
                    tgDetectorElement.push_back(new Acts::TGeoDetectorElement(Identifier(),sNode));
                    // create a layer out of the surfaces
                    //m_config.
                }
                
            }
            
        }
        
    }
    
    //return nullptr;
}

void Acts::TGeoGdmlCylinderGeometryBuilder::collectSensitive(TGeoVolume* tgVolume, TGeoNode* tNode, const std::string& sensitiveName, std::vector<TGeoNode*>& sensitiveNodes)const
{
    // if it's still the master volume
    if (tgVolume){
        auto daugthers = tgVolume->GetNodes();
        TIter iObj(daugthers);
        while (TObject* obj = iObj()) {
            // dynamic_cast to a node
            TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
            if (node)
                collectSensitive(nullptr, node, sensitiveName, sensitiveNodes);
        }
        
    } else {
        // get the name as a string and compare
        std::string tNodeName = tNode->GetName();
        MSG_VERBOSE("-- node :" << tNodeName );
        if (tNodeName.find(sensitiveName) != std::string::npos){
            // senstive volume found, collect it
            sensitiveNodes.push_back(tNode);
        }
        // get the children nodes from the
        auto daugthers = tNode->GetNodes();
        TIter iObj(daugthers);
        while (TObject* obj = iObj()) {
            // dynamic_cast to a node
            TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
            if (node)
                collectSensitive(nullptr, node, sensitiveName, sensitiveNodes);
            
        }
    }
}

