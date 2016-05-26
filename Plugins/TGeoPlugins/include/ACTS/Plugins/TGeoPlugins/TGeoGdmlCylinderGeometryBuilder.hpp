//
//  TGeoGdmlCylinderGeometryBuilder.hpp
//  ACTS-Development
//
//  Created by Andreas Salzburger on 25/05/16.
//
//

#ifndef ACTS_TGEOPLUGIN_TGEOGDMLCYLINDERGEOMETRYBUILDER_h
#define ACTS_TGEOPLUGIN_TGEOGDMLCYLINDERGEOMETRYBUILDER_h

// Geometry module
#include <memory>
#include <string>
#include <vector>
#include "ACTS/Tools/ITrackingGeometryBuilder.hpp"

class TGeoManager;
class TGeoVolume;
class TGeoNode;
class TObject;

namespace Acts {
    
    class ITrackingVolumeBuilder;
    class ITrackingVolumeHelper;
    class ILayerCreator;
    
    /** @ class TGeoGdmlCylinderGeometryBuilder 

     */
    class TGeoGdmlCylinderGeometryBuilder : virtual public Acts::ITrackingGeometryBuilder {
        
    public:
        
        /**  Helper config structs for volume parsin */
        struct LayerConfig {
            std::string                                 layerName;
            std::string                                 sensorName;
        };
        
        struct VolumeConfig {
            std::string                                 volumeName;
            std::vector<VolumeConfig>                   volumeConfigs;
            std::vector<LayerConfig>                    layerConfigs;
            
        };
        
        
        /** @struct Config
         Configuration for the GdmlCylinderGeometryBuilder */
        struct Config {
        public:
            std::shared_ptr<ITrackingVolumeBuilder>     volumeBuilder;  //!< building the contained sub detectors
            std::shared_ptr<ITrackingVolumeHelper>      volumeHelper;   //!< helper tool needed for volume building
            std::shared_ptr<ILayerCreator>              layerCreator;   //!< layer creator
            std::string                                 gdmlFile;       //!< file location of the GDML input file
            std::vector<VolumeConfig>                   volumeConfigs;  //!< the volume configurations
            
            Config() :
            volumeBuilder(nullptr),
            volumeHelper(nullptr),
            gdmlFile("")
            {}
        };
        
        /** Constructor */
        TGeoGdmlCylinderGeometryBuilder(const Config& dgbConfig);
        
        /** Destructor */
        virtual ~TGeoGdmlCylinderGeometryBuilder();
        
        /** setting the builders and helpers with the configuration object*/
        void setConfiguration(const Config& dgbConfig);
        
        /** get the configuration object */
        Config getConfiguration() const;
        
        /** TrackingGeometry Interface method */
        std::unique_ptr<TrackingGeometry> trackingGeometry() const final;
        
    private:
        
        void collectSensitive(TGeoVolume* tgVolume, TGeoNode* tgNode, const std::string& sensitiveName, std::vector<TGeoNode*>& collectedVolumes) const;
        
        /** configuration object */
        Config                       m_config;
        TGeoManager*                 m_tgManager;
        
    };
    
    inline TGeoGdmlCylinderGeometryBuilder::Config TGeoGdmlCylinderGeometryBuilder::getConfiguration() const
    { return m_config; }
    
}

#endif /* TGeoGdmlCylinderGeometryBuilder_h */
