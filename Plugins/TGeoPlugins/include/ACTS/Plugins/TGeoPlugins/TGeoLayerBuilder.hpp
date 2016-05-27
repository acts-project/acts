//
//  TGeoLayerBuilder.hpp
//  ACTS-Development
//
//  Created by Andreas Salzburger on 26/05/16.
//
//

#ifndef ACTS_TGEOPLUGINS_TGEOLAYERBUILDER_H
#define ACTS_TGEOPLUGINS_TGEOLAYERBUILDER_H

#include "ACTS/Tools/ILayerBuilder.hpp"

class TGeoVolume;
class TGeoNode;

namespace Acts {
    
    class ILayerCreator;
    class TGeoDetectorElement;

    /** @class TGeoLayerBuilder
         works on the gGeoManager, as this is filled from GDML */
    class TGeoLayerBuilder : public ILayerBuilder {
    public:
        /**  Helper config structs for volume parsin */
        struct LayerConfig {
        public:
            std::string                                 layerName;
            std::string                                 sensorName;
            size_t                                      binsLoc0;
            size_t                                      binsLoc1;
            
            LayerConfig():
            layerName(""),
            sensorName(""),
            binsLoc0(0),
            binsLoc1(0)
            {}
        };
        
        /** @struct Config
         nested configuration struct */
        class Config {
        public:
            std::string                                 configurationName;
            std::shared_ptr<ILayerCreator>              layerCreator;
            std::vector<LayerConfig>                    negativeLayerConfigs;
            std::vector<LayerConfig>                    centralLayerConfigs;
            std::vector<LayerConfig>                    positiveLayerConfigs;
            
            Config() :
            configurationName("Undefined"),
            layerCreator(nullptr),
            negativeLayerConfigs(),
            centralLayerConfigs(),
            positiveLayerConfigs()
            {}
            
        };
        
        /** Constructor */
        TGeoLayerBuilder(const Config& config);
        
        /** Destructor */
        ~TGeoLayerBuilder();
        
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector negativeLayers() const final;
        
        /** LayerBuilder interface method - returning the central layers */
        const LayerVector centralLayers() const final;
        
        /** LayerBuilder interface method - returning the layers at negative side */
        const LayerVector positiveLayers() const final;
        
        /** Name identification */
        const std::string& identification() const final;
        
        /** set the configuration object */
        void setConfiguration(const Config& glbConfig);

        /** get the configuration object */
        Config getConfiguration() const;
        
    private:
        Config  m_config; //!< configruation object

        mutable std::vector< std::shared_ptr<Acts::TGeoDetectorElement> > m_elementStore; //!< TODO make clear where the TGeoDetectorElement lives
        
        
        /** Private helper function to parse the geometry tree */
        void collectSensitive(TGeoVolume* tgVolume, TGeoNode* tgNode, const std::string& sensitiveName, std::vector<TGeoNode*>& collectedVolumes) const;
    
    };
    
    inline TGeoLayerBuilder::Config TGeoLayerBuilder::getConfiguration() const
    { return m_config; }
    
    inline const std::string& TGeoLayerBuilder::identification() const { return m_config.configurationName; }

}


#endif /* TGeoLayerBuilder_h */
