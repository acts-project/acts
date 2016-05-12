///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

// DD4hepPlugin
#include "ACTS/Plugins/DD4hepPlugins/IDD4hepGeometrySvc.hpp"
// Geometry module
#include "ACTS/Tools/ITrackingGeometryBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Layers/Layer.hpp"
// DD4hep
#include "DD4hep/Detector.hpp"

namespace Acts {
    class TrackingGeometry;
}
namespace Acts {
    
    class DD4hepLayerHelper;
    
    /** @ class DD4hepCylinderGeometryBuilder
     
    This class receives the DD4hep geometry from the given implementation of the IDD4hepGeomtrySvc, walks through the subvolumes and initiates their translation into the tracking geometry.
        It returns the world tracking geometry element.
     @TODO find replacement for Gaudi exeption and message stream
     
     @author julia.hrdinka@cern.ch, andreas.salzburger@cern.ch
     */
    
    
    class DD4hepCylinderGeometryBuilder : virtual public Acts::ITrackingGeometryBuilder {
        
    public:
        /** @struct Config
         Configuration for the DD4hepCylinderGeometryBuilder */
        struct Config {
        
            std::shared_ptr<ITrackingVolumeBuilder>     volumeBuilder; //!< building the contained sub detectors
            std::shared_ptr<ITrackingVolumeHelper>      volumeHelper; //!< helper tool needed for volume building
            std::shared_ptr<IDD4hepGeometrySvc>         DD4hepGeometrySvc; //!< service providing the DD4hep geometry
            
            Config() :
                volumeBuilder(nullptr),
                volumeHelper(nullptr),
                DD4hepGeometrySvc(nullptr)
            {}
        };
        
        /** Constructor */
        DD4hepCylinderGeometryBuilder(const Config dgbConfig);
        
        /** Destructor */
        virtual ~DD4hepCylinderGeometryBuilder();
        
        /** setting the builders and helpers with the configuration object*/
        void setConfiguration(const Config dgbConfig);
        
        /** get the configuration object */
        Config getConfiguration() const;
        
        /** TrackingGeometry Interface method */
        std::unique_ptr<TrackingGeometry> trackingGeometry() const override;
        
    private:
        /** configuration object */
        Config                                            m_config;
        /** Shared pointer class to build layers */
        DD4hepLayerHelper*                                m_layerHelper;
        
    };
    
    inline DD4hepCylinderGeometryBuilder::Config DD4hepCylinderGeometryBuilder::getConfiguration() const
    { return m_config; }
} //end of namespace

#endif //#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
