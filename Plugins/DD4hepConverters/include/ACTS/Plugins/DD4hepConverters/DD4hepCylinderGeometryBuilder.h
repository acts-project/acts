///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

// DD4hepPlugin
#include "ACTS/Plugins/DD4hepConverters/IDD4hepGeometrySvc.h"
// Geometry module
#include "ACTS/Tools/ITrackingGeometryBuilder.h"
#include "ACTS/Tools/ITrackingVolumeBuilder.h"
#include "ACTS/Tools/ITrackingVolumeHelper.h"
#include "ACTS/Layers/Layer.h"
// DD4hep
#include "DD4hep/Detector.h"

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
        /** Constructor */
        DD4hepCylinderGeometryBuilder();
        
        /** Destructor */
        virtual ~DD4hepCylinderGeometryBuilder();
        
        /** setting the builders and helpers */
        void setVolumeBuilder(std::shared_ptr<Acts::ITrackingVolumeBuilder> volumeBuilder) {
            m_volumeBuilder = std::move(volumeBuilder);
        }
        
        void setVolumeHelper(std::shared_ptr<Acts::ITrackingVolumeHelper> volumeHelper) {
            m_volumeHelper = std::move(volumeHelper);
        }
        
        void setDD4hepGeometrySvc(std::shared_ptr<Acts::IDD4hepGeometrySvc> DD4hepGeometrySvc) {
            m_DD4hepGeometrySvc = DD4hepGeometrySvc;
        }
        
        /** TrackingGeometry Interface method */
        std::unique_ptr<TrackingGeometry> trackingGeometry() const override;
        
    private:
        /** Shared pointer to the service providing the DD4hep geometry */
        std::shared_ptr<IDD4hepGeometrySvc>               m_DD4hepGeometrySvc;
        /** Shared pointer yo the volume building tool */
        std::shared_ptr<Acts::ITrackingVolumeBuilder>     m_volumeBuilder;
        /** Shared pointer to the tool helping to build the tracking volumes */
        std::shared_ptr<Acts::ITrackingVolumeHelper>      m_volumeHelper;
        /** Shared pointer class to build layers */
        DD4hepLayerHelper*                                m_layerHelper;
        
    };
} //end of namespace

#endif //#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
