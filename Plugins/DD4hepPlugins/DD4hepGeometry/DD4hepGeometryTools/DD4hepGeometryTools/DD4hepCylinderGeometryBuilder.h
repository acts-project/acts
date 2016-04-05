///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// DD4hepPlugin
#include "DD4hepGeometryInterfaces/IDD4hepGeometrySvc.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geometry module
#include "GeometryInterfaces/ITrackingGeometryBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeHelper.h"
#include "Detector/Layer.h"
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
     
     @author julia.hrdinka@cern.ch, andreas.salzburger@cern.ch
     */
    
    
    class DD4hepCylinderGeometryBuilder :  public Acts::AlgToolBase, virtual public Acts::ITrackingGeometryBuilder {
        
    public:
        /** Constructor */
        DD4hepCylinderGeometryBuilder(const std::string& t,const std::string& n,const IInterface* p);
        
        /** Destructor */
        virtual ~DD4hepCylinderGeometryBuilder();
        
        /** AlgTool initialize method */
        StatusCode initialize() override;
        
        /** AlgTool finalize method */
        StatusCode finalize() override;
        
        /** TrackingGeometry Interface method */
        Acts::TrackingGeometry* trackingGeometry() const override;
        
    private:
        /** Handle to the service providing the DD4hep geometry */
        ServiceHandle<IDD4hepGeometrySvc>            m_DD4hepGeometrySvc;
        /** Handle yo the volume building tool */
        ToolHandle<Acts::ITrackingVolumeBuilder>     m_volumeBuilder;
        /** Handle to the tool helping to build the tracking volumes */
        ToolHandle<Acts::ITrackingVolumeHelper>      m_volumeHelper;
        /** Helper class to build layers */
        DD4hepLayerHelper*                           m_layerHelper;
        
    };
}

#endif //#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
