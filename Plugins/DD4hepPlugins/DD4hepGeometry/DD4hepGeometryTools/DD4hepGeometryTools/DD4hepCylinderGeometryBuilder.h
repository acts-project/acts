///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ATS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// DD4hepPlugin
#include "DD4hepGeometryInterfaces/IDD4hepGeometrySvc.h"
#include "DD4hepGeometryInterfaces/IDD4hepLayerHelper.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geometry module
#include "GeometryInterfaces/ITrackingGeometryBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeHelper.h"
#include "Detector/Layer.h"
// DD4hep
#include "DD4hep/Detector.h"

namespace Ats {
    class TrackingGeometry;
}
namespace Add4hep {
    
    /** @ class DD4hepCylinderGeometryBuilder
     
    This class receives the DD4hep geometry from the given implementation of the IDD4hepGeomtrySvc, walks through the subvolumes and initiates their translation into the tracking geometry.
        It returns the world tracking geometry element.
     
     @author julia.hrdinka@cern.ch, andreas.salzburger@cern.ch
     */
    
    
  class DD4hepCylinderGeometryBuilder :  public Ats::AlgToolBase {
        
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
        Ats::TrackingGeometry* trackingGeometry() const override;
        
    private:
        /** Handle to the service providing the DD4hep geometry */
        ServiceHandle<IDD4hepGeometrySvc>           m_DD4hepGeometrySvc;
        /** Handle yo the volume building tool */
        ToolHandle<Ats::ITrackingVolumeBuilder>     m_volumeBuilder;
        /** Handle to the tool helping to build the tracking volumes */
        ToolHandle<Ats::ITrackingVolumeHelper>      m_trackingVolumeHelper;
        /** Handle to the tool helping to build the tracking layers */
        ToolHandle<IDD4hepLayerHelper>              m_DD4hepLayerHelper;
        
        
    };
}

#endif //#ifndef ATS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
