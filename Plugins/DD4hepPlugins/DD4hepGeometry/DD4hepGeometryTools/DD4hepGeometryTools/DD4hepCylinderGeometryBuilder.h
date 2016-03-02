///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ATS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

//DD4hepPlugin
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
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

namespace Ats {
    class TrackingGeometry;
}
namespace Add4hep {
    
    /** @ class DD4hepCylinderGeometryBuilder
     
     
     
     @author julia.hrdinka@cern.ch
     */
    
    
    class DD4hepCylinderGeometryBuilder : public Ats::AlgToolBase, virtual public Ats::ITrackingGeometryBuilder {
        
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
        
        ServiceHandle<IDD4hepGeometrySvc>           m_DD4hepGeometrySvc;
        ToolHandle<Ats::ITrackingVolumeBuilder>     m_volumeBuilder;
        ToolHandle<IDD4hepLayerHelper>              m_DD4hepLayerHelper;
        ToolHandle<Ats::ITrackingVolumeHelper>      m_trackingVolumeHelper;
        
    };
}

#endif //#ifndef ATS_ATS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
