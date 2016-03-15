//////////////////////////////////////////////////////////////////
// TrackingGeometryTest.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Examples module
#include "GeometryBuildingTest/TrackingGeometryTest.h"
// Geometry Module
#include "GeometryInterfaces/IGeometryProcessor.h"
#include "GeometryInterfaces/ITrackingGeometrySvc.h"
#include "Detector/TrackingGeometry.h"
#include "Detector/TrackingVolume.h"
#include "Detector/Layer.h"

// monitor memory usage
#ifdef TRKDETDESCR_MEMUSAGE   
#include <unistd.h>
#endif

Acts::TrackingGeometryTest::TrackingGeometryTest(const std::string& name, ISvcLocator* pSvcLocator) :
 Acts::AlgorithmBase(name, pSvcLocator),
#ifdef TRKDETDESCR_MEMUSAGE   
   m_memoryLogger(),
#endif
   m_executed(false),
   m_trackingGeometrySvc("TrackingGeometrySvc","AtlasTrackingGeometrySvc"),
   m_trackingGeometry(0),
   m_trackingGeometryName(""),
   m_trackingGeometryProcessors()
 {
     // get the service handle for the TrackingGeometry
     declareProperty("TrackingGeometrySvc"          , m_trackingGeometrySvc);
     // get the tools for display and recording
     declareProperty("TrackingGeometryProcessors"   , m_trackingGeometryProcessors);     
 }

StatusCode Acts::TrackingGeometryTest::initialize() 
{
    #ifdef TRKDETDESCR_MEMUSAGE   
        m_memoryLogger.refresh(getpid());
        MSG_INFO( "[ memory usage ] Memory monitoring activated through TRKDETDESCR_MEMUSAGE " );
        MSG_INFO( "[ memory usage ] initialize (start) : "                                     ); 
        MSG_INFO( m_memoryLogger );
    #endif   

    RETRIEVE_FATAL(m_trackingGeometrySvc);
    m_trackingGeometryName = m_trackingGeometrySvc->trackingGeometryName();
    
    // The Processors -------------------------------------------------------------
    RETRIEVE_NONEMPTY_FATAL( m_trackingGeometryProcessors );
    
    #ifdef TRKDETDESCR_MEMUSAGE   
        m_memoryLogger.refresh(getpid());
        MSG_INFO( "[ memory usage ] initialize (end) : " ); 
        MSG_INFO( m_memoryLogger );
    #endif  
    
    return StatusCode::SUCCESS;
}

StatusCode Acts::TrackingGeometryTest::finalize() 
{
    
    #ifdef TRKDETDESCR_MEMUSAGE   
        m_memoryLogger.refresh(getpid());
        MSG_INFO( "[ memory usage ] finalize (end) : " ); 
        MSG_INFO( m_memoryLogger );
    #endif  
    
    return StatusCode::SUCCESS;
}

StatusCode Acts::TrackingGeometryTest::execute()
{
    MSG_VERBOSE("Running the TrackingGeometryTest");
    
    // ------------------------------- get the trackingGeometry at first place
    if (!m_trackingGeometry) {
        m_trackingGeometry = m_trackingGeometrySvc->trackingGeometry();
    }
    // only run if it didn't already run before
    if (!m_executed && m_trackingGeometry){
        // push the geometry through the different processors
        for (  auto& tgp : m_trackingGeometryProcessors ){
            MSG_INFO("Parse geometry with processor " << tgp->name() );
            if ((tgp->process(*m_trackingGeometry)).isFailure()){
                MSG_FATAL("Could not process the TrackingGeometry with '" << tgp->name() <<"'. Aborting test.");
                return StatusCode::FAILURE;
            }
        }
        m_executed = true;    
    }
    return StatusCode::SUCCESS;
       
}
