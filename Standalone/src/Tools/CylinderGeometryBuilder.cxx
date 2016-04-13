///////////////////////////////////////////////////////////////////
// CylinderGeometryBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core include
#include "Algebra/AlgebraDefinitions.h"
// Geometry module
#include "GeometryTools/CylinderGeometryBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeHelper.h"
#include "Detector/TrackingGeometry.h"
#include "Detector/TrackingVolume.h"
#include "Volumes/CylinderVolumeBounds.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SystemOfUnits.h"

// STD/STL
#ifdef ACTS_GEOMETRY_MEMUSAGE                            
#include <unistd.h>
#endif


DECLARE_TOOL_FACTORY(Acts::CylinderGeometryBuilder)

// constructor
Acts::CylinderGeometryBuilder::CylinderGeometryBuilder(const std::string& t, const std::string& n, const IInterface* p)
: Acts::AlgToolBase(t,n,p)
//bug m_trackingVolumeBuilders(this)
//bug  m_beamPipeBuilder(this)
#ifndef ACTS_GAUDI
  ,m_trackingVolumeBuilders()
#else
  ,m_buildingTools()
#endif
//bug  m_trackingVolumeHelper(this)
#ifdef ACTS_GEOMETRY_MEMUSAGE
 , m_memoryLogger()
#endif
{
    declareInterface<ITrackingGeometryBuilder>(this);
    // retrieve the tools
#ifndef ACTS_GAUDI
    declareProperty("TrackingVolumeBuilders",  m_trackingVolumeBuilders);
#else
    declareProperty("TrackingVolumeBuilders", m_buildingTools);
#endif
    declareProperty("BeamPipeBuilder",         m_beamPipeBuilder);
    declareProperty("TrackingVolumeHelper",    m_trackingVolumeHelper);
}

// destructor
Acts::CylinderGeometryBuilder::~CylinderGeometryBuilder()
{}

// Athena standard methods
// initialize
StatusCode Acts::CylinderGeometryBuilder::initialize()
{
    //Tool needs to be initialized
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
    // screen output in debug mode
    MSG_DEBUG( "initialize()" );
    // BeamPipe & Helper ======================================================================
    RETRIEVE_NONEMPTY_FATAL(m_beamPipeBuilder);

    // Geometries =============================================================================
#ifndef ACTS_GAUDI
    RETRIEVE_FATAL(m_trackingVolumeBuilders);
#endif
    RETRIEVE_NONEMPTY_FATAL(m_trackingVolumeHelper);
    // return SUCCESS at this stage
    return StatusCode::SUCCESS;
}

// finalize
StatusCode Acts::CylinderGeometryBuilder::finalize()
{
    MSG_DEBUG( "finalize()" );
    return StatusCode::SUCCESS;
}

Acts::TrackingGeometry* Acts::CylinderGeometryBuilder::trackingGeometry() const
{
    // the return geometry -- and the highest volume 
    TrackingGeometry* trackingGeometry = nullptr;
    TrackingVolumePtr    highestVolume = nullptr;

#ifdef ACTS_GEOMETRY_MEMUSAGE       
    m_memoryLogger.refresh(getpid());
    MSG_INFO( "[ memory usage ] Start of TrackingGeometry building: "  );    
    MSG_INFO( m_memoryLogger );                     
#endif  

#ifndef ACTS_GAUDI
    // loop over the builders and wrap one around the other -----------------------------
    //bug
    for (auto& volumeBuilder : m_trackingVolumeBuilders) {
     // assign a new highest volume (and potentially wrap around the given highest volume so far)
     highestVolume = volumeBuilder->trackingVolume(highestVolume);
#ifdef TRKDETDESCR_MEMUSAGE
    m_memoryLogger.refresh(getpid());
    MSG_INFO( "[ memory usage ] After sub TrackingVolume building: "  );
    MSG_INFO( m_memoryLogger );
#endif
    } // --------------------------------------------------------------------------------
#else
    // temporary solution for Toolhandle array bug ******************************
    for (auto it : m_buildingTools) {
        ToolHandle<ITrackingVolumeBuilder> volumeBuilder(it);
        if (volumeBuilder.retrieve().isFailure()) MSG_FATAL("Could not retrieve volumeBuilder Tool");
        highestVolume = volumeBuilder->trackingVolume(highestVolume);
        
#ifdef ACTS_GEOMETRY_MEMUSAGE
        m_memoryLogger.refresh(getpid());
        MSG_INFO( "[ memory usage ] After sub TrackingVolume building: "  );
        MSG_INFO( m_memoryLogger );
#endif
    } // --------------------------------------------------------------------------------
#endif
#ifdef ACTS_GEOMETRY_MEMUSAGE
    m_memoryLogger.refresh(getpid());
    MSG_INFO( "[ memory usage ] End of TrackingGeometry building: "  );    
    MSG_INFO( m_memoryLogger );                     
#endif

    // if you have a highst volume, stuff it into a TrackingGeometry
    if (highestVolume) {
        // see if the beampipe needs to be wrapped
        if (!m_beamPipeBuilder.empty() && !m_trackingVolumeHelper.empty()){
            // some screen output
            MSG_DEBUG("BeamPipe is being built and inserted.");
            // cast to cylinder volume bounds
            const CylinderVolumeBounds* cvB = dynamic_cast<const CylinderVolumeBounds*>(&(highestVolume->volumeBounds()));
            if (cvB){
                // get the inner radius
                double innerR = cvB->innerRadius();
                double halfZ  = cvB->halflengthZ();
                // create bounds for the innermost Volume
                VolumeBoundsPtr   beamPipeBounds(new CylinderVolumeBounds(0.,innerR,halfZ));
                TrackingVolumePtr beamPipeVolume = m_beamPipeBuilder->trackingVolume(nullptr,beamPipeBounds);
                // update the highest volume with the beam pipe
                highestVolume = m_trackingVolumeHelper->createContainerTrackingVolume({beamPipeVolume,highestVolume});
            }
        }
        // create the TrackingGeoemtry
        trackingGeometry = new Acts::TrackingGeometry(highestVolume);      
    }
    // return the geometry to the service
    return trackingGeometry;
} 
