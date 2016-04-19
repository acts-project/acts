///////////////////////////////////////////////////////////////////
// CylinderGeometryBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core include
#include "ACTS/Utilities/AlgebraDefinitions.h"
// Geometry module
#include "ACTS/Tools/CylinderGeometryBuilder.h"
#include "ACTS/Tools/ITrackingVolumeBuilder.h"
#include "ACTS/Tools/ITrackingVolumeHelper.h"
#include "ACTS/TrackingGeometry.h"
#include "ACTS/Detector/TrackingVolume.h"
#include "ACTS/Volumes/CylinderVolumeBounds.h"

// STD/STL
#ifdef ACTS_GEOMETRY_MEMUSAGE                            
#include <unistd.h>
#endif


// constructor
Acts::CylinderGeometryBuilder::CylinderGeometryBuilder():
  Acts::ITrackingGeometryBuilder(),
#ifdef ACTS_GEOMETRY_MEMUSAGE
  m_memoryLogger(),
#endif
  m_beamPipeBuilder(),
  m_trackingVolumeBuilders(),
  m_trackingVolumeHelper()
{}

std::unique_ptr<Acts::TrackingGeometry> Acts::CylinderGeometryBuilder::trackingGeometry() const
{
  // the return geometry -- and the highest volume 
  std::unique_ptr<Acts::TrackingGeometry> trackingGeometry;
  TrackingVolumePtr    highestVolume = nullptr;

#ifdef ACTS_GEOMETRY_MEMUSAGE       
  m_memoryLogger.refresh(getpid());
    // MSG_INFO( "[ memory usage ] Start of TrackingGeometry building: "  );    
    // MSG_INFO( m_memoryLogger );                       
#endif  

  // loop over the builders and wrap one around the other -----------------------------
  //bug
  for (auto& volumeBuilder : m_trackingVolumeBuilders) {
    // assign a new highest volume (and potentially wrap around the given highest volume so far)
    highestVolume = volumeBuilder->trackingVolume(highestVolume);
#ifdef TRKDETDESCR_MEMUSAGE
    m_memoryLogger.refresh(getpid());
    // MSG_INFO( "[ memory usage ] After sub TrackingVolume building: "  );
    // MSG_INFO( m_memoryLogger );
#endif
  } // --------------------------------------------------------------------------------

#ifdef ACTS_GEOMETRY_MEMUSAGE
  m_memoryLogger.refresh(getpid());
  // MSG_INFO( "[ memory usage ] End of TrackingGeometry building: "  );    
  // MSG_INFO( m_memoryLogger );                     
#endif

  // if you have a highst volume, stuff it into a TrackingGeometry
  if (highestVolume) {
    // see if the beampipe needs to be wrapped
    if (m_beamPipeBuilder && m_trackingVolumeHelper){
      // some screen output
      // MSG_DEBUG("BeamPipe is being built and inserted.");
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
    trackingGeometry.reset(new Acts::TrackingGeometry(highestVolume));      
  }
  // return the geometry to the service
  return trackingGeometry;
} 
