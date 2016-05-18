///////////////////////////////////////////////////////////////////
// CylinderGeometryBuilder.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/CylinderGeometryBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

// constructor
Acts::CylinderGeometryBuilder::CylinderGeometryBuilder(const Acts::CylinderGeometryBuilder::Config& cgbConfig):
  m_config()
{    
  setConfiguration(cgbConfig);   
}

// configuration
void Acts::CylinderGeometryBuilder::setConfiguration(const Acts::CylinderGeometryBuilder::Config& cgbConfig)
{
  // @TODO check consistency    
  // copy the configuration 
  m_config = cgbConfig;
} 

std::unique_ptr<Acts::TrackingGeometry> Acts::CylinderGeometryBuilder::trackingGeometry() const
{
  // the return geometry -- and the highest volume
  std::unique_ptr<Acts::TrackingGeometry> trackingGeometry;
  TrackingVolumePtr    highestVolume = nullptr;
  // loop over the builders and wrap one around the other -----------------------------
  for (auto& volumeBuilder : m_config.trackingVolumeBuilders) {
    // assign a new highest volume (and potentially wrap around the given highest volume so far)
    highestVolume = volumeBuilder->trackingVolume(highestVolume);
  } // --------------------------------------------------------------------------------
  // if you have a highst volume, stuff it into a TrackingGeometry
  if (highestVolume) {
    // see if the beampipe needs to be wrapped
    if (m_config.beamPipeBuilder && m_config.trackingVolumeHelper){
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
        TrackingVolumePtr beamPipeVolume = m_config.beamPipeBuilder->trackingVolume(nullptr,beamPipeBounds);
        // update the highest volume with the beam pipe
        highestVolume = m_config.trackingVolumeHelper->createContainerTrackingVolume({beamPipeVolume,highestVolume});
      }
    }
    // create the TrackingGeoemtry
    trackingGeometry.reset(new Acts::TrackingGeometry(highestVolume));      
  }
  // return the geometry to the service
  return trackingGeometry;
} 
