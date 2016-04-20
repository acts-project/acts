 ///////////////////////////////////////////////////////////////////
// CylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_TRACKINGGEOMETRYBUILDER_H
#define ACTS_GEOMETRYTOOLS_TRACKINGGEOMETRYBUILDER_H 1

// STL include(s)
#include <memory>
#include <list>

// ATS include(s)
#include "ACTS/Utilities/AlgebraDefinitions.h"
#include "ACTS/Tools/ITrackingGeometryBuilder.h"
#include "ACTS/Tools/ITrackingVolumeBuilder.h"
#include "ACTS/Tools/ITrackingVolumeHelper.h"

#ifdef ACTS_GEOMETRY_MEMUSAGE   
#include "ACTS/Utilities/MemoryLogger.h"
#endif  

namespace Acts
{
  class TrackingGeometry;

  /** @class GeometryBuilder

      The Acts::TrackingGeometry Builder for volumes that wrap around another

      It retrieves ITrackingVolumeBuilders as configured and builds one
      detector around the other one.

      @TODO Julia: currently using work arround to use ToolHandleArray and private Tools because there is a bug in Gaudi. After update of Gaudi version go back to old usage -> currently marked with //bug

      @author Andreas.Salzburger@cern.ch   
  */

  class CylinderGeometryBuilder : public ITrackingGeometryBuilder
  {
  public:
    /** Constructor */
    CylinderGeometryBuilder();
        
    /** Destructor */
    virtual ~CylinderGeometryBuilder() = default;

    /** TrackingGeometry Interface method */
    virtual std::unique_ptr<TrackingGeometry> trackingGeometry() const override;

    void setBeamPipeBuilder(std::shared_ptr<ITrackingVolumeBuilder> beamPipe)
    {
      m_beamPipeBuilder = std::move(beamPipe);
    }

    void setVolumeBuilders(std::list<std::shared_ptr<ITrackingVolumeBuilder> > builders)
    {
      m_trackingVolumeBuilders = std::move(builders);
    }

    void setTrackingVolumeHelper(std::shared_ptr<ITrackingVolumeHelper> helper)
    {
      m_trackingVolumeHelper = std::move(helper);
    }

  private:

#ifdef ACTS_GEOMETRY_MEMUSAGE         
    MemoryLogger                              m_memoryLogger;           //!< in case the memory is logged
#endif                                            
    // -------------------------- Tools for geometry building ------------------------------------------------------ //
    std::shared_ptr<ITrackingVolumeBuilder>               m_beamPipeBuilder;        //!< a special builder for the beam pipe (for post-insertion)
    std::list<std::shared_ptr<ITrackingVolumeBuilder> >   m_trackingVolumeBuilders; //!< the sub detector TrackingVolume builder
    std::shared_ptr<ITrackingVolumeHelper>                m_trackingVolumeHelper;   //!< used for creating a container
  };
} // end of namespace

#endif // ACTS_GEOMETRYTOOLS_TRACKINGGEOMETRYBUILDER_H

