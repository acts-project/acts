 ///////////////////////////////////////////////////////////////////
// CylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_TRACKINGGEOMETRYBUILDER_H
#define ACTS_GEOMETRYTOOLS_TRACKINGGEOMETRYBUILDER_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "Algebra/AlgebraDefinitions.h"
// Geometry module
#include "GeometryInterfaces/ITrackingGeometryBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeBuilder.h"
#include "GeometryInterfaces/ITrackingVolumeHelper.h"
// Gaudi
#include "GaudiKernel/ToolHandle.h"

#ifdef ACTS_GEOMETRY_MEMUSAGE   
#include "GeometryUtils/MemoryLogger.h"
#endif  

namespace Acts {

    class TrackingGeometry;

    /** @class GeometryBuilder

      The Acts::TrackingGeometry Builder for volumes that wrap around another

      It retrieves ITrackingVolumeBuilders as configured and builds one
      detector around the other one.

      @TODO Julia: currently using work arround to use ToolHandleArray and private Tools because there is a bug in Gaudi. After update of Gaudi version go back to old usage -> currently marked with //bug

      @author Andreas.Salzburger@cern.ch   
     */

    class CylinderGeometryBuilder : public AlgToolBase, virtual public ITrackingGeometryBuilder {

      public:
        /** Constructor */
        CylinderGeometryBuilder(const std::string&,const std::string&,const IInterface*);
        
        /** Destructor */
        virtual ~CylinderGeometryBuilder();

        /** AlgTool initialize method */
        virtual StatusCode initialize() override;
        
        /** AlgTool finalize method */
        virtual StatusCode finalize() override;
        
        /** TrackingGeometry Interface method */
        TrackingGeometry* trackingGeometry() const override;

      private:

#ifdef ACTS_GEOMETRY_MEMUSAGE         
        MemoryLogger                              m_memoryLogger;           //!< in case the memory is logged
#endif                                            
        // -------------------------- Tools for geometry building ------------------------------------------------------ //
        ToolHandle<ITrackingVolumeBuilder>        m_beamPipeBuilder;        //!< a special builder for the beam pipe (for post-insertion)
#ifndef ACTS_GAUDI
        ToolHandleArray<ITrackingVolumeBuilder>   m_trackingVolumeBuilders; //!< the sub detector TrackingVolume builder
#else
        // list of tools to test
        typedef std::vector<std::string> ToolList;
        ToolList m_buildingTools;
#endif
        ToolHandle<ITrackingVolumeHelper>         m_trackingVolumeHelper;   //!< used for creating a container

    };

} // end of namespace

#endif // ACTS_GEOMETRYTOOLS_TRACKINGGEOMETRYBUILDER_H

