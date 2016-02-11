//////////////////////////////////////////////////////////////////
// TrackingGeometryTest.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXAMPLES_TRACKINGGEOMETRYTEST_H
#define ATS_EXAMPLES_TRACKINGGEOMETRYTEST_H

// Core
#include "CoreInterfaces/AlgorithmBase.h"
// Gauid Include
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// Trk includes
#include "GeometryInterfaces/ITrackingGeometrySvc.h"

#ifdef TRKDETDESCR_MEMUSAGE
#include "TrkDetDescrUtils/MemoryLogger.h"
#endif

namespace Ats {
     
    class IGeometryProcessor;
    class TrackingGeometry;
        
    /** @class TrackingGeometryTest
       
        Test to build the TrackingGeometry and process it with analizers,
        thosw are of type IGeometryProcessor
        
        @author Andreas.Salzburger@cern.ch       
      */
      
    class TrackingGeometryTest : public AlgorithmBase  {
     public:

       /** Standard Athena-Algorithm Constructor */
       TrackingGeometryTest(const std::string& name, ISvcLocator* pSvcLocator);

       /* initialize the test, i.e. retrieve the TrackingGeometry Svc */
       StatusCode initialize() override;
       
       /* exectue - and set executed to true */
       StatusCode execute() override;

       /* finalize  */
       StatusCode finalize() override;
              
       
     private:
#ifdef TRKDETDESCR_MEMUSAGE
       MemoryLogger                            m_memoryLogger;               //!< little memory logger to check the used memory
#endif                                              
       bool                                         m_executed;                   //!< Make sure it only runs once 
                                                    
       ServiceHandle<ITrackingGeometrySvc>          m_trackingGeometrySvc;        //!< Service handle for retrieving the TrackingGeometry
       mutable const TrackingGeometry*              m_trackingGeometry;           //!< The TrackingGeometry to be retrieved
       std::string                                  m_trackingGeometryName;       //!< The Name of the TrackingGeometry
       ToolHandleArray<IGeometryProcessor>          m_trackingGeometryProcessors; //!< Tool to write out a Display format for external viewers
                                    
   };
}

#endif
