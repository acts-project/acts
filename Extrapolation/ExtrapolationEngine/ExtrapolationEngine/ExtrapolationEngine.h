///////////////////////////////////////////////////////////////////
// ExtrapolationEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXTRAPOLATIONENGINE_EXTRAPOLATIONENGINE_H
#define ATS_EXTRAPOLATIONENGINE_EXTRAPOLATIONENGINE_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/GaudiException.h"
// Core module
#include "CoreInterfaces/ServiceBase.h"
// Extrapolation module
#include "ExtrapolationInterfaces/IExtrapolationEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationInterfaces/IPropagationEngine.h"
#include "ExtrapolationInterfaces/INavigationEngine.h"
// Geomtetry module
#include "GeometryInterfaces/ITrackingGeometrySvc.h"
// EventData module
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

namespace Ats {
  
  class TrackingGeometry;       

  /** @class ExtrapolationEngine 
      
      Master extrapolation engine for extrapolation through the TrackingGeometry,
      it delegates the extrapolation to optimised engines, handing over the ExtrapolationCell
      as internal cache.
  
      There are identical interfaces for charged and neutral track parameters.
      Providing a destination surface is optional, if no destination surface is given the extrapolation 
      process can be stopped by other directives, e.g. stopping at a certain path limit, material limit
      or with a change of detector signature.
  
      @author Andreas.Salzburger -at- cern.ch 
  */

class ExtrapolationEngine : public ServiceBase, virtual public IExtrapolationEngine {
      
      friend class NavigationInitTest;
      
      public:
        /** Constructor */
        ExtrapolationEngine(const std::string& name, ISvcLocator* svc);
        
        /** Destructor */
        ~ExtrapolationEngine();

        /** AlgTool initialize method */
        StatusCode initialize() final;
        
        /** AlgTool finalize method */
        StatusCode finalize() final;
        
        using IExtrapolationEngine::extrapolate;
        
        /** charged extrapolation - public interface */
        ExtrapolationCode extrapolate(ExCellCharged& ecCharged,
                                      const Surface* sf = 0,
                                      const BoundaryCheck& bcheck = true) const final;

        /** neutral extrapolation - public interface */
        ExtrapolationCode extrapolate(ExCellNeutral& ecNeutral,
                                      const Surface* sf = 0,
                                      const BoundaryCheck& bcheck = true) const final;
                         
                         
        /** define for which GeometrySignature this extrapolator is valid - this is GLOBAL */
        GeometryType geometryType() const final;                           
                         
     private:
        /** main loop extrapolation method */
        template <class T> ExtrapolationCode extrapolateT(ExtrapolationCell<T>& eCell,
                                                          const Surface* sf = 0,
                                                          PropDirection dir=alongMomentum,
                                                          const BoundaryCheck& bcheck = true) const;
            
        /** initialization method */                                      
        template <class T>  ExtrapolationCode initNavigation(ExtrapolationCell<T>& eCell,
                                                             const Surface* sf = 0,
                                                             PropDirection dir=alongMomentum) const throw (GaudiException);
                
                
        //!< retrieve TrackingGeometry
        StatusCode  updateTrackingGeometry() const; 

        //!< return and retrieve
        const TrackingGeometry& trackingGeometry() const throw (GaudiException);

        mutable const TrackingGeometry*                     m_trackingGeometry;          //!< the tracking geometry owned by the navigator
        ServiceHandle<ITrackingGeometrySvc>                 m_trackingGeometrySvc;       //!< ServiceHandle to the TrackingGeometrySvc
        std::string                                         m_trackingGeometryName;      //!< Name of the TrackingGeometry as given in Detector Store
        
        //!< the tool handle array for static / dense / detached
        ServiceHandleArray<IExtrapolationEngine>            m_extrapolationEngines;      //!< the extrapolation engines for retrieval
        ServiceHandle<IPropagationEngine>                   m_propagationEngine;         //!< the used propagation engine for navigation initialization
        std::vector<const IExtrapolationEngine*>            m_eeAccessor;                //!< the extrapolation engines for 

        ServiceHandle<INavigationEngine>                    m_navigationEngine;          //!< access to tracking geometry (unique?)

        //!< forces a global search for the initialization, allows to switch TrackingGeometries in one job
        bool                                                m_forceSearchInit; 
    
    };

  inline GeometryType  ExtrapolationEngine::geometryType() const 
      { return Ats::Master; }


  inline const Ats::TrackingGeometry& ExtrapolationEngine::trackingGeometry() const throw (GaudiException) {
      if (!m_trackingGeometry && updateTrackingGeometry().isFailure()){
          EX_MSG_FATAL("", "updateGeo", "", "Could not load TrackingGeometry with name '" << m_trackingGeometryName << "'. Aborting." );
          throw GaudiException("ExtrapolationEngine", "Problem with TrackingGeometry loading.", StatusCode::FAILURE);
      }
      return (*m_trackingGeometry);
  }
   
   

} // end of namespace

//!< define the templated function    
#include "ExtrapolationEngine.icc"  

#endif // ATS_EXTRAPOLATIONENGINE_EXTRAPOLATIONENGINE_H 

