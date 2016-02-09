///////////////////////////////////////////////////////////////////
// PropagationEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_TRKEXENGINE_PROPAGATIONENGINE_H
#define ATS_TRKEXENGINE_PROPAGATIONENGINE_H

// Ats
#include "CoreInterfaces/ServiceBase.h"
// Gaudi
#include "GaudiKernel/ToolHandle.h"
// Trk
#include "ExtrapolationInterfaces/IPropagationEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "EventDataUtils/PropDirection.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

#include "ExtrapolationInterfaces/IPropagator.h"
 
namespace Ats {
  
  class Surface;
    
  /** @class PropagationEngine 
  
      Wrapper around the IPropagator interface,
      the neutral propagation is solved by using the surface::straightLineIntersection
      mehtod.
  
      @author Andreas Salzburger -at - cern.ch 
  */
  class PropagationEngine : public Ats::ServiceBase, virtual public IPropagationEngine {

    public:

      /** Constructor */
      PropagationEngine(const std::string& name, ISvcLocator* svc);
      
      /** Destructor */
      ~PropagationEngine();

      /** AlgTool initialize method */
      StatusCode initialize() final;
      
      /** AlgTool finalize method */
      StatusCode finalize() final;

      /** Avoid shaddowing */
      using IPropagationEngine::propagate;

      /** resolve the boundary situation - for charged particles */
      ExtrapolationCode propagate(ExCellCharged& ecCell,
                                  const Surface& sf,
                                  PropDirection dir=alongMomentum,
                                  const BoundaryCheck& bcheck = true,
                                  bool returnCurvilinear = true) const final;  

      /** resolve the boundary situation - for neutral particles */
      ExtrapolationCode propagate(ExCellNeutral& enCell,
                                  const Surface& sf,
                                  PropDirection dir=alongMomentum,
                                  const BoundaryCheck& bcheck = true,
                                  bool returnCurvilinear = true) const final;
       
    protected:
        
      ToolHandle<IPropagator>             m_propagator;
      double                              m_pathLimitTolerance;


  };
      

} // end of namespace

#endif // ATS_TRKEXINTERFACES_PROPAGATIONENGINE_H

