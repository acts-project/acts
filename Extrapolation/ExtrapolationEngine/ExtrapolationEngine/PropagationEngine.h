///////////////////////////////////////////////////////////////////
// PropagationEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXTRAPOLATIONENGINE_PROPAGATIONENGINE_H
#define ATS_EXTRAPOLATIONENGINE_PROPAGATIONENGINE_H 1

// Gaudi
#include "GaudiKernel/ToolHandle.h"
// Core module
#include "CoreInterfaces/ServiceBase.h"
// Extrapolation module
#include "ExtrapolationInterfaces/IPropagationEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationInterfaces/IPropagator.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
// EventData module
#include "EventDataUtils/PropDirection.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"
 
namespace Ats {
  
  class Surface;

  /** @class PropagationEngine

      Wrapper around the IPropagator interface,
      the neutral propagation is solved by using the surface::straightLineIntersection
      mehtod.

      @author Andreas Salzburger -at - cern.ch
  */
  class PropagationEngine : public ServiceBase, virtual public IPropagationEngine {

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

      ToolHandle<IPropagator>     m_propagator;             // handle to propagor
      double                      m_pathLimitTolerance;     // path limit tolerance 


  };


} // end of namespace

#endif // ATS_EXTRAPOLATIONENGINE_PROPAGATIONENGINE_H

