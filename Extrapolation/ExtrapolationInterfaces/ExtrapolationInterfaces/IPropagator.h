///////////////////////////////////////////////////////////////////
// IPropagator.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_TRKEXINTERFACES_IPROPAGATOR_H
#define ATS_TRKEXINTERFACES_IPROPAGATOR_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Trk
#include "EventDataUtils/PropDirection.h"
#include "EventDataUtils/ParticleHypothesis.h"
#include "TrackParameters/TrackParameters.h"
#include "Surfaces/BoundaryCheck.h"

// STL
#include <utility>

namespace Ats {

  class TransportJacobian;
  class Surface;

  /** Interface ID for IPropagators*/
  static const InterfaceID IID_IPropagator("IPropagator", 1, 0);

  /** @class IPropagator

      Interface class IPropagators
      It inherits from IAlgTool.

    */
  class IPropagator : virtual public IAlgTool {
      public:

       /**Virtual destructor*/
       virtual ~IPropagator(){}

       /** AlgTool and IAlgTool interface methods */
       static const InterfaceID& interfaceID() { return IID_IPropagator; }

       /* The propagation method including the return of the TransportJacobian matrix. */
       virtual const TrackParameters* propagate( const TrackParameters& parm,
                                                 const Surface& sf,
                                                 PropDirection dir,
                                                 const BoundaryCheck& bcheck,
                                                 const MagneticFieldProperties& mprop,
                                                 TransportJacobian*&,
                                                 double& pathLength,
						                         bool returnCurv = false) const = 0;



 };

} // end of namespace

#endif // TRKEXINTERFACES_PROPAGATOR_H
