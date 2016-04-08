///////////////////////////////////////////////////////////////////
// IPropagator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_IPROPAGATOR_H
#define ACTS_EXTRAPOLATIONINTERFACES_IPROPAGATOR_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// EventData module
#include "EventDataUtils/PropDirection.h"
#include "EventDataUtils/ParticleHypothesis.h"
#include "TrackParameters/TrackParameters.h"
// Geometry module
#include "Surfaces/BoundaryCheck.h"

// STL
#include <utility>

namespace Acts {

  class TransportJacobian;
  class Surface;
  class MagneticFieldProperties;

  /** Interface ID for IPropagator*/
  static const InterfaceID IID_IPropagator("IPropagator", 1, 0);

  /** @class IPropagator

      Interface class IPropagator for simple factory type propagator
  
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

#endif // ACTS_EXTRAPOLATIONINTERFACES_IPROPAGATOR_H
