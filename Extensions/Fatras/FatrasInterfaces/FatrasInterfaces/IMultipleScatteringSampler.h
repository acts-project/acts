///////////////////////////////////////////////////////////////////
// IMultipleScatteringSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IMULTIPLESCATTERINGSAMPLER_H
#define ACTS_FATRASINTERFACES_IMULTIPLESCATTERINGSAMPLER_H 1

// Gaudi 
#include "GaudiKernel/IAlgTool.h"
// EventData module
#include "EventDataUtils/ParticleHypothesis.h"

namespace Acts {
  
  class MaterialProperties;
 
  /** Interface ID for IMultipleScatteringSampler*/  
  static const InterfaceID IID_IMultipleScatteringSampler("IMultipleScatteringSampler", 1, 0);
  
  /** @class IMultipleScatteringSampler
   * 
   * Interface class IMultipleScatteringSampler
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
  
  class IMultipleScatteringSampler : virtual public IAlgTool {

  public:
    /**Virtual destructor*/
    virtual ~IMultipleScatteringSampler(){}
    
    /** AlgTool and IAlgTool interface methods */
    static const InterfaceID& interfaceID() { return IID_IMultipleScatteringSampler; }
    
    virtual double simTheta (const MaterialProperties& mat,
			     double momentum,
			     double pathcorrection,
			     Acts::ParticleHypothesis particle = Acts::pion) const = 0;    
  };

}

#endif // ACTS_FATRASINTERFACES_IMULTIPLESCATTERINGSAMPLER_H