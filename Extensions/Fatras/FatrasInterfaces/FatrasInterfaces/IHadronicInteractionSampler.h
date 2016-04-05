///////////////////////////////////////////////////////////////////
// IHadronicInteractionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IHADRONICINTERACTIONSAMPLER_H
#define ACTS_FATRASINTERFACES_IHADRONICINTERACTIONSAMPLER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"
// EventData module
#include "EventDataUtils/ParticleHypothesis.h"

namespace Acts {
  
  class InteractionVertex;

  /** Interface ID for IHadronicInteractionSampler */
  static const InterfaceID IID_IHadronicInteractionSampler("IHadronicInteractionSampler", 1, 0);

  /**
   * @class IHadronicInteractionSampler
   * Interface definition for the handling of nuclear/hadronic interactions,
   * to be used by the MC based material effects updater
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace   <Noemi.Calace@cern.ch>
   * 
  */
    class IHadronicInteractionSampler : virtual public IAlgTool {

  public:
    
    /** Virtual destructor */    
    virtual ~IHadronicInteractionSampler() {}
    
    /** AlgTool interface methods */
    static const InterfaceID& interfaceID() { return IID_IHadronicInteractionSampler; }
    
    /** interface for processing of the presampled nuclear interactions on layer*/
    virtual std::vector<Acts::InteractionVertex> doHadronicInteraction(double time, 
								       const Acts::Vector3D& position, 
								       const Acts::Vector3D& momentum, 
								       Acts::ParticleHypothesis particle=Acts::pion) const = 0;

  };

}

#endif // ACTS_FATRASINTERFACES_IHADRONICINTERACTIONSAMPLER_H