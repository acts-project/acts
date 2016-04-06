///////////////////////////////////////////////////////////////////
// IPhotonConversionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IPHOTONCONVERSIONSAMPLER_H
#define ACTS_FATRASINTERFACES_IPHOTONCONVERSIONSAMPLER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"

namespace Acts {
  
  class InteractionVertex;

  /** Interface ID for IPhotonConversionSampler */
  static const InterfaceID IID_IPhotonConversionSampler("IPhotonConversionSampler", 1, 0);

  /**
   * @class IPhotonConversionSampler
   * Interface definition for the handling of photon conversion,
   * to be used by the FatrasMaterialEffectsEngine
   * 
   * @author Sarka Todorova <Sarka.Todorova@cern.ch>
   * @author Noemi Calace   <Noemi.Calace@cern.ch>
   * 
  */
    class IPhotonConversionSampler : virtual public IAlgTool {

  public:
    
    /** Virtual destructor */    
    virtual ~IPhotonConversionSampler() {}
    
    /** AlgTool interface methods */
    static const InterfaceID& interfaceID() { return IID_IPhotonConversionSampler; }
    
    /** interface for processing of the presampled conversion on layer*/
    virtual std::vector<Acts::InteractionVertex> doConversion(double time, 
							      const Vector3D& position , 
							      const Vector3D& momentum) const = 0;

  };

}

#endif // ACTS_FATRASINTERFACES_IPHOTONCONVERSIONSAMPLER_H