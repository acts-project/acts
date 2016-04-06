///////////////////////////////////////////////////////////////////
// IEnergyLossSampler.h, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IENERGYLOSSSAMPLER_H
#define ACTS_FATRASINTERFACES_IENERGYLOSSSAMPLER_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
//EventData module
#include "EventDataUtils/ParticleHypothesis.h"
#include "EventDataUtils/PropDirection.h"

namespace Acts {
  
  class MaterialProperties;
  class EnergyLoss;

  /** Interface ID for IEnergyLossSampler*/  
  static const InterfaceID IID_IEnergyLossSampler("IEnergyLossSampler", 1, 0);
  
  /** @class IEnergyLossSampler
   *
   * Interface class IEnergyLossSampler
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   */

  class IEnergyLossSampler : virtual public IAlgTool {
    
  public:
    /**Virtual destructor*/
    virtual ~IEnergyLossSampler(){}
    
    /** AlgTool and IAlgTool interface methods */
    static const InterfaceID& interfaceID() { return IID_IEnergyLossSampler; } 
    
    /** deltaE calculation
     * using dEdX and integrating along pathlength,
     * assuming constant dEdX during for the path.
     * - The sign depends on the given propagation direction 
     * - Units: [MeV]
    */
    virtual Acts::EnergyLoss* energyLoss(const Acts::MaterialProperties& mat,
				 	 double momentum,
				 	 double pathcorrection,
				 	 Acts::PropDirection dir=Acts::alongMomentum,
				 	 Acts::ParticleHypothesis particle=Acts::pion) const = 0;  
					 
    /** dEdX calculation when providing MaterialProperties,
     * a momentum, and a ParicleHypothesis. 
     * - Units: [Mev/mm]
     */
    virtual double dEdX(const Acts::MaterialProperties& mat,
			double momentum,
			Acts::ParticleHypothesis particle=Acts::pion) const = 0;
   
  };

} // end of namespace

#endif // ACTS_FATRASINTERFACES_IENERGYLOSSSAMPLER_H