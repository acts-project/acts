///////////////////////////////////////////////////////////////////
// EnergyLossSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H
#define ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "CoreInterfaces/IRandomNumberSvc.h"
#include "Algebra/AlgebraDefinitions.h"
// EventData module
#include "EventDataUtils/ParticleHypothesis.h"
// Extrapolation module
#include "ExtrapolationUtils/MaterialInteraction.h"
// Geometry module
#include "Material/MaterialProperties.h"
// Fatras module
#include "FatrasInterfaces/IEnergyLossSampler.h"

namespace Acts{
  
  class EnergyLoss;

  /** @class EnergyLossSampler
   * 
   * Sampler for a eloss of a track
   * It uses the Bethe-Bloch calculation
   * it extends the IEnergyLossSampler interface 
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   */

  class EnergyLossSampler : public AlgToolBase, virtual public IEnergyLossSampler {
   
  public:
   
    /** Constructor with AlgTool parameters */
    EnergyLossSampler( const std::string&, const std::string&, const IInterface* );
   
    /** Destructor */
    ~EnergyLossSampler();
    
    /** AlgTool interface methods */
    static const InterfaceID& interfaceID() { return IID_IEnergyLossSampler; }
       
    /** AlgTool initialise method */
    StatusCode initialize() final;
   
    /** AlgTool finalise method */
    StatusCode finalize() final;
    
    /** IEnergyLossSampler public method to compute dEdX */
    double dEdX( const Acts::MaterialProperties& materialProperties,
		 double momentum,
		 Acts::ParticleHypothesis particleHypothesis = Acts::pion ) const final;
   
    /** IEnergyLossSampler public method to compute the mean and variance of the energy loss */
    Acts::EnergyLoss* energyLoss( const Acts::MaterialProperties& mat,
				  double momentum,
				  double pathcorrection,
				  Acts::PropDirection dir=Acts::alongMomentum,
				  Acts::ParticleHypothesis particle=Acts::pion) const final;
   
       
  private:
    
    double dEdXBetheBloch(const MaterialProperties& mat,
			  double gamma,
			  double beta,
			  ParticleHypothesis particle=pion) const;
			  
    double dEdXBetheHeitler(const MaterialProperties& mat,
			    double initialE,
			    ParticleHypothesis particle=pion) const;
   
    /** Random Generator service  */
    ServiceHandle<Acts::IRandomNumberSvc>      m_rndGenSvc;  
    
    /** the formulas for energy loss evaluation */
    MaterialInteraction                        m_interactionFormulae;     
    
    /** struct of Particle masses  */
    static ParticleMasses                      s_particleMasses; 
    
    /** KOverA factor in Bethe-Bloch equation [MeV*cm2/gram] */
    static double                              s_ka_BetheBloch;          
};

} 
#endif //  ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H